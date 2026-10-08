"""Recheck previous constructions against the current geometry and rules.

Warm starts are optional lower-bound proposals, not generator metadata or
optimality certificates. No previous audit or old graph node ID is trusted.
"""
from dataclasses import asdict
import json
import math
from pathlib import Path
import re
import time

import numpy as np
import shapely
from scipy.spatial import cKDTree
from shapely.geometry import LineString
from shapely.ops import unary_union

from curved_centerline import CurveSettings
from island_router import _attachment,betti,disk
from ordered_pad_fanout import rebuild_ordered_outer
from electrode_region import region_for_navigation
from solver_time_policy import deadline_after, deadline_expired


def discover_incumbents(runs,input_sha256,support_layer,rules):
    reports=[]
    for status in Path(runs).glob('*/status.json'):
        try:
            with status.open(encoding='utf-8') as stream:header=stream.read(4096)
            if (not re.search(r'"status"\s*:\s*"complete"',header) or
                    not re.search(r'"input_sha256"\s*:\s*"'+re.escape(input_sha256)+'"',header)):
                continue
            path=status.parent/'summary.json'
            report=json.loads(path.read_text(encoding='utf-8'))
            if (report.get('method')!='four_side_pads' or
                    report['selected_support_layer']!=list(support_layer) or
                    report['routing']['rules']!=asdict(rules) or
                    report['routing'].get('outer_exit_policy',{}).get('revision') not in
                        ('global_radial_outer_ports_v1','directional_convex_envelope_outer_ports_v2') or
                    not report['capacity_interval'].get('integer_polygon_lower_verified')):
                continue
            reports.append({'job_id':status.parent.name,'report':report})
        except (OSError,ValueError,KeyError,TypeError):continue
    return sorted(reports,key=lambda r:(-r['report']['routing']['retained_routes'],
        r['report']['routing'].get('center_preference',{}).get('selected',{}).get('mean_radius_um') or math.inf))


def recertify_incumbent(records,nav,rules,settings,pads,frame,policy,*,minimum_count=0,
                       progress=None,time_limit_s=None):
    start=time.perf_counter();deadline=deadline_after(time_limit_s,started=start)
    if not records or not nav.outlets:return [],frame,{'status':'no_previous_complete_construction'}
    outlets=list(nav.outlets);tree=cKDTree([nav.graph.nodes[n]['xy_um'] for n in outlets])
    radius=rules.electrode_diameter_um/2+rules.margin_um+settings.numeric_guard_um
    curves=CurveSettings();support=nav.support;shapely.prepare(support)
    electrode_region=region_for_navigation(nav,rules)
    tried=[];chosen=[];chosen_frame=frame;chosen_record=None
    for record in records:
        report=record['report'];routing=report['routing']
        if routing['retained_routes']<max(minimum_count,len(chosen)):break
        if deadline_expired(deadline):break
        if (report['geometry']['sha256']!=nav.summary['sha256'] or routing['rules']!=asdict(rules) or
                report['selected_support_layer']!=nav.summary['support_layer']):continue
        if progress:progress(f'按当前几何重新认证已有完整方案：{routing["retained_routes"]} 条',.57)
        pool=[];failure=None
        for index,old in enumerate(routing['routes']):
            if deadline_expired(deadline):failure='incumbent_recertification_budget';break
            xy=np.asarray(old['source_um']);line=LineString(old['points_um'])
            if not electrode_region.contains(xy,export_safe=True):
                failure='incumbent_outside_electrode_region';break
            island,reason=_attachment(support,xy,radius,settings.pad_gap_um)
            if island is None:failure=reason;break
            if (np.linalg.norm(np.asarray(line.coords[0])-xy)>policy.grid_um or
                    any(np.linalg.norm(xy-c['source_um'])<rules.minimum_center_spacing_um+settings.numeric_guard_um
                        for c in pool)):
                failure='incumbent_source_spacing_or_endpoint';break
            wire=line.buffer((rules.wire_width_um/2+curves.max_chord_error_um)/math.cos(math.pi/64),quad_segs=16)
            if not support.covers(wire.buffer(rules.margin_um+settings.numeric_guard_um,quad_segs=32)):
                failure='incumbent_original_support_margin';break
            electrode=disk(xy,rules.electrode_diameter_um/2);metal=wire.union(electrode)
            if any(metal.distance(c['metal'])<rules.spacing_um+settings.numeric_guard_um or
                   island.distance(c['island'])<settings.pad_gap_um for c in pool):
                failure='incumbent_metal_or_island_spacing';break
            _,nearest=tree.query(np.asarray(old.get('declared_terminal_um',old['outlet_um'])))
            outlet=outlets[int(nearest)]
            pool.append({**old,'source_node':('recertified_incumbent',record['job_id'],index),
                'source_um':xy.tolist(),'outlet_node':outlet,'navigation_view':0,
                'declared_terminal_um':list(nav.graph.nodes[outlet]['xy_um']),
                'points_um':np.asarray(line.coords),'island':island,'electrode':electrode,
                'wire':wire,'central_wire':wire,'metal':metal,
                'route_shape_strategy':'recertified_previous_complete_construction'})
        if failure is None and betti(support.union(unary_union([c['island'] for c in pool])))!=betti(support):
            failure='incumbent_joint_island_topology'
        if failure is not None:
            tried.append({'job_id':record['job_id'],'reason':failure,'checked_before_failure':len(pool)});continue
        trial_frame=frame
        rebuilt,diagnostic=rebuild_ordered_outer(nav,pool,rules,settings,pads,trial_frame,policy)
        if rebuilt is None and diagnostic.get('reason')=='ordered_bank_has_too_few_slots':
            from pad_router import make_pad_frame
            trial_frame=make_pad_frame(support,pads,minimum_side_um=frame['square_side_um'],spacing_um=rules.spacing_um,
                margin_um=rules.margin_um,wire_width_um=rules.wire_width_um,
                exit_policy=policy,required_slots_per_side=diagnostic['required'])
            rebuilt,diagnostic=rebuild_ordered_outer(nav,pool,rules,settings,pads,trial_frame,policy)
        if rebuilt is None:
            tried.append({'job_id':record['job_id'],'reason':diagnostic});continue
        cost=sum(float(np.sum((np.asarray(c['source_um'])-policy.origin)**2)) for c in rebuilt)
        oldcost=sum(float(np.sum((np.asarray(c['source_um'])-policy.origin)**2)) for c in chosen)
        if len(rebuilt)>len(chosen) or (len(rebuilt)==len(chosen) and cost<oldcost):
            chosen=rebuilt;chosen_frame=trial_frame;chosen_record=record['job_id']
        tried.append({'job_id':record['job_id'],'certified_complete_count':len(rebuilt)})
    return chosen,chosen_frame,{'status':'previous_constructions_rechecked_in_current_geometry',
        'selected_previous_job':chosen_record,'certified_complete_count':len(chosen),
        'records_checked':tried,'seconds':time.perf_counter()-start,
        'wall_clock_limit_enabled':time_limit_s is not None,'time_limit_s':time_limit_s,
        'uses_previous_generator_graph':False,'optimality_proven':False}
