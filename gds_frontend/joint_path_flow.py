"""Joint paths on a geometry-only navigation graph, certified in real metal.

Node-disjoint flow is a conservative proposal generator. Unit graph capacity
is never a fabrication capacity or an upper bound: later residual searches
may add tracks. The common source and sink have no prescribed net pairing.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import replace
import math
import time

import numpy as np
import shapely
from scipy.optimize import linprog
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import maximum_flow
from shapely.geometry import GeometryCollection,LineString,Point,Polygon
from shapely.ops import unary_union,substring

from curved_centerline import CurveSettings,UncertifiableCurve,smooth_centerline
from island_router import _attachment,_path_line,_shorten,_select,betti,disk
from navigation_simplification import check_deadline,SearchDeadlineExceeded
from ordered_pad_fanout import rebuild_ordered_outer
from electrode_region import region_for_navigation
from solver_time_policy import deadline_after, deadline_expired, solver_options


def joint_path_routes(nav,rules,settings,pad_settings,frame,policy,progress=None,*,time_limit_s=None):
    start=time.perf_counter();deadline=deadline_after(time_limit_s,started=start)
    source_deadline=deadline_after(time_limit_s/3 if time_limit_s is not None else None,started=start)
    geometry_deadline=deadline_after(3*time_limit_s if time_limit_s is not None else None,started=start)
    timings={}
    if not hasattr(nav,'graph') or not hasattr(nav,'candidates'):
        return [],{'status':'navigation_has_no_joint_flow_graph','global_maximum_proven':False}
    reject=Counter();source_points={};cells=defaultdict(list)
    gap=max(rules.minimum_center_spacing_um,rules.electrode_diameter_um+2*rules.margin_um+settings.pad_gap_um)
    origin=policy.origin;support=nav.support;shapely.prepare(support)
    electrode_region=region_for_navigation(nav,rules)
    grid=policy.grid_um;island_radius=rules.electrode_diameter_um/2+rules.margin_um+settings.numeric_guard_um
    candidate_ids=sorted(nav.candidates,key=lambda n:float(np.sum((np.asarray(nav.graph.nodes[n]['xy_um'])-origin)**2)))
    # Source spacing is a constructive independent set, not an optimal packing.
    for node in candidate_ids:
        if deadline_expired(source_deadline):break
        xy=np.asarray(nav.graph.nodes[node]['xy_um']);key=tuple(np.floor(xy/gap).astype(int))
        if not electrode_region.contains(xy,export_safe=True):
            reject['outside_electrode_region']+=1;continue
        if any(np.linalg.norm(xy-other)<gap+settings.numeric_guard_um for i in range(-1,2) for j in range(-1,2)
               for other in cells[(key[0]+i,key[1]+j)]):continue
        island,reason=_attachment(support,xy,island_radius,settings.pad_gap_um)
        if island is None:reject[reason]+=1;continue
        source_points[node]=(xy,island);cells[key].append(xy)
    timings['sources_seconds']=time.perf_counter()-start
    nodes=list(nav.graph);indices={n:i for i,n in enumerate(nodes)};count=len(nodes)
    source=2*count;sink=source+1;size=sink+1
    rows=[];cols=[];cost=[]
    def add(u,v,weight=0):rows.append(u);cols.append(v);cost.append(weight)
    for n in nodes:add(2*indices[n],2*indices[n]+1)
    max_radius=max(1,policy.radius_ticks*grid)
    for u,v,data in nav.graph.edges(data=True):
        weight=data['weight']/max_radius*1e-3
        add(2*indices[u]+1,2*indices[v],weight);add(2*indices[v]+1,2*indices[u],weight)
    for n,(xy,island) in source_points.items():
        add(source,2*indices[n],float(np.sum((xy-origin)**2))/max_radius**2)
    for n in nav.outlets:add(2*indices[n]+1,sink)
    if not nav.outlets or not source_points:return [],{'status':'no_joint_flow_terminals_or_sources','seconds':time.perf_counter()-start}
    rows=np.asarray(rows);cols=np.asarray(cols);m=len(rows)
    capacity=coo_matrix((np.ones(m,dtype=np.int32),(rows,cols)),shape=(size,size)).tocsr()
    if progress:progress('联合分配电极、内部路径与外边界出口',.43)
    maximum=maximum_flow(capacity,source,sink)
    flow_count=int(maximum.flow_value)
    check_deadline(deadline,'joint maximum flow')
    timings['flow_seconds']=time.perf_counter()-start
    chosen_flow=maximum.flow.tocsr();mincost_status='maximum_flow_incumbent'
    # Two stages implement lexicographic count/centrality on this guide only.
    incidence=coo_matrix((np.r_[np.ones(m),-np.ones(m)],
                         (np.r_[rows,cols],np.r_[np.arange(m),np.arange(m)])),shape=(size,m)).tocsr()
    supply=np.zeros(size);supply[source]=flow_count;supply[sink]=-flow_count
    remaining=(max(.01,min(time_limit_s/3,deadline-time.perf_counter()))
               if deadline is not None else None)
    result=linprog(cost,A_eq=incidence,b_eq=supply,bounds=(0,1),method='highs',
                   options=solver_options(remaining))
    if result.x is not None and np.max(np.abs(result.x-np.rint(result.x)))<1e-6:
        integral=np.rint(result.x).astype(np.int32)
        if np.max(np.abs(incidence@integral-supply))<1e-6:
            chosen_flow=coo_matrix((integral,(rows,cols)),shape=(size,size)).tocsr()
            mincost_status='optimal_integral_min_cost_flow' if result.success else 'feasible_integral_flow'
    timings['cost_seconds']=time.perf_counter()-start
    # Decompose the positive integral flow into electrode-to-portal paths.
    paths=[]
    for n in source_points:
        if chosen_flow[source,2*indices[n]]<=0:continue
        node=2*indices[n];path=[];seen=set()
        while node!=sink:
            if node in seen:raise RuntimeError('Joint integral flow contains a source-connected cycle')
            seen.add(node)
            if node<2*count and node%2==0:path.append(nodes[node//2])
            row=chosen_flow.getrow(node)
            targets=row.indices[row.data>0]
            if len(targets)!=1:raise RuntimeError('Joint integral flow violates node capacity')
            node=int(targets[0])
        paths.append((n,path))
    curves=CurveSettings();proposals=[];greedy=[];greedy_ids=[];geometry_topology=betti(support)
    original_anchor_count=0;suffix_contact_count=0
    # A wire suffix stays in its certified corridor. Electrode occupancy must
    # nevertheless be certified again at each contact. Keep the original
    # central option; alternatives are part of the same joint finite model.
    contact_tail_lengths=sorted({
        rules.electrode_diameter_um/2+island_radius+rules.spacing_um+settings.numeric_guard_um,
        max(rules.minimum_center_spacing_um,gap)})
    # Geometry acceptance, not the guide flow, determines the actual count.
    for ordinal,(source_node,path) in enumerate(sorted(paths,key=lambda p:float(np.linalg.norm(source_points[p[0]][0]-origin))),1):
        if deadline_expired(geometry_deadline):
            reject['construction_budget_reached']+=len(paths)-ordinal+1;break
        if progress and ordinal%10==0:
            progress(f'认证联合路径：{ordinal}/{len(paths)}，已有 {len(proposals)} 条几何提案',.48)
        xy,island=source_points[source_node];outlet=path[-1]
        raw,corridors=_path_line(nav.graph,path)
        try:
            line=_shorten(raw.coords,nav.navigation_domain,rules.wire_width_um/2)
            line,curve=smooth_centerline(line.coords,nav.center_domain,curves)
        except UncertifiableCurve:reject['uncertifiable_joint_curve']+=1;continue
        wire=line.buffer((rules.wire_width_um/2+curves.max_chord_error_um)/math.cos(math.pi/64),quad_segs=16)
        if not support.covers(wire.buffer(rules.margin_um+settings.numeric_guard_um,quad_segs=32)):
            reject['joint_wire_margin']+=1;continue
        electrode=disk(xy,rules.electrode_diameter_um/2);metal=wire.union(electrode)
        proposal={'source_node':(0,source_node),'source_um':xy.tolist(),
            'outlet_node':outlet,'outlet_um':list(line.coords[-1]),'declared_terminal_um':list(nav.graph.nodes[outlet]['xy_um']),
            'navigation_view':0,'terminal_kind':'external_boundary','points_um':np.asarray(line.coords),
            'corridor_ids':corridors,'lane_offset_um':0,'line_length_um':line.length,'depth_um':line.length,
            'curve':curve,'route_shape_strategy':'joint_integral_geometry_only_flow',
            'island':island,'electrode':electrode,'wire':wire,'central_wire':wire,'metal':metal}
        if all(metal.distance(c['metal'])>=rules.spacing_um+settings.numeric_guard_um for c in greedy):
            greedy_ids.append(len(proposals));greedy.append(proposal)
        proposals.append(proposal)
        original_anchor_count+=1
        for tail_length in contact_tail_lengths:
            if tail_length>=line.length-settings.numeric_guard_um:continue
            trimmed=substring(line,line.length-tail_length,line.length)
            contact=np.asarray(trimmed.coords[0])
            if not electrode_region.contains(contact,export_safe=True):
                reject['suffix_outside_electrode_region']+=1;continue
            new_island,reason=_attachment(support,contact,island_radius,settings.pad_gap_um)
            if new_island is None:reject['suffix_contact_'+reason]+=1;continue
            try:trimmed,new_curve=smooth_centerline(trimmed.coords,nav.center_domain,curves)
            except UncertifiableCurve:reject['uncertifiable_suffix_curve']+=1;continue
            new_wire=trimmed.buffer((rules.wire_width_um/2+curves.max_chord_error_um)/math.cos(math.pi/64),quad_segs=16)
            if not support.covers(new_wire.buffer(rules.margin_um+settings.numeric_guard_um,quad_segs=32)):
                reject['suffix_wire_margin']+=1;continue
            new_electrode=disk(contact,rules.electrode_diameter_um/2)
            proposals.append({**proposal,'source_um':contact.tolist(),'island':new_island,
                'electrode':new_electrode,'wire':new_wire,'central_wire':new_wire,
                'metal':new_wire.union(new_electrode),'points_um':np.asarray(trimmed.coords),
                'outlet_um':list(trimmed.coords[-1]),'line_length_um':trimmed.length,'depth_um':trimmed.length,
                'curve':new_curve,'route_shape_strategy':'joint_integral_flow_suffix_contact',
                'contact_tail_length_um':tail_length,
                'contact_variant_scope':'one electrode per shared guide-route family; all spacing and island rules remain joint constraints'})
            suffix_contact_count+=1
    # Greedy deletion can waste independent routes in a conflict graph. Solve
    # cardinality first over all certified proposals, retaining the greedy
    # subset as a nondegrading incumbent. This is still only a finite library.
    accepted,finite_selection=_select(proposals,replace(rules,spacing_um=rules.spacing_um+settings.numeric_guard_um),
        settings,progress or (lambda *_:None),baseline_candidate_ids=greedy_ids,center_origin=origin)
    if len(accepted)<len(greedy):accepted=greedy
    reject['joint_metal_spacing']=len(proposals)-len(accepted)
    while accepted and betti(support.union(unary_union([c['island'] for c in accepted])))!=geometry_topology:
        reject['joint_islands_topology']+=1
        accepted.pop(max(range(len(accepted)),key=lambda i:np.linalg.norm(np.asarray(accepted[i]['source_um'])-origin)))
    # One geometric failure cannot erase the rest of a complete flow proposal.
    complete=[];outer_reject=Counter();pool=list(accepted);outer_cache={}
    timings['internal_geometry_seconds']=time.perf_counter()-start
    while pool:
        if progress:progress(f'联合安排外部连续 Pad：{len(pool)} 条内部路径',.57)
        rebuilt,diagnostic=rebuild_ordered_outer(nav,pool,rules,settings,pad_settings,frame,policy,
                                                prepared_cache=outer_cache)
        if rebuilt is None and diagnostic.get('reason')=='ordered_bank_has_too_few_slots':
            from pad_router import make_pad_frame
            required=diagnostic['required']
            frame=make_pad_frame(support,pad_settings,minimum_side_um=frame['square_side_um'],spacing_um=rules.spacing_um,
                                 margin_um=rules.margin_um,wire_width_um=rules.wire_width_um,
                                 exit_policy=policy,required_slots_per_side=required)
            rebuilt,diagnostic=rebuild_ordered_outer(nav,pool,rules,settings,pad_settings,frame,policy,
                                                    prepared_cache=outer_cache)
        if rebuilt is not None:
            complete=rebuilt;break
        outer_reject[diagnostic.get('reason','unknown')]+=1
        # Remove only the implicated proposal; keep the rest and jointly
        # rematch. This does not treat a failed finite proposal as impossible.
        assignment,_=__import__('ordered_pad_fanout').ordered_pad_assignment(
            [nav.graph.nodes[c['outlet_node']]['xy_um'] for c in pool],frame)
        bad_pads=set(diagnostic.get('pad_ids',()))
        if diagnostic.get('pad_id'):bad_pads.add(diagnostic['pad_id'])
        implicated=diagnostic.get('route_indices') or [i for i,c in enumerate(pool) if c['outlet_node']==diagnostic.get('outlet_node') or
                    (assignment and assignment[i][0]['pad_id'] in bad_pads)]
        remove=max(implicated or range(len(pool)),key=lambda i:np.linalg.norm(np.asarray(pool[i]['source_um'])-origin))
        pool.pop(remove)
    return complete,{'status':'certified_joint_complete_path_proposals','guide_flow_count':flow_count,
        'guide_cost_status':mincost_status,'spaced_source_candidates':len(source_points),
        'certified_internal_count':len(accepted),'certified_complete_count':len(complete),
        'individually_certified_path_proposals':len(proposals),'greedy_internal_baseline_count':len(greedy),
        'original_anchor_proposals':original_anchor_count,'suffix_contact_proposals':suffix_contact_count,
        'suffix_contact_tail_lengths_um':contact_tail_lengths,
        'joint_metal_conflict_selection':finite_selection,
        'geometry_rejections':dict(reject),'outer_rejections':dict(outer_reject),'seconds':time.perf_counter()-start,
        'phase_cumulative_seconds':timings,
        'wall_clock_limit_enabled':time_limit_s is not None,'time_limit_s':time_limit_s,
        '_frame':frame,'guide_capacity_scope':'unit node capacity is only a conservative proposal policy; it is not a physical upper bound and does not forbid later multitrack residual routes',
        'global_maximum_proven':False,'uses_generator_source_or_sites_csv':False}
