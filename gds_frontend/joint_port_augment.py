"""Generate new interior routes while jointly rematching all outside Pads.

The incumbent complete layout is retained unless one full-width, spacing-
checked network is added. Internal residual paths are generated from the
continuous geometry, not from the original nearest-outlet library. Pad
assignments and outside paths are reopened jointly for each proposal.
"""
from collections import Counter
from dataclasses import dataclass
import time

import numpy as np
from shapely.geometry import Point,Polygon
from shapely.ops import unary_union

from center_compaction import _try_replacement
from island_router import betti
from ordered_pad_fanout import rebuild_ordered_outer
from port_coverage import coverage_record,outer_port_groups
from navigation_simplification import SearchDeadlineExceeded
from solver_time_policy import deadline_after, deadline_expired


@dataclass(frozen=True)
class AugmentSettings:
    time_limit_s: float | None = None
    candidate_budget_per_port: int = 32
    max_sweeps: int = 2


def _interior_seed(nav,outlet,rules):
    point=list(nav.graph.nodes[outlet]['xy_um'])
    # This is an unconnected construction proposal. Empty Pad and outside
    # placeholders cannot be exported: successful proposals are rebuilt and
    # checked with the whole physical Pad bank before being committed.
    contact=Point(point).buffer(rules.wire_width_um/2+.01,quad_segs=16)
    return {'source_um':point,'outlet_node':outlet,'outlet_um':point,
            'declared_terminal_um':point,'points_um':[point,point],
            'terminal_kind':'external_boundary','navigation_view':0,
            'pad_id':f'proposal-{outlet}','pad_index':0,'pad_side':'unassigned',
            'pad':Polygon(),'join_wire':contact,'outer_wire':contact,
            'outer_length_um':0.0,'line_length_um':0.0}


def augment_joint_ports(nav,chosen,rules,settings,pad_settings,frame,policy,
                        *,progress=None,search=AugmentSettings(),capacity_upper=None,
                        navigation_views=None):
    if (search.candidate_budget_per_port<1 or
            search.max_sweeps<1):raise ValueError('Augmentation computational limits must be positive')
    mapping,groups=outer_port_groups(nav,policy)
    current=list(chosen);before=coverage_record(current,mapping,groups)
    started=time.perf_counter();attempts=[];accepted=[];budget_hit=False;cache={}
    deadline=deadline_after(search.time_limit_s,started=started)
    if capacity_upper is not None and len(current)>capacity_upper:
        raise ValueError('Incumbent exceeds certified continuous capacity upper bound')
    for sweep in range(search.max_sweeps):
        if capacity_upper is not None and len(current)>=capacity_upper:break
        changed=False
        used={mapping[c['outlet_node']] for c in current if c.get('outlet_node') in mapping}
        occupied=unary_union([c['metal'] for c in current]) if current else Polygon()
        # Explore unused windows first, then spare geometry in used windows.
        # A window has no artificial unit-capacity constraint: another track
        # is permitted whenever its real metal and source island fit.
        candidates=[]
        for group in groups:
            if not group['outlets']:continue
            middle=np.asarray(group['representative_um'])
            for outlet in group['outlets']:
                point=np.asarray(nav.graph.nodes[outlet]['xy_um'])
                endpoint_gap=occupied.distance(Point(point))
                if endpoint_gap<rules.wire_width_um/2+rules.spacing_um+settings.numeric_guard_um:
                    continue
                candidates.append((group['group'] in used,
                                   float(np.linalg.norm(point-middle)),group['group'],outlet))
        # Interleave the best sample of every window before taking its second
        # sample, so a large cap cannot consume all of a finite search budget.
        rank_in_group=Counter();ranked=[]
        for item in sorted(candidates,key=lambda c:(c[0],c[2],c[1])):
            rank=rank_in_group[item[2]];rank_in_group[item[2]]+=1
            ranked.append((item[0],rank,item[2],item[3]))
        for ordinal,(_,_,group_id,outlet) in enumerate(sorted(ranked),1):
            if capacity_upper is not None and len(current)>=capacity_upper:break
            if deadline_expired(deadline):
                budget_hit=True;break
            if progress:
                progress(f'出口增广与 Pad 联合重排：{len(current)} 条网络，窗口 {ordinal}/{len(ranked)}',.78)
            point=Point(nav.graph.nodes[outlet]['xy_um'])
            occupied=unary_union([c['metal'] for c in current]) if current else Polygon()
            if occupied.distance(point)<rules.wire_width_um/2+rules.spacing_um+settings.numeric_guard_um:
                continue
            seed=_interior_seed(nav,outlet,rules)
            try:
                internal,detail=_try_replacement(nav.support,seed,current,rules,settings,
                    policy.origin,policy.grid_um,search.candidate_budget_per_port,
                    inward_only=False,deadline=deadline)
            except SearchDeadlineExceeded as exc:
                attempts.append({'status':'search_deadline_reached','stage':str(exc),
                                 'outlet':outlet,'group':group_id})
                budget_hit=True;break
            record={'sweep':sweep+1,'group':group_id,'outlet':outlet,
                    'incumbent_count':len(current),'internal':detail}
            if internal is None:
                record['status']='no_certified_internal_proposal';attempts.append(record);continue
            rebuilt,outer_diagnostic=rebuild_ordered_outer(nav,[*current,internal],rules,settings,
                pad_settings,frame,policy,prepared_cache=cache,navigation_views=navigation_views)
            proposal_frame=frame
            if rebuilt is None and outer_diagnostic.get('reason')=='ordered_bank_has_too_few_slots':
                from pad_router import make_pad_frame
                required=outer_diagnostic['required']
                proposal_frame=make_pad_frame(nav.support,pad_settings,minimum_side_um=frame['square_side_um'],
                    spacing_um=rules.spacing_um,margin_um=rules.margin_um,wire_width_um=rules.wire_width_um,exit_policy=policy,
                    required_slots_per_side=required)
                rebuilt,outer_diagnostic=rebuild_ordered_outer(nav,[*current,internal],rules,settings,
                    pad_settings,proposal_frame,policy,prepared_cache=cache,navigation_views=navigation_views)
                record['frame_growth_proposed_um']=proposal_frame['square_side_um']
            record['outer']=outer_diagnostic
            if rebuilt is None:
                record['status']='joint_outer_proposal_rejected';attempts.append(record);continue
            final_islands=unary_union([c['island'] for c in rebuilt])
            if betti(nav.support.union(final_islands))!=betti(nav.support):
                record['status']='joint_islands_change_source_topology';attempts.append(record);continue
            if len(rebuilt)!=len(current)+1:raise RuntimeError('Joint augmentation did not add exactly one full network')
            current=rebuilt;frame=proposal_frame;changed=True
            record['status']='accepted_additional_complete_network'
            record['new_count']=len(current);attempts.append(record);accepted.append(record)
        if budget_hit or not changed:break
    after=coverage_record(current,mapping,groups)
    if len(current)<len(chosen):raise RuntimeError('Joint augmentation lost existing networks')
    return current,{'method':'residual_interior_route_generation_with_joint_ordered_pad_rematching',
        'before':before,'after':after,'initial_connected_count':len(chosen),
        'final_connected_count':len(current),'accepted_updates':accepted,'attempts':attempts,
        'seconds':time.perf_counter()-started,'time_budget_hit':budget_hit,
        'wall_clock_limit_enabled':search.time_limit_s is not None,
        'stopping_policy':'certified count bound reached, no accepted addition, configured sweep limit, or explicitly requested offline deadline; candidate/cell limits still apply',
        'count_nondecreasing':True,'one_track_per_window_imposed':False,
        'center_spacing_um':rules.minimum_center_spacing_um,
        'search_limits':{'time_limit_s':search.time_limit_s,
                         'candidate_budget_per_port':search.candidate_budget_per_port,
                         'max_sweeps':search.max_sweeps},
        'global_maximum_proven':False,
        '_frame':frame,
        'scope':'constructive complete-network improvements with all current interior nets fixed; outside Pad assignments are changed jointly; failure is not a continuous infeasibility or maximum proof'}
