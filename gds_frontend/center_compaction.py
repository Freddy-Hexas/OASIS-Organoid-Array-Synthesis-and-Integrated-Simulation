"""Move electrode endpoints inward in continuous, occupied-aware geometry.

The Pad, bridge, exit and all other nets stay fixed in each coordinate update.
The moving net is removed before computing its residual configuration space.
Only strictly inward, fully checked replacements are committed. Search
budgets limit construction effort, never prove continuous optimality.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import math
import time

import numpy as np
import shapely
from shapely.geometry import GeometryCollection, LineString, Point, Polygon
from shapely.ops import nearest_points, unary_union

from curved_centerline import CurveSettings, UncertifiableCurve, smooth_centerline
from frontend import _polygonal
from island_router import betti, disk, _shorten
from residual_vector_router import ResidualVectorRouter
from navigation_simplification import SearchDeadlineExceeded,check_deadline
from electrode_region import region_for
from attachment_placement import AttachmentAwarePlacementSearch, PLACEMENT_SEARCH_REVISION
from solver_time_policy import deadline_after, deadline_expired


@dataclass(frozen=True)
class CompactionSettings:
    # Computational limits, deliberately not geometric/process thresholds.
    max_sweeps: int = 3
    candidate_budget_per_net: int = 32
    time_limit_s: float | None = None


def radial_statistics(routes,origin,source_radius):
    radii=[float(np.linalg.norm(np.asarray(c['source_um'])-origin)) for c in routes]
    return {'count':len(radii),
            'within_half_source_radius':sum(r<=source_radius/2+1e-7 for r in radii),
            'mean_radius_um':float(np.mean(radii)) if radii else None,
            'maximum_radius_um':max(radii) if radii else None,
            'sum_squared_radius_um2':float(sum(r*r for r in radii))}


def _placement_proposals(domain,origin,grid_um,budget):
    """Historical finite proposals, retained for inspection/replay comparisons.

    Production compaction and residual port augmentation now use the
    attachment-aware adaptive cell search below instead of this vertex list.
    """
    reference=Point(origin)
    ranked=[]
    for piece in _polygonal(domain):
        inset=piece.buffer(-2*grid_um,quad_segs=16)
        for interior in _polygonal(inset):
            closest=nearest_points(reference,interior)[1]
            ranked.append((closest.distance(reference),interior,closest))
    ranked.sort(key=lambda item:item[0])
    proposals=[];seen=set()
    def add(point):
        xy=np.asarray(point.coords[0],dtype=float)
        key=tuple(int(round(v/grid_um)) for v in xy)
        if key not in seen:
            seen.add(key);proposals.append(xy)
    for _,piece,closest in ranked:
        add(closest)
        if len(proposals)>=budget:break
    # Geometric interior witnesses provide alternatives if the optimal
    # projection attaches across two different substrate branches.
    for _,piece,closest in ranked:
        add(piece.representative_point())
        vertices=np.asarray(piece.exterior.coords)
        order=np.argsort(np.sum((vertices-origin)**2,axis=1))
        for index in order[:budget]:
            vertex=Point(vertices[index])
            add(vertex)
            if len(proposals)>=2*budget:break
        if len(proposals)>=2*budget:break
    proposals.sort(key=lambda xy:float(np.sum((xy-origin)**2)))
    return proposals[:budget],{'feasible_components':len(ranked),
                               'minimum_projected_radius_um':ranked[0][0] if ranked else None}


def _try_replacement(support,old,others,rules,settings,origin,grid_um,budget,
                     *,inward_only=True,deadline=None):
    check_deadline(deadline,'center replacement')
    rejects=Counter();curve_settings=CurveSettings()
    guard=settings.numeric_guard_um
    # Reserve chord approximation and native-grid export error in the
    # residual obstacles as well as in the original-support erosion.
    geometric_guard=guard+curve_settings.max_chord_error_um+4*grid_um
    metal_radius=(rules.wire_width_um/2+curve_settings.max_chord_error_um)/math.cos(math.pi/64)
    residual_guard=geometric_guard+(metal_radius-rules.wire_width_um/2)
    island_radius=rules.electrode_diameter_um/2+rules.margin_um+guard
    neighborhood=support.buffer(max(island_radius+settings.pad_gap_um,
                                    rules.electrode_diameter_um/2+rules.spacing_um)+
                                 geometric_guard,quad_segs=16)
    other_metals=[c['metal'].intersection(neighborhood) for c in others]
    other_metals=[g for g in other_metals if not g.is_empty]
    occupied=unary_union(other_metals) if other_metals else Polygon()
    target=np.asarray(old['points_um'][-1],dtype=float)
    router=ResidualVectorRouter(support,other_metals,
        wire_width_um=rules.wire_width_um,spacing_um=rules.spacing_um,
        margin_um=rules.margin_um,numeric_guard_um=residual_guard,
        target_points=[target],deadline=deadline)
    detail={'pad_id':old['pad_id'],'old_source_um':list(old['source_um']),
            'residual_topology':betti(router.center_domain),
            'component_filter':router.component_filter_certificate,
            'residual_triangles':len(router.triangles)}
    if router.center_domain.is_empty:
        detail.update(status='fixed_exit_has_no_guarded_residual_component')
        return None,detail
    electrode_region=region_for(support,rules,grid_um,origin=origin)
    placement=electrode_region.clip(router.center_domain)
    detail['electrode_region']=electrode_region.record()
    if not occupied.is_empty:
        placement=placement.difference(occupied.buffer(
            rules.electrode_diameter_um/2+rules.spacing_um+geometric_guard,
            quad_segs=32))
    for c in others:
        placement=placement.difference(c['island'].buffer(
            island_radius+settings.pad_gap_um+geometric_guard,quad_segs=32))
        if rules.minimum_center_spacing_um:
            placement=placement.difference(Point(c['source_um']).buffer(
                (rules.minimum_center_spacing_um+geometric_guard)/math.cos(math.pi/128),
                quad_segs=32))
    old_radius=float(np.linalg.norm(np.asarray(old['source_um'])-origin))
    islands_fixed=unary_union([c['island'] for c in others]) if others else Polygon()
    source_topology=betti(support)
    placement_search=AttachmentAwarePlacementSearch(placement,support,origin,
        island_radius_um=island_radius,minimum_substrate_gap_um=settings.pad_gap_um,
        grid_um=grid_um,resolution_um=max(4*grid_um,rules.wire_width_um/4),
        candidate_budget=budget,
        maximum_radius_um=(old_radius-2*grid_um if inward_only else math.inf),
        deadline=deadline)
    best=None;best_audit=None;deadline_hit=False
    def reject(reason):
        rejects[reason]+=1
        placement_search.record_result(False,reason=reason)
    def certify(xy,island):
        if not electrode_region.contains(xy,export_safe=True):
            return None,'outside_electrode_region',None
        if inward_only and float(np.linalg.norm(xy-origin))>=old_radius-2*grid_um:
            return None,'not_strictly_inward',None
        if betti(support.union(islands_fixed).union(island))!=source_topology:
            return None,'combined_islands_change_source_topology',None
        electrode=disk(xy,rules.electrode_diameter_um/2)
        if not support.union(island).covers(electrode.buffer(rules.margin_um,quad_segs=32)):
            return None,'electrode_support_margin',None
        routed=router.route(xy,[target])
        if routed is None:
            return None,'no_certified_residual_route',None
        line=routed[0]
        line=_shorten(line.coords,router.center_domain,rules.wire_width_um/2,
                      visibility_first=True)
        try:
            line,curve=smooth_centerline(line.coords,router.center_domain,curve_settings)
        except UncertifiableCurve:
            return None,'uncertifiable_residual_bend',None
        wire=line.buffer(metal_radius,quad_segs=16)
        if not support.covers(wire.buffer(rules.margin_um+guard,quad_segs=32)):
            return None,'residual_wire_support_margin',None
        central=wire.union(electrode)
        if not support.union(island).covers(central.buffer(rules.margin_um,quad_segs=32)):
            return None,'complete_central_margin',None
        metal=GeometryCollection([central,old['join_wire'],old['outer_wire'],old['pad']])
        if any(metal.distance(c['metal'])<rules.spacing_um+guard for c in others):
            return None,'complete_net_metal_spacing',None
        if any(island.distance(c['island'])<settings.pad_gap_um+guard for c in others):
            return None,'island_spacing',None
        if any(np.linalg.norm(xy-np.asarray(c['source_um']))<
               rules.minimum_center_spacing_um+guard for c in others):
            return None,'center_spacing',None
        if not central.intersects(old['join_wire']):
            return None,'fixed_suffix_disconnected',None
        result={**old,'source_um':xy.tolist(),
                'source_node':('residual_center',old['pad_id'],*map(float,xy)),
                'points_um':np.asarray(line.coords),'line_length_um':line.length,
                'depth_um':line.length,'curve':curve,
                'corridor_ids':[], 'lane_offset_um':None,
                'route_shape_strategy':('continuous_residual_center_compaction' if inward_only else
                                        'continuous_residual_port_augmentation'),
                'island':island,'electrode':electrode,'central_wire':wire,
                'wire':wire,'metal':metal,
                'complete_length_um':line.length+old['outer_length_um'],
                ('center_compaction' if inward_only else 'port_augmentation'):detail}
        return result,None,routed[2]
    try:
        while True:
            proposal=placement_search.next_proposal()
            if proposal is None:break
            xy,island=proposal
            check_deadline(deadline,'electrode placement and curve certification')
            result,reason,audit=certify(xy,island)
            if result is None:
                reject(reason);continue
            best=result;best_audit=audit
            placement_search.record_result(True,radius_um=float(np.linalg.norm(xy-origin)))
    except SearchDeadlineExceeded:
        deadline_hit=True
        # A later search timeout cannot discard an already fully checked move.
        # When no move exists the caller keeps the original complete incumbent.
    search_record=placement_search.certificate(deadline_hit=deadline_hit)
    rejects.update(placement_search.attachment_rejections)
    detail.update(placement_search=search_record,
                  feasible_components=search_record['feasible_components'],
                  minimum_projected_radius_um=search_record['minimum_projected_radius_um'],
                  candidates_tested=search_record['attachment_valid_proposals'],
                  rejections=dict(rejects),search_deadline_hit=deadline_hit)
    if best is not None:
        detail.update(status=('moved_inward' if inward_only else 'certified_internal_proposal'),
                      new_source_um=best['source_um'],old_radius_um=old_radius,
                      new_radius_um=float(np.linalg.norm(np.asarray(best['source_um'])-origin)),
                      path_audit=best_audit)
        return best,detail
    detail.update(status=('no_accepted_inward_replacement' if inward_only else 'no_certified_internal_proposal'))
    return None,detail


def compact_centers(support,chosen,rules,settings,origin,*,grid_um=.001,
                    progress=None,search=CompactionSettings()):
    """Monotone feasible coordinate descent with fixed physical Pad identities."""
    origin=np.asarray(origin,dtype=float);current=list(chosen)
    if grid_um<=0 or search.max_sweeps<1 or search.candidate_budget_per_net<1:
        raise ValueError('Compaction precision, sweeps and candidate limits must be positive')
    source_radius=float(np.linalg.norm(shapely.get_coordinates(support)-origin,axis=1).max())
    before=radial_statistics(current,origin,source_radius)
    started=time.perf_counter();updates=[];attempts=[];completed_sweeps=0;budget_hit=False
    deadline=deadline_after(search.time_limit_s,started=started)
    for sweep in range(search.max_sweeps):
        changed=False
        order=sorted(range(len(current)),key=lambda i:
                     -np.linalg.norm(np.asarray(current[i]['source_um'])-origin))
        for ordinal,index in enumerate(order,1):
            if deadline_expired(deadline):
                budget_hit=True;break
            if progress:
                progress(f'残余空间向中心重布线：第 {sweep+1} 轮，电极 {ordinal}/{len(order)}',.79)
            try:
                replacement,detail=_try_replacement(support,current[index],
                    [c for j,c in enumerate(current) if j!=index],rules,settings,
                    origin,grid_um,search.candidate_budget_per_net,
                    deadline=deadline)
            except SearchDeadlineExceeded as exc:
                attempts.append({'status':'search_deadline_reached','stage':str(exc),
                                 'sweep':sweep+1,'network_index':index+1})
                budget_hit=True;break
            detail={**detail,'sweep':sweep+1,'network_index':index+1}
            attempts.append(detail)
            if replacement is not None:
                current[index]=replacement;updates.append(detail);changed=True
        completed_sweeps+=1
        if budget_hit or not changed:break
    after=radial_statistics(current,origin,source_radius)
    if after['count']!=before['count'] or after['sum_squared_radius_um2']>before['sum_squared_radius_um2']+1e-6:
        raise RuntimeError('Center compaction violated its monotonicity invariant')
    placement_records=[attempt['placement_search'] for attempt in attempts
                       if attempt.get('placement_search')]
    return current,{'method':'attachment_aware_residual_configuration_space_coordinate_descent',
        'placement_search_revision':PLACEMENT_SEARCH_REVISION,
        'placement_search_summary':{
            'net_updates_searched':len(placement_records),
            'visited_cells':sum(r['visited_cells'] for r in placement_records),
            'attachment_valid_proposals':sum(r['attachment_valid_proposals'] for r in placement_records),
            'whole_cell_exclusion_count':sum(r['whole_cell_exclusion_count'] for r in placement_records),
            'unresolved_cell_count':sum(r['unresolved_cell_count'] for r in placement_records),
            'unresolved_cells_are_not_infeasibility_proofs':True},
        'before':before,'after':after,'accepted_updates':updates,
        'attempts':attempts,'completed_sweeps':completed_sweeps,
        'time_budget_hit':budget_hit,'seconds':time.perf_counter()-started,
        'wall_clock_limit_enabled':search.time_limit_s is not None,
        'stopping_policy':'no accepted move in a sweep, configured sweep limit, or explicitly requested offline deadline; candidate/cell limits still apply',
        'search_limits':{'max_sweeps':search.max_sweeps,
                         'candidate_budget_per_net':search.candidate_budget_per_net,
                         'time_limit_s':search.time_limit_s},
        'connected_count_preserved':True,'each_electrode_radius_nonincreasing':True,
        'pad_bridge_exit_and_other_nets_fixed_per_update':True,
        'global_center_optimality_proven':False,
        'scope':'certified feasible improvement relative to the finite-library incumbent; no exhausted search is an impossibility or maximum certificate'}
