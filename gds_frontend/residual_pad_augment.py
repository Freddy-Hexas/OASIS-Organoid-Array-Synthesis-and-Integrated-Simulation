"""Try one additional electrode-to-Pad network around an audited layout.

The existing output GDS is authoritative. The search reconstructs metal,
electrode centers, islands and Pad occupancy from that file plus its report,
then triangulates residual original-support space. Any found route is written
as a new GDS network and independently audited. Failure to find one is only
a search failure, never an upper-bound or maximality certificate.
"""
from __future__ import annotations

from argparse import ArgumentParser
from dataclasses import asdict
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
from tempfile import TemporaryDirectory
import json
import math

import gdstk
import numpy as np
import shapely
from shapely.geometry import LineString, Point, Polygon
from shapely.ops import unary_union

from exact_gds_audit import audit_exported_gds
from exterior_bridge_router import (DirectExteriorPadRouter,
                                    ExteriorBridgeRouter)
from frontend import read_support, _polygonal
from island_router import (IslandSettings, _attachment, betti, build_navigation,
                           disk)
from pad_router import (PadSettings, _add_polygon, _make_outer_option,
                        _nearest_exterior_escape, _pad_choices,
                        _prepare_escape, make_pad_frame)
from process_geometry import ProcessRules
from residual_vector_router import ResidualVectorRouter
from outer_exit_policy import make_outer_exit_policy, policy_from_record


def _grouped_polygons(gds_path):
    with TemporaryDirectory(prefix='residual_readback_') as temporary:
        local=Path(temporary)/'readback.gds'
        local.write_bytes(Path(gds_path).read_bytes())
        lib=gdstk.read_gds(str(local),unit=1e-6)
    groups={}
    for top in lib.top_level():
        for polygon in top.get_polygons():
            groups.setdefault((polygon.layer,polygon.datatype),[]).append(
                shapely.make_valid(Polygon(polygon.points)))
    return {key:unary_union(value) for key,value in groups.items()}


def _adjacent_slots(frame,occupancy):
    result=[]
    for side in ('top','right','bottom','left'):
        bank=occupancy[side]['slot_indices']
        available=[slot for slot in frame['slots'] if slot['side']==side]
        wanted=({bank[0]-1,bank[-1]+1} if bank else
                {available[len(available)//2]['index']})
        result.extend(slot for slot in available if slot['index'] in wanted)
    return result


def _readback_pad_occupancy(groups,frame,layers,count,grid_um):
    """Identify every occupied physical Pad slot from exported marker shapes."""
    occupied={side:[] for side in ('top','right','bottom','left')}
    slots=frame['slots']
    for network in range(1,count+1):
        marker=groups.get((layers['pad_contact_layer'],network))
        if marker is None or marker.is_empty:
            raise ValueError(f'Network {network} has no exported Pad marker')
        matches=[]
        for slot in slots:
            if max(abs(a-b) for a,b in zip(marker.bounds,
                                          slot['polygon'].bounds))<=2*grid_um+.002 and (
                    marker.symmetric_difference(slot['polygon']).area<=
                    2*grid_um*(marker.length+slot['polygon'].length)):
                matches.append(slot)
        if len(matches)!=1:
            raise ValueError(f'Network {network} does not match one physical Pad slot')
        slot=matches[0]
        occupied[slot['side']].append(slot['index'])
    for side,indices in occupied.items():
        indices.sort()
        if len(indices)!=len(set(indices)) or (indices and
                indices[-1]-indices[0]+1!=len(indices)):
            raise ValueError(f'{side} exported Pad bank has a duplicate or vacancy')
    return occupied


def _occupancy_record(indices_by_side,pitch,width):
    record={}
    for side,indices in indices_by_side.items():
        indices=sorted(indices)
        record[side]={'connected_pad_count':len(indices),
                      'slot_indices':indices,
                      'first_slot':indices[0] if indices else None,
                      'last_slot':indices[-1] if indices else None,
                      'internal_empty_slots':0,
                      'adjacent_center_pitch_um':pitch if len(indices)>1 else None,
                      'adjacent_edge_gap_um':pitch-width if len(indices)>1 else None}
    return record


def _existing_geometry(report,groups):
    layers=report['audit']['output_layers']
    count=report['connected_pad_count']
    metal_layer=layers['metal_layer']
    marker_layer=layers['electrode_marker_layer']
    island_layer=layers['island_marker_layer']
    shell=groups[(layers['shell_marker_layer'],0)]
    metals=[groups[(metal_layer,i)] for i in range(1,count+1)]
    islands=[groups[(island_layer,i)] for i in range(1,count+1)]
    centers=[]
    for i in range(1,count+1):
        marker=groups[(marker_layer,i)]
        centers.append(np.asarray(marker.centroid.coords[0]))
    return layers,shell,metals,islands,np.asarray(centers)


def _shell_detour_option(nav,frame,rules,prepared,slot,shell_router,occupied):
    """Route from an original-support exit through the residual outer shell."""
    bridge_line=prepared.get('bridge_line')
    bridge_vertices=(list(bridge_line.coords)[1:] if bridge_line is not None
                     else [prepared['shell_point']])
    vertices=[*prepared['inside_line'].coords,*bridge_vertices]
    points=[]
    for vertex in vertices:
        xy=np.asarray(vertex,dtype=float)
        if not points or np.linalg.norm(xy-points[-1])>1e-8:
            points.append(xy)
    if len(points)<2:
        return None,'shell_detour_has_no_bridge_line'
    bridge_line=LineString(points)
    bridge_wire=bridge_line.buffer(rules.wire_width_um/2+.01,
                                   quad_segs=32,cap_style='round',join_style='round')
    if bridge_wire.distance(occupied)<rules.spacing_um+.05:
        return None,'shell_detour_bridge_metal_gap'
    supported=unary_union([nav.support,prepared['bridge'],frame['shell']])
    if not supported.covers(bridge_wire.buffer(rules.margin_um+.005,
                                               quad_segs=32)):
        return None,'shell_detour_bridge_margin'
    routed=shell_router.route(prepared['shell_point'],[slot['target_um']])
    if routed is None:
        return None,'shell_detour_no_shell_path'
    shell_line=routed[0]
    line=LineString([*bridge_line.coords,*list(shell_line.coords)[1:]])
    wire=line.buffer(rules.wire_width_um/2+.01,quad_segs=32,
                     cap_style='round',join_style='round')
    outer_metal=unary_union([wire,slot['polygon']])
    if not wire.intersects(slot['polygon']):
        return None,'shell_detour_pad_contact'
    if outer_metal.distance(occupied)<rules.spacing_um+.05:
        return None,'shell_detour_metal_gap'
    if not supported.covers(outer_metal.buffer(rules.margin_um+.005,
                                               quad_segs=32)):
        return None,'shell_detour_support_margin'
    return {'pad_id':slot['pad_id'],'pad_side':slot['side'],
            'pad_index':slot['index'],'pad':slot['polygon'],
            'pad_target_um':slot['target_um'],
            'bridge':prepared['bridge'],'outer_wire':wire,
            'outer_line':line,
            'outer_route_strategy':'residual_shell_detour',
            'outer_curve':{'method':'residual_shell_constrained_triangles',
                           'path_audit':routed[2]},
            'outer_length_um':line.length,
            'portal_escape_um':prepared['exit_point'].tolist()},None


def _fanout_escape_candidates(nav,rules,pad_settings,frame,outlet,
                              *,numeric_guard_um=.05):
    """Propose bridge bends from physical pitch and local exit geometry.

    These are construction candidates only.  The exact exported-GDS audit is
    still required, and exhausting them says nothing about optimality.
    """
    start=np.asarray(nav.graph.nodes[outlet]['xy_um'],dtype=float)
    exit_point,inside=_nearest_exterior_escape(nav,start)
    if exit_point is None:
        return []
    origin=np.asarray(frame['origin_um'],dtype=float)
    radial=exit_point-origin
    radius=float(np.linalg.norm(radial))
    outward=exit_point-start
    normal_length=float(np.linalg.norm(outward))
    if radius<1 or normal_length<1e-6:
        return []
    outward/=normal_length
    base_angle=math.atan2(radial[1],radial[0])
    shell_radius=(frame['inner_radius_um']+
                  pad_settings.bridge_shell_overlap_um)
    # A displacement of one metal pitch at the launch radius gives one
    # angular step.  No structure name or prescribed number of tracks enters.
    pitch=rules.wire_width_um+rules.spacing_um+2*numeric_guard_um
    angle_step=math.atan2(pitch,max(radius,pitch))
    offsets=(0,.5,-.5,1,-1,2,-2,4,-4)
    candidates=[]
    seen=set()
    source_holes=[Polygon(ring) for component in _polygonal(nav.support)
                  for ring in component.interiors]
    stored=getattr(nav,'summary',{}).get('outer_exit_policy')
    policy=(policy_from_record(stored) if stored else make_outer_exit_policy(nav.support,
        wire_width_um=rules.wire_width_um,margin_um=rules.margin_um,
        bridge_width_um=pad_settings.bridge_width_um,numeric_guard_um=numeric_guard_um))
    for factor in offsets:
        angle=base_angle+factor*angle_step
        shell_point=origin+shell_radius*np.asarray((math.cos(angle),
                                                     math.sin(angle)))
        # The short outward step is derived from the substrate bridge width.
        # It lets a bridge leave a concave boundary before turning.
        for local_step in (0,pad_settings.bridge_width_um):
            vertices=[exit_point]
            if local_step:
                vertices.append(exit_point+local_step*outward)
            vertices.append(shell_point)
            line=LineString(vertices)
            signature=tuple((round(float(x),6),round(float(y),6))
                            for x,y in line.coords)
            if signature in seen:
                continue
            seen.add(signature)
            if line.intersection(nav.support).length>.002:
                continue
            bridge=line.buffer(pad_settings.bridge_width_um/2,
                               quad_segs=32,cap_style='round',
                               join_style='round')
            if (bridge.intersection(nav.support).area<=1 or
                    bridge.intersection(frame['shell']).area<=1):
                continue
            added=bridge.difference(nav.support)
            if any(added.intersects(hole) for hole in source_holes):
                continue
            allowed,_=policy.check_bridge(nav.support,bridge,exit_point=exit_point)
            if not allowed:continue
            candidates.append({'start':start,'inside_line':inside,
                               'exit_point':exit_point,
                               'shell_point':shell_point,'bridge':bridge,
                               'bridge_line':line,
                               'fanout_pitch_factor':factor,
                               'fanout_local_step_um':local_step})
    return candidates


def _outer_options(nav,frame,pad_settings,rules,slots,existing_metals,
                   *,numeric_guard_um=.05,max_outlets=None,
                   direct_exterior_only=False):
    by_outlet={}
    rejected={}
    diagnostics={}
    occupied=unary_union(existing_metals)
    shell_router=None
    exterior_router=None
    direct_exterior_router=None
    policy=policy_from_record(nav.summary['outer_exit_policy'])
    # A complete independent Pad must itself clear every old metal net.
    # Checking the exported Pad polygon first avoids expensive route searches
    # that cannot make a blocked adjacent slot usable.
    blocked=[]
    usable=[]
    for slot in slots:
        blockers=[network for network,metal in enumerate(existing_metals,1)
                  if slot['polygon'].distance(metal)<rules.spacing_um+
                  numeric_guard_um]
        if blockers:
            blocked.append({'pad_id':slot['pad_id'],
                            'blocking_existing_networks':blockers})
        else:
            usable.append(slot)
    diagnostics['adjacent_pad_slot_blockers']=blocked
    diagnostics['physically_usable_adjacent_pad_slots']=len(usable)
    if not usable:
        diagnostics['outlets_sampled']=0
        return by_outlet,rejected,diagnostics
    slot_ids={slot['pad_id'] for slot in usable}
    outlets=list(nav.outlets)
    if max_outlets is not None and len(outlets)>max_outlets:
        if max_outlets<1:
            raise ValueError('max_outlets must be positive')
        selected={round(i*(len(outlets)-1)/(max_outlets-1))
                  for i in range(max_outlets)} if max_outlets>1 else {len(outlets)//2}
        outlets=[outlets[i] for i in sorted(selected)]
    diagnostics['outlets_sampled']=len(outlets)
    for outlet in outlets:
        point=nav.graph.nodes[outlet]['xy_um']
        eligible=[slot for slot in _pad_choices(point,frame,pad_settings)
                  if slot['pad_id'] in slot_ids]
        if not eligible:
            continue
        if direct_exterior_only:
            if direct_exterior_router is None:
                direct_exterior_router=DirectExteriorPadRouter(
                    nav.support,frame,existing_metals,
                    bridge_width_um=pad_settings.bridge_width_um,
                    bridge_shell_overlap_um=
                        pad_settings.bridge_shell_overlap_um,
                    wire_width_um=rules.wire_width_um,
                    spacing_um=rules.spacing_um,
                    margin_um=rules.margin_um,
                    numeric_guard_um=numeric_guard_um,exit_policy=policy)
            start=np.asarray(point,dtype=float)
            exit_point,inside=_nearest_exterior_escape(nav,start)
            if exit_point is not None:
                accepted=list(direct_exterior_router.propose(
                    exit_point,start,inside,eligible))
                if accepted:
                    by_outlet[outlet]=sorted(
                        accepted,key=lambda option:option['outer_length_um'])
            continue
        prepared,reason=_prepare_escape(nav,rules,pad_settings,frame,outlet)
        if prepared is None:
            rejected[reason]=rejected.get(reason,0)+1
        accepted=[]
        if prepared is not None:
            for slot in eligible:
                option,reason=_make_outer_option(
                    nav,rules,pad_settings,frame,outlet,slot,prepared)
                if option is None:
                    rejected[reason]=rejected.get(reason,0)+1
                    continue
                outer_metal=unary_union([option['outer_wire'],option['pad']])
                if outer_metal.distance(occupied)<rules.spacing_um+.05:
                    rejected['outer_metal_existing_gap']=(
                        rejected.get('outer_metal_existing_gap',0)+1)
                    continue
                accepted.append(option)
        if not accepted:
            if shell_router is None:
                shell_neighborhood=frame['shell'].buffer(
                    rules.wire_width_um+rules.spacing_um+rules.margin_um+1)
                shell_metals=[part.intersection(shell_neighborhood)
                              for part in existing_metals
                              if part.intersects(shell_neighborhood)]
                shell_router=ResidualVectorRouter(
                    frame['shell'],shell_metals,
                    wire_width_um=rules.wire_width_um,
                    spacing_um=rules.spacing_um,margin_um=rules.margin_um,
                    numeric_guard_um=numeric_guard_um)
            if prepared is not None:
                for slot in eligible:
                    option,reason=_shell_detour_option(
                        nav,frame,rules,prepared,slot,shell_router,occupied)
                    if option is None:
                        rejected[reason]=rejected.get(reason,0)+1
                    else:
                        accepted.append(option)
            if not accepted:
                if direct_exterior_router is None:
                    direct_exterior_router=DirectExteriorPadRouter(
                        nav.support,frame,existing_metals,
                        bridge_width_um=pad_settings.bridge_width_um,
                        bridge_shell_overlap_um=
                            pad_settings.bridge_shell_overlap_um,
                        wire_width_um=rules.wire_width_um,
                        spacing_um=rules.spacing_um,
                        margin_um=rules.margin_um,
                        numeric_guard_um=numeric_guard_um,exit_policy=policy)
                start=np.asarray(point,dtype=float)
                exit_point,inside=_nearest_exterior_escape(nav,start)
                if exit_point is not None:
                    accepted.extend(direct_exterior_router.propose(
                        exit_point,start,inside,eligible))
            if not accepted:
                for bent in _fanout_escape_candidates(
                        nav,rules,pad_settings,frame,outlet,
                        numeric_guard_um=numeric_guard_um):
                    bridge_vertices=[*bent['inside_line'].coords,
                                     *list(bent['bridge_line'].coords)[1:]]
                    bridge_wire=LineString(bridge_vertices).buffer(
                        rules.wire_width_um/2+.01,quad_segs=16,
                        cap_style='round',join_style='round')
                    if bridge_wire.distance(occupied)<rules.spacing_um+.05:
                        rejected['fanout_bridge_metal_gap']=(
                            rejected.get('fanout_bridge_metal_gap',0)+1)
                        continue
                    for slot in eligible:
                        option,reason=_shell_detour_option(
                            nav,frame,rules,bent,slot,shell_router,occupied)
                        if option is None:
                            rejected[reason]=rejected.get(reason,0)+1
                        else:
                            option['outer_route_strategy']='pitch_adaptive_bent_bridge_and_shell'
                            accepted.append(option)
                    if accepted:
                        break
            if not accepted:
                if exterior_router is None:
                    exterior_router=ExteriorBridgeRouter(
                        nav.support,frame,existing_metals,
                        bridge_width_um=pad_settings.bridge_width_um,
                        bridge_shell_overlap_um=
                            pad_settings.bridge_shell_overlap_um,
                        wire_width_um=rules.wire_width_um,
                        spacing_um=rules.spacing_um,
                        numeric_guard_um=numeric_guard_um,exit_policy=policy)
                start=np.asarray(point,dtype=float)
                exit_point,inside=_nearest_exterior_escape(nav,start)
                if exit_point is not None:
                    for curved in exterior_router.propose(
                            exit_point,start,inside):
                        bridge_vertices=[*inside.coords,
                                         *list(curved['bridge_line'].coords)[1:]]
                        bridge_wire=LineString(bridge_vertices).buffer(
                            rules.wire_width_um/2+.01,quad_segs=16,
                            cap_style='round',join_style='round')
                        if bridge_wire.distance(occupied)<rules.spacing_um+.05:
                            rejected['exterior_bridge_metal_gap']=(
                                rejected.get('exterior_bridge_metal_gap',0)+1)
                            continue
                        for slot in eligible:
                            option,reason=_shell_detour_option(
                                nav,frame,rules,curved,slot,shell_router,
                                occupied)
                            if option is None:
                                rejected[reason]=rejected.get(reason,0)+1
                            else:
                                option['outer_route_strategy']=curved[
                                    'external_strategy']
                                option['exterior_path_audit']=curved[
                                    'external_path_audit']
                                accepted.append(option)
                        if accepted:
                            break
        if accepted:
            permitted=[]
            for option in accepted:
                allowed,reason=policy.check_bridge(nav.support,option['bridge'],
                    exit_point=option['portal_escape_um'])
                if allowed:allowed,reason=policy.check_exterior_line(nav.support,option['outer_line'])
                if allowed:permitted.append(option)
                else:rejected[reason]=rejected.get(reason,0)+1
            if permitted:by_outlet[outlet]=sorted(permitted,key=lambda option:option['outer_length_um'])
            else:by_outlet.pop(outlet,None)
    if direct_exterior_router is not None:
        diagnostics['joint_exterior_to_pad']=direct_exterior_router.diagnostics
    return by_outlet,rejected,diagnostics


def find_augmenting_route(original_gds,report,*,candidate_step_um=600,
                          external_portal_pitch_um=40,max_sources=200,
                          external_portal_mode='radial_envelope',
                          exclude_source_um=None,
                          navigation_cache=None,
                          anchor_order='farthest_existing',
                          target_trials=1,max_outlets=None,
                          direct_exterior_only=False):
    original=Path(original_gds).resolve(strict=True)
    routed=Path(report['output']).resolve(strict=True)
    if (sha256(original.read_bytes()).hexdigest()!=report['input_sha256'] or
            sha256(routed.read_bytes()).hexdigest()!=report['output_sha256']):
        raise ValueError('Original or routed GDS changed since the audit')
    if not report.get('integer_polygon_audit',{}).get('outer_exit_policy_verified'):
        raise ValueError('Existing Pad layout lacks the outer exit certificate; rerun the base Pad task under the current policy before augmentation')
    rules=ProcessRules(**report['rules'])
    prior_settings=report.get('island_settings',{})
    settings=IslandSettings(
        candidate_step_um=candidate_step_um,
        external_portal_pitch_um=external_portal_pitch_um,
        external_portal_mode=external_portal_mode,
        numeric_guard_um=prior_settings.get('numeric_guard_um',.05),
        pad_gap_um=prior_settings.get('pad_gap_um',4.0),
        navigation_component_exclusion_distance_um=0.0,
        navigation_prune_leaf_length_um=0.0)
    pad_settings=PadSettings(**report['pad_settings'])
    source_support,_,source_meta=read_support(original,*report['support_layer'])
    policy=make_outer_exit_policy(source_support,wire_width_um=rules.wire_width_um,
        margin_um=rules.margin_um,numeric_guard_um=settings.numeric_guard_um,
        bridge_width_um=pad_settings.bridge_width_um,
        grid_um=source_meta['gds_native_precision_m']*1e6)
    cache_key=(report['input_sha256'],tuple(report['support_layer']),
               rules,settings,pad_settings.bridge_width_um)
    if navigation_cache is not None and navigation_cache.get('key')==cache_key:
        nav=navigation_cache['navigation']
    else:
        nav=build_navigation(original,rules,
                             layer=report['support_layer'][0],
                             datatype=report['support_layer'][1],
                             outlet_mode='external_boundary',settings=settings,
                             outer_exit_policy=policy)
        if navigation_cache is not None:
            navigation_cache.clear()
            navigation_cache.update({'key':cache_key,'navigation':nav})
    groups=_grouped_polygons(routed)
    layers,shell,metals,islands,centers=_existing_geometry(report,groups)
    side=float(shell.bounds[2]-shell.bounds[0])
    frame=make_pad_frame(nav.support,pad_settings,minimum_side_um=side,
                         spacing_um=rules.spacing_um,margin_um=rules.margin_um,
                         exit_policy=policy)
    actual_occupancy=_readback_pad_occupancy(
        groups,frame,layers,report['connected_pad_count'],
        report['integer_polygon_audit']['native_grid_um'])
    claimed={side:report['pad_bank_occupancy'][side]['slot_indices']
             for side in actual_occupancy}
    if actual_occupancy!=claimed:
        raise ValueError('Report Pad bank differs from the exported GDS')
    slots=_adjacent_slots(frame,report['pad_bank_occupancy'])
    options,rejected,outer_diagnostics=_outer_options(
        nav,frame,pad_settings,rules,slots,metals,
        numeric_guard_um=settings.numeric_guard_um,
        max_outlets=max_outlets,
        direct_exterior_only=direct_exterior_only)
    diagnostics={'candidate_anchors':len(nav.candidates),
                 'external_boundary_outlets':len(nav.outlets),
                 'external_boundary_sampling':
                     nav.summary.get('external_boundary_sampling'),
                 'outlets_sampled':outer_diagnostics['outlets_sampled'],
                 'adjacent_pad_slots':len(slots),
                 'physically_usable_adjacent_pad_slots':
                     outer_diagnostics['physically_usable_adjacent_pad_slots'],
                 'adjacent_pad_slot_blockers':
                     outer_diagnostics['adjacent_pad_slot_blockers'],
                 'exported_pad_bank_reverified':True,
                 'outlets_with_outer_pad_options':len(options),
                 'outer_option_rejections':rejected}
    if 'joint_exterior_to_pad' in outer_diagnostics:
        diagnostics['joint_exterior_to_pad']=outer_diagnostics[
            'joint_exterior_to_pad']
    if not options:
        return None,diagnostics
    # Only metal near the original support can block an internal route.
    near_support=nav.support.buffer(
        rules.wire_width_um+rules.spacing_um+rules.margin_um+1)
    core_metals=[part.intersection(near_support) for part in metals
                 if part.intersects(near_support)]
    router=ResidualVectorRouter(
        nav.support,core_metals,wire_width_um=rules.wire_width_um,
        spacing_um=rules.spacing_um,margin_um=rules.margin_um,
        numeric_guard_um=settings.numeric_guard_um)
    diagnostics['residual_triangles']=len(router.triangles)
    if not router.triangles:
        return None,diagnostics
    option_outlets=list(options)
    targets=[nav.graph.nodes[outlet]['xy_um'] for outlet in option_outlets]
    used_islands=unary_union(islands)
    occupied=unary_union(metals)
    island_radius=(rules.electrode_diameter_um/2+rules.margin_um+
                   settings.numeric_guard_um)
    required_gap=max(rules.minimum_center_spacing_um,
                     rules.electrode_diameter_um+rules.spacing_um,
                     2*island_radius+settings.pad_gap_um)
    if anchor_order not in ('farthest_existing','closest_eligible_exit') or target_trials<1:
        raise ValueError('Unknown geometric anchor order')
    # Both orders are geometry-only heuristics. Neither contributes to U.
    candidates=[]
    seen=set()
    target_coordinates=np.asarray(targets,dtype=float)
    for node in nav.candidates:
        xy=tuple(map(float,nav.graph.nodes[node]['xy_um']))
        if xy in seen:
            continue
        seen.add(xy)
        distance=min(np.linalg.norm(centers-np.asarray(xy),axis=1))
        if distance+1e-7<required_gap:
            continue
        nearest_target=min(np.linalg.norm(target_coordinates-
                                          np.asarray(xy),axis=1))
        priority=((-distance,nearest_target) if
                  anchor_order=='farthest_existing' else
                  (nearest_target,-distance))
        candidates.append((*priority,node,xy))
    candidates.sort()
    diagnostics['anchor_order']=anchor_order
    if exclude_source_um is not None:
        excluded=np.asarray(exclude_source_um,dtype=float)
        if excluded.shape==(2,):
            excluded=excluded.reshape(1,2)
        if (excluded.ndim!=2 or excluded.shape[1]!=2 or
                len(excluded)<1 or not np.isfinite(excluded).all()):
            raise ValueError('Excluded sources must be finite XY coordinates')
        candidates=[item for item in candidates if
                    np.min(np.linalg.norm(excluded-np.asarray(item[3]),
                                          axis=1))>=required_gap-1e-7]
        diagnostics['excluded_source_centers_um']=excluded.tolist()
        if len(excluded)==1:
            diagnostics['excluded_source_center_um']=excluded[0].tolist()
        diagnostics['excluded_source_radius_um']=required_gap
    checked=0
    rejected_sources={}
    def reject(reason):
        rejected_sources[reason]=rejected_sources.get(reason,0)+1
    for _,_,node,xy in candidates[:max_sources]:
        checked+=1
        island,reason=_attachment(nav.support,xy,island_radius,
                                   settings.pad_gap_um)
        if island is None:
            reject(reason)
            continue
        if island.distance(used_islands)<settings.pad_gap_um+.05:
            reject('island_existing_gap')
            continue
        electrode=disk(xy,rules.electrode_diameter_um/2)
        if electrode.distance(occupied)<rules.spacing_um+.05:
            reject('electrode_existing_metal_gap')
            continue
        remaining=list(range(len(targets)))
        for trial in range(min(target_trials,len(remaining))):
            result=router.route(xy,[targets[i] for i in remaining])
            if result is None:
                reject('residual_route_not_found')
                break
            line,local_id,path_audit=result
            target_id=remaining.pop(local_id)
            outlet=option_outlets[target_id]
            central_wire=line.buffer(
                rules.wire_width_um/2+.01,quad_segs=32,
                cap_style='round',join_style='round')
            for option in options[outlet]:
                metal=unary_union([electrode,central_wire,
                                   option['outer_wire'],option['pad']])
                if betti(metal)['components']!=1:
                    reject('complete_metal_disconnected')
                    continue
                if metal.distance(occupied)<rules.spacing_um+.05:
                    reject('complete_metal_existing_gap')
                    continue
                final_support=unary_union([groups[tuple(layers['support_layer'])],
                                           island,option['bridge']])
                if not final_support.covers(
                        metal.buffer(rules.margin_um+.005,quad_segs=32)):
                    reject('complete_metal_support_margin')
                    continue
                if betti(unary_union([nav.support,used_islands,island]))!=betti(
                        unary_union([nav.support,used_islands])):
                    reject('island_changes_original_topology')
                    continue
                route={**option,'source_node':node,'source_um':list(xy),
                       'island':island,'electrode':electrode,'wire':central_wire,
                       'central_wire':central_wire,'metal':metal,
                       'points_um':np.asarray(line.coords),
                       'declared_terminal_um':targets[target_id],
                       'residual_path_audit':path_audit,
                       'outlet_node':outlet,
                       'route_shape_strategy':'residual_constrained_triangles'}
                diagnostics['anchors_checked']=checked
                diagnostics['selected_outlet']=outlet
                diagnostics['selected_pad']=option['pad_id']
                diagnostics['selected_outer_route_strategy']=option.get(
                    'outer_route_strategy','direct_smoothed_outer_line')
                diagnostics['selected_source_um']=list(xy)
                diagnostics['residual_path_audit']=path_audit
                diagnostics['selected_target_trial']=trial+1
                diagnostics['source_rejections']=rejected_sources
                return route,diagnostics
    diagnostics['anchors_checked']=checked
    diagnostics['source_rejections']=rejected_sources
    return None,diagnostics


def export_augmented_gds(original_gds,report,route,output_path):
    output_path=Path(output_path)
    old=Path(report['output'])
    layers=report['audit']['output_layers']
    count=report['connected_pad_count']
    grid_um=report['integer_polygon_audit']['native_grid_um']
    with TemporaryDirectory(prefix='residual_export_') as temporary:
        local=Path(temporary)/'old.gds'
        local.write_bytes(old.read_bytes())
        lib=gdstk.read_gds(str(local),unit=1e-6)
        tops=lib.top_level()
        if len(tops)!=1:
            raise ValueError('Audited Pad GDS must have one routed top cell')
        cell=tops[0]
        network=count+1
        support_layer,support_datatype=layers['support_layer']
        for geom in (route['bridge'],route['island']):
            _add_polygon(cell,geom,support_layer,support_datatype)
        for geom,layer,datatype in [
                (route['bridge'],layers['bridge_marker_layer'],network),
                (route['island'],layers['island_marker_layer'],network),
                (route['metal'],layers['metal_layer'],network),
                (route['electrode'],layers['metal_layer'],network),
                (route['pad'],layers['metal_layer'],network),
                (route['electrode'],layers['electrode_marker_layer'],network),
                (route['pad'],layers['pad_contact_layer'],network)]:
            _add_polygon(cell,geom,layer,datatype)
        local_out=Path(temporary)/'new.gds'
        lib.write_gds(str(local_out),max_points=4000)
        result=local_out.read_bytes()
    output_path.parent.mkdir(parents=True,exist_ok=True)
    output_path.write_bytes(result)
    old_witness=json.loads(Path(layers['wire_width_witness_path']).read_text(
        encoding='utf-8'))
    if old_witness['output_sha256']!=report['output_sha256']:
        raise ValueError('Previous width witness is not bound to routed GDS')
    vertices=[*route['points_um'],*route['outer_line'].coords]
    points=[]
    for vertex in vertices:
        xy=[int(round(float(value)/grid_um)) for value in vertex]
        if not points or points[-1]!=xy:
            points.append(xy)
    old_witness['networks'].append({'network':network,
                                    'points_grid_ticks':points})
    old_witness['output_sha256']=sha256(result).hexdigest()
    new_witness=output_path.with_suffix('.width_witness.json')
    new_witness.write_text(json.dumps(old_witness,ensure_ascii=False,
                                      separators=(',',':')),encoding='utf-8')
    return {**layers,'wire_width_witness_path':str(new_witness.resolve())}


def _augment_once(report,source,destination,*,candidate_step_um,
                  external_portal_pitch_um,max_sources,navigation_cache=None,
                  external_portal_mode='radial_envelope',
                  exclude_source_um=None,
                  anchor_order='farthest_existing',target_trials=1,
                  max_outlets=None,direct_exterior_only=False):
    if destination.exists() and any(destination.iterdir()):
        raise FileExistsError(f'Augmentation evidence directory is not empty: {destination}')
    route,diagnostics=find_augmenting_route(
        source,report,candidate_step_um=candidate_step_um,
        external_portal_pitch_um=external_portal_pitch_um,
        external_portal_mode=external_portal_mode,
        exclude_source_um=exclude_source_um,
        max_sources=max_sources,navigation_cache=navigation_cache,
        anchor_order=anchor_order,target_trials=target_trials,
        max_outlets=max_outlets,
        direct_exterior_only=direct_exterior_only)
    destination.mkdir(parents=True,exist_ok=True)
    result={'input':str(source),'input_sha256':report['input_sha256'],
            'previous_output':report['output'],
            'previous_output_sha256':report['output_sha256'],
            'previous_connected_pad_count':report['connected_pad_count'],
            'search_diagnostics':diagnostics}
    if route is None:
        result.update({'status':'no_augmenting_route_found',
                       'optimality_proven':False})
    else:
        output=destination/'residual_augmented.gds'
        layers=export_augmented_gds(source,report,route,output)
        rules=ProcessRules(**report['rules'])
        prior_settings=report.get('island_settings',{})
        settings=IslandSettings(
            numeric_guard_um=prior_settings.get('numeric_guard_um',.05),
            pad_gap_um=prior_settings.get('pad_gap_um',4.0))
        pads=PadSettings(**report['pad_settings'])
        island_radius=(rules.electrode_diameter_um/2+rules.margin_um+
                       settings.numeric_guard_um)
        required_gap=max(rules.minimum_center_spacing_um,
                         rules.electrode_diameter_um+rules.spacing_um,
                         2*island_radius+settings.pad_gap_um)
        audit=audit_exported_gds(
            source,output,layers,wire_spacing_um=rules.spacing_um,
            metal_support_margin_um=rules.margin_um,
            expected_nets=report['connected_pad_count']+1,
            minimum_electrode_diameter_um=rules.electrode_diameter_um,
            minimum_electrode_center_distance_um=required_gap,
            minimum_substrate_disk_radius_um=island_radius,
            maximum_substrate_disk_radius_um=island_radius+.01,
            minimum_island_spacing_um=settings.pad_gap_um,
            minimum_pad_short_side_um=pads.pad_width_um,
            minimum_pad_long_side_um=pads.pad_length_um,
            minimum_wire_width_um=rules.wire_width_um)
        if audit['passed'] is not True:
            raise RuntimeError('Augmented GDS did not pass the integer audit')
        groups=_grouped_polygons(output)
        shell=groups[(layers['shell_marker_layer'],0)]
        nav_support=read_support(source,layer=report['support_layer'][0],
                                 datatype=report['support_layer'][1])[0]
        side=float(shell.bounds[2]-shell.bounds[0])
        frame=make_pad_frame(nav_support,pads,minimum_side_um=side,
                             spacing_um=rules.spacing_um,margin_um=rules.margin_um)
        bank=_readback_pad_occupancy(groups,frame,layers,
                                     report['connected_pad_count']+1,
                                     audit['native_grid_um'])
        if route['pad_index'] not in bank[route['pad_side']]:
            raise RuntimeError('New Pad is absent from exported Pad bank')
        occupancy=_occupancy_record(bank,pads.pad_pitch_um,
                                     pads.pad_width_um)
        output_hash=sha256(output.read_bytes()).hexdigest()
        next_report={
            'created_at':datetime.now(timezone.utc).isoformat(),
            'input':str(source),'input_sha256':report['input_sha256'],
            'support_layer':report['support_layer'],'rules':report['rules'],
            'pad_settings':report['pad_settings'],
            'island_settings':{'candidate_step_um':candidate_step_um,
                               'external_portal_pitch_um':external_portal_pitch_um,
                               'external_portal_mode':external_portal_mode,
                               'max_sources':max_sources,
                               'anchor_order':anchor_order,
                               'target_trials':target_trials,
                               'numeric_guard_um':settings.numeric_guard_um,
                               'pad_gap_um':settings.pad_gap_um,
                               'navigation_component_exclusion_distance_um':0.0,
                               'navigation_prune_leaf_length_um':0.0,
                               'method':'residual_constrained_triangles'},
            'connected_pad_count':report['connected_pad_count']+1,
            'continuous_upper_bound':report.get('continuous_upper_bound'),
            'output':str(output),'output_sha256':output_hash,
            'pad_bank_occupancy':occupancy,
            'audit':{'passed':True,'pad_connection_verified':True,
                     'connected_pad_count':report['connected_pad_count']+1,
                     'output_layers':layers,'pad_bank_occupancy':occupancy,
                     'audit_scope':'exported GDS integer polygon audit'},
            'integer_polygon_audit':audit,
            'optimization':{'status':'constructive_incremental_search',
                            'finite_upper_bound':None,
                            'continuous_optimality_proven':False},
            'construction_provenance':{
                'previous_output_sha256':report['output_sha256'],
                'previous_connected_pad_count':report['connected_pad_count'],
                'new_pad_id':route['pad_id'],
                'new_source_um':list(route['source_um']),
                'method':'residual constrained triangulation and exact readback'}}
        report_path=destination/'report.json'
        report_path.write_text(json.dumps(next_report,ensure_ascii=False,indent=2),
                               encoding='utf-8')
        result.update({'status':'audited_additional_network',
                       'optimality_proven':False,
                       'certified_lower_bound':report['connected_pad_count']+1,
                       'connected_pad_count':report['connected_pad_count']+1,
                       'output':str(output),'output_sha256':output_hash,
                       'output_layers':layers,'pad_bank_occupancy':occupancy,
                       'integer_polygon_audit':audit,'report':str(report_path),
                       'rules':report['rules']})
    path=destination/'residual_search_result.json'
    path.write_text(json.dumps(result,ensure_ascii=False,indent=2),
                    encoding='utf-8')
    print(json.dumps({key:result[key] for key in
                      ('status','previous_connected_pad_count',
                       'connected_pad_count','optimality_proven')
                      if key in result},ensure_ascii=False),flush=True)
    return result


def main():
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--report',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--candidate-step-um',type=float,default=600)
    parser.add_argument('--external-portal-pitch-um',type=float,default=40)
    parser.add_argument('--external-portal-mode',
                        choices=('radial_envelope','full_boundary'),
                        default='radial_envelope')
    parser.add_argument('--exclude-source-um',type=float,nargs=2,action='append',
                        metavar=('X','Y'),
                        help='Repeat to exclude each former electrode center within the rule-derived separation')
    parser.add_argument('--max-sources',type=int,default=200)
    parser.add_argument('--iterations',type=int,default=1,
                        help='Repeatedly search and audit one new network at a time')
    parser.add_argument('--anchor-order',choices=('farthest_existing',
                        'closest_eligible_exit'),default='farthest_existing')
    parser.add_argument('--target-trials',type=int,default=1,
                        help='Try alternate geometrically reachable outlets if the nearest cannot form a complete Pad net')
    parser.add_argument('--max-outlets',type=int,
                        help='Evenly sample this many exits for a bounded exploratory run')
    parser.add_argument('--direct-exterior-only',action='store_true',
                        help='Search every selected exit with one joint exterior-to-Pad path')
    parser.add_argument('--capacity-audit-on-completion',action='store_true',
                        help='Recompute a candidate-independent bound for the best exported GDS')
    args=parser.parse_args()
    report_path=args.report.resolve(strict=True)
    report=json.loads(report_path.read_text(encoding='utf-8'))
    source=Path(report['input']).resolve(strict=True)
    destination=args.output_dir.resolve()
    if (args.iterations<1 or args.max_sources<1 or args.target_trials<1 or
            (args.max_outlets is not None and args.max_outlets<1)):
        parser.error('iterations, source/outlet limits and target trials must be positive')
    navigation_cache={}
    best_result=None
    for iteration in range(args.iterations):
        folder=(destination if args.iterations==1 else
                destination/f'iteration_{iteration+1:03d}')
        result=_augment_once(
            report,source,folder,candidate_step_um=args.candidate_step_um,
            external_portal_pitch_um=args.external_portal_pitch_um,
            external_portal_mode=args.external_portal_mode,
            exclude_source_um=args.exclude_source_um,
            max_sources=args.max_sources,navigation_cache=navigation_cache,
            anchor_order=args.anchor_order,target_trials=args.target_trials,
            max_outlets=args.max_outlets,
            direct_exterior_only=args.direct_exterior_only)
        if result['status']!='audited_additional_network':
            break
        best_result=result
        report=json.loads(Path(result['report']).read_text(encoding='utf-8'))
    if args.capacity_audit_on_completion and best_result is not None:
        from proof_portfolio import _physical_pad_policy
        from dataclasses import asdict
        if (report['rules']!=asdict(ProcessRules(
                minimum_center_spacing_um=
                    report['rules']['minimum_center_spacing_um'])) or
            _physical_pad_policy(report['pad_settings'])!=
                _physical_pad_policy(asdict(PadSettings()))):
            raise ValueError('Automatic capacity audit requires default physical rules')
        import subprocess
        import sys
        capacity=destination/'best_capacity_audit.json'
        subprocess.run([
            sys.executable,'-B',
            str(Path(__file__).with_name('run_capacity_audit.py')),
            '--input-dir',str(source.parent),'--input-names',source.name,
            '--output',str(capacity),'--layer',str(report['support_layer'][0]),
            '--datatype',str(report['support_layer'][1]),
            '--center-spacing-um',
            str(report['rules']['minimum_center_spacing_um']),
            '--reference-reports',best_result['report']],check=True)
        record=json.loads(capacity.read_text(encoding='utf-8'))['records'][0]
        if record['strict_integer_electrode_lower_bound']!=best_result['connected_pad_count']:
            raise RuntimeError('Independent capacity audit did not confirm best GDS')
        summary={'input':str(source),'input_sha256':report['input_sha256'],
                 'selected_report':best_result['report'],
                 'selected_output_sha256':best_result['output_sha256'],
                 'certified_lower_bound':best_result['connected_pad_count'],
                 'continuous_upper_bound':record['continuous_geometric_upper_bound'],
                 'declared_geometric_model_optimality_proven':
                     record['declared_geometric_model_optimality_proven'],
                 'capacity_certificate':str(capacity),
                 'capacity_certificate_sha256':sha256(capacity.read_bytes()).hexdigest()}
        (destination/'best_proof_summary.json').write_text(
            json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
        print(json.dumps(summary,ensure_ascii=False,indent=2),flush=True)


if __name__=='__main__':
    main()
