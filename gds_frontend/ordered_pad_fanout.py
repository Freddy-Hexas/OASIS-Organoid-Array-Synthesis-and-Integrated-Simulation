"""Jointly assign contiguous Pad banks and route an ordered annular fanout.

The annular isotopy preserves the cyclic order of source exits. This removes
spurious crossings caused by choosing an unrelated shortest path per Pad.
Finite-width feasibility is always checked on actual buffered metal. No
window is assigned unit capacity and no source-shape label is consulted.
"""
import math

import numpy as np
from shapely.geometry import GeometryCollection,LineString
from shapely.ops import unary_union

from curved_centerline import CurveSettings,UncertifiableCurve,smooth_centerline


SIDES=('top','right','bottom','left')


def _clockwise(point,origin):
    x,y=np.asarray(point)-origin
    return (math.pi/2-math.atan2(y,x))%(2*math.pi)


def _wrap(angle):
    return (angle+math.pi)%(2*math.pi)-math.pi


def ordered_pad_assignment(exits,frame):
    """Exact best consecutive block for each nearest-side order partition.

    The nearest-side partition is one declared constructive policy, not an
    exhaustive joint Pad-assignment search. Within it all possible contiguous
    blocks in the current frame are evaluated with squared angular cost.
    """
    origin=np.asarray(frame['origin_um']);banks={side:[] for side in SIDES}
    for i,exit_point in enumerate(exits):
        angle=_clockwise(exit_point,origin)
        # Boundary ties are resolved clockwise; no geometry-specific choice.
        side_index=int(math.floor((angle+math.pi/4)/(math.pi/2)))%4
        banks[SIDES[side_index]].append((side_index*math.pi/2+_wrap(angle-side_index*math.pi/2),i))
    # Size a uniformly expandable square for every bank, not just the first
    # overflowing side. Otherwise one expansion can leave a later bank full
    # and unnecessarily remove a perfectly feasible internal route.
    overflow=[(len(banks[side]),side) for side in SIDES
              if len(banks[side])>sum(slot['side']==side for slot in frame['slots'])]
    if overflow:
        required,side=max(overflow,key=lambda item:item[0])
        return None,{'reason':'ordered_bank_has_too_few_slots','side':side,'required':required,
                     'required_slots_by_side':{side:len(banks[side]) for side in SIDES}}
    assignment={};bank_record={}
    for side_index,side in enumerate(SIDES):
        sources=sorted(banks[side])
        slots=sorted([s for s in frame['slots'] if s['side']==side],key=lambda s:s['index'])
        n=len(sources)
        if not n:
            bank_record[side]={'count':0,'slot_indices':[]};continue
        center=side_index*math.pi/2
        angles=[center+_wrap(_clockwise(s['target_um'],origin)-center) for s in slots]
        alternatives=[]
        for start in range(len(slots)-n+1):
            cost=sum((a-angles[start+j])**2 for j,(a,_) in enumerate(sources))
            alternatives.append((cost,start))
        cost,start=min(alternatives)
        for j,(angle,i) in enumerate(sources):
            assignment[i]=(slots[start+j],angle,angles[start+j])
        bank_record[side]={'count':n,'slot_indices':[s['index'] for s in slots[start:start+n]],
                           'squared_angular_cost':cost}
    ordered=sorted((alpha,beta,i) for i,(_,alpha,beta) in assignment.items())
    if len(ordered)>1:
        gaps=[(ordered[(j+1)%len(ordered)][0]+(2*math.pi if j+1==len(ordered) else 0)-alpha,
               ordered[(j+1)%len(ordered)][1]+(2*math.pi if j+1==len(ordered) else 0)-beta)
              for j,(alpha,beta,_) in enumerate(ordered)]
        if any(a<=1e-10 or b<=1e-10 for a,b in gaps):
            return None,{'reason':'source_or_target_cyclic_order_is_degenerate'}
        minimum=min(min(a,b) for a,b in gaps)
    else:minimum=None
    return assignment,{'method':'cyclic_order_preserving_nearest_side_contiguous_bank_assignment',
        'banks':bank_record,'minimum_endpoint_angular_gap_rad':minimum,
        'all_contiguous_blocks_evaluated_within_partition':True,
        'globally_optimal_pad_assignment_proven':False}


def polar_fanout(alpha,beta,inner,outer,origin,chord_error):
    """G1 radial-endpoint spiral with a conservative second-derivative bound."""
    delta=beta-alpha;length=outer-inner
    bound=3*length*abs(delta)+6*outer*abs(delta)+2.25*outer*delta*delta
    segments=max(1,math.ceil(math.sqrt(bound/max(8*chord_error,1e-15))))
    t=np.linspace(0,1,segments+1)
    angle=alpha+delta*(3*t*t-2*t*t*t)
    radius=inner+length*t
    # Clockwise angles are measured from positive y, hence (sin,cos).
    points=origin+np.column_stack((radius*np.sin(angle),radius*np.cos(angle)))
    return points,{'method':'radial_cubic_angular_isotopy',
        'angle_start_clockwise_rad':alpha,'angle_end_clockwise_rad':beta,
        'radius_start_um':inner,'radius_end_um':outer,
        'second_derivative_norm_upper_um':bound,'polyline_segments':segments,
        'maximum_chord_error_um':bound/(8*segments*segments),
        'zero_width_disjointness':'strict cyclic order is preserved by the same monotone radial parameter; finite-width geometry is checked separately'}


def rebuild_ordered_outer(nav,routes,rules,settings,pad_settings,frame,policy,*,prepared_cache=None,
                          navigation_views=None):
    # Late import avoids a module cycle with the public Pad solver.
    from pad_router import _prepare_escape,_source_escape_cache
    from pad_sizing import pad_settings_for_frame
    pad_settings=pad_settings_for_frame(pad_settings,frame)
    cache=prepared_cache if prepared_cache is not None else {}
    views=tuple(navigation_views or (nav,))
    frame_key=(frame['square_side_um'],frame['inner_radius_um'],tuple(frame['origin_um']),
               pad_settings.pad_width_um,pad_settings.pad_length_um,pad_settings.pad_pitch_um)
    if cache.get('frame_key')!=frame_key:
        cache.clear();cache['frame_key']=frame_key
    prepared=[]
    for route in routes:
        outlet=route['outlet_node'];view_id=route.get('navigation_view',0)
        view=views[view_id]
        source_key=('source',view_id);key=(view_id,outlet)
        if source_key not in cache:cache[source_key]=_source_escape_cache(view,policy)
        if key not in cache:
            cache[key]=_prepare_escape(view,rules,pad_settings,frame,outlet,source_cache=cache[source_key])
        option,reason=cache[key]
        if option is None:return None,{'reason':reason,'outlet_node':outlet}
        prepared.append(option)
    assignment,diagnostic=ordered_pad_assignment([p['exit_point'] for p in prepared],frame)
    if assignment is None:return None,diagnostic
    guard=settings.numeric_guard_um;curve_settings=CurveSettings()
    inner=frame['inner_radius_um']+pad_settings.bridge_shell_overlap_um
    # Keep the entire interpolation circle strictly inside every Pad's inner
    # edge and inside the square. Final radial spokes meet only their Pad.
    outer=frame['square_side_um']/2-pad_settings.pad_outer_setback_um-pad_settings.pad_length_um-(
        rules.wire_width_um/2+rules.margin_um+guard+4*policy.grid_um)
    if outer<=inner:return None,{'reason':'no_common_annular_fanout_room'}
    result=[]
    for i,(old,escape) in enumerate(zip(routes,prepared)):
        slot,alpha,beta=assignment[i]
        points,analytic=polar_fanout(alpha,beta,inner,outer,policy.origin,
                                    curve_settings.max_chord_error_um)
        # The radial bridge reaches this same inner-circle point. All paths
        # then interpolate together to the same outer circle.
        vertices=[*escape['inside_line'].coords,*points,slot['target_um']]
        clean=[]
        for point in vertices:
            if not clean or np.linalg.norm(np.asarray(point)-clean[-1])>1e-7:clean.append(np.asarray(point))
        raw=LineString(clean)
        window=raw.buffer(pad_settings.bridge_width_um+rules.margin_um+rules.wire_width_um,quad_segs=4).union(
            slot['polygon'].buffer(rules.margin_um+rules.wire_width_um,quad_segs=4))
        support=unary_union([nav.support.intersection(window),frame['shell'].intersection(window),escape['bridge']])
        domain=support.buffer(-(rules.wire_width_um/2+rules.margin_um+guard),quad_segs=32)
        if not domain.covers(raw):return None,{'reason':'ordered_fanout_center_outside_support','pad_id':slot['pad_id']}
        try:line,curve=smooth_centerline(raw.coords,domain,curve_settings)
        except UncertifiableCurve:return None,{'reason':'ordered_fanout_uncertifiable_curve','pad_id':slot['pad_id']}
        allowed,reason=policy.check_exterior_line(nav.support,line)
        if not allowed:return None,{'reason':reason,'pad_id':slot['pad_id']}
        wire=line.buffer((rules.wire_width_um/2+curve_settings.max_chord_error_um)/math.cos(math.pi/64),quad_segs=16)
        if not support.covers(wire.buffer(rules.margin_um,quad_segs=32)):
            return None,{'reason':'ordered_fanout_wire_margin','pad_id':slot['pad_id']}
        view=views[old.get('navigation_view',0)]
        join=LineString([old['outlet_um'],view.graph.nodes[old['outlet_node']]['xy_um']]).buffer(
            rules.wire_width_um/2+.01,quad_segs=16)
        if not nav.support.covers(join.buffer(rules.margin_um,quad_segs=16)):
            return None,{'reason':'ordered_fanout_join_outside_source'}
        # Seed routes used during augmentation have no conductor until the
        # residual interior search provides it; only their outside is checked.
        central=old.get('central_wire')
        if central is not None:
            central_metal=central.union(old['electrode'])
            if not (central_metal.intersects(join) and join.intersects(wire)):
                return None,{'reason':'ordered_fanout_join_disconnected'}
            metal=GeometryCollection([central_metal,join,wire,slot['polygon']])
        else:metal=GeometryCollection([join,wire,slot['polygon']])
        if not wire.intersects(slot['polygon']):
            return None,{'reason':'ordered_fanout_pad_disconnected'}
        curve['ordered_fanout_analytic']=analytic
        result.append({**old,'pad_id':slot['pad_id'],'pad_side':slot['side'],'pad_index':slot['index'],
            'pad_target_um':slot['target_um'],'pad':slot['polygon'],'bridge':escape['bridge'],
            'portal_escape_um':escape['exit_point'].tolist(),'outer_line':line,'outer_wire':wire,
            'outer_curve':curve,'outer_length_um':line.length,'join_wire':join,'metal':metal,
            'outer_route_strategy':'joint_cyclic_order_annular_fanout',
            'complete_length_um':old.get('line_length_um',0)+line.length})
    # Full metal checks also cover a different network's Pad rectangle.
    for i,route in enumerate(result):
        for j,other in enumerate(result[:i]):
            if route['metal'].distance(other['metal'])<rules.spacing_um+guard:
                return None,{'reason':'ordered_fanout_full_metal_conflict',
                             'route_indices':[i,j],
                             'pad_ids':[route['pad_id'],other['pad_id']]}
    diagnostic.update({'status':'geometry_checked','route_count':len(result),
                       'inner_radius_um':inner,'outer_radius_um':outer,
                       'finite_width_checked':True})
    return result,diagnostic
