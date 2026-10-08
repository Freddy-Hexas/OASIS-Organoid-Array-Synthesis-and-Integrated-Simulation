"""Four-side pad routing on geometry-derived portals.

Every accepted column contains an electrode, its original-support route, a
certified escape to the exterior, an added substrate bridge, and one real pad.
The initial finite library is heuristic. Continuous residual additions reopen
the outside Pad assignment jointly. Only exported, re-read nets are a
constructive lower bound. No input filename or generator metadata is used.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, replace
from hashlib import sha256
from pathlib import Path
from tempfile import TemporaryDirectory
import json
import math

import gdstk
import networkx as nx
import numpy as np
import shapely
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra
from shapely.geometry import GeometryCollection, LineString, Point, Polygon, box
from shapely.ops import nearest_points, unary_union
from shapely.strtree import STRtree

from curved_centerline import CurveSettings, UncertifiableCurve, smooth_centerline
from capacity_bounds import (gds_convex_hull_cut_upper_bound,
                             gds_radial_cut_upper_bound,
                             verify_pads_outside_hull_cut,
                             verify_pads_outside_radial_cut)
from frontend import _polygonal, read_support
from island_router import (continuous_center_upper_bound, _route_columns, _select, _shorten,
                           betti, disk, build_navigation, _external_boundary_samples)
from outer_exit_policy import (make_outer_exit_policy, policy_from_record,
                               outside_circle_edge_intervals)
from center_compaction import compact_centers, radial_statistics
from joint_port_augment import augment_joint_ports
from port_coverage import coverage_record, outer_port_groups
from joint_path_flow import joint_path_routes
from navigation_simplification import SearchDeadlineExceeded
from route_incumbent import recertify_incumbent
from pad_sizing import pad_settings_for_frame,size_pad_bank,validate_pad_sizing
from electrode_region import region_for_navigation

PAD_LAYOUT_REVISION='center_preferred_pad_routing_v14_no_default_time_limits'


@dataclass(frozen=True)
class PadSettings:
    # Reference minimum footprint and dimensions; the actual side and the
    # number of slots are derived from geometry and may grow during routing.
    square_side_um: float = 32000.0
    minimum_inner_radius_um: float = 4000.0
    source_to_shell_gap_um: float = 200.0
    pad_width_um: float = 500.0
    pad_length_um: float = 3000.0
    pad_pitch_um: float = 1000.0
    pad_outer_setback_um: float = 100.0
    bridge_width_um: float = 40.0
    bridge_shell_overlap_um: float = 100.0
    columns_per_pad: int = 12
    deep_extra_columns_per_pad: int = 0
    # In addition to the shortest routes, retain distinct anchors nearest to
    # the original support's minimum-enclosing-circle center. This preserves the old
    # finite library while making central placements available to the MILP.
    center_extra_columns_per_pad: int = 12
    # None tries every slot on a geometrically eligible side.  A small local
    # shortlist can make a contiguous bank impossible before optimization.
    pad_choices_per_portal: int | None = None
    sizing_mode: str = 'expand_frame'
    # None preserves the user's preferred size as its lower limit.
    minimum_pad_width_um: float | None = None
    minimum_pad_length_um: float | None = None


def _side_point(side, u, v, origin, half):
    x, y = origin
    return {'top': (x + u, y + half - v),
            'right': (x + half - v, y - u),
            'bottom': (x - u, y - half + v),
            'left': (x - half + v, y + u)}[side]


def make_pad_frame(support, settings=PadSettings(), *, minimum_side_um=None,
                   spacing_um=4.0, margin_um=4.0, wire_width_um=5.0,
                   exit_policy=None,required_slots_per_side=None):
    """Derive every legal Pad slot from the current square and corner keepout.

    The 32 mm square and 0.5 x 3 mm Pad are starting dimensions, not a fixed
    bank of 28 positions.  Larger requested sides create more slots at the
    same process pitch; orthogonal banks have a geometric corner clearance.
    """
    if (settings.pad_width_um<=0 or settings.pad_length_um<=0 or
            settings.pad_pitch_um<=0 or settings.pad_outer_setback_um<0 or
            settings.square_side_um<=0 or spacing_um<0 or margin_um<0):
        raise ValueError('Pad dimensions, pitch, and process clearances are invalid')
    if settings.pad_pitch_um<settings.pad_width_um+spacing_um:
        raise ValueError('Pad pitch cannot provide the required inter-Pad spacing')
    if settings.pad_outer_setback_um<margin_um:
        raise ValueError('Pad outer setback is smaller than the support edge margin')
    policy=exit_policy or make_outer_exit_policy(support)
    origin=policy.origin
    radius=policy.radius_ticks*policy.grid_um
    opening=max(settings.minimum_inner_radius_um,
                radius+settings.source_to_shell_gap_um)
    validate_pad_sizing(settings,wire_width_um=wire_width_um,margin_um=margin_um)
    # Quantize even an initially unsaturated bank. Half-grid rectangle edges
    # must not round below the requested dimensions, and a reference pitch at
    # the process limit must also fit the same guard used by metal checking.
    settings,minimum_side_um,sizing=size_pad_bank(settings,opening,
        required_slots_per_side if required_slots_per_side is not None else 1,
        max(settings.square_side_um,float(minimum_side_um or 0)),
        spacing_um=spacing_um,margin_um=margin_um,grid_um=policy.grid_um,
        wire_width_um=wire_width_um,
        clearance_guard_um=policy.numeric_guard_um+2*policy.grid_um)
    sizing['sizing_trigger']='bank_demand' if required_slots_per_side is not None else 'initial_quantization'
    corner_keepout=settings.pad_outer_setback_um+settings.pad_length_um+spacing_um
    side=max(settings.square_side_um,
             float(minimum_side_um or 0),
             2*(opening+settings.pad_length_um+max(1000.0,margin_um)),
             2*corner_keepout+settings.pad_width_um)
    half=side/2
    available=side-2*corner_keepout
    slots_per_side=1+math.floor((available-settings.pad_width_um+1e-8)/settings.pad_pitch_um)
    if slots_per_side<1:
        raise ValueError('Square cannot contain a Pad after corner clearance')
    if settings.bridge_width_um<wire_width_um+2*margin_um:
        raise ValueError('Bridge is narrower than the wire with its required edge margins')
    shell=box(origin[0]-half,origin[1]-half,origin[0]+half,origin[1]+half).difference(
        Point(origin).buffer(opening,quad_segs=128))
    slots=[]
    for side_name in ('top','right','bottom','left'):
        for index in range(slots_per_side):
            u=(index-(slots_per_side-1)/2)*settings.pad_pitch_um
            a=settings.pad_width_um/2
            v0=settings.pad_outer_setback_um
            v1=v0+settings.pad_length_um
            pad=Polygon([_side_point(side_name,u-a,v0,origin,half),
                         _side_point(side_name,u+a,v0,origin,half),
                         _side_point(side_name,u+a,v1,origin,half),
                         _side_point(side_name,u-a,v1,origin,half)])
            target=_side_point(side_name,u,v1-min(100.0,settings.pad_length_um/4),origin,half)
            slots.append({'pad_id':f'{side_name}-{index+1:02d}','side':side_name,
                          'index':index+1,'u_um':u,'polygon':pad,'target_um':target})
    if any(not shell.buffer(-margin_um).covers(s['polygon']) for s in slots):
        raise ValueError('Pad layout violates the substrate edge margin')
    return {'origin_um':origin.tolist(),'square_side_um':side,
            'inner_radius_um':opening,'source_max_radius_um':radius,
            'pad_slots_per_side':slots_per_side,'corner_keepout_um':corner_keepout,
            'pad_pitch_um':settings.pad_pitch_um,
            'effective_pad_dimensions':{key:getattr(settings,key) for key in
                ('pad_width_um','pad_length_um','pad_pitch_um')},'pad_sizing':sizing,
            'shell':shell,'slots':slots}


def _visible_escape(nav, point, source_cache=None):
    """Propose a source-contained continuation to a permitted outer window.

    The complete wire and margin are certified by the outer option checker.
    No proposal is repaired by adding substrate across an original hole.
    """
    # This direct continuation is a proposal in the original support; it can
    # never select an inner component's local radius or a hard-coded band.
    return _nearest_exterior_escape(nav,point,source_cache)


def _source_escape_cache(nav,exit_policy=None):
    components=_polygonal(nav.support)
    holes=[Polygon(ring) for component in components for ring in component.interiors]
    stored=getattr(nav,'summary',{}).get('outer_exit_policy')
    policy=(exit_policy or (policy_from_record(stored) if stored else
                           make_outer_exit_policy(nav.support)))
    edges=[]
    for component in components:
        coords=np.asarray(component.exterior.coords)
        for i in range(len(coords)-1):
            for a,b in policy.edge_intervals(coords[i],coords[i+1],extra_guard_um=4*policy.grid_um):
                if np.linalg.norm(b-a)>1e-9:edges.append(LineString([a,b]))
    return {'exteriors':[LineString(component.exterior.coords) for component in components],
            'support_guard':nav.support.buffer(.002),
            'hole_tree':STRtree(holes) if holes else None,
            'exit_policy':policy,'allowed_edges':edges,'allowed_edge_tree':STRtree(edges)}


def _nearest_exterior_escape(nav,point,source_cache=None):
    """Join a sampled exterior wire-center portal to its source-support edge."""
    cache=source_cache if source_cache is not None else _source_escape_cache(nav)
    p=Point(point)
    candidates=[]
    tree=cache['allowed_edge_tree']
    if not cache['allowed_edges']:return None,None
    reach=cache['exit_policy'].bridge_width_um+2*cache['exit_policy'].center_inset_um
    indices=tree.query(p.buffer(reach))
    if not len(indices):indices=[int(tree.nearest(p))]
    for index in indices:
        outer=cache['allowed_edges'][int(index)]
        edge=nearest_points(p,outer)[1]
        candidates.append((p.distance(edge),edge))
    for _,edge in sorted(candidates,key=lambda item:item[0]):
        line=LineString([point,(edge.x,edge.y)])
        if line.length>1e-6 and cache['support_guard'].covers(line) and cache['exit_policy'].allows_exit((edge.x,edge.y)):
            # An exterior arc can face back into a strut. Search another
            # nearby cap rather than discard the entire portal after choosing
            # the first geometrically close, but inward-facing source edge.
            direction=np.asarray((edge.x,edge.y))-cache['exit_policy'].origin
            length=np.linalg.norm(direction)
            if length<1e-9:continue
            far=cache['exit_policy'].origin+direction/length*(
                cache['exit_policy'].radius_ticks*cache['exit_policy'].grid_um+2*reach)
            if LineString([(edge.x,edge.y),far]).intersection(nav.support).length>.002:
                continue
            return np.asarray((edge.x,edge.y)),line
    return None,None


class _CollectorEscapeGraph:
    """Multi-source shortest paths through the exact polygonal collector.

    The old graph stops at the wide region.  This auxiliary dual navigates
    that region and is built identically for every GDS.  It is only a path
    proposal; subsequent continuous containment checks certify each route.
    """
    def __init__(self,nav,source_cache=None):
        domain=nav.center_domain.intersection(nav.collector)
        triangles=list(shapely.constrained_delaunay_triangles(domain).geoms)
        self.domain=domain;self.triangles=triangles;self.tree=STRtree(triangles)
        self.center_check=nav.center_domain.buffer(.002)
        shapely.prepare(self.center_check)
        self.centers=np.asarray([np.asarray(t.exterior.coords)[:3].mean(axis=0) for t in triangles])
        ownership={};rows=[];cols=[];weights=[];portals={}
        for i,tri in enumerate(triangles):
            verts=np.asarray(tri.exterior.coords)[:3]
            for j in range(3):
                a=tuple(verts[j]);b=tuple(verts[(j+1)%3]);key=tuple(sorted((a,b)))
                if key in ownership:
                    other=ownership.pop(key)
                    midpoint=(np.asarray(a)+np.asarray(b))/2
                    length=np.linalg.norm(self.centers[i]-midpoint)+np.linalg.norm(self.centers[other]-midpoint)
                    rows.extend((i,other));cols.extend((other,i));weights.extend((length,length))
                    portals[(i,other)]=portals[(other,i)]=midpoint
                else:ownership[key]=i
        self.portals=portals
        cache=source_cache or _source_escape_cache(nav)
        policy=cache['exit_policy'];self.terminal_exits={}
        samples,_=_external_boundary_samples(domain,ownership,policy.origin,
                    policy.wire_width_um+4.0,exit_policy=policy)
        for owner,point in samples:
            exit_point,inside=_nearest_exterior_escape(nav,point,cache)
            if exit_point is None:continue
            length=float(np.linalg.norm(self.centers[owner]-point))+inside.length
            if owner in self.terminal_exits and length>=self.terminal_exits[owner][2]:continue
            self.terminal_exits[owner]=(point,exit_point,length)
        terminals=list(self.terminal_exits)
        if not len(terminals):raise ValueError('Collector does not reach the exterior')
        super_id=len(triangles)
        rows.extend([super_id]*len(terminals));cols.extend(terminals)
        weights.extend(self.terminal_exits[i][2] for i in terminals)
        matrix=coo_matrix((weights,(rows,cols)),shape=(super_id+1,super_id+1)).tocsr()
        _,predecessors=dijkstra(matrix,directed=True,indices=super_id,return_predecessors=True)
        self.predecessors=predecessors;self.super_id=super_id

    def escape(self,nav,point):
        i=int(self.tree.nearest(Point(point)))
        if self.predecessors[i]<0:return None,None
        chain=[i]
        while chain[-1]!=self.super_id and len(chain)<len(self.triangles)+1:
            chain.append(int(self.predecessors[chain[-1]]))
        if chain[-1]!=self.super_id:return None,None
        chain=chain[:-1]
        coords=[np.asarray(point),self.centers[chain[0]]]
        for a,b in zip(chain[:-1],chain[1:]):
            coords.extend((self.portals[(a,b)],self.centers[b]))
        line=LineString(coords)
        if not self.center_check.covers(line):return None,None
        point,exit_point,_=self.terminal_exits[chain[-1]]
        coords.append(point)
        result=_shorten(coords,self.center_check,5,visibility_first=True)
        if not self.center_check.covers(result):return None,None
        result=LineString([*result.coords,exit_point])
        if not nav.support.buffer(.002).covers(result):return None,None
        return np.asarray(exit_point),result


def _pad_choices(exit_point,frame,settings):
    d=np.asarray(exit_point)-np.asarray(frame['origin_um'])
    half=frame['square_side_um']/2
    sides=[]
    if d[1]>=0 and abs(d[1])>=.7*abs(d[0]):sides.append(('top',d[0]/max(abs(d[1]),1e-9)))
    if d[0]>=0 and abs(d[0])>=.7*abs(d[1]):sides.append(('right',-d[1]/max(abs(d[0]),1e-9)))
    if d[1]<0 and abs(d[1])>=.7*abs(d[0]):sides.append(('bottom',-d[0]/max(abs(d[1]),1e-9)))
    if d[0]<0 and abs(d[0])>=.7*abs(d[1]):sides.append(('left',d[1]/max(abs(d[0]),1e-9)))
    if not sides:
        sides=[('top' if d[1]>=0 else 'bottom',0)]
    inward=half-settings.pad_outer_setback_um-settings.pad_length_um
    selected=[]
    for side,ratio in sides:
        expected=ratio*inward
        eligible=[slot for slot in frame['slots'] if slot['side']==side]
        ordered=sorted(eligible,key=lambda slot:abs(slot['u_um']-expected))
        selected.extend(ordered if settings.pad_choices_per_portal is None else
                        ordered[:settings.pad_choices_per_portal])
    return selected


def _pad_occupancy(chosen,settings):
    """Certificate that each side's used Pad slots have no internal vacancy."""
    by_side={side:[] for side in ('top','right','bottom','left')}
    for route in chosen:by_side[route['pad_side']].append(route['pad_index'])
    result={}
    for side,indices in by_side.items():
        indices.sort()
        consecutive=(len(indices)==len(set(indices)) and
                     (not indices or indices[-1]-indices[0]+1==len(indices)))
        if not consecutive:
            raise RuntimeError(f'{side} Pad bank has unoccupied slots between connected Pads: {indices}')
        result[side]={'connected_pad_count':len(indices),'slot_indices':indices,
                      'first_slot':indices[0] if indices else None,
                      'last_slot':indices[-1] if indices else None,
                      'internal_empty_slots':0,
                      'adjacent_center_pitch_um':settings.pad_pitch_um if len(indices)>1 else None,
                      'adjacent_edge_gap_um':settings.pad_pitch_um-settings.pad_width_um if len(indices)>1 else None}
    return result


def _prepare_escape(nav,rules,settings,frame,outlet,escape_graph=None,
                    source_cache=None):
    stored=getattr(nav,'summary',{}).get('outer_exit_policy')
    cache=source_cache if source_cache is not None else _source_escape_cache(nav,
        policy_from_record(stored) if stored else
        make_outer_exit_policy(nav.support,wire_width_um=rules.wire_width_um,
            margin_um=rules.margin_um,bridge_width_um=settings.bridge_width_um,
            grid_um=getattr(nav,'summary',{}).get('gds_native_precision_m',1e-9)*1e6))
    policy=cache['exit_policy']
    start=np.asarray(nav.graph.nodes[outlet]['xy_um'],dtype=float)
    external=nav.graph.nodes[outlet].get('terminal_kind')=='external_boundary'
    exit_point,inside=(_nearest_exterior_escape(nav,start,cache) if external else
                       escape_graph.escape(nav,start) if escape_graph is not None else
                       _visible_escape(nav,start,cache))
    if exit_point is None:return None,'no_certified_exterior_escape'
    if not policy.allows_exit(exit_point):return None,'exit_not_in_outer_extremity_window'
    origin=policy.origin
    direction=exit_point-origin
    radius=np.linalg.norm(direction)
    if radius<1:return None,'undefined_outward_direction'
    shell_point=origin+direction/radius*(frame['inner_radius_um']+settings.bridge_shell_overlap_um)
    if LineString([exit_point,shell_point]).intersection(nav.support).length>.002:
        return None,'bridge_reenters_original_support'
    bridge=LineString([exit_point,shell_point]).buffer(settings.bridge_width_um/2,cap_style='round')
    if bridge.intersection(nav.support).area<=1 or bridge.intersection(frame['shell']).area<=1:
        return None,'bridge_does_not_join_both_supports'
    added=bridge.difference(nav.support)
    if cache['hole_tree'] is not None and len(cache['hole_tree'].query(
            added,predicate='intersects')):
        return None,'bridge_enters_protected_source_hole'
    allowed,reason=policy.check_bridge(nav.support,bridge,exit_point=exit_point)
    if not allowed:return None,reason
    return {'start':start,'inside_line':inside,'exit_point':exit_point,'shell_point':shell_point,
            'bridge':bridge,'exit_policy':policy},None


def _make_outer_option(nav,rules,settings,frame,outlet,slot,prepared):
    start=prepared['start'];exit_point=prepared['exit_point']
    shell_point=prepared['shell_point'];bridge=prepared['bridge']
    target=np.asarray(slot['target_um'],dtype=float)
    vertices=[*prepared['inside_line'].coords,shell_point,target]
    points=[np.asarray(vertices[0])]
    for point in vertices[1:]:
        point=np.asarray(point)
        if np.linalg.norm(point-points[-1])>1e-6:points.append(point)
    raw=LineString(points)
    window=unary_union([raw.buffer(100,quad_segs=4),slot['polygon'].buffer(100,quad_segs=4)])
    support=unary_union([nav.support.intersection(window),
                         frame['shell'].intersection(window),bridge])
    inset=rules.wire_width_um/2+rules.margin_um+.05
    center=support.buffer(-inset,quad_segs=16)
    if not center.buffer(.002).covers(raw):
        return None,'outer_line_exits_supported_wire_domain'
    try:
        line,curve=smooth_centerline(raw.coords,center,CurveSettings())
    except UncertifiableCurve:
        return None,'outer_curve_uncertifiable'
    allowed,reason=prepared['exit_policy'].check_exterior_line(nav.support,line)
    if not allowed:return None,reason
    radius_wire=rules.wire_width_um/2+CurveSettings().max_chord_error_um
    wire=line.buffer(radius_wire/math.cos(math.pi/64),quad_segs=16)
    if not support.covers(wire.buffer(rules.margin_um,quad_segs=16)):
        return None,'outer_wire_margin'
    if not wire.intersects(slot['polygon']) or not support.buffer(-rules.margin_um+.01).covers(slot['polygon']):
        return None,'pad_contact_or_margin'
    return {'pad_id':slot['pad_id'],'pad_side':slot['side'],'pad_index':slot['index'],
            'pad':slot['polygon'],'pad_target_um':slot['target_um'],
            'bridge':bridge,'outer_wire':wire,'outer_line':line,
            'outer_curve':curve,'outer_length_um':line.length,
            'portal_escape_um':exit_point.tolist()},None


def solve_with_pads(nav,rules,*,settings,pad_settings=PadSettings(),progress=None,
                    additional_navigations=(),incumbent_records=()):
    output=progress or (lambda *args:None)
    last_progress=0.0
    def notify(stage,value):
        nonlocal last_progress
        last_progress=max(last_progress,float(value))
        output(stage,last_progress)
    if pad_settings.pad_choices_per_portal is not None and pad_settings.pad_choices_per_portal<1:
        raise ValueError('pad_choices_per_portal must be positive or None')
    validate_pad_sizing(pad_settings,wire_width_um=rules.wire_width_um,margin_um=rules.margin_um)
    policy=make_outer_exit_policy(nav.support,wire_width_um=rules.wire_width_um,
        margin_um=rules.margin_um,numeric_guard_um=settings.numeric_guard_um,
        bridge_width_um=pad_settings.bridge_width_um,
        grid_um=nav.summary.get('gds_native_precision_m',1e-9)*1e6)
    # Pad routing always has a full-support geometric frontend. Legacy
    # caller views are rebuilt from the GDS, never from generator metadata.
    if (nav.summary.get('input_gds_path') and
            (nav.summary.get('outlet_mode')!='external_boundary' or
             nav.summary.get('outer_exit_policy')!=policy.record())):
        notify('按原始 GDS 最外端窗口重建 Pad 导航',.12)
        nav=build_navigation(nav.summary['input_gds_path'],rules,
            layer=nav.summary['support_layer'][0],datatype=nav.summary['support_layer'][1],
            outlet_mode='external_boundary',settings=settings,outer_exit_policy=policy,
            progress=lambda stage,value:notify(stage,.12+.2*value))
    island_radius=(rules.electrode_diameter_um/2+rules.margin_um+
                   settings.numeric_guard_um)
    required_center_gap=max(rules.minimum_center_spacing_um,
                            rules.electrode_diameter_um+rules.spacing_um,
                            2*island_radius+settings.pad_gap_um if settings.pad_gap_um>0 else 0)
    center_bound=continuous_center_upper_bound(nav,required_center_gap,rules=rules)
    navigations=(nav,*additional_navigations)
    if any(other.support.symmetric_difference(nav.support).area>1e-5
           for other in navigations[1:]):
        raise ValueError('All navigation views must use the same original support')
    frame=make_pad_frame(nav.support,pad_settings,spacing_um=rules.spacing_um,
                         margin_um=rules.margin_um,wire_width_um=rules.wire_width_um,exit_policy=policy)
    constrained=[]
    for view in navigations:
        if not hasattr(view,'graph'):
            constrained.append(view);continue
        cache=_source_escape_cache(view,policy)
        try:
            escape_graph=_CollectorEscapeGraph(view,cache) if not view.collector.is_empty else None
        except ValueError:
            escape_graph=None
        approved=[];rejections=Counter();exit_points=[]
        for outlet in view.outlets:
            prepared,reason=_prepare_escape(view,rules,pad_settings,frame,outlet,
                                            escape_graph,cache)
            if prepared is None:rejections[reason]+=1
            else:
                approved.append(outlet);exit_points.append(prepared['exit_point'].tolist())
        certificate={'candidate_terminal_count':len(view.outlets),
                     'approved_terminal_count':len(approved),
                     'rejections':dict(rejections),'approved_source_exits_um':exit_points,
                     'applied_before_shortest_path_search':True,
                     'empty_scope':'no permitted exit found in this finite navigation view; not a zero-capacity proof'}
        constrained.append(replace(view,outlets=approved,
            summary={**view.summary,'outer_exit_policy':policy.record(),
                     'pad_exit_filter':certificate}))
    navigations=tuple(constrained);nav=navigations[0]
    view_columns=[];view_info=[]
    joint_routes=[];joint_diagnostic=None
    try:
        joint_routes,joint_diagnostic=joint_path_routes(nav,rules,settings,pad_settings,frame,policy,notify)
        frame=joint_diagnostic.pop('_frame',frame)
    except SearchDeadlineExceeded as exc:
        joint_diagnostic={'status':'joint_proposal_deadline','stage':str(exc),'global_maximum_proven':False}
    incumbent_routes,incumbent_frame,incumbent_diagnostic=recertify_incumbent(
        incumbent_records,nav,rules,settings,pad_settings,frame,policy,
        minimum_count=len(joint_routes),progress=notify)
    if incumbent_routes:
        def construction_key(pool):
            return (len(pool),-sum(float(np.sum((np.asarray(c['source_um'])-policy.origin)**2)) for c in pool))
        if construction_key(incumbent_routes)>construction_key(joint_routes):
            joint_routes=incumbent_routes;frame=incumbent_frame
    for view_id,view in enumerate(navigations):
        def view_progress(stage,value):
            if .4<=value<=.52:
                value=.4+.12*(view_id+(value-.4)/.12)/len(navigations)
            notify(stage,value)
        # Joint certified proposals supply a diverse incumbent. For every
        # input the same independent-column/residual construction stages remain
        # available; no structure label determines the solver.
        if joint_routes:
            columns=[];info={'status':'joint_geometry_checked_incumbent','joint_diagnostic':joint_diagnostic}
        else:
            columns,info=_route_columns(view,rules,settings,view_progress)
        columns=[{**c,'navigation_view':view_id,
                  'source_node':(view_id,c['source_node'])} for c in columns]
        view_columns.append(columns)
        view_info.append({'view_id':view_id,'outlet_mode':view.summary['outlet_mode'],
                          'pad_exit_filter':view.summary.get('pad_exit_filter'),
                          'candidate_library':info})
    source_count=max(len(joint_routes),len({column['source_node'] for columns in view_columns for column in columns}))
    cut_bound=None
    hull_cut_bound=None
    input_path=nav.summary.get('input_gds_path')
    if input_path:
        try:
            candidate=gds_radial_cut_upper_bound(input_path,*nav.summary['support_layer'],
                                                   rules.wire_width_um,rules.spacing_um)
            if candidate['input_sha256']!=nav.summary['sha256']:
                raise ValueError('Input GDS changed after navigation construction')
            cut_origin=Point(candidate['circle_center_um'])
            cut_radius=candidate['circle_radius_um']
            pad_clearance=min(slot['polygon'].distance(cut_origin)-cut_radius
                              for slot in frame['slots'])
            if verify_pads_outside_radial_cut(
                    candidate,[slot['polygon'].bounds for slot in frame['slots']]):
                candidate['pad_outside_cut_verified']=True
                candidate['minimum_pad_to_cut_clearance_um']=pad_clearance
                candidate['pad_exclusion_method']='exact_integer_circle_vs_outward_rounded_pad_boxes'
                cut_bound=candidate
            hull=gds_convex_hull_cut_upper_bound(input_path,*nav.summary['support_layer'],
                                                  rules.wire_width_um,rules.spacing_um)
            if hull['input_sha256']!=nav.summary['sha256']:
                raise ValueError('Input GDS changed after navigation construction')
            exclusion_center=Point(hull['pad_exclusion_circle_center_um'])
            exclusion_radius=hull['pad_exclusion_circle_radius_um']
            hull_pad_clearance=min(slot['polygon'].distance(exclusion_center)-exclusion_radius
                                   for slot in frame['slots'])
            if verify_pads_outside_hull_cut(
                    hull,[slot['polygon'].bounds for slot in frame['slots']]):
                hull['pad_outside_cut_verified']=True
                hull['minimum_pad_to_exclusion_circle_um']=hull_pad_clearance
                hull['pad_exclusion_method']='exact_integer_hull_vs_outward_rounded_pad_boxes'
                hull_cut_bound=hull
        except ValueError as exc:
            if 'changed after navigation' in str(exc):
                raise
    bounds=[b['value'] for b in (center_bound,cut_bound,hull_cut_bound) if b is not None]
    capacity_bound={'value':min(bounds)} if bounds else None
    best=None;trials=[];previous_count=None
    while True:
        if cut_bound:
            if not verify_pads_outside_radial_cut(
                    cut_bound,[slot['polygon'].bounds for slot in frame['slots']]):
                raise RuntimeError('Expanded Pad frame violates the certified radial cut')
        if hull_cut_bound:
            if not verify_pads_outside_hull_cut(
                    hull_cut_bound,[slot['polygon'].bounds for slot in frame['slots']]):
                raise RuntimeError('Expanded Pad frame violates the certified convex-hull cut')
        result=_solve_on_pad_frame(navigations,rules,settings,pad_settings,frame,
                                   view_columns,view_info,capacity_bound,notify,
                                   **({'preselected':joint_routes} if joint_routes else {}))
        count=result['retained_routes']
        trials.append({'square_side_um':frame['square_side_um'],
                       'available_slots_per_side':frame['pad_slots_per_side'],
                       'selected_pad_count':count,
                       'saturated_sides':[side for side,bank in result['pad_bank_occupancy'].items()
                                          if bank['connected_pad_count']==frame['pad_slots_per_side']]})
        if best is None or count>best['retained_routes']:best=result
        # Expand only when the current bank itself is demonstrably full.  The
        # finite source library and any certified packing bound stop expansion.
        saturated=bool(trials[-1]['saturated_sides'])
        if (not saturated or count>=source_count or
                (capacity_bound and count>=capacity_bound['value']) or
                (previous_count is not None and count<=previous_count)):
            break
        previous_count=count
        next_slots=min(source_count,max(frame['pad_slots_per_side']+1,
                                        math.ceil(frame['pad_slots_per_side']*1.5)))
        enlarged=make_pad_frame(nav.support,pad_settings,minimum_side_um=frame['square_side_um'],
                                spacing_um=rules.spacing_um,margin_um=rules.margin_um,
                                wire_width_um=rules.wire_width_um,exit_policy=policy,
                                required_slots_per_side=next_slots)
        if enlarged['pad_slots_per_side']<=frame['pad_slots_per_side']:
            break
        frame=enlarged
    if best['_chosen']:
        before=best['center_preference']['selected']
        compacted,certificate=compact_centers(nav.support,best['_chosen'],rules,settings,
            policy.origin,grid_um=policy.grid_um,progress=notify)
        best['_chosen']=compacted
        best['center_preference']['finite_library_selected']=before
        best['center_preference']['selected']=certificate['after']
        best['center_compaction']=certificate
        best['optimization']['secondary_result_scope']='finite candidate selection before continuous residual refinement'
        best['optimization']['continuous_center_optimality_proven']=False
    if hasattr(nav,'graph'):
        # Count takes priority over centrality. Reopen exterior assignment for
        # every added route instead of freezing a locally chosen Pad slot.
        # The finite-library dual bound has no authority over this new space.
        active_frame=make_pad_frame(nav.support,pad_settings_for_frame(pad_settings,best['frame']),
            minimum_side_um=best['frame']['square_side_um'],spacing_um=rules.spacing_um,
            margin_um=rules.margin_um,wire_width_um=rules.wire_width_um,exit_policy=policy)
        active_frame['pad_sizing']=best['frame'].get('pad_sizing')
        initial_optimization=dict(best['optimization'])
        initial_optimization['scope']='only the initial pruned finite column library; does not bound dynamically generated routes or the continuous problem'
        augmented,augmentation=augment_joint_ports(nav,best['_chosen'],rules,settings,
            pad_settings,active_frame,policy,progress=notify,
            capacity_upper=capacity_bound['value'] if capacity_bound else None,
            navigation_views=navigations)
        active_frame=augmentation.pop('_frame')
        for bound,verify in ((cut_bound,verify_pads_outside_radial_cut),
                             (hull_cut_bound,verify_pads_outside_hull_cut)):
            if bound and not verify(bound,[s['polygon'].bounds for s in active_frame['slots']]):
                raise RuntimeError('Augmented Pad frame violates a certified continuous cut')
        best['initial_finite_library_optimization']=initial_optimization
        best['port_augmentation']=augmentation
        best['_chosen']=augmented
        if augmentation['accepted_updates']:
            best['initial_center_compaction']=best.get('center_compaction')
            best['center_preference']['pre_augmentation_selected']=best['center_preference']['selected']
            best['_chosen'],best['center_compaction']=compact_centers(nav.support,augmented,
                rules,settings,policy.origin,grid_um=policy.grid_um,progress=notify)
        best['_shell']=active_frame['shell']
        best['frame']={k:v for k,v in active_frame.items() if k not in ('shell','slots')}
        best['pad_slot_count']=len(active_frame['slots'])
        best['current_frame_pad_slot_upper_bound']=len(active_frame['slots'])
        best['optimization']={
            'status':'constructive_dynamic_route_generation',
            'finite_lower_bound':len(best['_chosen']),'finite_upper_bound':None,
            'finite_model_proven_optimal':False,'continuous_center_optimality_proven':False,
            'secondary_objective':'minimize_sum_squared_geometric_center_distance',
            'secondary_proven_optimal':False,
            'pad_contiguous_per_side':True,'shared_corridors_permitted':True,
            'unit_edge_or_node_capacity_imposed':False,
            'scope':'initial finite optimization followed by continuous residual route generation and joint exterior rematching; only complete constructions are retained; no exhaustive pricing or continuous optimum certificate'}
        mapping,groups=outer_port_groups(nav,policy)
        best['outer_port_coverage']=coverage_record(best['_chosen'],mapping,groups)
        if active_frame['square_side_um']!=trials[-1]['square_side_um']:
            trials.append({'phase':'joint_port_augmentation',
                'square_side_um':active_frame['square_side_um'],
                'available_slots_per_side':active_frame['pad_slots_per_side'],
                'selected_pad_count':len(best['_chosen'])})
    # Rebuild all public counts and geometry from the final complete layout.
    # This is also necessary when only centers changed, without an added net.
    if hasattr(nav,'graph'):
        chosen=best['_chosen']
        best['retained_routes']=best['selected_pad_count']=len(chosen)
        effective_pads=pad_settings_for_frame(pad_settings,active_frame)
        best['pad_settings']=asdict(effective_pads)
        best['pad_sizing']=active_frame.get('pad_sizing')
        best['status']='awaiting_export_audit' if chosen else 'no_finite_complete_pad_route'
        best['pad_bank_occupancy']=_pad_occupancy(chosen,effective_pads)
        best['selected_navigation_views']=dict(Counter(c['navigation_view'] for c in chosen))
        radius=float(np.linalg.norm(shapely.get_coordinates(nav.support)-policy.origin,axis=1).max())
        best['center_preference']['selected']=radial_statistics(chosen,policy.origin,radius)
        islands=unary_union([c['island'] for c in chosen]) if chosen else Polygon()
        source_with_islands=nav.support.union(islands)
        if betti(source_with_islands)!=betti(nav.support):
            raise RuntimeError('Refined electrode islands changed source topology')
        best['_final_support']=unary_union([source_with_islands,best['_shell'],*[c['bridge'] for c in chosen]])
        best['routes']=[{k:v for k,v in c.items() if k not in
            ('metal','wire','electrode','island','central_wire','join_wire','outer_wire',
             'outer_line','bridge','pad')} for c in chosen]
        best['substrate_topology_after']=betti(best['_final_support'])
    best['pad_frame_trials']=trials
    best['joint_path_flow']=joint_diagnostic
    best['rechecked_incumbent']=incumbent_diagnostic
    best['center_distance_upper_bound']=center_bound
    best['radial_routing_cut_upper_bound']=cut_bound
    best['convex_hull_routing_cut_upper_bound']=hull_cut_bound
    best['continuous_upper_bound']=capacity_bound['value'] if capacity_bound else None
    if best['continuous_upper_bound'] is not None and best['retained_routes']>best['continuous_upper_bound']:
        raise RuntimeError('Constructed Pad layout exceeds a continuous capacity certificate')
    best['pad_frame_policy']='reference minimum square; geometry-derived slots; attempt growth when a bank is saturated; keep best audited candidate if growth stalls'
    best['finite_frame_search_exhaustive']=False
    best['electrode_region']=region_for_navigation(nav,rules).audit(best['_chosen'])
    if not best['electrode_region']['passed']:
        raise RuntimeError('Complete Pad solution has an electrode outside its allowed region')
    return best


def _diverse_pad_columns(choices, short_limit, deep_extra):
    """Keep the complete old short list, then add distinct deep anchors."""
    if short_limit < 1 or deep_extra < 0:
        raise ValueError('columns_per_pad must be positive')
    short=sorted(choices,key=lambda c:(c['complete_length_um'],-c['depth_um']))
    deep=sorted(choices,key=lambda c:(-c['depth_um'],c['complete_length_um']))
    retained=[];seen=set()
    def add(pool,target):
        for candidate in pool:
            if len(retained)>=target:
                break
            source=candidate['source_node']
            if source in seen:
                continue
            seen.add(source);retained.append(candidate)
    add(short,short_limit)
    add(deep,short_limit+deep_extra)
    return retained


class _PadColumnReservoir:
    """Stream exact short, graph-deep and geometrically central choices."""
    def __init__(self,short_limit,deep_extra,center_extra=0,center_origin=None):
        if short_limit<1 or deep_extra<0 or center_extra<0:
            raise ValueError('Pad candidate budgets must be nonnegative, with at least one short route')
        if center_extra and center_origin is None:
            raise ValueError('A geometric center is required for central Pad candidates')
        self.short_limit=short_limit
        self.deep_extra=deep_extra
        self.center_extra=center_extra
        self.deep_limit=short_limit+deep_extra
        self.center_limit=short_limit+deep_extra+center_extra
        self.center_origin=(np.asarray(center_origin,dtype=float)
                            if center_origin is not None else None)
        self.short={}
        self.deep={}
        self.center={}
        self.count=0

    @staticmethod
    def _keep(table,limit,source,rank,column):
        if limit<=0:
            return
        prior=table.get(source)
        if prior is not None:
            if rank<prior[0]:
                table[source]=(rank,column)
            return
        if len(table)<limit:
            table[source]=(rank,column)
            return
        worst=max(table,key=lambda key:table[key][0])
        if rank<table[worst][0]:
            del table[worst]
            table[source]=(rank,column)

    def add(self,column):
        serial=self.count
        self.count+=1
        source=column['source_node']
        length=column['complete_length_um']
        depth=column['depth_um']
        self._keep(self.short,self.short_limit,source,
                   (length,-depth,serial),column)
        self._keep(self.deep,self.deep_limit,source,
                   (-depth,length,serial),column)
        if self.center_extra:
            point=np.asarray(column['source_um'],dtype=float)
            radius_squared=float(np.sum((point-self.center_origin)**2))
            self._keep(self.center,self.center_limit,source,
                       (radius_squared,length,serial),column)

    def selected(self):
        pool={id(column):(rank[-1],column) for rank,column in
              [*self.short.values(),*self.deep.values(),*self.center.values()]}
        ordered=[column for _,column in sorted(pool.values(),
                                              key=lambda item:item[0])]
        retained=_diverse_pad_columns(ordered,self.short_limit,self.deep_extra)
        if self.center_extra:
            seen={column['source_node'] for column in retained}
            for _,column in sorted(self.center.values(),key=lambda item:item[0]):
                if len(retained)>=self.short_limit+self.deep_extra+self.center_extra:
                    break
                if column['source_node'] not in seen:
                    retained.append(column)
                    seen.add(column['source_node'])
        return retained


def _solve_on_pad_frame(navigations,rules,settings,pad_settings,frame,view_columns,view_info,
                        capacity_bound,notify,preselected=()):
    pad_settings=pad_settings_for_frame(pad_settings,frame)
    nav=navigations[0]
    policy=(policy_from_record(nav.summary['outer_exit_policy']) if nav.summary.get('outer_exit_policy') else
            make_outer_exit_policy(nav.support,wire_width_um=rules.wire_width_um,
                margin_um=rules.margin_um,numeric_guard_um=settings.numeric_guard_um,
                bridge_width_um=pad_settings.bridge_width_um))
    center_origin=policy.origin
    source_radius=float(np.linalg.norm(shapely.get_coordinates(nav.support)-center_origin,
                                       axis=1).max())
    option_cache={};rejected=Counter();by_pad={}
    total_outlets=sum(len({c['outlet_node'] for c in columns}) for columns in view_columns)
    processed=0
    for view_id,(view,central) in enumerate(zip(navigations,view_columns) if not preselected else ()):
        used_outlets={c['outlet_node'] for c in central}
        source_cache=_source_escape_cache(view,policy) if used_outlets else None
        try:
            escape_graph=(_CollectorEscapeGraph(view,source_cache) if used_outlets and
                          not view.collector.is_empty else None)
        except ValueError:
            escape_graph=None
        for outlet in sorted(used_outlets):
            point=view.graph.nodes[outlet]['xy_um']
            prepared,reason=_prepare_escape(view,rules,pad_settings,frame,outlet,
                                            escape_graph,source_cache)
            if prepared is None:
                rejected[reason]+=1;continue
            for slot in _pad_choices(point,frame,pad_settings):
                option,reason=_make_outer_option(view,rules,pad_settings,frame,outlet,slot,prepared)
                if option is None:rejected[reason]+=1
                else:option_cache[(view_id,outlet,slot['pad_id'])]=option
            processed+=1
            if processed%50==0:
                notify('认证门户、承载桥与 Pad 接入',.54+.06*processed/max(1,total_outlets))
        for c in central:
            outlet=c['outlet_node']
            for slot in _pad_choices(view.graph.nodes[outlet]['xy_um'],frame,pad_settings):
                option=option_cache.get((view_id,outlet,slot['pad_id']))
                if option is None:continue
                join=LineString([c['outlet_um'],view.graph.nodes[outlet]['xy_um']])
                join_wire=join.buffer(rules.wire_width_um/2+.01,quad_segs=16)
                if not view.support.covers(join_wire.buffer(rules.margin_um,quad_segs=16)):
                    rejected['shifted_lane_cannot_join_portal']+=1;continue
                # Each member is connected.  These three intersections prove
                # their union is connected without overlaying the complex
                # GDS-derived polygons for every source/Pad combination.
                # A GeometryCollection has exactly the same point set as its
                # union for STRtree dwithin conflict tests and for GDS export.
                if not (c['metal'].intersects(join_wire) and
                        join_wire.intersects(option['outer_wire']) and
                        option['outer_wire'].intersects(option['pad'])):
                    rejected['complete_metal_disconnected']+=1;continue
                metal=GeometryCollection(
                    [c['metal'],join_wire,option['outer_wire'],option['pad']])
                full={**c,**option,'join_wire':join_wire,'metal':metal,
                      'central_wire':c['wire'],'depth_um':c['depth_um'],
                      'complete_length_um':c['line_length_um']+option['outer_length_um']}
                key=(option['pad_id'],view_id)
                if key not in by_pad:
                    by_pad[key]=_PadColumnReservoir(
                        pad_settings.columns_per_pad,
                        pad_settings.deep_extra_columns_per_pad,
                        pad_settings.center_extra_columns_per_pad,
                        center_origin)
                by_pad[key].add(full)
    # Prune the finite library per *physical* pad. Keep distinct source
    # anchors and genuinely reserve places for both short and deep routes.
    # This search heuristic is never treated as a continuous upper bound.
    columns=list(preselected)
    baseline_column_ids=list(range(len(columns)))
    for pad_id in dict.fromkeys(key[0] for key in by_pad):
        retained=[]
        for view_id in range(len(navigations)):
            reservoir=by_pad.get((pad_id,view_id))
            selected=reservoir.selected() if reservoir is not None else []
            start=len(columns)+len(retained)
            baseline_column_ids.extend(range(start,start+min(len(selected),pad_settings.columns_per_pad)))
            retained.extend(selected)
        columns.extend(retained)
    notify('联合选择电极、门户、桥与专属 Pad',.65)
    pad_order={side:[slot['pad_id'] for slot in frame['slots'] if slot['side']==side]
               for side in ('top','right','bottom','left')}
    chosen,optimization=_select(columns,rules,settings,notify,
                                capacity_bound['value'] if capacity_bound else None,
                                pad_order=pad_order,
                                 baseline_candidate_ids=baseline_column_ids,
                                 center_origin=center_origin)
    origin=center_origin
    radius=source_radius
    def radial_summary(pool):
        distances=[float(np.linalg.norm(np.asarray(c['source_um'])-origin))
                   for c in pool]
        return {'count':len(distances),
                'within_half_source_radius':sum(d<=radius/2+1e-7 for d in distances),
                'mean_radius_um':(float(np.mean(distances)) if distances else None),
                'maximum_radius_um':(max(distances) if distances else None)}
    center_preference={'reference_um':origin.tolist(),
                       'reference_definition':'minimum enclosing circle center of original GDS support, rounded to its native grid; same reference as outer ports and Pad frame',
                       'source_max_radius_um':radius,
                       'priority':'maximize complete connected nets, then minimize sum of squared electrode radii',
                       'candidate_after_pruning':radial_summary(columns),
                       'selected':radial_summary(chosen)}
    occupancy=_pad_occupancy(chosen,pad_settings)
    pad_ids=[c['pad_id'] for c in chosen]
    if len(set(pad_ids))!=len(pad_ids):raise RuntimeError('Pad assignment is not injective')
    islands=unary_union([c['island'] for c in chosen]) if chosen else Polygon()
    center_with_islands=nav.support.union(islands)
    if betti(center_with_islands)!=betti(nav.support):
        raise RuntimeError('Selected electrode islands changed the protected source topology')
    final=unary_union([center_with_islands,frame['shell'],*[c['bridge'] for c in chosen]])
    records=[]
    for c in chosen:
        records.append({k:v for k,v in c.items() if k not in
                        ('metal','wire','electrode','island','central_wire','join_wire','outer_wire',
                         'outer_line','bridge','pad')})
    # A bank is allowed to grow.  Current slot count therefore cannot be a
    # continuous capacity upper bound; only the independent electrode-center
    # packing argument applies to all later square sizes.
    upper=capacity_bound['value'] if capacity_bound else None
    return {'status':'awaiting_export_audit' if chosen else 'no_finite_complete_pad_route',
            'model':'joint_electrode_portal_bridge_four_side_pad',
            'pad_layout_revision':PAD_LAYOUT_REVISION,
            'outer_exit_policy':policy.record(),
            'rules':asdict(rules),'settings':asdict(settings),'pad_settings':asdict(pad_settings),
            'retained_routes':len(chosen),'selected_pad_count':len(chosen),
            'pad_slot_count':len(frame['slots']),'pad_connection_verified':False,
            'pad_bank_occupancy':occupancy,
            'selected_navigation_views':dict(Counter(c['navigation_view'] for c in chosen)),
            'current_frame_pad_slot_upper_bound':len(frame['slots']),
            'current_frame_upper_bound_scope':'applies only if this square and its Pad pitch are fixed; expandable-frame capacity remains governed by electrode-center geometry',
            'route_completion_scope':'electrode_to_physical_pad_pending_gds_audit',
            'candidate_library':{'navigation_views':view_info,
                                 'complete_pad_columns_before_pruning':sum(
                                     reservoir.count for reservoir in by_pad.values()),
                                 'complete_pad_columns_after_pruning':len(columns),
                                 'geometric_outer_rejections':dict(rejected)},
             'center_preference':center_preference,
            'optimization':optimization,'routes':records,'_chosen':chosen,
            '_final_support':final,'_shell':frame['shell'],
            'frame':{k:v for k,v in frame.items() if k not in ('shell','slots')},
            'substrate_topology_before':betti(nav.support),
            'substrate_topology_after':betti(final),
            'continuous_upper_bound':upper,
            'center_distance_upper_bound':capacity_bound,
            'capacity_scope':'GDS-audited complete Pad routes give a constructive lower bound; exact-grid center packing and the Pad-exterior radial cut give continuous upper bounds independent of finite candidates and expandable Pad frame'}


def _add_polygon(cell,geom,layer,datatype=0):
    for polygon in _polygonal(geom):
        shell=gdstk.Polygon(np.asarray(polygon.exterior.coords),layer=layer,datatype=datatype)
        holes=[gdstk.Polygon(np.asarray(ring.coords)) for ring in polygon.interiors]
        pieces=gdstk.boolean(shell,holes,'not',precision=.001,layer=layer,datatype=datatype) if holes else [shell]
        for piece in pieces:cell.add(*piece.fracture(max_points=4000,precision=.001))


def write_pad_gds(input_path,output_path,routing,support_spec):
    with TemporaryDirectory(prefix='pad_gds_') as temp:
        source=Path(temp)/'input.gds';output=Path(temp)/'output.gds'
        source.write_bytes(Path(input_path).read_bytes())
        lib=gdstk.read_gds(str(source),unit=1e-6)
        tops=lib.top_level();names={cell.name for cell in lib.cells}
        name='FOUR_SIDE_PAD_ROUTING'
        while name in names:name+='X'
        cell=lib.new_cell(name)
        for top in tops:cell.add(gdstk.Reference(top))
        used={p.layer for top in tops for p in top.get_polygons()}
        def allocate(preferred):
            layer=preferred if preferred not in used else next(k for k in range(1020,65000) if k not in used)
            used.add(layer);return layer
        metal=allocate(20);electrode_marker=allocate(30);pad_marker=allocate(31);bridge_marker=allocate(32)
        island_marker=allocate(33);shell_marker=allocate(34)
        _add_polygon(cell,routing['_shell'],support_spec[0],support_spec[1])
        _add_polygon(cell,routing['_shell'],shell_marker,0)
        for i,route in enumerate(routing['_chosen'],1):
            _add_polygon(cell,route['bridge'],support_spec[0],support_spec[1])
            _add_polygon(cell,route['island'],support_spec[0],support_spec[1])
            _add_polygon(cell,route['bridge'],bridge_marker,i)
            _add_polygon(cell,route['island'],island_marker,i)
            _add_polygon(cell,route['metal'],metal,i)
            # Include these exact exported polygons in the conductor layer as
            # well.  Independent rounding of a union can otherwise leave a
            # 1 nm sliver of the electrode/Pad contact marker outside metal.
            _add_polygon(cell,route['electrode'],metal,i)
            _add_polygon(cell,route['pad'],metal,i)
            _add_polygon(cell,route['electrode'],electrode_marker,i)
            _add_polygon(cell,route['pad'],pad_marker,i)
        # Preserve the input cell polygons exactly on their native grid.  The
        # gdstk default (199 points) re-fractures long source polygons and can
        # remove tiny slivers even though the input cell is only referenced.
        lib.write_gds(str(output),max_points=4000)
        Path(output_path).write_bytes(output.read_bytes())
    # A separate rational-grid path witness lets the independent GDS reader
    # check that a full-width metal disk can travel from each electrode to its
    # Pad.  The witness carries the output hash and is never trusted without
    # geometric checks against the exported metal polygons.
    witness_path=Path(output_path).with_suffix('.width_witness.json')
    grid_um=lib.precision*1e6
    paths=[]
    for i,route in enumerate(routing['_chosen'],1):
        vertices=[*route['points_um'],route['declared_terminal_um'],
                  *route['outer_line'].coords]
        grid_points=[]
        for point in vertices:
            xy=[int(round(float(value)/grid_um)) for value in point]
            if not grid_points or xy!=grid_points[-1]:grid_points.append(xy)
        paths.append({'network':i,'points_grid_ticks':grid_points})
    witness_path.write_text(json.dumps({
        'format':'exact_grid_wire_centerline_v1',
        'output_sha256':sha256(Path(output_path).read_bytes()).hexdigest(),
        'native_grid_um':grid_um,
        'outer_exit_policy':routing.get('outer_exit_policy'),
        'networks':paths},ensure_ascii=False,separators=(',',':')),
        encoding='utf-8')
    return {'support_layer':list(support_spec),'metal_layer':metal,
            'electrode_marker_layer':electrode_marker,'pad_contact_layer':pad_marker,
            'bridge_marker_layer':bridge_marker,
            'island_marker_layer':island_marker,'shell_marker_layer':shell_marker,
            'outer_exit_policy':routing.get('outer_exit_policy'),
            'electrode_region':routing.get('electrode_region'),
            'wire_width_witness_path':str(witness_path.resolve())}


def audit_pad_gds(path,nav,routing,layers,*,progress=None):
    occupancy=_pad_occupancy(routing['_chosen'],PadSettings(**routing['pad_settings']))
    policy=policy_from_record(routing['outer_exit_policy'])
    with TemporaryDirectory(prefix='pad_audit_') as temp:
        copy=Path(temp)/'audit.gds';copy.write_bytes(Path(path).read_bytes())
        lib=gdstk.read_gds(str(copy),unit=1e-6)
    grouped=defaultdict(list)
    for top in lib.top_level():
        for p in top.get_polygons():grouped[(p.layer,p.datatype)].append(p)
    def layer_geom(spec):
        records=grouped[spec]
        if not records:return Polygon()
        merged=gdstk.boolean(records,[],'or',precision=.001)
        return unary_union([part for p in merged for part in _polygonal(shapely.make_valid(Polygon(p.points)))])
    # A single large shell with a circular void can trigger a topology error
    # in gdstk.boolean(records, [], 'or') even though its GDS polygons are
    # individually valid. Read the exported support as the vector union of
    # its actual polygons; the separate integer audit still checks containment.
    actual_support,_,_=read_support(Path(path),*layers['support_layer'])
    expected=routing['_final_support']
    support_diff=actual_support.symmetric_difference(expected).area
    generated_perimeter=(routing['_shell'].boundary.length+
                         sum(c['bridge'].boundary.length+c['island'].boundary.length for c in routing['_chosen']))
    tolerance=.004*generated_perimeter+.1
    support_ok=(support_diff<=tolerance and expected.buffer(.015).covers(actual_support) and
                actual_support.buffer(.015).covers(expected))
    support_provenance=None
    if not support_ok and support_diff<=tolerance and nav.summary.get('input_gds_path'):
        # Large coincident-boundary lattices can make a floating union/buffer
        # miss its own tiny slivers. Raw native-grid polygon identity proves
        # unchanged original support and exactly marked additions, without
        # relaxing the physical metal checks or independent integer audit.
        from exact_gds_audit import audit_support_polygon_provenance
        support_provenance=audit_support_polygon_provenance(
            nav.summary['input_gds_path'],path,layers,
            expected_precision_um=nav.summary.get('gds_native_precision_m',1e-9)*1e6)
        if support_provenance['input_sha256']!=nav.summary['sha256']:
            raise ValueError('Input GDS changed since navigation construction')
        support_ok=support_provenance['passed']
    # These source geometries are identical for every exported network.
    # Buffering a million-vertex lattice once per net made the audit dominate
    # runtime without adding any independent evidence.
    source_tolerance=nav.support.buffer(.003)
    shapely.prepare(source_tolerance);shapely.prepare(actual_support)
    actual_boundary=actual_support.boundary
    nets=[];metals=[];rules=routing['rules'];pad_ids=[];marker_centers=[]
    for i,c in enumerate(routing['_chosen'],1):
        if progress and (i==1 or i%5==0 or i==len(routing['_chosen'])):
            progress(f'GDS 回读检查：网络 {i}/{len(routing["_chosen"])}',.86+.025*i/len(routing['_chosen']))
        metal=layer_geom((layers['metal_layer'],i))
        electrode=layer_geom((layers['electrode_marker_layer'],i))
        marker_centers.append({'source_um':list(electrode.centroid.coords[0])})
        pad=layer_geom((layers['pad_contact_layer'],i))
        bridge=layer_geom((layers['bridge_marker_layer'],i))
        margin=metal.distance(actual_boundary) if actual_support.covers(metal) else -1
        electrode_ok=electrode.buffer(.003).covers(c['electrode'])
        pad_ok=(pad.buffer(.003).covers(c['pad']) and c['pad'].buffer(.003).covers(pad) and
                metal.buffer(.003).covers(pad) and metal.covers(Point(c['pad_target_um'])))
        bridge_ok=(bridge.buffer(.003).covers(c['bridge']) and
                   bridge.intersection(nav.support).area>1 and bridge.intersection(routing['_shell']).area>1)
        gate_ok,gate_reason=policy.check_bridge(nav.support,bridge,
                                               exit_point=c['portal_escape_um'])
        line_ok,line_reason=policy.check_exterior_line(nav.support,c['outer_line'])
        connected=betti(metal)=={'components':1,'holes':0} or betti(metal)['components']==1
        central_ok=source_tolerance.covers(c['central_wire'])
        geometry_error=metal.symmetric_difference(c['metal']).area/max(1,metal.area)
        record={'net':i,'pad_id':c['pad_id'],'connected':connected,
                'electrode_contact':electrode_ok and metal.covers(Point(c['source_um'])),
                'pad_contact':pad_ok,'bridge_join':bridge_ok,
                'outer_exit_window':gate_ok and line_ok,
                'outer_exit_failure_reason':gate_reason or line_reason,
                'central_wire_inside_original_support':central_ok,
                'support_containment':actual_support.covers(metal),
                'margin_um':margin,'metal_relative_difference':geometry_error}
        nets.append(record);metals.append(metal);pad_ids.append(c['pad_id'])
    gaps=[]
    for i,a in enumerate(metals):
        if progress and (i%5==0 or i==len(metals)-1):
            progress(f'GDS 回读检查：异网金属间距 {i+1}/{len(metals)}',
                     .885+.005*(i+1)/len(metals))
        gaps.extend(a.distance(b) for b in metals[:i])
    center_distances=[math.dist(a['source_um'],b['source_um']) for i,a in enumerate(routing['_chosen'])
                      for b in routing['_chosen'][:i]]
    min_center=min(center_distances,default=None)
    region_audit=region_for_navigation(nav,rules).audit(marker_centers)
    passed=(support_ok and region_audit['passed'] and len(set(pad_ids))==len(pad_ids) and
            all(c['connected'] and c['electrode_contact'] and c['pad_contact'] and c['bridge_join']
                and c['outer_exit_window']
                and c['central_wire_inside_original_support'] and c['support_containment']
                and c['margin_um']>=rules['margin_um']-.02 and c['metal_relative_difference']<.001
                for c in nets) and (not gaps or min(gaps)>=rules['spacing_um']-.02) and
            (min_center is None or min_center>=rules['minimum_center_spacing_um']-.003))
    if not passed:
        raise RuntimeError(f'Pad GDS roundtrip failed: support={support_ok}, diff={support_diff}, '
                           f'gap={min(gaps,default=None)}, center={min_center}, nets={nets[:8]}')
    routing['pad_connection_verified']=True
    routing['route_completion_scope']='electrode_to_dedicated_physical_pad'
    routing['status']='checked_complete_pad_lower_bound'
    return {'passed':True,'pad_connection_verified':True,'exported_nets':len(nets),
            'electrode_region_audit':region_audit,
            'outer_exit_policy_verified':True,'outer_exit_policy':policy.record(),
            'exact_support_provenance_when_float_overlay_unreliable':support_provenance,
            'connected_pad_count':len(nets),'audit_scope':'electrode_to_dedicated_physical_pad',
            'pad_bank_occupancy':occupancy,
            'output_layers':layers,'support_symmetric_difference_um2':support_diff,
            'support_area_tolerance_um2':tolerance,
            'source_support_topology':betti(nav.support),
            'output_support_topology':betti(actual_support),
            'minimum_metal_to_support_boundary_um':min((c['margin_um'] for c in nets),default=None),
            'minimum_inter_net_gap_um':min(gaps,default=None),
            'minimum_electrode_center_distance_um':min_center,
            'pad_ids':pad_ids,'net_checks':nets}
