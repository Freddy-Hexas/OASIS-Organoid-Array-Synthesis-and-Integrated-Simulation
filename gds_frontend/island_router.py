"""GDS-only circular support attachments and geometry-certified route columns.

The same vector erosion, constrained triangulation, attachment test and MILP
apply to every input. Filenames, generators and A/B/C labels are not consulted.
The continuous problem is approximated by an explicit finite candidate library.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, asdict
from hashlib import sha256
from heapq import heappop, heappush
from itertools import combinations
from pathlib import Path
from tempfile import TemporaryDirectory
import math
import time
import gzip
import json

import gdstk
import networkx as nx
import numpy as np
import shapely
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.spatial import cKDTree
from scipy.sparse import coo_array
from shapely.geometry import LineString, Point, Polygon
from shapely.ops import unary_union, substring
from shapely.strtree import STRtree

from frontend import read_support, _polygonal
from process_geometry import ProcessRules
from curved_centerline import CurveSettings, UncertifiableCurve, smooth_centerline
from capacity_bounds import (gds_cell_cover_upper_bound,
                             gds_enclosing_disk_angular_upper_bound,
                             verify_gds_angular_certificate)
from small_cell_conflicts import (gds_small_cell_conflict_upper_bound,
                                  verify_small_cell_conflict_certificate)
from outer_exit_policy import outside_circle_edge_intervals
from navigation_simplification import simplify_inside
from electrode_region import region_for, region_for_navigation, add_region_candidates
from solver_time_policy import deadline_after, deadline_expired, solver_options


@dataclass(frozen=True)
class IslandSettings:
    candidate_step_um: float = 600.0
    routes_per_candidate: int = 1
    lane_half_steps: int = 1
    # Search resolution only; never used in a capacity certificate.
    external_portal_pitch_um: float | None = None
    # The full boundary mode samples concave exterior arcs as well as the
    # radial envelope.  Both modes only construct lower-bound candidates.
    external_portal_mode: str = 'radial_envelope'
    # None derives the navigation tolerance from the existing numerical guard.
    # An explicit value is only an upper bound; topology is still certified.
    navigation_inset_um: float | None = None
    # None preserves the guarded finite-column search. Zero retains every
    # center-domain component for searches that permit a short island-to-exit
    # connection; the terminal guard is not a fabrication rule.
    navigation_component_exclusion_distance_um: float | None = None
    navigation_prune_leaf_length_um: float | None = None
    # Search restriction only. None keeps the historical 23.05 um guard;
    # zero permits an island-supported electrode immediately next to an exit.
    anchor_terminal_search_gap_um: float | None = None
    include_guarded_outlet_alternative: bool = False
    numeric_guard_um: float = .05
    pad_gap_um: float = 4.0
    milp_time_limit_s: float | None = None


@dataclass
class VectorNavigation:
    support: object
    center_domain: object
    navigation_domain: object
    collector: object
    graph: nx.Graph
    summary: dict
    candidates: list
    outlets: list


def betti(geom):
    pieces=_polygonal(geom)
    return {'components':len(pieces), 'holes':sum(len(p.interiors) for p in pieces)}


class NavigationGeometryError(ValueError):
    def __init__(self,message,diagnostics):
        super().__init__(message)
        self.diagnostics=diagnostics


def _navigation_domain(center,rules,settings,grid_um):
    """Apply an optional terminal-distance search guard without losing topology.

    A wire centerline cannot move between positively separated components of
    the exact eroded support.  The bounding-box diagonal bounds the distance
    between any two points of a component.  The default finite-column search
    excludes a component only when its diagonal is smaller than that search's
    anchor-terminal distance.  Zero disables the exclusion for constructions
    that permit short island-supported connections to the exterior.
    """
    if (settings.numeric_guard_um<0 or settings.pad_gap_um<0 or
            (settings.navigation_inset_um is not None and settings.navigation_inset_um<0) or
            (settings.navigation_component_exclusion_distance_um is not None and
             settings.navigation_component_exclusion_distance_um<0) or
            (settings.navigation_prune_leaf_length_um is not None and
             settings.navigation_prune_leaf_length_um<0) or
            (settings.anchor_terminal_search_gap_um is not None and
             settings.anchor_terminal_search_gap_um<0)):
        raise ValueError('Navigation settings and numerical clearances must be nonnegative')
    pieces=sorted(_polygonal(center),key=lambda p:(*p.bounds,p.area))
    tree=STRtree(pieces)
    guarded_separation=(rules.electrode_diameter_um/2+rules.margin_um+
                        settings.numeric_guard_um+settings.pad_gap_um)
    separation=(guarded_separation if
                settings.navigation_component_exclusion_distance_um is None else
                settings.navigation_component_exclusion_distance_um)
    target=settings.numeric_guard_um if settings.navigation_inset_um is None else settings.navigation_inset_um
    scale=max(1.0,*[abs(v) for v in center.bounds])
    precision=max(float(grid_um),8*math.ulp(scale))
    kept=[];excluded=[]
    for index,piece in enumerate(pieces):
        x0,y0,x1,y1=piece.bounds
        diameter_bound=math.hypot(x1-x0,y1-y0)
        isolated=not any(j!=index and piece.intersects(pieces[j]) for j in tree.query(piece))
        if isolated and diameter_bound+2*precision<separation:
            sample=piece.representative_point()
            excluded.append({'component_index':index,'area_um2':float(piece.area),
                             'bounds_um':[float(v) for v in piece.bounds],
                             'representative_point_um':[float(sample.x),float(sample.y)],
                             'component_wkb_sha256':sha256(piece.wkb).hexdigest(),
                             'bbox_diameter_upper_um':diameter_bound,
                             'precision_allowance_um':2*precision,
                             'required_center_to_terminal_distance_um':separation,
                             'isolated_from_other_center_components':True,
                             'reason':'no_anchor_terminal_pair_can_meet_clearance',
                             'proof':'Every centerline stays in this isolated exact wire-center component; every pair of its points is at most the bounding-box diagonal apart, below the required electrode-center/terminal separation.'})
        else:
            kept.append((index,piece))
    diagnostics={'method':'componentwise_topology_certified_navigation',
                 'wire_center_domain_topology':betti(center),
                 'gds_precision_um':float(grid_um),
                 'required_center_to_terminal_distance_um':separation,
                 'target_inset_upper_bound_um':target,
                 'target_source':('existing_numeric_guard' if settings.navigation_inset_um is None else 'explicit_upper_bound'),
                 'excluded_component_count':len(excluded),
                 'excluded_area_um2':sum(p['area_um2'] for p in excluded),
                 'excluded_components':excluded,
                 'exclusion_scope':(
                     'no components excluded: short island-to-exit connections retained'
                     if separation==0 else
                     'guarded finite-column electrode-to-terminal search only; '
                     'not a physical impossibility claim for shorter connections')}
    if not kept:
        raise NavigationGeometryError('No wire-center component can meet electrode-to-terminal clearance',diagnostics)

    def inset_piece(piece,inset):
        if inset<=0:return piece
        return piece.buffer(-inset,quad_segs=16).simplify(inset/2,preserve_topology=True).intersection(piece)

    navigation_pieces=[];inset_records=[]
    for index,piece in kept:
        expected={'components':1,'holes':len(piece.interiors)}
        selected=0.0;nav_piece=piece
        if target>0:
            trial=inset_piece(piece,target)
            if trial.is_valid and betti(trial)==expected:
                selected=target;nav_piece=trial
            else:
                # The GDS grid defines the search resolution.  No fixed retry
                # count or structure label determines the local inset.
                low=0.0;high=target
                while high-low>precision:
                    mid=(low+high)/2
                    trial=inset_piece(piece,mid)
                    if trial.is_valid and betti(trial)==expected:
                        low=mid;selected=mid;nav_piece=trial
                    else:
                        high=mid
        navigation_pieces.append(nav_piece)
        inset_records.append({'component_index':index,'area_um2':float(piece.area),
                              'holes':len(piece.interiors),'selected_inset_um':selected,
                              'navigation_topology':betti(nav_piece)})
    retained=unary_union([p for _,p in kept])
    navigation=unary_union(navigation_pieces)
    expected={'components':len(kept),'holes':sum(len(p.interiors) for _,p in kept)}
    diagnostics['retained_center_topology']=betti(retained)
    diagnostics['navigation_topology']=betti(navigation)
    diagnostics['retained_center_area_um2']=float(retained.area)
    diagnostics['navigation_area_um2']=float(navigation.area)
    diagnostics['component_insets']=inset_records
    diagnostics['zero_inset_component_count']=sum(p['selected_inset_um']==0 for p in inset_records)
    if diagnostics['retained_center_topology']!=expected or diagnostics['navigation_topology']!=expected or not retained.covers(navigation):
        raise NavigationGeometryError('Adaptive navigation failed to preserve retained wire-domain topology',diagnostics)
    navigation,coarse_certificate=simplify_inside(navigation,rules.wire_width_um/4,grid_um)
    diagnostics['inner_simplification']=coarse_certificate
    diagnostics['navigation_area_um2']=float(navigation.area)
    return navigation,diagnostics


def disk(center,radius,quad=64):
    # Circumscribe the ideal disk with 2 nm extra radial headroom. Exporting
    # vertices to the native 1 nm GDS grid can otherwise shave roughly 1 nm
    # from the actual polygon inradius even though the floating polygon passed.
    native_grid_guard_um=.002
    return Point(center).buffer((radius+native_grid_guard_um)/
                                math.cos(math.pi/(4*quad)),quad_segs=quad)


def _shorten(points, domain, tolerance, *, visibility_first=True):
    """Construct a certified shortcut with either of two uniform path policies.

    ``visibility_first`` removes every directly visible graph detour.  The
    guide-preserving policy retains offsets larger than ``tolerance``; in a
    tight curve it can leave more room for a certified smooth bend.  Both are
    tested for every structure, then the same MILP chooses among them.
    """
    pts=np.asarray(points,dtype=float)
    keep=np.r_[True,np.linalg.norm(np.diff(pts,axis=0),axis=1)>1e-10]
    pts=pts[keep]
    if len(pts)<2:raise ValueError('Degenerate path')
    shapely.prepare(domain)
    selected={0,len(pts)-1};stack=[(0,len(pts)-1)]
    while stack:
        a,b=stack.pop()
        if b-a<=1:continue
        # The graph dual is only a topology guide.  Its triangle-center
        # detours should not be kept merely because they are farther than a
        # fixed tolerance from a visible direct chord.
        if visibility_first and domain.covers(LineString([pts[a],pts[b]])):
            continue
        delta=pts[b]-pts[a];norm=float(np.dot(delta,delta))
        inner=pts[a+1:b]
        t=np.clip(((inner-pts[a])@delta)/norm,0,1) if norm>1e-18 else np.zeros(len(inner))
        distances=np.linalg.norm(inner-(pts[a]+t[:,None]*delta),axis=1)
        split=a+1+int(np.argmax(distances))
        if float(distances.max())<=tolerance:
            if not visibility_first and domain.covers(LineString([pts[a],pts[b]])):
                continue
            split=(a+b)//2
        selected.add(split);stack.extend([(a,split),(split,b)])
    result=LineString(pts[sorted(selected)])
    if not domain.covers(result):
        outside=result.difference(domain)
        coordinates=shapely.get_coordinates(outside)
        numerical_boundary_only=(outside.length<=1e-7 and
            (not len(coordinates) or bool(np.all(shapely.distance(shapely.points(coordinates),domain)<=1e-7))))
        if not numerical_boundary_only:raise RuntimeError('Certified path shortcut left its domain')
    return result


def _orient_points(points, source):
    pts=np.asarray(points)
    xy=np.asarray(source)
    return pts[::-1] if np.linalg.norm(pts[-1]-xy)<np.linalg.norm(pts[0]-xy) else pts


def _prune_short_leaves(graph,protected,limit):
    graph=graph.copy()
    changed=True
    while changed:
        changed=False
        for leaf in list(graph.nodes):
            if leaf not in graph or leaf in protected or graph.degree(leaf)!=1:
                continue
            chain=[leaf];previous=None;current=leaf;length=0.0
            while True:
                nxt=[n for n in graph[current] if n!=previous]
                if len(nxt)!=1:break
                node=nxt[0];length+=graph[current][node]['weight']
                previous,current=current,node
                if current in protected or graph.degree(current)!=2:break
                chain.append(current)
            if length<limit and current not in chain:
                graph.remove_nodes_from(chain);changed=True
    return graph


def _compress_dual(graph,xy,protected,domain,tolerance):
    """Retain edge geometry, loops and ports while collapsing degree-two chains."""
    shapely.prepare(domain)
    terminals={n for n in graph if graph.degree(n)!=2 or n in protected}
    for comp in nx.connected_components(graph):
        if not terminals.intersection(comp):terminals.add(min(comp))
    result=nx.MultiGraph()
    for node in terminals:
        result.add_node(node,xy_um=xy[node].tolist(),kind='outlet' if node in protected else
                        'leaf' if graph.degree(node)==1 else 'junction')
    visited=set()
    for start in sorted(terminals):
        for nxt in sorted(graph[start]):
            edge=tuple(sorted((start,nxt)))
            if edge in visited:continue
            points=[xy[start]];previous=start;current=nxt
            while True:
                visited.add(tuple(sorted((previous,current))))
                points.append(graph[previous][current]['portal'])
                points.append(xy[current])
                if current in terminals:break
                following=next(n for n in graph[current] if n!=previous)
                previous,current=current,following
            line=_shorten(points,domain,tolerance)
            result.add_edge(start,current,points_um=np.asarray(line.coords),length_um=line.length)
    # Subdivide every parallel edge and self-loop into an explicit unique edge
    # node so shortest paths cannot silently discard a topological alternative.
    simple=nx.Graph()
    simple.add_nodes_from(result.nodes(data=True))
    next_id=max(graph.nodes,default=-1)+1
    for edge_id,(u,v,k,data) in enumerate(result.edges(keys=True,data=True)):
        line=LineString(_orient_points(data['points_um'],result.nodes[u]['xy_um']))
        a=substring(line,0,line.length/2);b=substring(line,line.length/2,line.length)
        mid=next_id;next_id+=1
        simple.add_node(mid,xy_um=list(a.coords[-1]),kind='corridor')
        if u==v:
            # A self-loop requires two distinct intermediate nodes in a simple graph.
            mid2=next_id;next_id+=1
            third=line.length/3
            pieces=[substring(line,0,third),substring(line,third,2*third),substring(line,2*third,line.length)]
            simple.nodes[mid]['xy_um']=list(pieces[0].coords[-1])
            simple.add_node(mid2,xy_um=list(pieces[1].coords[-1]),kind='corridor')
            for x,y,p in zip([u,mid,mid2],[mid,mid2,u],pieces):
                simple.add_edge(x,y,points_um=np.asarray(p.coords),weight=p.length,corridor_id=edge_id)
        else:
            for x,y,p in ((u,mid,a),(mid,v,b)):
                simple.add_edge(x,y,points_um=np.asarray(p.coords),weight=p.length,corridor_id=edge_id)
    return simple


def _external_boundary_samples(domain,ownership,global_origin,pitch,
                               mode='radial_envelope',exit_policy=None):
    """Sample each component's exterior without assuming a structure class.

    The terminal is a real edge of the certified wire-center domain, never a
    width-class boundary or a boundary of a protected hole.  The historical
    radial mode is small; full_boundary also samples inward-facing exterior
    arcs.  Neither finite proposal set is a continuous upper bound.
    """
    if mode not in ('radial_envelope','full_boundary') or pitch<=0:
        raise ValueError('Unknown exterior sampling mode or nonpositive pitch')
    if exit_policy is not None:
        # Sample every connected finite outer cap, including caps much shorter
        # than pitch.  Angular bins and the component's own centroid must not
        # manufacture an exit in an internal concavity or a small inner island.
        samples=[];arcs=0;matched=0;unmatched=0;rejected=0
        for part in _polygonal(domain):
            coords=np.asarray(part.exterior.coords)
            groups=[]
            for a,b in zip(coords[:-1],coords[1:]):
                owner=ownership.get(tuple(sorted((tuple(a),tuple(b)))))
                if owner is None:
                    unmatched+=1;continue
                matched+=1
                intervals=exit_policy.edge_intervals(a,b,search=True)
                if not intervals:rejected+=1
                for aa,bb in intervals:
                    length=float(np.linalg.norm(bb-aa))
                    if length<1e-9:continue
                    item=(owner,aa,bb,length)
                    if groups and np.linalg.norm(groups[-1][-1][2]-aa)<1e-8:
                        groups[-1].append(item)
                    else:groups.append([item])
            if len(groups)>1 and np.linalg.norm(groups[-1][-1][2]-groups[0][0][1])<1e-8:
                groups[0]=groups[-1]+groups[0];groups.pop()
            for group in groups:
                arcs+=1
                total=sum(edge[3] for edge in group)
                count=max(1,math.ceil(total/pitch));index=0;before=0.0
                for j in range(count):
                    target=total*(j+.5)/count
                    while index+1<len(group) and before+group[index][3]<target:
                        before+=group[index][3];index+=1
                    owner,a,b,length=group[index]
                    fraction=min(1.0,max(0.0,(target-before)/length))
                    samples.append((owner,a+fraction*(b-a)))
        return samples,{'sampling_mode':'process_sized_outer_extremity_caps',
                        'sampling_pitch_um':pitch,'outer_cap_count':arcs,
                        'sampled_exterior_portals':len(samples),
                        'exterior_mesh_edges_matched':matched,
                        'exterior_mesh_edges_unmatched':unmatched,
                        'inner_boundary_edges_excluded':rejected,
                        'sampling_scope':'all connected directional outer-envelope caps; candidate superset followed by mandatory launch/whole-bridge checks'}
    samples=[];rings=0;matched=0;unmatched=0
    for part in _polygonal(domain):
        rings+=1
        local_origin=np.asarray(part.centroid.coords[0],dtype=float)
        bins=max(12,math.ceil(part.convex_hull.length/max(pitch,1.0)))
        best={}
        coords=np.asarray(part.exterior.coords)
        matching_edges=[]
        for a,b in zip(coords[:-1],coords[1:]):
            key=tuple(sorted((tuple(a),tuple(b))))
            owner=ownership.get(key)
            if owner is None:
                unmatched+=1;continue
            matched+=1
            if mode=='full_boundary':
                length=float(np.linalg.norm(b-a))
                if length>0:
                    matching_edges.append((owner,a,b,length))
                continue
            point=(a+b)/2
            angle=math.atan2(*(point-local_origin)[::-1])
            sector=min(bins-1,int((angle+math.pi)/(2*math.pi)*bins))
            radius=float(np.linalg.norm(point-global_origin))
            if sector not in best or radius>best[sector][0]:
                best[sector]=(radius,owner,point)
        if mode=='full_boundary' and matching_edges:
            total=sum(edge[3] for edge in matching_edges)
            count=max(1,math.ceil(total/pitch))
            index=0
            before=0.0
            for sample_index in range(count):
                target=total*(sample_index+.5)/count
                while (index+1<len(matching_edges) and
                       before+matching_edges[index][3]<target):
                    before+=matching_edges[index][3]
                    index+=1
                owner,a,b,length=matching_edges[index]
                fraction=min(1.0,max(0.0,(target-before)/length))
                samples.append((owner,a+fraction*(b-a)))
        elif mode=='radial_envelope' and best:
            samples.extend((owner,point) for _,owner,point in best.values())
        else:
            raise RuntimeError('Exterior boundary could not be matched to constrained mesh edges')
    return samples,{'exterior_component_rings':rings,
                    'exterior_mesh_edges_matched':matched,
                    'exterior_mesh_edges_unmatched':unmatched,
                    'sampled_exterior_portals':len(samples),
                    'sampling_pitch_um':pitch,
                    'sampling_scope':('full wire-center exterior arclength'
                                      if mode=='full_boundary' else
                                      'finite angular envelope of every wire-center exterior ring'),
                    'sampling_mode':mode}


def build_navigation(path,rules:ProcessRules,*,layer=10,datatype=0,
                     outlet_mode='collector',settings=IslandSettings(),progress=None,
                     diagnostics_callback=None,outer_exit_policy=None):
    notify=progress or (lambda *a,**k:None)
    rules.validate()
    if outlet_mode not in ('collector','open_tips','auto_geometry','external_boundary'):
        raise ValueError('Unknown outlet mode')
    if settings.external_portal_pitch_um is not None and settings.external_portal_pitch_um<=0:
        raise ValueError('External portal sampling pitch must be positive')
    if settings.external_portal_mode not in ('radial_envelope','full_boundary'):
        raise ValueError('Unknown exterior sampling mode')
    start=time.perf_counter()
    support,_,meta=read_support(Path(path),layer,datatype)
    radius=rules.wire_width_um/2+rules.margin_um+settings.numeric_guard_um
    center=support.buffer(-radius,quad_segs=32)
    if center.is_empty:raise ValueError('No manufacturable wire center domain')
    try:
        nav,navigation_certificate=_navigation_domain(center,rules,settings,
            meta['gds_native_precision_m']*1e6)
    except NavigationGeometryError as exc:
        exc.diagnostics['input_sha256']=meta['sha256']
        exc.diagnostics['support_layer']=[layer,datatype]
        raise
    navigation_certificate['input_sha256']=meta['sha256']
    navigation_certificate['support_layer']=[layer,datatype]
    if diagnostics_callback is not None:
        diagnostics_callback(navigation_certificate)
    collector=Polygon()
    if outlet_mode in ('collector','auto_geometry'):
        core=support.buffer(-rules.collector_clearance_um,quad_segs=32)
        minimum_area=4*math.pi*rules.collector_clearance_um**2
        broad=[p for p in _polygonal(core) if p.area>=minimum_area]
        collector=unary_union([p.buffer(rules.collector_clearance_um,quad_segs=32) for p in broad]).intersection(support)
    routing_domain=nav.difference(collector)
    notify('约束三角剖分与拓扑认证',.18)
    triangles=list(shapely.constrained_delaunay_triangles(routing_domain).geoms)
    if not triangles:raise ValueError('No triangulatable routing domain')
    mesh_union=shapely.coverage_union_all(triangles)
    mesh_error=mesh_union.symmetric_difference(routing_domain).area
    if mesh_error/max(1,routing_domain.area)>1e-8:
        raise RuntimeError('Constrained-triangulation coverage differs from vector routing domain')
    xy=np.asarray([np.asarray(t.exterior.coords)[:3].mean(axis=0) for t in triangles])
    dual=nx.Graph();dual.add_nodes_from(range(len(triangles)))
    ownership={};boundary=[]
    for i,tri in enumerate(triangles):
        vertices=np.asarray(tri.exterior.coords)[:3]
        for j in range(3):
            a=tuple(vertices[j]);b=tuple(vertices[(j+1)%3]);key=tuple(sorted((a,b)))
            if key in ownership:
                other=ownership.pop(key);portal=(np.asarray(a)+np.asarray(b))/2
                dual.add_edge(i,other,portal=portal,
                              weight=float(np.linalg.norm(xy[i]-portal)+np.linalg.norm(xy[other]-portal)))
            else:ownership[key]=i
    topo=betti(routing_domain)
    dual_beta={'components':nx.number_connected_components(dual),
               'holes':dual.number_of_edges()-dual.number_of_nodes()+nx.number_connected_components(dual)}
    if dual_beta!=topo:raise RuntimeError(f'Triangulation dual topology mismatch: {dual_beta} != {topo}')
    protected=set();external_outlets=set();external_sampling=None
    if outlet_mode=='external_boundary':
        origin=(outer_exit_policy.origin if outer_exit_policy is not None else
                np.asarray([(support.bounds[0]+support.bounds[2])/2,
                            (support.bounds[1]+support.bounds[3])/2]))
        portal_pitch=(settings.external_portal_pitch_um if
                      settings.external_portal_pitch_um is not None else
                      max(100.0,settings.candidate_step_um/2))
        if outer_exit_policy is not None:
            portal_pitch=min(portal_pitch,rules.wire_width_um+rules.spacing_um)
        samples,external_sampling=_external_boundary_samples(
            routing_domain,ownership,origin,portal_pitch,
            settings.external_portal_mode,outer_exit_policy)
        for owner,point in samples:
            terminal=len(xy)+len(boundary)
            boundary.append({'triangle':owner,'midpoint':point.tolist(),
                             'kind':'external_boundary'})
            dual.add_node(terminal)
            dual.add_edge(owner,terminal,portal=point,
                          weight=float(np.linalg.norm(xy[owner]-point)))
            protected.add(terminal);external_outlets.add(terminal)
        if boundary:
            xy=np.vstack([xy,np.asarray([p['midpoint'] for p in boundary])])
    if outlet_mode in ('collector','auto_geometry') and not collector.is_empty:
        for (a,b),owner in ownership.items():
            midpoint=(np.asarray(a)+np.asarray(b))/2
            if collector.distance(Point(midpoint))<1e-6:
                # The terminal is the actual boundary portal. A triangle
                # centroid can be hundreds of micrometers away in a long
                # triangle and must never be reported as an outlet.
                terminal=len(xy)+len(boundary)
                boundary.append({'triangle':owner,'a':list(a),'b':list(b),'midpoint':midpoint.tolist()})
                dual.add_node(terminal)
                dual.add_edge(owner,terminal,portal=midpoint,
                              weight=float(np.linalg.norm(xy[owner]-midpoint)))
                protected.add(terminal)
        if boundary:
            xy=np.vstack([xy,np.asarray([p['midpoint'] for p in boundary])])
    notify('压缩导航图与生成圆岛候选',.3)
    leaf_limit=(rules.electrode_diameter_um+2*rules.margin_um if
                settings.navigation_prune_leaf_length_um is None else
                settings.navigation_prune_leaf_length_um)
    pruned=_prune_short_leaves(dual,protected,leaf_limit)
    collector_outlets=protected.intersection(pruned.nodes)-external_outlets
    collector_terminal_count=len(collector_outlets)
    if outlet_mode in ('open_tips','auto_geometry'):
        protected.update(n for n in pruned if pruned.degree(n)==1)
    open_tip_terminal_count=len(protected-collector_outlets-external_outlets)
    selected_insets=[entry['selected_inset_um'] for entry in navigation_certificate['component_insets']]
    graph=_compress_dual(pruned,xy,protected,routing_domain,max(.1,max(selected_insets,default=0)))
    outlets=[n for n in sorted(protected) if n in graph]
    for outlet in outlets:
        graph.nodes[outlet]['terminal_kind']=(
            'external_boundary' if outlet in external_outlets else
            'collector_interface' if outlet in collector_outlets else 'open_tip')
    candidates=[]
    # Uniform dyadic subdivision provides nested candidate sets when step halves.
    next_id=max(graph.nodes,default=-1)+1
    corridor_endpoints={n for n,d in graph.nodes(data=True) if d['kind'] in ('junction','leaf') and n not in protected}
    for u,v,data in list(graph.edges(data=True)):
        line=LineString(_orient_points(data['points_um'],graph.nodes[u]['xy_um']));length=line.length
        if length<1e-7:continue
        divisions=max(1,2**int(math.ceil(math.log2(max(1,length/settings.candidate_step_um)))))
        distances=[length*j/divisions for j in range(1,divisions)]
        # Include a midpoint even for a short edge.
        if not distances:distances=[length/2]
        chain=[u];cuts=[0.0,*distances,length];graph.remove_edge(u,v)
        for distance in distances:
            c=list(line.interpolate(distance).coords[0]);node=next_id;next_id+=1
            graph.add_node(node,xy_um=c,kind='candidate')
            candidates.append(node);chain.append(node)
        chain.append(v)
        for a,b,x,y in zip(chain[:-1],chain[1:],cuts[:-1],cuts[1:]):
            piece=substring(line,x,y)
            graph.add_edge(a,b,points_um=np.asarray(piece.coords),weight=piece.length,
                           corridor_id=data['corridor_id'])
    candidates+=sorted(corridor_endpoints)
    electrode_region=region_for(support,rules,meta['gds_native_precision_m']*1e6,
        origin=outer_exit_policy.origin if outer_exit_policy is not None else None)
    unrestricted_candidate_count=len(candidates)
    candidates,region_candidates_added=add_region_candidates(graph,candidates,electrode_region)
    # Cyclic orders record the embedding; feasibility is checked on actual paths.
    for n,d in graph.nodes(data=True):
        ports=[]
        origin=np.asarray(d['xy_um'])
        for other,edge in graph[n].items():
            p=edge['points_um']
            if np.linalg.norm(p[-1]-origin)<np.linalg.norm(p[0]-origin):p=p[::-1]
            tangent=p[min(2,len(p)-1)]-origin
            ports.append({'other':other,'angle_rad':float(math.atan2(tangent[1],tangent[0]))})
        d['ports_ccw']=sorted(ports,key=lambda p:p['angle_rad'])
    summary={**meta,'input_gds_path':str(Path(path).resolve()),
             'extractor_uses_generator_source':False,
             'method':'vector_erosion_constrained_triangulation_dual',
             'wire_center_domain_topology':betti(center),'navigation_domain_topology':betti(nav),
             'navigation_geometry_certificate':navigation_certificate,
             'routing_domain_topology':topo,'triangle_dual_topology':dual_beta,
             'triangle_count':len(triangles),'graph_nodes':graph.number_of_nodes(),
             'triangle_mesh_symmetric_difference_um2':float(mesh_error),
             'graph_edges':graph.number_of_edges(),'candidate_anchors':len(candidates),
             'electrode_region':electrode_region.record(),
             'electrode_region_candidate_filter':{'unrestricted_candidates':unrestricted_candidate_count,
                 'inside_candidates':len(candidates),'disk_intersection_candidates_added':region_candidates_added,
                 'wire_graph_clipped':False},
             'collector_interface_windows':collector_terminal_count,
             'external_boundary_terminals':len(external_outlets),
             'external_boundary_sampling':external_sampling,
             'outer_exit_policy':(outer_exit_policy.record() if outer_exit_policy is not None else None),
             'open_tip_terminals':open_tip_terminal_count,
             'geometry_candidate_terminals':len(outlets),'outlet_mode':outlet_mode,
             'selected_pitch_um':None,'full_graph_matches_vector_topology':True,
             'topology_certificate_scope':'vector routing domain and constrained-triangulation dual',
             'navigation_seconds':time.perf_counter()-start,
             'collector_clearance_threshold_um':(
                 rules.collector_clearance_um if outlet_mode in ('collector','auto_geometry') else None),
             'collector_minimum_core_area_um2':(
                 4*math.pi*rules.collector_clearance_um**2
                 if outlet_mode in ('collector','auto_geometry') else None),
             'outlet_semantics':'geometry-derived candidate terminals, not assigned electrical pads'}
    return VectorNavigation(support,center,routing_domain,collector,graph,summary,candidates,outlets)


def _path_line(graph,nodes):
    points=[graph.nodes[nodes[0]]['xy_um']];corridors=[]
    for a,b in zip(nodes[:-1],nodes[1:]):
        edge=graph[a][b];p=edge['points_um'];origin=np.asarray(points[-1])
        if np.linalg.norm(p[-1]-origin)<np.linalg.norm(p[0]-origin):p=p[::-1]
        points.extend(p[1:].tolist());corridors.append(edge['corridor_id'])
    return LineString(points),sorted(set(corridors))


def _attachment(support,c,radius,minimum_gap):
    island=disk(c,radius)
    overlap=island.intersection(support)
    if betti(overlap)!={'components':1,'holes':0} or overlap.area<1:
        return None,'noncontractible_or_disconnected_attachment'
    # Check distinct nearby substrate branches for a sub-limit remaining gap.
    neighborhood=disk(c,radius+minimum_gap+.1).intersection(support)
    nearby=[p for p in _polygonal(neighborhood) if not p.intersects(overlap)]
    if any(island.distance(p)<minimum_gap-1e-6 for p in nearby):
        return None,'substrate_gap'
    return island,None


def _route_columns(nav,rules,settings,notify):
    g=nav.graph;support=nav.support;shapely.prepare(support)
    electrode_region=region_for_navigation(nav,rules)
    shapely.prepare(nav.center_domain)
    curve_settings=CurveSettings()
    if not nav.outlets:return [],{'reason':'no_declared_geometry_terminal'}
    terminal=('super','terminal')
    search=g.copy()
    for outlet in nav.outlets:search.add_edge(terminal,outlet,weight=0)
    distances,all_paths=nx.single_source_dijkstra(search,terminal,weight='weight')
    wire_radius=rules.wire_width_um/2
    pad_radius=rules.electrode_diameter_um/2+rules.margin_um+settings.numeric_guard_um
    terminal_separation=(pad_radius+settings.pad_gap_um if
                         settings.anchor_terminal_search_gap_um is None else
                         settings.anchor_terminal_search_gap_um)
    outlet_set=set(nav.outlets)
    def eligible_outlet_paths(source,origin,limit,min_separation):
        """Shortest valid portals, without treating a too-near portal as fatal.

        A single Dijkstra traversal is cheaper than re-solving once for each
        blocked terminal when the external boundary is sampled densely.
        """
        queue=[(0.0,source)]
        distance={source:0.0}
        parent={}
        found=[]
        while queue and len(found)<limit:
            length,node=heappop(queue)
            if length>distance[node]+1e-9:
                continue
            if (node in outlet_set and length>=min_separation and
                    np.linalg.norm(np.asarray(origin)-
                                   np.asarray(g.nodes[node]['xy_um']))>=min_separation):
                path=[node]
                while path[-1]!=source:
                    path.append(parent[path[-1]])
                path.reverse()
                found.append((length,node,path))
            for neighbor,data in g[node].items():
                candidate=length+data['weight']
                if candidate<distance.get(neighbor,math.inf)-1e-9:
                    distance[neighbor]=candidate
                    parent[neighbor]=node
                    heappush(queue,(candidate,neighbor))
        return found
    offsets=[0.0,*[sign*i*(rules.wire_width_um+rules.spacing_um)/2
                   for i in range(1,settings.lane_half_steps+1) for sign in (-1,1)]]
    columns=[];reject=Counter();valid_sources=set()
    for index,source in enumerate(nav.candidates):
        c=g.nodes[source]['xy_um']
        if not electrode_region.contains(c,export_safe=True):
            reject['outside_electrode_region']+=1;continue
        if source not in distances:reject['unreachable_anchor']+=1;continue
        island,reason=_attachment(support,c,pad_radius,settings.pad_gap_um)
        if island is None:reject[reason]+=1;continue
        electrode=disk(c,rules.electrode_diameter_um/2)
        local_support=support.union(island);shapely.prepare(local_support)
        base_path=all_paths[source][1:]
        if not base_path:reject['terminal_guard']+=1;continue
        closest=base_path[0]
        closest_valid=(distances[source]>=terminal_separation and
                       np.linalg.norm(np.asarray(c)-
                                      np.asarray(g.nodes[closest]['xy_um']))>=terminal_separation)
        choices=([(distances[source],closest,base_path[::-1])]
                 if closest_valid and settings.routes_per_candidate==1 else
                 eligible_outlet_paths(source,c,settings.routes_per_candidate,
                                       terminal_separation))
        if settings.include_guarded_outlet_alternative:
            guarded_gap=pad_radius+settings.pad_gap_um
            used={outlet for _,outlet,_ in choices}
            choices.extend((length,outlet,path) for length,outlet,path in
                           eligible_outlet_paths(source,c,
                                                 settings.routes_per_candidate,
                                                 guarded_gap)
                           if outlet not in used)
        if not choices:reject['terminal_guard']+=1;continue
        for depth,outlet,nodes in choices:
            # This geometric requirement also certifies components excluded
            # before triangulation: graph detours cannot create real clearance.
            if np.linalg.norm(np.asarray(c)-np.asarray(g.nodes[outlet]['xy_um']))<terminal_separation:
                reject['terminal_clearance']+=1;continue
            raw,corridors=_path_line(g,nodes)
            base_variants=[]
            for policy,visibility_first in (('visibility_shortcut',True),('guide_preserving',False)):
                base=_shorten(raw.coords,nav.navigation_domain,rules.wire_width_um/2,
                              visibility_first=visibility_first)
                if any(base.equals_exact(prior,1e-7) for _,prior in base_variants):
                    continue
                base_variants.append((policy,base))
            for policy,base in base_variants:
                for offset in offsets:
                    shifted=base if offset==0 else base.offset_curve(offset,quad_segs=16,join_style='round')
                    if shifted.is_empty or shifted.geom_type!='LineString':continue
                    # Start stays at the elected anchor; no new routing shortcut through the island.
                    line=LineString([c,*shifted.coords])
                    # The center curve is first screened against the exact eroded
                    # vector domain. Expensive full-buffer certification follows
                    # only for paths that survive this necessary condition.
                    if not nav.center_domain.covers(line):
                        reject['wire_outside_original_support']+=1;continue
                    try:
                        line,curve=smooth_centerline(line.coords,nav.center_domain,curve_settings)
                    except UncertifiableCurve:
                        reject['uncertifiable_smooth_bend']+=1;continue
                    if np.linalg.norm(np.asarray(c)-np.asarray(line.coords[-1]))<terminal_separation:
                        reject['terminal_clearance']+=1;continue
                    conservative_wire_radius=wire_radius+curve_settings.max_chord_error_um
                    if not support.covers(line.buffer(conservative_wire_radius+rules.margin_um+settings.numeric_guard_um,quad_segs=16)):
                        reject['wire_outside_original_support']+=1;continue
                    wire=line.buffer(conservative_wire_radius/math.cos(math.pi/64),quad_segs=16)
                    metal=wire.union(electrode)
                    if not local_support.covers(metal.buffer(rules.margin_um,quad_segs=32)):
                        reject['metal_margin']+=1;continue
                    valid_sources.add(source)
                    columns.append({'source_node':source,'outlet_node':outlet,'source_um':c,
                                    'outlet_um':list(line.coords[-1]),'points_um':np.asarray(line.coords),
                                    'declared_terminal_um':g.nodes[outlet]['xy_um'],
                                    'terminal_kind':g.nodes[outlet]['terminal_kind'],
                                    'corridor_ids':corridors,'lane_offset_um':offset,
                                    'line_length_um':line.length,'depth_um':depth,
                                    'curve':curve,'route_shape_strategy':policy,
                                    'island':island,'electrode':electrode,'wire':wire,'metal':metal})
        if index%100==0:notify('生成并认证电极—轨道—出口路径列',.4+.12*index/max(1,len(nav.candidates)))
    return columns,{'anchors_checked':len(nav.candidates),'valid_source_anchors':len(valid_sources),
                    'path_columns':len(columns),'rejected':dict(reject),
                    'centerline_model':'piecewise_quadratic_bezier_G1',
                    'maximum_chord_error_um':curve_settings.max_chord_error_um,
                    'route_shape_strategies':['visibility_shortcut','guide_preserving'],
                    'lane_offsets_um':offsets,'candidate_step_um':settings.candidate_step_um,
                    'outlet_search':'shortest_clearance_eligible_geometry_terminal',
                    'anchor_terminal_search_gap_um':terminal_separation,
                    'include_guarded_outlet_alternative':
                        settings.include_guarded_outlet_alternative,
                    'routes_per_candidate':settings.routes_per_candidate}


def _physical_conflicts(shapes,gap,source_codes,pad_codes,chunk_size=128):
    """Find close or intersecting pairs, omitting group-redundant pairs."""
    if gap<0:raise ValueError('Physical gap must be nonnegative')
    # Zero process gap still forbids touching/overlapping independent nets.
    # Scale the equality tolerance for small positive gaps instead of turning
    # dwithin's distance negative and silently dropping every conflict.
    conflict_distance=gap-min(1e-6,gap*1e-6)
    # Pad route candidates can be stored as collections of already-certified
    # connected pieces.  Index each piece separately: a collection bounding
    # box stretches from the electrode to an outside Pad and otherwise makes
    # nearly every pair a costly broad-phase candidate.  A pair of unions is
    # closer than gap iff at least one pair of their constituent pieces is.
    if any(shape.geom_type=='GeometryCollection' for shape in shapes):
        parts=[];owners=[]
        for owner,shape in enumerate(shapes):
            pieces=(shape.geoms if shape.geom_type=='GeometryCollection' else (shape,))
            for piece in pieces:
                if not piece.is_empty:
                    parts.append(piece);owners.append(owner)
        if not parts:return set()
        owners=np.asarray(owners,dtype=np.int32)
        tree=STRtree(parts)
        conflicts=set()
        for start in range(0,len(parts),chunk_size):
            pair=tree.query(parts[start:start+chunk_size],
                            predicate='dwithin',distance=conflict_distance)
            if pair.shape[1]==0:continue
            left=owners[pair[0]+start]
            right=owners[pair[1]]
            eligible=(right<left)&(source_codes[left]!=source_codes[right])
            eligible &= ((pad_codes[left]<0)|(pad_codes[left]!=pad_codes[right]))
            conflicts.update((int(j),int(i)) for i,j in
                             zip(left[eligible],right[eligible]))
        return conflicts
    tree=STRtree(shapes)
    conflicts=set()
    for start in range(0,len(shapes),chunk_size):
        pair=tree.query(shapes[start:start+chunk_size],
                        predicate='dwithin',distance=conflict_distance)
        if pair.shape[1]==0:continue
        left=pair[0]+start
        right=pair[1]
        eligible=(right<left)&(source_codes[left]!=source_codes[right])
        eligible &= ((pad_codes[left]<0)|(pad_codes[left]!=pad_codes[right]))
        conflicts.update((int(j),int(i)) for i,j in
                         zip(left[eligible],right[eligible]))
    return conflicts


def _select(columns,rules,settings,notify,center_upper=None,pad_order=None,
            baseline_candidate_ids=None,center_origin=None):
    n=len(columns)
    if not n:return [],{'status':'empty_finite_library','finite_lower_bound':0,'finite_upper_bound':0}
    center_squared=(np.asarray([float(np.sum((np.asarray(c['source_um'],dtype=float)-
                                                np.asarray(center_origin,dtype=float))**2))
                                for c in columns]) if center_origin is not None else None)
    secondary_name=('minimize_sum_squared_radius_after_maximum_count'
                    if center_squared is not None else 'maximize_graph_depth_after_maximum_count')
    if center_upper==1:
        # Every pair in the continuous anchor domain is incompatible. Since
        # every column is already certified, the finite optimum is exactly 1.
        chosen=(columns[int(np.argmin(center_squared))] if center_squared is not None else
                max(columns,key=lambda c:c['depth_um']))
        return [chosen],{'status':0,'message':'continuous center-distance bound proves at most one',
                         'finite_lower_bound':1,'finite_upper_bound':1,
                         'finite_model_proven_optimal':True,'secondary_status':0,
                         'secondary_objective':secondary_name,
                         'secondary_proven_optimal':True,
                         'conflict_pairs':0,'center_spacing_conflict_pairs':0,
                         'constraint_rows':1,'variables':n,
                         'route_candidate_variables':n,
                         'pad_contiguous_per_side':pad_order is not None,
                         'pad_interval_start_variables':0,
                         'shared_corridors_permitted':True,'unit_edge_or_node_capacity_imposed':False}
    notify('构造实体金属与衬底圆岛冲突约束',.6)
    metals=[c['metal'] for c in columns];islands=[c['island'] for c in columns]
    source_index={key:i for i,key in enumerate(dict.fromkeys(c['source_node'] for c in columns))}
    source_codes=np.asarray([source_index[c['source_node']] for c in columns],dtype=np.int32)
    pad_keys=[key for key in dict.fromkeys(c.get('pad_id') for c in columns)
              if key is not None]
    pad_index={key:i for i,key in enumerate(pad_keys)}
    pad_codes=np.asarray([pad_index.get(c.get('pad_id'),-1) for c in columns],dtype=np.int32)
    conflicts=set()
    for shapes,gap in ((metals,rules.spacing_um),(islands,settings.pad_gap_um)):
        conflicts.update(_physical_conflicts(shapes,gap,source_codes,pad_codes))
    center_conflicts=set()
    if rules.minimum_center_spacing_um>0:
        centers=np.asarray([c['source_um'] for c in columns],dtype=float)
        # A center-distance rule is independent of metal and island clearance.
        # The tiny tolerance only avoids classifying an exact d pair as a conflict.
        center_conflicts={tuple(sorted((int(i),int(j)))) for i,j in
                          cKDTree(centers).query_pairs(rules.minimum_center_spacing_um-1e-7)}
        conflicts.update(center_conflicts)
    groups=defaultdict(list)
    for i,c in enumerate(columns):groups[c['source_node']].append(i)
    pad_groups=defaultdict(list)
    for i,c in enumerate(columns):
        if 'pad_id' in c:pad_groups[c['pad_id']].append(i)
    constructive_ids=None
    if center_upper == 4 and pad_order is not None:
        # Exact four-side witness search over already geometry-certified
        # columns. One Pad per side is automatically contiguous; bit masks
        # enforce every pairwise metal, island and center conflict plus source
        # exclusivity. This constructs a lower bound, never an upper bound.
        side_indices=defaultdict(list)
        for i,c in enumerate(columns):side_indices[c.get('pad_side')].append(i)
        if len(side_indices)>=3:
            invalid=[1<<i for i in range(n)]
            for i,j in conflicts:
                invalid[i]|=1<<j
                invalid[j]|=1<<i
            for indices in groups.values():
                if len(indices)>1:
                    mask=sum(1<<i for i in indices)
                    for i in indices:invalid[i]|=mask
            side_masks={side:sum(1<<i for i in indices)
                        for side,indices in side_indices.items()}
            ordered_sides=sorted(side_masks,key=lambda side:len(side_indices[side]))
            deadline=deadline_after(settings.milp_time_limit_s)
            def find_distinct_sides(sides,depth,available,chosen):
                if deadline_expired(deadline):return None
                if depth==len(sides):return chosen
                options=side_masks[sides[depth]]&available
                candidates=[]
                while options:
                    low=options&-options
                    i=low.bit_length()-1
                    options-=low
                    candidates.append(i)
                if center_squared is not None:
                    candidates.sort(key=lambda i:center_squared[i])
                for i in candidates:
                    result=find_distinct_sides(sides,depth+1,
                                               available&~invalid[i],chosen+[i])
                    if result is not None:return result
                return None
            if len(ordered_sides)>=4:
                constructive_ids=find_distinct_sides(ordered_sides[:4],0,
                                                     (1<<n)-1,[])
                if constructive_ids is not None:
                    return [columns[i] for i in constructive_ids],{
                        'status':0,
                        'message':'four distinct Pad-side routes attain the exact continuous center upper bound',
                        'finite_lower_bound':4,'finite_upper_bound':4,
                        'finite_model_proven_optimal':True,
                        'secondary_status':None,
                        'secondary_objective':secondary_name,
                        'secondary_proven_optimal':False,
                        'conflict_pairs':len(conflicts),
                        'center_spacing_conflict_pairs':len(center_conflicts),
                        'constraint_rows':None,'variables':n,
                        'route_candidate_variables':n,
                        'pad_contiguous_per_side':True,
                        'pad_interval_start_variables':0,
                        'shared_corridors_permitted':True,
                        'unit_edge_or_node_capacity_imposed':False,
                        'constructive_witness_ids':constructive_ids}
            if not deadline_expired(deadline):
                for sides in combinations(ordered_sides,3):
                    constructive_ids=find_distinct_sides(sides,0,(1<<n)-1,[])
                    if constructive_ids is not None:break
    if center_upper == 3 and pad_order is not None:
        # A constructive three-net witness closes the *continuous* upper
        # bound. Use distinct sides so each occupied Pad bank is automatically
        # contiguous; exact pairwise metal, island and center conflicts above
        # still apply. This avoids asking a huge MILP to rediscover a tiny
        # feasible witness before its time limit.
        side_indices=defaultdict(list)
        for i,c in enumerate(columns):
            side_indices[c.get('pad_side')].append(i)
        if len(side_indices)>=3:
            invalid=[1<<i for i in range(n)]
            for i,j in conflicts:
                invalid[i] |= 1<<j
                invalid[j] |= 1<<i
            for indices in groups.values():
                if len(indices)>1:
                    mask=sum(1<<i for i in indices)
                    for i in indices:invalid[i] |= mask
            side_masks={side:sum(1<<i for i in indices)
                        for side,indices in side_indices.items()}
            side_names=sorted(side_masks)
            for side_a,side_b,side_c in combinations(side_names,3):
                for i in sorted(side_indices[side_a],
                                key=(lambda index:center_squared[index])
                                    if center_squared is not None else
                                    (lambda index:-columns[index]['depth_um'])):
                    available_b=side_masks[side_b] & ~invalid[i]
                    while available_b:
                        low=available_b & -available_b
                        j=low.bit_length()-1
                        available_b-=low
                        available_c=side_masks[side_c] & ~invalid[i] & ~invalid[j]
                        if available_c:
                            k=(available_c & -available_c).bit_length()-1
                            return [columns[index] for index in (i,j,k)],{
                                'status':0,
                                'message':'three distinct Pad-side routes attain the exact continuous center-distance upper bound',
                                'finite_lower_bound':3,'finite_upper_bound':3,
                                'finite_model_proven_optimal':True,
                                'secondary_status':None,
                                    'secondary_objective':secondary_name,
                                    'secondary_proven_optimal':False,
                                'conflict_pairs':len(conflicts),
                                'center_spacing_conflict_pairs':len(center_conflicts),
                                'constraint_rows':None,'variables':n,
                                'route_candidate_variables':n,
                                'pad_contiguous_per_side':True,
                                'pad_interval_start_variables':0,
                                'shared_corridors_permitted':True,
                                'unit_edge_or_node_capacity_imposed':False,
                                'constructive_witness_ids':[i,j,k]}
    rows=[];cols=[];values=[];rhs=[]
    for indices in list(groups.values())+list(pad_groups.values()):
        if len(indices)>1:
            row=len(rhs);rows.extend([row]*len(indices));cols.extend(indices)
            values.extend([1.0]*len(indices));rhs.append(1.0)
    for i,j in sorted(conflicts):
        row=len(rhs);rows.extend([row,row]);cols.extend([i,j])
        values.extend((1.0,1.0));rhs.append(1.0)
    if center_upper is not None:
        row=len(rhs);rows.extend([row]*n);cols.extend(range(n))
        values.extend([1.0]*n);rhs.append(float(center_upper))
    # In Pad mode, occupied slots on each side form one interval.  A binary
    # start variable is required at every 0 -> 1 transition; at most one start
    # is allowed per side.  The slot occupancy itself is the sum of the
    # already binary route columns assigned to that physical Pad.
    starts={}
    if pad_order is not None:
        declared=[pad_id for ids in pad_order.values() for pad_id in ids]
        if len(declared)!=len(set(declared)) or set(pad_groups)-set(declared):
            raise ValueError('Pad order must uniquely include every candidate Pad')
        for side,ids in pad_order.items():
            if not any(pad_groups.get(pad_id) for pad_id in ids):continue
            side_starts=[]
            for index,pad_id in enumerate(ids):
                start=n+len(starts)
                starts[(side,pad_id)]=start;side_starts.append(start)
                row=len(rhs)
                current=pad_groups.get(pad_id,())
                previous=pad_groups.get(ids[index-1],()) if index else ()
                rows.extend([row]*(len(current)+len(previous)+1))
                cols.extend([*current,*previous,start])
                values.extend([1.0]*len(current)+[-1.0]*(len(previous)+1))
                rhs.append(0.0)
            row=len(rhs);rows.extend([row]*len(side_starts));cols.extend(side_starts)
            values.extend([1.0]*len(side_starts))
            rhs.append(1.0)
    variable_count=n+len(starts)
    matrix=coo_array((np.asarray(values),(np.asarray(rows,dtype=np.int32),np.asarray(cols,dtype=np.int32))),
                     shape=(len(rhs),variable_count)).tocsc()
    constraint=LinearConstraint(matrix,-np.inf,np.asarray(rhs))
    notify('联合选择电极、几何轨道、节点转向和出口',.72)
    objective=np.r_[-np.ones(n),np.zeros(len(starts))]
    integrality=np.ones(variable_count)
    bounds=Bounds(np.zeros(variable_count),np.ones(variable_count))
    baseline=None
    if baseline_candidate_ids is not None and len(baseline_candidate_ids)<n:
        allowed=np.zeros(variable_count)
        allowed[n:]=1
        allowed[np.asarray(baseline_candidate_ids,dtype=int)]=1
        baseline=milp(objective,integrality=integrality,
                      bounds=Bounds(np.zeros(variable_count),allowed),
                      constraints=constraint,
                      options=solver_options(settings.milp_time_limit_s,mip_rel_gap=0))
    opt=milp(objective,integrality=integrality,bounds=bounds,
             constraints=constraint,options=solver_options(settings.milp_time_limit_s,mip_rel_gap=0))
    incumbents=[result for result in (baseline,opt) if result is not None and result.x is not None]
    if not incumbents:
        ids=constructive_ids if constructive_ids is not None else [
            int(np.argmin(center_squared)) if center_squared is not None else
            max(range(n),key=lambda i:columns[i]['depth_um'])]
        return [columns[i] for i in ids],{
            'status':int(opt.status),
            'message':f'No finite MILP incumbent; retained a directly checked constructive witness. {opt.message}',
            'finite_lower_bound':len(ids),
            'finite_upper_bound':center_upper if center_upper is not None else n,
            'finite_model_proven_optimal':False,
            'secondary_status':None,
            'secondary_objective':secondary_name,
            'secondary_proven_optimal':False,
            'conflict_pairs':len(conflicts),
            'center_spacing_conflict_pairs':len(center_conflicts),
            'constraint_rows':len(rhs),'variables':variable_count,
            'route_candidate_variables':n,
            'physical_pad_identity_groups':len(pad_groups),
            'pad_contiguous_per_side':pad_order is not None,
            'pad_interval_start_variables':len(starts),
            'shared_corridors_permitted':True,
            'unit_edge_or_node_capacity_imposed':False,
            'constructive_witness_ids':ids}
    incumbent=max(incumbents,key=lambda result:np.count_nonzero(result.x[:n]>.5))
    ids=np.flatnonzero(incumbent.x[:n]>.5);count=len(ids)
    # Preserve the first-stage cardinality exactly.  In Pad mode the second
    # objective is geometric centrality, measured independently of graph path
    # depth.  The original non-Pad mode retains its graph-depth preference.
    if center_squared is not None:
        secondary_cost=center_squared/max(1.0,float(center_squared.max()))
    else:
        depth=np.asarray([c['depth_um'] for c in columns],dtype=float)
        secondary_cost=-depth/max(1.0,float(depth.max()))
    fixed=LinearConstraint(np.r_[np.ones(n),np.zeros(len(starts))][None,:],count,count)
    improve=milp(np.r_[secondary_cost,np.zeros(len(starts))],integrality=integrality,bounds=bounds,
                 constraints=[constraint,fixed],options=solver_options(settings.milp_time_limit_s,mip_rel_gap=0))
    secondary_accepted=False
    if improve.x is not None:
        proposed=np.flatnonzero(improve.x[:n]>.5)
        # A time-limited second solve can return a worse feasible incumbent
        # than the first one. Keep the first unless count and centrality
        # are both nondegrading; a solver timeout is not improvement.
        if (len(proposed)==count and
                secondary_cost[proposed].sum()<=secondary_cost[ids].sum()+1e-9):
            ids=proposed;secondary_accepted=True
    chosen=[columns[int(i)] for i in ids]
    upper=math.floor(-float(opt.mip_dual_bound)+1e-6) if getattr(opt,'mip_dual_bound',None) is not None else None
    return chosen,{'status':int(opt.status),'message':opt.message,'finite_lower_bound':count,
                   'finite_upper_bound':upper,'finite_model_proven_optimal':opt.status==0,
                   'baseline_finite_lower_bound':(int(np.count_nonzero(baseline.x[:n]>.5))
                                                   if baseline is not None and baseline.x is not None else None),
                   'baseline_status':(int(baseline.status) if baseline is not None else None),
                   'baseline_candidate_variables':(len(baseline_candidate_ids)
                                                    if baseline_candidate_ids is not None else None),
                    'secondary_status':int(improve.status),
                    'secondary_objective':secondary_name,
                    'secondary_proven_optimal':improve.status==0,
                    'secondary_incumbent_accepted':secondary_accepted,
                    'conflict_pairs':len(conflicts),
                   'center_spacing_conflict_pairs':len(center_conflicts),
                   'constraint_rows':len(rhs),'variables':variable_count,
                   'route_candidate_variables':n,
                   'physical_pad_identity_groups':len(pad_groups),
                   'pad_contiguous_per_side':pad_order is not None,
                   'pad_interval_start_variables':len(starts),
                   'shared_corridors_permitted':True,'unit_edge_or_node_capacity_imposed':False}


def _center_distance_upper_bound(domain,distance):
    """Conservative continuous bound for centers in any polygonal domain.

    The polygon's bounding vertices give an enclosing disk; the convex hull
    dilation gives a disjoint-disk area bound. Neither uses sampled anchors.
    """
    if distance<=0 or domain.is_empty:return None
    coords=shapely.get_coordinates(domain)
    left,bottom,right,top=domain.bounds
    origin=np.asarray([(left+right)/2,(bottom+top)/2])
    radius=float(np.linalg.norm(coords-origin,axis=1).max())+.01
    angular=None
    if distance>2*radius+.01:
        angular=1
    elif distance>radius+.01:
        for impossible_count in range(2,16):
            pair_limit=max(radius,2*radius*math.sin(math.pi/impossible_count))
            if distance>pair_limit+.01:
                angular=impossible_count-1
                break
    disk_radius=distance/2
    outer=disk_radius/math.cos(math.pi/128)
    enclosing=domain.convex_hull.buffer(outer,quad_segs=32)
    area_bound=math.ceil((enclosing.area*(1+1e-9)+1e-6)/(math.pi*disk_radius**2))
    return {'value':min(area_bound,angular) if angular is not None else area_bound,
            'minimum_center_separation_um':distance,
            'angular_upper_bound':angular,'area_upper_bound':area_bound,
            'enclosing_disk_center_um':origin.tolist(),'enclosing_disk_radius_um':radius,
            'domain':'supplied continuous polygonal center domain, independent of finite candidate sampling',
            'proof':'enclosing-disk angular pigeonhole and disjoint center disks in convex-hull dilation',
            'uses_interval_exact_arithmetic':False}


def continuous_center_upper_bound(nav, distance, *, rules=None):
    """Bound all admissible centers on the ORIGINAL support, across views.

    Navigation insets, collector removal, pruned components and alternative
    navigation views are search devices. None may narrow a continuous upper
    bound. The exact-grid cell cover is preferred whenever available.
    """
    if distance <= 0:
        return None
    path = nav.summary.get('input_gds_path')
    if not path:
        return None
    try:
        cell = gds_cell_cover_upper_bound(path, *nav.summary['support_layer'], distance)
    except ValueError:
        # Never promote a floating geometry estimate to a strict certificate.
        return None
    if cell['input_sha256'] != nav.summary['sha256']:
        raise ValueError('Input GDS changed after navigation construction')
    angular = gds_enclosing_disk_angular_upper_bound(
        path, *nav.summary['support_layer'], distance)
    if angular is not None and angular['input_sha256'] != cell['input_sha256']:
        raise ValueError('Input GDS changed while calculating angular bound')
    if angular is not None and not verify_gds_angular_certificate(path,angular):
        raise ValueError('Angular upper-bound certificate failed independent verification')
    small_cells=(gds_small_cell_conflict_upper_bound(
        path,*nav.summary['support_layer'],distance)
        if cell['cell_count']<=18 else None)
    if small_cells is not None and not verify_small_cell_conflict_certificate(path,small_cells):
        raise ValueError('Small-cell conflict certificate failed independent verification')
    values=[cell['value']]
    if angular:values.append(angular['value'])
    if small_cells:values.append(small_cells['value'])
    region_bound=region_for_navigation(nav,rules).packing_upper_bound(distance) if rules is not None else None
    if region_bound:values.append(region_bound['value'])
    return {'value':min(values),
            'minimum_center_separation_um': distance,
            'domain': 'original GDS support; superset of electrode centers in every navigation view',
            'grid_cell_cover': cell,
            'enclosing_disk_angular_certificate': angular,
            'small_cell_conflict_certificate':small_cells,
            'electrode_region_packing_certificate':region_bound,
            'upper_bound_method': ('exact_rational_electrode_disk_packing'
                                   if region_bound and region_bound['value']==min(values) else
                                   'exact_integer_enclosing_disk_angular_pigeonhole'
                                   if angular and angular['value']==min(values) else
                                   'exact_rational_cell_intersection_conflict_graph'
                                   if small_cells and small_cells['value']==min(values) else
                                   'exact_integer_square_cover_and_expanded_area'),
            'scope_note': 'No navigation inset, collector subtraction, component exclusion, route sampling, or fixed Pad frame is used.',
            'uses_interval_exact_arithmetic': False,
            'uses_exact_integer_grid_predicates': True}


def solve(nav,rules,*,settings=IslandSettings(),progress=None):
    notify=progress or (lambda *a,**k:None)
    island_radius=(rules.electrode_diameter_um/2+rules.margin_um+
                   settings.numeric_guard_um)
    required_center_gap=max(rules.minimum_center_spacing_um,
                            rules.electrode_diameter_um+rules.spacing_um,
                            2*island_radius+settings.pad_gap_um if settings.pad_gap_um>0 else 0)
    center_bound=continuous_center_upper_bound(nav,required_center_gap,rules=rules)
    columns,library=_route_columns(nav,rules,settings,notify)
    chosen,optimization=_select(columns,rules,settings,notify,
                                center_bound['value'] if center_bound else None)
    region_audit=region_for_navigation(nav,rules).audit(chosen)
    if not region_audit['passed']:
        raise RuntimeError('Selected electrode center outside permitted disk')
    centers=np.asarray([c['source_um'] for c in chosen],dtype=float)
    center_distances=[float(np.linalg.norm(centers[i]-centers[j]))
                      for i in range(len(chosen)) for j in range(i)]
    minimum_selected_center_distance=min(center_distances,default=None)
    if (minimum_selected_center_distance is not None and
            minimum_selected_center_distance<rules.minimum_center_spacing_um-1e-6):
        raise RuntimeError('Selected electrodes violate minimum center spacing')
    final=nav.support.union(unary_union([c['island'] for c in chosen])) if chosen else nav.support
    if betti(final)!=betti(nav.support):raise RuntimeError('Selected islands changed original substrate topology')
    records=[{k:v for k,v in c.items() if k not in ('metal','wire','electrode','island')} for c in chosen]
    terminal_counts=dict(Counter(c['terminal_kind'] for c in chosen))
    gaps=[a['metal'].distance(b['metal']) for i,a in enumerate(chosen) for b in chosen[:i]]
    # Each substrate island contains an ideal disk of radius R. Disjoint
    # islands at gap g imply disjoint center disks of radius rho=R+g/2.
    # All centers lie in S0, so the rho disks lie in S0 (+) B_rho.
    # A circumscribed polygonal dilation and upward integer rounding keep
    # this area-packing bound conservative; it is valid but usually loose.
    rho=rules.electrode_diameter_um/2+rules.margin_um+settings.numeric_guard_um+settings.pad_gap_um/2
    outer_radius=rho/math.cos(math.pi/128)
    enclosing=nav.support.buffer(outer_radius,quad_segs=32)
    area_upper=enclosing.area*(1+1e-9)+1e-6
    packing_upper=math.ceil(area_upper/(math.pi*rho*rho))
    upper=center_bound['value'] if center_bound else None
    return {'status':'checked_geometry_lower_bound' if chosen else 'no_finite_compatible_route',
            'model':'topology_preserving_attached_disks','rules':asdict(rules),'settings':asdict(settings),
            'electrode_region':region_audit,
            'outlet_mode':nav.summary['outlet_mode'],'retained_routes':len(chosen),
             'terminal_counts':terminal_counts,
             'route_completion_scope':'electrode_to_geometry_candidate_terminal',
             'pad_connection_verified':False,
            'candidate_library':library,'optimization':optimization,
            'substrate_topology_before':betti(nav.support),'substrate_topology_after':betti(final),
            'added_substrate_area_um2':final.area-nav.support.area,
            'substrate_island_diameter_um':rules.electrode_diameter_um+2*rules.margin_um+2*settings.numeric_guard_um,
            'minimum_inter_net_gap_um':min(gaps,default=None),'routes':records,'_chosen':chosen,
            'minimum_selected_center_distance_um':minimum_selected_center_distance,
            '_final_support':final,
            'continuous_packing_upper_bound':{'value':packing_upper,'packing_radius_um':rho,
                'enclosing_area_um2':area_upper,'proof':'disjoint equal disks contained in original-support dilation',
                'scope':'centers in original support; fixed circular islands with stated island gap',
                'tight_routing_cut_bound':False,'uses_interval_exact_arithmetic':False,
                'usable_as_strict_certificate':False},
            'center_distance_upper_bound':center_bound,'continuous_upper_bound':upper,
            'capacity_scope':'geometry-audited finite-library lower bound; exact-grid continuous center upper bound on original support when available; floating area estimate is diagnostic only'}


def write_gds(input_path,output_path,routing,support_spec):
    with TemporaryDirectory(prefix='island_gds_') as temp:
        source=Path(temp)/'input.gds';output=Path(temp)/'output.gds'
        source.write_bytes(Path(input_path).read_bytes());lib=gdstk.read_gds(str(source),unit=1e-6)
        tops=lib.top_level();names={c.name for c in lib.cells}
        name='TOPOLOGICAL_ISLAND_ROUTING'
        while name in names:name+='X'
        cell=lib.new_cell(name)
        for top in tops:cell.add(gdstk.Reference(top))
        used={p.layer for top in tops for p in top.get_polygons()}
        metal=20 if 20 not in used else next(n for n in range(1020,65000) if n not in used)
        marker=30 if 30 not in used and metal!=30 else next(n for n in range(1030,65000) if n not in used and n!=metal)
        for i,route in enumerate(routing['_chosen'],1):
            for geom,layer,datatype in ((route['island'],*support_spec),(route['metal'],metal,i),(route['electrode'],marker,i)):
                for polygon in _polygonal(geom):
                    shell=gdstk.Polygon(np.asarray(polygon.exterior.coords),layer=layer,datatype=datatype)
                    holes=[gdstk.Polygon(np.asarray(h.coords)) for h in polygon.interiors]
                    pieces=gdstk.boolean(shell,holes,'not',precision=.001,layer=layer,datatype=datatype) if holes else [shell]
                    for piece in pieces:cell.add(*piece.fracture(max_points=4000,precision=.001))
            # Keep the exact electrode marker polygon in the metal layer.
            for polygon in _polygonal(route['electrode']):
                cell.add(gdstk.Polygon(np.asarray(polygon.exterior.coords),layer=metal,datatype=i))
        # Avoid re-fracturing source polygons merely because they are brought
        # into a new wrapper cell; the default 199-point limit changes their
        # union by native-grid slivers on real curved inputs.
        lib.write_gds(str(output),max_points=4000)
        Path(output_path).write_bytes(output.read_bytes())
    return {'support_layer':list(support_spec),'metal_layer':metal,'electrode_marker_layer':marker}


def audit_gds(path,nav,routing,layers):
    with TemporaryDirectory(prefix='island_audit_') as temp:
        p=Path(temp)/'audit.gds';p.write_bytes(Path(path).read_bytes());lib=gdstk.read_gds(str(p),unit=1e-6)
    grouped=defaultdict(list)
    for top in lib.top_level():
        for p in top.get_polygons():grouped[(p.layer,p.datatype)].append(p)
    def union_layer(spec):
        # Join GDS fracture pieces on their actual 1 nm grid before computing
        # vector predicates. This avoids an unnecessarily huge GEOS overlay.
        records=grouped[spec]
        if not records:return Polygon()
        merged=gdstk.boolean(records,[],'or',precision=.001)
        return unary_union([part for p in merged for part in _polygonal(shapely.make_valid(Polygon(p.points)))])
    support=unary_union([part for p in grouped[tuple(layers['support_layer'])]
                         for part in _polygonal(shapely.make_valid(Polygon(p.points)))])
    expected=routing['_final_support']
    support_diff=support.symmetric_difference(expected).area
    # Newly generated island vertices and Boolean fractures are rounded to a
    # 1 nm GDS grid. Check a local 10 nm envelope as well as an area budget
    # proportional to the boundary that was actually written on that grid.
    island_perimeter=sum(c['island'].boundary.length for c in routing['_chosen'])
    support_area_tolerance=.003*island_perimeter+.01
    support_local_grid_check=(expected.buffer(.01).covers(support) and
                              support.buffer(.01).covers(expected))
    metal=[];checks=[];marker_centers=[];rules=routing['rules']
    for i,c in enumerate(routing['_chosen'],1):
        m=union_layer((layers['metal_layer'],i));marker=union_layer((layers['electrode_marker_layer'],i))
        covered=support.covers(m);margin=m.distance(support.boundary) if covered else -1
        wire_original=nav.support.buffer(.002).covers(c['wire'])
        outside=m.difference(nav.support)
        outside_only_electrode=(outside.is_empty or c['electrode'].buffer(.002).covers(outside))
        connected=betti(m)['components']==1
        endpoints=m.covers(Point(c['source_um'])) and m.covers(Point(c['outlet_um']))
        terminal_error=float(np.linalg.norm(np.asarray(c['outlet_um'])-np.asarray(c['declared_terminal_um'])))
        terminal_ok=terminal_error<=abs(c['lane_offset_um'])+.01
        disk_ok=marker.buffer(.002).covers(disk(c['source_um'],rules['electrode_diameter_um']/2))
        marker_center=np.asarray(marker.centroid.coords[0]) if not marker.is_empty else np.asarray([math.inf,math.inf])
        center_error=float(np.linalg.norm(marker_center-np.asarray(c['source_um'])))
        marker_centers.append(marker_center)
        difference=m.symmetric_difference(c['metal']).area/max(1,c['metal'].area)
        checks.append({'net':i,'connected':connected,'endpoints':endpoints,'support_containment':covered,
                       'margin_um':margin,'wire_inside_original_support':wire_original,
                       'outside_original_only_own_electrode':outside_only_electrode,
                       'electrode_diameter_check':disk_ok,'metal_relative_difference':difference,
                       'electrode_center_error_um':center_error,
                       'terminal_position_error_um':terminal_error,'declared_terminal_window_check':terminal_ok})
        metal.append(m)
    gaps=[a.distance(b) for i,a in enumerate(metal) for b in metal[:i]]
    center_distances=[float(np.linalg.norm(marker_centers[i]-marker_centers[j]))
                      for i in range(len(marker_centers)) for j in range(i)]
    minimum_center_distance=min(center_distances,default=None)
    center_spacing_ok=(minimum_center_distance is None or
                       minimum_center_distance>=rules['minimum_center_spacing_um']-.002)
    region_audit=region_for_navigation(nav,rules).audit([{'source_um':p} for p in marker_centers])
    separated=all(not a.intersects(b) for i,a in enumerate(metal) for b in metal[:i])
    passed=(support_diff<=support_area_tolerance and support_local_grid_check and
            betti(support)==betti(nav.support) and separated and center_spacing_ok and region_audit['passed'] and
            all(c['connected'] and c['endpoints'] and c['support_containment'] and c['margin_um']>=rules['margin_um']-.01
                and c['wire_inside_original_support'] and c['outside_original_only_own_electrode']
                and c['electrode_diameter_check'] and c['electrode_center_error_um']<=.002
                and c['declared_terminal_window_check']
                and c['metal_relative_difference']<1e-3 for c in checks) and
            (not gaps or min(gaps)>=rules['spacing_um']-.01))
    if not passed:raise RuntimeError(f'Exported island GDS audit failed: {checks}, gap {min(gaps,default=None)}, center gap {minimum_center_distance}, '
                                    f'support_difference={support_diff}, area_tolerance={support_area_tolerance}, '
                                    f'local_grid_check={support_local_grid_check}, topology={betti(support)}')
    return {'passed':True,'exported_nets':len(metal),'output_layers':layers,
            'electrode_region_audit':region_audit,
            'audit_scope':'electrode_to_declared_geometry_terminal',
            'pad_connection_verified':False,
            'support_symmetric_difference_um2':support_diff,
            'support_area_tolerance_um2':support_area_tolerance,
            'support_local_grid_check':support_local_grid_check,
            'topology_preserved':True,
            'minimum_metal_to_support_boundary_um':min((c['margin_um'] for c in checks),default=None),
            'minimum_inter_net_gap_um':min(gaps,default=None),'all_nets_connected_and_endpoints_included':True,
            'minimum_electrode_center_distance_um':minimum_center_distance,
            'minimum_center_spacing_rule_um':rules['minimum_center_spacing_um'],
            'center_spacing_passed':center_spacing_ok,
            'all_wires_inside_original_support':True,'net_checks':checks,
            'verification_tolerance_um':.01,'independent_geometry_kernel_used':False}


def save_navigation_graph(path,nav):
    record={'input_sha256':nav.summary['sha256'],'method':nav.summary['method'],
            'nodes':[{'id':int(n),**d} for n,d in nav.graph.nodes(data=True)],
            'corridors':[{'start':int(a),'end':int(b),
                          'points_um':e['points_um'].tolist(),
                          'length_um':float(e['weight']),
                          'corridor_id':int(e['corridor_id'])}
                         for a,b,e in nav.graph.edges(data=True)],
            'outlets':nav.outlets,'candidate_anchors':nav.candidates,
            'topology_certificate':nav.summary['triangle_dual_topology']}
    with gzip.open(path,'wt',encoding='utf-8') as file:
        json.dump(record,file,ensure_ascii=False,separators=(',',':'))


def preview_navigation(path,nav,routing=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection,PolyCollection
    fig,ax=plt.subplots(figsize=(9,9))
    fills=[];colors=[]
    def fill_geometry(geometry,color):
        for p in _polygonal(geometry):
            fills.append(np.asarray(p.exterior.coords));colors.append(color)
            for hole in p.interiors:
                fills.append(np.asarray(hole.coords));colors.append('white')
    fill_geometry(nav.support,'#d8e5e1')
    if routing is not None and '_shell' in routing:
        fill_geometry(routing['_shell'],'#e9efeb')
    if routing is None:
        lines=[e['points_um'] for _,_,e in nav.graph.edges(data=True)]
        ax.add_collection(LineCollection(lines,colors='#315f71',linewidths=.5,alpha=.7))
        if nav.outlets:
            points=np.asarray([nav.graph.nodes[n]['xy_um'] for n in nav.outlets])
            ax.scatter(points[:,0],points[:,1],s=5,color='#d26342',zorder=3)
    else:
        for i,c in enumerate(routing['_chosen']):
            color=plt.cm.turbo((i+.5)/max(1,len(routing['_chosen'])))
            if 'bridge' in c:
                fill_geometry(c['bridge'],'#a9c8bc')
            fill_geometry(c['island'],'#a678d3')
            fill_geometry(c['metal'],color)
    # Batch polygon paths: adding a Patch per million-vertex polygon made
    # Axes repeatedly calculate Bezier bounds, even for straight GDS edges.
    ax.add_collection(PolyCollection(fills,facecolors=colors,edgecolors='none',linewidths=0))
    ax.autoscale_view()
    ax.set_aspect('equal');ax.axis('off');fig.tight_layout(pad=.1)
    fig.savefig(path,dpi=180,facecolor='white');plt.close(fig)
