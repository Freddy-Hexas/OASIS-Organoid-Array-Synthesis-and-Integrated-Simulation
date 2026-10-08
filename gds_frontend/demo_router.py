"""One-track GDS-only routing baseline on a recovered support graph.

This is deliberately a lower-bound witness, not a multi-track or optimal
continuous router. Broad collector interfaces are geometry-derived candidates;
they are not electrical pads or assigned pins.
"""
from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory
import math

import gdstk
import networkx as nx
import numpy as np
import shapely
from scipy.optimize import Bounds, LinearConstraint, milp
from shapely.geometry import LineString, Point
from shapely.ops import unary_union

from frontend import FrontendResult, _polygonal
from candidate_sampling import add_corridor_candidates
from electrode_region import region_for


def route(result: FrontendResult, *, electrode_diameter_um=6.0,
          wire_width_um=2.0, spacing_um=2.0, margin_um=1.0,
          electrode_region_radius_um=3000.0,
          outlet_mode='collector', candidate_strategy='legacy'):
    g = result.graph.copy()
    if outlet_mode not in ('collector','open_tips'):
        raise ValueError('outlet_mode must be collector or open_tips')
    if outlet_mode=='collector':
        initial_outlets=[n for n,d in g.nodes(data=True) if d['kind']=='collector_interface']
    else:
        initial_outlets=[n for n,d in g.nodes(data=True) if g.degree(n)==1
                         and d['kind']=='leaf']
    if not initial_outlets:
        return {'status': 'no_geometry_derived_outlet', 'outlet_mode':outlet_mode,
                'routes': [], 'reason': 'No candidate outlet of the selected type was found.'}
    support = result.support
    placement_region=region_for(support,{'electrode_region_radius_um':electrode_region_radius_um},
        result.summary.get('gds_native_precision_m',1e-9)*1e6)
    shapely.prepare(support)
    if candidate_strategy == 'all_corridors':
        g, sampling = add_corridor_candidates(
            g, electrode_diameter_um=electrode_diameter_um,
            spacing_um=spacing_um)
        midpoint_count = sampling['candidates_added']
    elif candidate_strategy == 'legacy':
        # Preserve the previous baseline for an apples-to-apples comparison.
        simple=nx.Graph(g)
        simple.remove_edges_from(nx.selfloop_edges(simple))
        cycle_core=set(nx.k_core(simple,2).nodes)
        next_node=max(g.nodes)+1
        midpoint_count=0
        for u,v,k,data in list(g.edges(keys=True,data=True)):
            if u not in cycle_core or v not in cycle_core:
                continue
            points=np.asarray(data['points_um'])
            if len(points)<4:
                continue
            if np.linalg.norm(points[0]-np.asarray(g.nodes[u]['xy_um'])) > np.linalg.norm(points[-1]-np.asarray(g.nodes[u]['xy_um'])):
                points=points[::-1]
            middle=len(points)//2
            mid=next_node;next_node+=1;midpoint_count+=1
            g.remove_edge(u,v,k)
            g.add_node(mid,xy_um=points[middle].tolist(),kind='corridor_candidate',
                       collector_ids=[],window_radius_um=data['min_raster_clearance_um'])
            for a,b,segment in ((u,mid,points[:middle+1]),(mid,v,points[middle:])):
                attrs=dict(data)
                attrs['points_um']=segment
                attrs['length_um']=float(np.linalg.norm(np.diff(segment,axis=0),axis=1).sum())
                g.add_edge(a,b,**attrs)
        sampling={'candidates_added':midpoint_count,'scope':'cycle-bearing core only'}
    else:
        raise ValueError('candidate_strategy must be legacy or all_corridors')
    # Exact vector checks gate every corridor and proposed electrode centre.
    legal_edges = {}
    for eid, (u, v, k, data) in enumerate(g.edges(keys=True, data=True)):
        line = LineString(data['points_um']).simplify(result.pitch * .25,
                                                      preserve_topology=False)
        if support.covers(line.buffer(wire_width_um / 2 + margin_um, quad_segs=8)):
            legal_edges[eid] = (u, v, k, data)
    sources = [n for n, d in g.nodes(data=True)
               if d['kind'] in ('junction','corridor_candidate') and
               placement_region.contains(d['xy_um'],export_safe=True) and support.covers(
                   Point(d['xy_um']).buffer(electrode_diameter_um / 2 + margin_um,
                                            quad_segs=32))]
    outlets=initial_outlets
    if not sources:
        return {'status':'no_sampled_legal_electrode_candidate',
                'outlet_mode':outlet_mode,'candidate_strategy':candidate_strategy,
                'candidate_sampling':sampling,'sources_checked':0,
                'routes':[],
                'reason':'No legal electrode disk among sampled source positions; this does not prove geometric infeasibility.'}
    H = nx.DiGraph()
    start, finish = ('special', 'start'), ('special', 'finish')
    for n in g.nodes:
        H.add_edge(('node_in', n), ('node_out', n), capacity=1, weight=0)
    edge_token_to_id = {}
    for eid, (u, v, _, data) in legal_edges.items():
        entry, leave = ('corridor_in', eid), ('corridor_out', eid)
        edge_token_to_id[entry] = eid
        H.add_edge(entry, leave, capacity=1, weight=0)
        weight = max(1, round(data['length_um']))
        H.add_edge(('node_out', u), entry, capacity=1, weight=weight)
        H.add_edge(('node_out', v), entry, capacity=1, weight=weight)
        H.add_edge(leave, ('node_in', u), capacity=1, weight=0)
        H.add_edge(leave, ('node_in', v), capacity=1, weight=0)
    center = np.asarray(result.support.centroid.coords[0])
    for n in sources:
        radius = np.linalg.norm(np.asarray(g.nodes[n]['xy_um']) - center)
        H.add_edge(start, ('node_in', n), capacity=1, weight=int(radius * 10))
    for n in outlets:
        H.add_edge(('node_out', n), finish, capacity=1, weight=0)
    flow = nx.max_flow_min_cost(H, start, finish)
    proposals = []
    for source in sources:
        if flow[start].get(('node_in', source), 0) != 1:
            continue
        current = ('node_in', source)
        chosen_edges = []
        nodes = [source]
        seen = set()
        while current != finish:
            if current in seen:
                raise RuntimeError('Flow contains a cycle')
            seen.add(current)
            choices = [v for v, count in flow[current].items() if count > 0]
            if len(choices) != 1:
                raise RuntimeError('Flow is not a single path')
            nxt = choices[0]
            if nxt in edge_token_to_id:
                chosen_edges.append(edge_token_to_id[nxt])
            if nxt[0] == 'node_in' and isinstance(nxt[1], int) and nxt[1] != nodes[-1]:
                nodes.append(nxt[1])
            current = nxt
        if len(chosen_edges) != len(nodes)-1:
            raise RuntimeError('Flow corridor/node sequence mismatch')
        points = [g.nodes[source]['xy_um']]
        for a, b, eid in zip(nodes[:-1], nodes[1:], chosen_edges):
            _, _, _, data = legal_edges[eid]
            aa=np.asarray(g.nodes[a]['xy_um'])
            segment = (data['points_um'] if np.linalg.norm(data['points_um'][0]-aa)
                       <= np.linalg.norm(data['points_um'][-1]-aa)
                       else data['points_um'][::-1])
            points.extend(segment.tolist())
            points.append(g.nodes[b]['xy_um'])
        line = LineString(points).simplify(result.pitch * .25, preserve_topology=False)
        electrode = Point(g.nodes[source]['xy_um']).buffer(electrode_diameter_um / 2,
                                                           quad_segs=32)
        metal = unary_union([line.buffer(wire_width_um / 2, quad_segs=16), electrode])
        margin = support.boundary.distance(metal) if support.covers(metal) else -1.0
        valid = support.covers(metal.buffer(margin_um, quad_segs=8))
        proposals.append({'source_node': source, 'outlet_node': nodes[-1],
                          'source_um': g.nodes[source]['xy_um'],
                          'outlet_um': g.nodes[nodes[-1]]['xy_um'],
                          'corridor_ids': chosen_edges, 'graph_nodes': nodes,
                          'line_length_um': float(line.length), 'metal': metal,
                          'margin_um': float(margin), 'individual_valid': bool(valid)})
    if not proposals:
        return {'status':'no_sampled_route_to_selected_outlet',
                'outlet_mode':outlet_mode,'candidate_strategy':candidate_strategy,
                'candidate_sampling':sampling,'sources_checked':len(sources),
                'vector_legal_corridors':len(legal_edges), 'routes':[],
                'reason':'No path in the sampled one-track graph; this does not prove continuous geometric infeasibility.'}
    conflicts=[]
    for i, a in enumerate(proposals):
        for j in range(i):
            if (a['metal'].intersects(proposals[j]['metal']) or
                    a['metal'].distance(proposals[j]['metal']) < spacing_um - 1e-6):
                conflicts.append((i,j))
    n=len(proposals)
    if n:
        upper = np.asarray([int(r['individual_valid']) for r in proposals], dtype=float)
        constraints=[]
        if conflicts:
            rows=np.zeros((len(conflicts),n))
            for k,(i,j) in enumerate(conflicts):
                rows[k,i]=rows[k,j]=1
            constraints=[LinearConstraint(rows,-np.inf,1)]
        opt=milp(-np.ones(n),integrality=np.ones(n),bounds=Bounds(np.zeros(n),upper),
                 constraints=constraints)
        if opt.x is None or opt.status != 0:
            raise RuntimeError(f'Fixed-route conflict MILP did not prove optimality: status={opt.status}, message={opt.message}')
        chosen=[r for i,r in enumerate(proposals) if opt.x[i]>.5]
    else:
        chosen=[]
    gaps=[a['metal'].distance(b['metal']) for i,a in enumerate(chosen) for b in chosen[:i]]
    record_routes=[]
    for r in chosen:
        record_routes.append({k:v for k,v in r.items() if k!='metal'})
    return {'status': 'checked_geometry_lower_bound',
            'candidate_strategy':candidate_strategy,
            'candidate_sampling':sampling,
            'outlet_mode':outlet_mode,
            'outlet_semantics':('broad support collector interface, not a pad' if outlet_mode=='collector'
                                else 'experimental open support tips, not assigned pads or fanout'),
            'sources_checked': len(sources),
            'collector_interfaces':sum(d['kind']=='collector_interface' for _,d in g.nodes(data=True)),
            'outlet_candidates':len(outlets),
            'corridor_midpoint_candidates_added': midpoint_count,
            'vector_legal_corridors': len(legal_edges), 'graph_corridors': g.number_of_edges(),
            'initial_flow_paths': len(proposals),
            'individual_invalid': sum(not r['individual_valid'] for r in proposals),
            'pair_conflicts': len(conflicts), 'retained_routes': len(chosen),
            'fixed_proposal_milp_optimal': True,
            'minimum_metal_to_support_boundary_um': min((r['margin_um'] for r in chosen),default=None),
            'minimum_inter_net_gap_um': min(gaps,default=None),
            'rules': {'electrode_diameter_um': electrode_diameter_um,
                      'electrode_region_radius_um':electrode_region_radius_um,
                      'wire_width_um': wire_width_um, 'spacing_um': spacing_um,
                      'support_margin_um': margin_um},
            'routes': record_routes, '_chosen': chosen}


def write_demo_gds(input_path: Path, output_path: Path, routing: dict):
    if routing['status']!='checked_geometry_lower_bound':
        return
    raw=Path(input_path).read_bytes()
    with TemporaryDirectory(prefix='gds_routing_') as tmp:
        tmp=Path(tmp)
        ascii_input=tmp/'source.gds';ascii_output=tmp/'routed.gds'
        ascii_input.write_bytes(raw)
        lib=gdstk.read_gds(str(ascii_input),unit=1e-6)
        tops=lib.top_level()
        if not tops:
            raise ValueError('GDS has no top-level cell')
        if len(tops)==1:
            cell=tops[0]
        else:
            # Preserve every source top cell under one output wrapper, so the
            # added metal shares the same coordinate system as the extracted union.
            names={c.name for c in lib.cells}
            wrapper='GDS_WORKBENCH_ROUTED'
            while wrapper in names:
                wrapper+='X'
            cell=lib.new_cell(wrapper)
            for top in tops:
                cell.add(gdstk.Reference(top))
        used_layers={p.layer for top in tops for p in top.get_polygons()}
        metal_layer=20 if 20 not in used_layers else next(x for x in range(1020,65000) if x not in used_layers)
        marker_layer=30 if 30 not in used_layers and metal_layer!=30 else next(x for x in range(1030,65000) if x not in used_layers and x!=metal_layer)
        for idx,r in enumerate(routing['_chosen'],1):
            for polygon in _polygonal(r['metal']):
                exterior=gdstk.Polygon(np.asarray(polygon.exterior.coords),layer=metal_layer,datatype=idx)
                holes=[gdstk.Polygon(np.asarray(h.coords)) for h in polygon.interiors]
                pieces=gdstk.boolean([exterior],holes,'not',precision=.001,
                                     layer=metal_layer,datatype=idx) if holes else [exterior]
                for piece in pieces:
                    cell.add(*piece.fracture(max_points=4000,precision=.001))
            # Layer 30 is a location marker; passivation openings are not designed.
            cell.add(gdstk.ellipse(r['source_um'],routing['rules']['electrode_diameter_um']/2,
                                   tolerance=.002,layer=marker_layer,datatype=idx))
        lib.write_gds(str(ascii_output))
        Path(output_path).write_bytes(ascii_output.read_bytes())
    return {'metal_layer':metal_layer,'electrode_marker_layer':marker_layer}


def audit_demo_gds(output_path: Path, result: FrontendResult, routing: dict,
                   output_layers:dict):
    """Re-open exported GDS and check every network against its support layer."""
    raw=Path(output_path).read_bytes()
    with TemporaryDirectory(prefix='gds_audit_') as tmp:
        path=Path(tmp)/'routed.gds'
        path.write_bytes(raw)
        lib=gdstk.read_gds(str(path),unit=1e-6)
    grouped={}
    for cell in lib.top_level():
        for item in cell.get_polygons():
            grouped.setdefault((item.layer,item.datatype),[]).append(item.points)
    def union(spec):
        shapes=[]
        for points in grouped.get(spec,[]):
            p=shapely.make_valid(shapely.geometry.Polygon(points))
            shapes.extend(_polygonal(p))
        return unary_union(shapes)
    support=union(tuple(result.summary['support_layer']))
    if support.is_empty:
        raise RuntimeError('Exported GDS lacks support layer')
    support_difference_um2=float(support.symmetric_difference(result.support).area)
    support_relative_difference=(support_difference_um2/result.support.area)
    networks=[];marker_centers=[]
    for idx,r in enumerate(routing['_chosen'],1):
        marker=union((output_layers['electrode_marker_layer'],idx))
        if not marker.is_empty:marker_centers.append({'source_um':list(marker.centroid.coords[0])})
        metal=union((output_layers['metal_layer'],idx))
        if metal.is_empty:
            raise RuntimeError(f'Exported net {idx} missing')
        connected=len(_polygonal(metal))==1
        endpoint_ok=metal.covers(Point(r['source_um'])) and metal.covers(Point(r['outlet_um']))
        containment=support.covers(metal)
        margin=float(metal.distance(support.boundary)) if containment else -1.0
        marker_count=len(grouped.get((output_layers['electrode_marker_layer'],idx),[]))
        networks.append({'net':idx,'connected':connected,'endpoint_ok':endpoint_ok,
                         'support_containment':containment,'margin_um':margin,
                         'electrode_marker_count':marker_count,'metal':metal})
    pairs=[(a['metal'],b['metal']) for i,a in enumerate(networks) for b in networks[:i]]
    gaps=[a.distance(b) for a,b in pairs]
    minimum_gap=min(gaps,default=None)
    required_margin=routing['rules']['support_margin_um']
    required_gap=routing['rules']['spacing_um']
    # A zero spacing rule still means electrically separate networks.
    # Distance alone cannot distinguish a legal small gap from contact.
    all_nets_separate=all(not a.intersects(b) for a,b in pairs)
    region_audit=region_for(result.support,routing['rules'],
        result.summary.get('gds_native_precision_m',1e-9)*1e6).audit(marker_centers)
    passed=(support_relative_difference < 1e-6 and region_audit['passed'] and
            all(n['connected'] and n['endpoint_ok'] and n['support_containment']
                and n['margin_um'] >= required_margin-0.01 and n['electrode_marker_count']==1
                for n in networks) and
            all_nets_separate and
            (minimum_gap is None or minimum_gap >= required_gap-0.01))
    if not passed:
        raise RuntimeError(f'Exported GDS geometry audit failed: relative support difference '
                           f'{support_relative_difference:.3g}, net margins '
                           f'{[round(n["margin_um"],3) for n in networks]}, '
                           f'minimum gap {minimum_gap}')
    return {'passed':True,'exported_nets':len(networks),
            'electrode_region_audit':region_audit,
            'output_layers':output_layers,
            'support_symmetric_difference_um2':support_difference_um2,
            'support_relative_symmetric_difference':support_relative_difference,
            'minimum_metal_to_support_boundary_um':min((n['margin_um'] for n in networks),default=None),
            'minimum_inter_net_gap_um':minimum_gap,
            'all_nets_connected_and_endpoints_included':True,
            'electrode_location_markers_present':True}
