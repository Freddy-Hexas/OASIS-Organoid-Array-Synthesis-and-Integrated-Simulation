"""Source-blind GDS support graph extraction.

The input is a GDS support layer and geometric rule values. No generator source,
pre-labeled ports, structure name, or centerline is read here. A raster skeleton
is used for topology; the original vector GDS remains the geometry authority.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from collections import defaultdict
import hashlib
import math

import cv2
import gdstk
import networkx as nx
import numpy as np
import shapely
from scipy import ndimage
from shapely.geometry import Polygon
from shapely.ops import unary_union
from skimage.measure import euler_number
from skimage.morphology import skeletonize


@dataclass
class FrontendResult:
    support: object
    graph: nx.MultiGraph
    summary: dict
    mask: np.ndarray
    skeleton: np.ndarray
    wide: np.ndarray
    origin: tuple[float, float]  # left x, top y in micrometres
    pitch: float


def _polygonal(geom):
    if geom.is_empty:
        return []
    if isinstance(geom, Polygon):
        return [geom]
    return [p for item in getattr(geom, 'geoms', ()) for p in _polygonal(item)]


def read_support(gds_path: Path, layer: int = 10, datatype: int = 0):
    """Read support as vector union, including GDS paths flattened by gdstk."""
    gds_path = Path(gds_path)
    raw = gds_path.read_bytes()
    # gdstk on this Windows host cannot open a Unicode path directly.
    with TemporaryDirectory(prefix='gds_geometry_') as tmp:
        path = Path(tmp) / 'input.gds'
        path.write_bytes(raw)
        lib = gdstk.read_gds(str(path), unit=1e-6)
    tops = lib.top_level()
    if not tops:
        raise ValueError('GDS has no top-level cell')
    pieces = []
    for cell in tops:
        for item in cell.get_polygons():
            if (item.layer, item.datatype) != (layer, datatype):
                continue
            shape = Polygon(item.points)
            if not shape.is_valid:
                shape = shapely.make_valid(shape)
            pieces.extend(_polygonal(shape))
    if not pieces:
        raise ValueError(f'GDS has no support polygons on {layer}/{datatype}')
    support = unary_union(pieces)
    return support, pieces, {
        'sha256': hashlib.sha256(raw).hexdigest(),
        'gds_native_unit_m': lib.unit,
        'gds_native_precision_m': lib.precision,
        'top_cells': [cell.name for cell in tops],
        'support_layer': [layer, datatype],
        'support_polygon_records': len(pieces),
        'support_components': len(_polygonal(support)),
        'support_holes': sum(len(p.interiors) for p in _polygonal(support)),
        'bbox_um': list(support.bounds),
        'support_area_um2': support.area,
    }


def _rasterize(pieces, bounds, pitch):
    minx, miny, maxx, maxy = bounds
    left = math.floor(minx / pitch) * pitch - 10 * pitch
    top = math.ceil(maxy / pitch) * pitch + 10 * pitch
    right = math.ceil(maxx / pitch) * pitch + 10 * pitch
    bottom = math.floor(miny / pitch) * pitch - 10 * pitch
    width = int(round((right - left) / pitch)) + 1
    height = int(round((top - bottom) / pitch)) + 1
    if width*height>80_000_000:
        raise RuntimeError(f'Raster requires {width*height:,} pixels; tiling or a larger pitch is required')
    mask = np.zeros((height, width), dtype=np.uint8)
    for p in pieces:
        q = np.rint(np.column_stack(((np.asarray(p.exterior.xy[0]) - left) / pitch,
                                    (top - np.asarray(p.exterior.xy[1])) / pitch))).astype(np.int32)
        cv2.fillPoly(mask, [q], 1)
        for hole in p.interiors:
            q = np.rint(np.column_stack(((np.asarray(hole.xy[0]) - left) / pitch,
                                        (top - np.asarray(hole.xy[1])) / pitch))).astype(np.int32)
            cv2.fillPoly(mask, [q], 0)
    return mask, (left, top)


def _neighbors(y, x, skeleton):
    h, w = skeleton.shape
    for dy, dx in ((-1, 0), (1, 0), (0, -1), (0, 1),
                   (-1, -1), (-1, 1), (1, -1), (1, 1)):
        yy, xx = y + dy, x + dx
        if yy < 0 or xx < 0 or yy >= h or xx >= w or not skeleton[yy, xx]:
            continue
        # Drop diagonal shortcuts around a filled pixel corner. This keeps
        # raster branches from acquiring false triangular cycles.
        if dy and dx and (skeleton[y, xx] or skeleton[yy, x]):
            continue
        yield yy, xx


def _pixel_graph(skeleton):
    ys, xs = np.nonzero(skeleton)
    coords = list(zip(ys.tolist(), xs.tolist()))
    ids = np.full(skeleton.shape, -1, dtype=np.int32)
    ids[ys, xs] = np.arange(len(coords), dtype=np.int32)
    adjacency = [[] for _ in coords]
    for i, (y, x) in enumerate(coords):
        for yy, xx in _neighbors(y, x, skeleton):
            j = int(ids[yy, xx])
            if j > i:
                adjacency[i].append(j)
                adjacency[j].append(i)
    return coords, ids, adjacency


def _compress(skeleton, wide, wide_labels, origin, pitch, clearance):
    narrow = skeleton & ~wide
    coords, ids, adjacency = _pixel_graph(narrow)
    terminals = {i for i, adj in enumerate(adjacency) if len(adj) != 2}
    exit_to_collector = {}
    for i, (y, x) in enumerate(coords):
        neighbors = list(_neighbors(y, x, skeleton))
        attached = [int(wide_labels[yy, xx]) for yy, xx in neighbors if wide[yy, xx]]
        if attached:
            terminals.add(i)
            exit_to_collector[i] = min(attached)
    # A closed cycle without a junction still needs an anchor node.
    seen = set()
    for i in range(len(coords)):
        if i in seen:
            continue
        stack = [i]
        seen.add(i)
        component = []
        while stack:
            u = stack.pop()
            component.append(u)
            for v in adjacency[u]:
                if v not in seen:
                    seen.add(v)
                    stack.append(v)
        if not terminals.intersection(component):
            terminals.add(min(component))
    # Adjacent terminal pixels are one junction window, not separate nodes.
    terminal_node = {}
    node_members = []
    for i in sorted(terminals):
        if i in terminal_node:
            continue
        node_id = len(node_members)
        stack = [i]
        terminal_node[i] = node_id
        members = []
        while stack:
            u = stack.pop()
            members.append(u)
            for v in adjacency[u]:
                if v in terminals and v not in terminal_node:
                    terminal_node[v] = node_id
                    stack.append(v)
        node_members.append(members)
    g = nx.MultiGraph()
    left, top = origin
    def xy(i):
        y, x = coords[i]
        return (left + x * pitch, top - y * pitch)
    for node_id, members in enumerate(node_members):
        points = np.asarray([xy(i) for i in members])
        collectors = sorted({exit_to_collector[i] for i in members if i in exit_to_collector})
        kinds = ('collector_interface' if collectors else 'junction' if len(members)>1 or max(len(adjacency[i]) for i in members)>2 else 'leaf' if min(len(adjacency[i]) for i in members)<2 else 'cycle_anchor')
        g.add_node(node_id, xy_um=points.mean(axis=0).tolist(), kind=kinds,
                   collector_ids=collectors, window_radius_um=float(max(clearance[coords[i]] for i in members) * pitch),
                   raster_pixel_count=len(members))
    visited = set()
    def key(a, b):
        return (a, b) if a < b else (b, a)
    for start, members in enumerate(node_members):
        for u in members:
            for v in adjacency[u]:
                if v in terminal_node and terminal_node[v] == start:
                    visited.add(key(u, v))
                    continue
                if key(u, v) in visited:
                    continue
                path = [u]
                prev, current = u, v
                visited.add(key(u, v))
                while True:
                    path.append(current)
                    if current in terminal_node:
                        break
                    nxt = [q for q in adjacency[current] if q != prev]
                    if len(nxt) != 1:
                        raise RuntimeError('Nonterminal skeleton pixel is not degree two')
                    prev, current = current, nxt[0]
                    visited.add(key(prev, current))
                end = terminal_node[current]
                if start == end and len(path) <= 2:
                    continue
                points = np.asarray([xy(i) for i in path], dtype=np.float64)
                local = np.asarray([clearance[coords[i]] * pitch for i in path])
                g.add_edge(start, end, points_um=points,
                           length_um=float(np.linalg.norm(np.diff(points, axis=0), axis=1).sum()),
                           min_raster_clearance_um=float(local.min()),
                           median_raster_clearance_um=float(np.median(local)))
    expected = sum(len(a) for a in adjacency) // 2
    if len(visited) != expected:
        raise RuntimeError(f'Skeleton compression missed {expected-len(visited)} pixel edges')
    return g


def _absorb_junction_spurs(graph,maximum_length_um):
    """Move short terminal skeleton spurs into their surrounding node window.

    Geometry is not deleted: the original support remains the material region,
    and each absorbed branch is exported as a record. This is a graph
    simplification for routing, not a claim that the branch cannot host an electrode.
    """
    g=graph.copy()
    absorbed=[]
    def oriented(edge,node):
        p=edge['points_um']
        center=np.asarray(g.nodes[node]['xy_um'])
        return p if np.linalg.norm(p[0]-center)<=np.linalg.norm(p[-1]-center) else p[::-1]
    changed=True
    while changed:
        changed=False
        degree_before=dict(g.degree())
        for node in list(g.nodes):
            if degree_before[node]!=1 or g.nodes[node]['kind']=='collector_interface':
                continue
            u,v,k,data=next(iter(g.edges(node,keys=True,data=True)))
            neighbor=v if u==node else u
            if data['length_um']>maximum_length_um or degree_before[neighbor]<3:
                continue
            absorbed.append({'absorbed_leaf_xy_um':g.nodes[node]['xy_um'],
                             'junction_xy_um':g.nodes[neighbor]['xy_um'],
                             'length_um':data['length_um'],
                             'points_um':data['points_um'].tolist()})
            g.remove_node(node)
            changed=True
        # Merge degree-two sample nodes created by the above absorption.
        for node in list(g.nodes):
            if g.degree(node)!=2 or g.nodes[node]['kind']=='collector_interface':
                continue
            incident=list(g.edges(node,keys=True,data=True))
            if len(incident)!=2:
                continue  # One self-loop has degree two and remains a loop.
            (u1,v1,k1,e1),(u2,v2,k2,e2)=incident
            a=v1 if u1==node else u1
            b=v2 if u2==node else u2
            if a==node or b==node:
                continue
            p1=oriented(e1,node)[::-1];p2=oriented(e2,node)
            p=np.vstack([p1,np.asarray(g.nodes[node]['xy_um'])[None,:],p2])
            g.remove_node(node)
            g.add_edge(a,b,points_um=p,
                       length_um=float(np.linalg.norm(np.diff(p,axis=0),axis=1).sum()),
                       min_raster_clearance_um=min(e1['min_raster_clearance_um'],e2['min_raster_clearance_um']),
                       median_raster_clearance_um=min(e1['median_raster_clearance_um'],e2['median_raster_clearance_um']))
            changed=True
    g.graph['absorbed_junction_spurs']=absorbed
    return g


def extract(gds_path: Path, *, layer=10, datatype=0, first_pitch=2.0,
            minimum_pitch=0.5, collector_clearance_um=15.0,junction_spur_length_um=12.0):
    support, pieces, meta = read_support(gds_path, layer, datatype)
    expected_euler = meta['support_components'] - meta['support_holes']
    pitches = []
    pitch = first_pitch
    raster_repair = None
    while pitch >= minimum_pitch - 1e-9:
        mask, origin = _rasterize(_polygonal(support), support.bounds, pitch)
        raster_euler = euler_number(mask, connectivity=2)
        cc, _ = ndimage.label(mask, np.ones((3,3), dtype=np.uint8))
        raster_components = int(cc.max())
        pitches.append({'pitch_um': pitch, 'raster_euler': int(raster_euler),
                        'raster_components': raster_components,
                        'vector_euler': expected_euler})
        if raster_euler == expected_euler and raster_components == meta['support_components']:
            break
        # Polygon quantization can leave one-pixel false bridges or islands.
        # A one-pixel inward repair is permitted only if the recovered skeleton
        # has the exact vector topology and every skeleton point lies in the
        # original vector support. Later route geometry still uses that vector.
        if pitch <= 0.5:
            repaired = cv2.erode(mask, cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(3,3)))
            repaired_euler = euler_number(repaired, connectivity=2)
            repaired_components = ndimage.label(repaired, np.ones((3,3), dtype=np.uint8))[1]
            if repaired_euler == expected_euler and repaired_components == meta['support_components']:
                trial_skeleton = skeletonize(repaired.astype(bool))
                yy,xx = np.nonzero(trial_skeleton)
                inside = shapely.contains_xy(support, origin[0]+xx*pitch, origin[1]-yy*pitch)
                if bool(np.all(inside)):
                    mask = repaired
                    raster_repair={'method':'one_pixel_inward_topology_repair',
                                   'radius_um':pitch,'raw_euler':int(raster_euler),
                                   'repaired_euler':int(repaired_euler),
                                   'skeleton_points_inside_vector':int(len(xx))}
                    break
        pitch /= 2
    else:
        raise RuntimeError(f'Raster topology did not match vector GDS: {pitches}')
    skeleton = skeletonize(mask.astype(bool))
    clearance = cv2.distanceTransform(mask, cv2.DIST_L2, cv2.DIST_MASK_PRECISE)
    full_graph = _compress(skeleton, np.zeros_like(skeleton), np.zeros_like(mask),
                           origin, pitch, clearance)
    full_cycle_rank = (full_graph.number_of_edges() - full_graph.number_of_nodes()
                       + nx.number_connected_components(full_graph))
    if (nx.number_connected_components(full_graph) != meta['support_components']
            or full_cycle_rank != meta['support_holes']):
        raise RuntimeError('Recovered full skeleton graph disagrees with vector support topology: '
                           f'components={nx.number_connected_components(full_graph)} '
                           f'cycles={full_cycle_rank}; expected {meta["support_components"]} '
                           f'and {meta["support_holes"]}')
    wide = skeleton & (clearance * pitch >= collector_clearance_um)
    wide_labels, count = ndimage.label(wide, np.ones((3,3), dtype=np.uint8))
    # Components shorter than five raster pixels are isolated junction bulges,
    # not a broad collector region.
    sizes = np.bincount(wide_labels[wide], minlength=count+1)
    for label in range(1, count+1):
        if sizes[label] < 5:
            wide[wide_labels == label] = False
    wide_labels, count = ndimage.label(wide, np.ones((3,3), dtype=np.uint8))
    graph = _compress(skeleton, wide, wide_labels, origin, pitch, clearance)
    raw_n,raw_e=graph.number_of_nodes(),graph.number_of_edges()
    graph=_absorb_junction_spurs(graph,junction_spur_length_um)
    if graph.number_of_nodes() == 0:
        raise RuntimeError('No narrow support graph recovered')
    graph_cycles = graph.number_of_edges() - graph.number_of_nodes() + nx.number_connected_components(graph)
    # Every broad collector component that contains a cycle (the outside frame)
    # removes that cycle from the narrow graph; do not claim equality of raw ranks.
    interface_nodes = [n for n, d in graph.nodes(data=True) if d['kind']=='collector_interface']
    summary = {**meta,
        'raster_trials': pitches,
        'raster_repair': raster_repair,
        'selected_pitch_um': pitch,
        'skeleton_pixels': int(skeleton.sum()),
        'full_graph_cycle_rank': full_cycle_rank,
        'full_graph_matches_vector_topology': True,
        'collector_clearance_threshold_um': collector_clearance_um,
        'broad_collector_components': int(count),
        'collector_interface_windows': len(interface_nodes),
        'graph_nodes': graph.number_of_nodes(),
        'graph_edges': graph.number_of_edges(),
        'raw_graph_nodes':raw_n,
        'raw_graph_edges':raw_e,
        'absorbed_junction_spurs':len(graph.graph['absorbed_junction_spurs']),
        'junction_spur_length_threshold_um':junction_spur_length_um,
        'graph_components': nx.number_connected_components(graph),
        'graph_cycle_rank_after_collector_removal': graph_cycles,
        'graph_node_kinds': dict(__import__('collections').Counter(d['kind'] for _, d in graph.nodes(data=True))),
        'extractor_uses_generator_source': False,
        'collector_semantics': 'geometry-only candidate broad region, not a confirmed electrical pad',
    }
    return FrontendResult(support, graph, summary, mask, skeleton, wide, origin, pitch)
