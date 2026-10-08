"""Exact continuous center upper bound for a small occupied-cell cover.

Each cell contains at most one center.  Exact rational clipping of source GDS
polygons to a cell gives a conservative bounding rectangle for every possible
center in it.  Two cells conflict only if *all* pairs in those rectangles are
strictly closer than the required center distance.  The maximum independent
set of this necessary-condition graph is therefore a global upper bound.
"""
from __future__ import annotations

from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations
from pathlib import Path
from tempfile import TemporaryDirectory
import json
import math

import gdstk

from capacity_bounds import (gds_cell_cover_upper_bound,_grid_vertices,
                             _validated_native_grid_um,_simple_integer_polygon)


SHIFTS=((0,0),(0,.5),(.5,0),(.5,.5))


def _clip(vertices,rectangle):
    points=[(Q(x),Q(y)) for x,y in vertices]
    x0,y0,x1,y1=rectangle
    for axis,edge,keep_greater in ((0,x0,True),(0,x1,False),
                                   (1,y0,True),(1,y1,False)):
        if not points:break
        result=[]
        previous=points[-1]
        def inside(point):
            return point[axis]>=edge if keep_greater else point[axis]<=edge
        for current in points:
            pi,ci=inside(previous),inside(current)
            if pi!=ci:
                delta=current[axis]-previous[axis]
                if delta==0:raise RuntimeError('Parallel edge changed clipping side')
                t=(Q(edge)-previous[axis])/delta
                result.append((previous[0]+t*(current[0]-previous[0]),
                               previous[1]+t*(current[1]-previous[1])))
            if ci:result.append(current)
            previous=current
        points=result
    return points


def _qpair(value):
    return [value.numerator,value.denominator]


def _parse_qpair(pair):
    if (not isinstance(pair,list) or len(pair)!=2 or
            not all(isinstance(v,int) for v in pair) or pair[1]<=0):
        raise ValueError('Invalid rational coordinate')
    return Q(*pair)


def _source_polygons(path,layer,datatype):
    raw=Path(path).read_bytes()
    with TemporaryDirectory(prefix='cell_conflict_source_') as temporary:
        local=Path(temporary)/'source.gds'
        local.write_bytes(raw)
        lib=gdstk.read_gds(str(local),unit=1e-6)
    grid=_validated_native_grid_um(lib)
    polygons=[]
    for top in lib.top_level():
        for polygon in top.get_polygons():
            if (polygon.layer,polygon.datatype)!=(layer,datatype):continue
            vertices=_grid_vertices(polygon.points,grid)
            xs=[p[0] for p in vertices];ys=[p[1] for p in vertices]
            bounds=(min(xs),min(ys),max(xs),max(ys))
            valid=_simple_integer_polygon(vertices)
            polygons.append((vertices,bounds,valid))
    if not polygons:raise ValueError('Selected GDS support is empty')
    return polygons,sha256(raw).hexdigest(),grid


def _cell_rectangles(occupied,cover,polygons):
    hx,hy=cover['cell_side_grid_ticks']
    ox,oy=cover['cell_origin_shift_grid_ticks']
    cells=[]
    for ix,iy in sorted(occupied):
        box=(ox+ix*hx,oy+iy*hy,ox+(ix+1)*hx,oy+(iy+1)*hy)
        points=[]
        for vertices,bounds,valid in polygons:
            if (bounds[2]<box[0] or bounds[0]>box[2] or
                    bounds[3]<box[1] or bounds[1]>box[3]):continue
            if not valid:
                # The cell cover already treats invalid polygons by their
                # bounding boxes; a full cell is a safe outer rectangle.
                points.extend(((Q(box[0]),Q(box[1])),
                               (Q(box[2]),Q(box[3]))))
                continue
            points.extend(_clip(vertices,box))
        if not points:
            # Closed-cell occupation can be caused by a boundary-only touch;
            # the full cell remains a safe overestimate.
            points=[(Q(box[0]),Q(box[1])),(Q(box[2]),Q(box[3]))]
        xs=[p[0] for p in points];ys=[p[1] for p in points]
        bounds=(min(xs),min(ys),max(xs),max(ys))
        cells.append({'index':[ix,iy],
                      'possible_center_bbox_ticks':[_qpair(v) for v in bounds]})
    return cells


def _conflicts(cells,distance_ticks):
    boxes=[tuple(_parse_qpair(pair) for pair in c['possible_center_bbox_ticks'])
           for c in cells]
    edges=[]
    for i,j in combinations(range(len(boxes)),2):
        a,b=boxes[i],boxes[j]
        dx=max(abs(a[0]-b[2]),abs(a[2]-b[0]))
        dy=max(abs(a[1]-b[3]),abs(a[3]-b[1]))
        if dx*dx+dy*dy<distance_ticks*distance_ticks:
            edges.append([i,j])
    return edges


def _independent_number(n,edges):
    adjacency=[0]*n
    for i,j in edges:
        adjacency[i]|=1<<j
        adjacency[j]|=1<<i
    best=[]
    def search(remaining,chosen):
        nonlocal best
        if len(chosen)+remaining.bit_count()<=len(best):return
        if not remaining:
            best=list(chosen);return
        vertex=max((v for v in range(n) if remaining>>v&1),
                   key=lambda v:(adjacency[v]&remaining).bit_count())
        without=remaining&~(1<<vertex)
        search(without&~adjacency[vertex],chosen+[vertex])
        search(without,chosen)
    search((1<<n)-1,[])
    return len(best),best


def gds_small_cell_conflict_upper_bound(gds_path,layer,datatype,
                                        center_distance_um,max_cells=18):
    polygons,digest,grid=_source_polygons(gds_path,layer,datatype)
    candidates=[]
    for shift in SHIFTS:
        cover,occupied=gds_cell_cover_upper_bound(
            gds_path,layer,datatype,center_distance_um,
            offset_fraction=shift,_return_occupied=True)
        if len(occupied)>max_cells:continue
        if cover['input_sha256']!=digest or cover['native_grid_um']!=grid:
            raise ValueError('GDS changed during small-cell certification')
        cells=_cell_rectangles(occupied,cover,polygons)
        edges=_conflicts(cells,cover['conservative_distance_grid_ticks'])
        value,witness=_independent_number(len(cells),edges)
        candidates.append({'value':value,
                           'method':'exact_rational_cell_intersection_conflict_graph',
                           'input_sha256':digest,
                           'support_layer':[layer,datatype],
                           'native_grid_um':grid,
                           'required_center_distance_um':center_distance_um,
                           'conservative_distance_grid_ticks':cover['conservative_distance_grid_ticks'],
                           'cell_shift_fraction':list(shift),
                           'cell_cover':cover,
                           'cells':cells,
                           'conflict_edges':edges,
                           'maximum_independent_cell_indices':witness,
                           'proof':'each occupied closed cell has diameter below d; rational GDS-clipped rectangles contain every possible center in that cell; every listed conflict pair has maximum possible squared distance strictly below d_grid squared; exhaustive independent-set size bounds all center placements'})
    return min(candidates,key=lambda item:item['value']) if candidates else None


def verify_small_cell_conflict_certificate(gds_path,certificate):
    if certificate is None:return False
    try:
        layer,datatype=certificate['support_layer']
        polygons,digest,grid=_source_polygons(gds_path,layer,datatype)
        if digest!=certificate['input_sha256'] or grid!=certificate['native_grid_um']:
            return False
        cover=certificate['cell_cover']
        shift=tuple(certificate['cell_shift_fraction'])
        if shift not in SHIFTS:return False
        fresh,occupied=gds_cell_cover_upper_bound(
            gds_path,layer,datatype,certificate['required_center_distance_um'],
            offset_fraction=shift,_return_occupied=True)
        if fresh!=cover:return False
        cells=_cell_rectangles(occupied,cover,polygons)
        if cells!=certificate['cells']:return False
        edges=_conflicts(cells,cover['conservative_distance_grid_ticks'])
        if edges!=certificate['conflict_edges']:return False
        # Independent verifier enumerates every subset instead of reusing
        # the production branch-and-bound search.
        n=len(cells)
        if n>20:return False
        edge_masks=[(1<<i)|(1<<j) for i,j in edges]
        best=max((mask.bit_count() for mask in range(1<<n)
                  if all(mask&edge !=edge for edge in edge_masks)),default=0)
        return best==certificate['value']
    except (KeyError,TypeError,ValueError,ZeroDivisionError):
        return False
