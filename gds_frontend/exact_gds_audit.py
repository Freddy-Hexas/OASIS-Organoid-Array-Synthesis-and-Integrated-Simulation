"""Independent integer-grid audit of exported GDS metal and contacts.

This checks polygon connectivity, subset relations, boundary clearances and
inter-net spacing. It is deliberately separate from the Shapely checks used
while constructing routes. A passing result covers only the listed polygonal
GDS rules; it does not certify that the candidate library is exhaustive or
that every support addition obeys a future process rule.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import cached_property
from hashlib import sha256
from pathlib import Path
from tempfile import TemporaryDirectory
import json
import math

import gdstk
import numpy as np
import pyclipper
import shapely
from electrode_region import geometric_origin, exported_center_in_disk
from shapely.geometry import Polygon as ShapelyPolygon
from shapely.ops import unary_union
from shapely.strtree import STRtree


def _signed_area2(path):
    return sum(x0*y1-x1*y0 for (x0,y0),(x1,y1) in
               zip(path,path[1:]+path[:1]))


def _read_layer_paths(path, expected_precision_um):
    raw=Path(path).read_bytes()
    with TemporaryDirectory(prefix='integer_gds_audit_') as temporary:
        local=Path(temporary)/'layout.gds'
        local.write_bytes(raw)
        lib=gdstk.read_gds(str(local),unit=1e-6)
    tick_um=lib.precision*1e6
    if not math.isfinite(tick_um) or tick_um<=0 or abs(tick_um-expected_precision_um)>1e-10:
        raise ValueError('Unexpected GDS native grid')
    layers=defaultdict(list)
    for top in lib.top_level():
        for polygon in top.get_polygons():
            scaled=np.asarray(polygon.points,dtype=float)/tick_um
            rounded=np.rint(scaled)
            if (not np.isfinite(scaled).all() or np.max(np.abs(scaled))>=2**52 or
                    np.max(np.abs(scaled-rounded))>1e-3):
                raise ValueError('Flattened polygon has a vertex off the native grid')
            points=[(int(x),int(y)) for x,y in rounded]
            if len(points)>1 and points[0]==points[-1]:points.pop()
            if len(set(points))<3:continue
            if _signed_area2(points)==0:continue
            if _signed_area2(points)<0:points.reverse()
            layers[(polygon.layer,polygon.datatype)].append(points)
    return layers,sha256(raw).hexdigest(),tick_um


def _union(paths):
    if not paths:return []
    clip=pyclipper.Pyclipper()
    clip.AddPaths(paths,pyclipper.PT_SUBJECT,True)
    tree=clip.Execute2(pyclipper.CT_UNION,pyclipper.PFT_NONZERO,
                       pyclipper.PFT_NONZERO)
    return tree


def _canonical_path(path):
    points=list(path)
    smallest=min(points)
    starts=(index for index,point in enumerate(points) if point==smallest)
    return min(tuple(points[index:]+points[:index]) for index in starts)


def _contours(tree):
    if not tree:return []
    paths=[]
    def visit(node):
        for child in node.Childs:
            if len(child.Contour)>=3:paths.append((list(map(tuple,child.Contour)),child.IsHole))
            visit(child)
    visit(tree)
    return paths


def _paths(tree):
    return [path for path,_ in _contours(tree)]


def _outside(subject,container):
    if not subject:return False
    clip=pyclipper.Pyclipper()
    clip.AddPaths(_paths(subject),pyclipper.PT_SUBJECT,True)
    clip.AddPaths(_paths(container),pyclipper.PT_CLIP,True)
    difference=clip.Execute2(pyclipper.CT_DIFFERENCE,
                             pyclipper.PFT_NONZERO,pyclipper.PFT_NONZERO)
    return bool(_contours(difference))


def _filled_overlap(first,second):
    if not first or not second:return False
    clip=pyclipper.Pyclipper()
    clip.AddPaths(_paths(first),pyclipper.PT_SUBJECT,True)
    clip.AddPaths(_paths(second),pyclipper.PT_CLIP,True)
    return bool(_contours(clip.Execute2(pyclipper.CT_INTERSECTION,
                                        pyclipper.PFT_NONZERO,
                                        pyclipper.PFT_NONZERO)))


def _filled_intersection(first,second):
    clip=pyclipper.Pyclipper()
    clip.AddPaths(_paths(first),pyclipper.PT_SUBJECT,True)
    clip.AddPaths(_paths(second),pyclipper.PT_CLIP,True)
    return clip.Execute2(pyclipper.CT_INTERSECTION,
                         pyclipper.PFT_NONZERO,pyclipper.PFT_NONZERO)


def _audit_outer_policy_record(record,raw_source_paths,tick_um,margin_um,
                                minimum_wire_width_um):
    """Recompute the reference from original GDS vertices, not route metadata."""
    if record.get('revision') not in ('global_radial_outer_ports_v1', 'directional_convex_envelope_outer_ports_v2'):
        raise ValueError('Unrecognized outer exit policy')
    center=record.get('center_grid_ticks')
    if (not isinstance(center,list) or len(center)!=2 or
            any(type(v) is not int for v in center)):
        raise ValueError('Invalid native-grid outer reference center')
    if abs(record.get('native_grid_um',0)-tick_um)>1e-10:
        raise ValueError('Outer exit policy grid does not match GDS')
    points=np.asarray([p for path in raw_source_paths for p in path],dtype=float)
    circle=shapely.minimum_bounding_circle(shapely.multipoints(points).convex_hull)
    expected_center=[int(round(v)) for v in circle.centroid.coords[0]]
    if center!=expected_center:
        raise ValueError('Outer reference is not the grid-rounded GDS enclosing-circle center')
    maximum=max((x-center[0])**2+(y-center[1])**2
                for path in raw_source_paths for x,y in path)
    radius=math.isqrt(maximum)+(math.isqrt(maximum)**2<maximum)
    if (record.get('enclosing_radius_grid_ticks')!=radius or
            record.get('maximum_vertex_radius_squared_ticks')!=maximum):
        raise ValueError('Outer enclosing radius is not certified by original GDS vertices')
    width=record.get('wire_width_um');guard=record.get('numeric_guard_um')
    bridge_width=record.get('bridge_width_um')
    if (any(not isinstance(v,(int,float)) or not math.isfinite(v) for v in
            (width,guard,bridge_width)) or width<=0 or guard<0 or bridge_width<=0 or
            record.get('margin_um')!=margin_um or
            (minimum_wire_width_um is not None and width!=minimum_wire_width_um)):
        raise ValueError('Outer window fabrication rules do not match the audit')
    launch_depth=math.ceil((width/2+margin_um+guard)/tick_um-1e-9)+2
    contact_depth=launch_depth+math.ceil(bridge_width/(2*tick_um)-1e-9)+2
    if (record.get('launch_depth_grid_ticks')!=launch_depth or
            record.get('contact_depth_grid_ticks')!=contact_depth):
        raise ValueError('Outer windows are not derived from fabrication dimensions')
    if record['revision']=='directional_convex_envelope_outer_ports_v2':
        hull=_integer_hull([point for path in raw_source_paths for point in path])
        provided=record.get('envelope_vertices_grid_ticks',[])
        if len(provided)!=len(hull) or any(len(p)!=2 or any(type(v) is not int for v in p) for p in provided):
            raise ValueError('Invalid original-GDS convex envelope')
        expected_edges={frozenset((a,b)) for a,b in zip(hull,hull[1:]+hull[:1])}
        provided=[tuple(p) for p in provided]
        if {frozenset((a,b)) for a,b in zip(provided,provided[1:]+provided[:1])}!=expected_edges:
            raise ValueError('Outer envelope differs from the exact integer GDS hull')
        return tuple(center),_HullGate(tuple(hull),launch_depth-2),_HullGate(tuple(hull),contact_depth-2)
    return tuple(center),max(0,radius-launch_depth),max(0,radius-contact_depth)


def _integer_hull(points):
    points=sorted(set(map(tuple,points)))
    def cross(o,a,b):return (a[0]-o[0])*(b[1]-o[1])-(a[1]-o[1])*(b[0]-o[0])
    lower=[];upper=[]
    for point in points:
        while len(lower)>=2 and cross(lower[-2],lower[-1],point)<=0:lower.pop()
        lower.append(point)
    for point in reversed(points):
        while len(upper)>=2 and cross(upper[-2],upper[-1],point)<=0:upper.pop()
        upper.append(point)
    return lower[:-1]+upper[:-1]


@dataclass(frozen=True)
class _HullGate:
    hull: tuple
    depth: int

    @cached_property
    def planes(self):
        # floor(|edge|) gives a larger forbidden core than the ideal erosion.
        # Thus this integer check is conservative, including GDS rounding.
        result=[]
        for a,b in zip(self.hull,self.hull[1:]+self.hull[:1]):
            dx=b[0]-a[0];dy=b[1]-a[1]
            result.append((-dy,dx,-dy*a[0]+dx*a[1]+max(0,self.depth)*math.isqrt(dx*dx+dy*dy)))
        return tuple(result)

    def inequalities(self):
        return iter(self.planes)

    @cached_property
    def core_reference(self):
        # The enclosing-circle center can lie outside an eroded acute hull.
        # Clip the integer hull with exact rational half-planes to obtain a
        # point belonging to the actual forbidden core, including degeneracy.
        polygon=[tuple(map(Fraction,p)) for p in self.hull]
        for nx,ny,bound in self.planes:
            if not polygon:return None
            clipped=[]
            for a,b in zip(polygon,polygon[1:]+polygon[:1]):
                va=nx*a[0]+ny*a[1]-bound;vb=nx*b[0]+ny*b[1]-bound
                if va>=0:clipped.append(a)
                if (va>=0)!=(vb>=0):
                    fraction=va/(va-vb)
                    clipped.append((a[0]+fraction*(b[0]-a[0]),a[1]+fraction*(b[1]-a[1])))
            polygon=clipped
        if not polygon:return None
        return (sum(p[0] for p in polygon)/len(polygon),
                sum(p[1] for p in polygon)/len(polygon))

    def segment_enters(self,a,b):
        low=Fraction(0);high=Fraction(1)
        for nx,ny,bound in self.inequalities():
            start=nx*a[0]+ny*a[1]-bound
            change=nx*(b[0]-a[0])+ny*(b[1]-a[1])
            if change==0:
                if start<0:return False
                continue
            cut=Fraction(-start,change)
            if change>0:low=max(low,cut)
            else:high=min(high,cut)
            if low>high:return False
        return low<=high


def _rational_point_in_filled_tree(point,tree):
    x,y=point;crossings=0
    for path,_ in _contours(tree):
        for a,b in zip(path,path[1:]+path[:1]):
            cross=(b[0]-a[0])*(y-a[1])-(b[1]-a[1])*(x-a[0])
            if cross==0 and min(a[0],b[0])<=x<=max(a[0],b[0]) and min(a[1],b[1])<=y<=max(a[1],b[1]):
                return True
            if (a[1]>y)!=(b[1]>y):
                intersection=Fraction(a[0])+Fraction(b[0]-a[0],b[1]-a[1])*(y-a[1])
                crossings+=intersection>x
    return crossings%2==1


def _audit_bridge_outer_window(bridge,source,center,core_ticks):
    if isinstance(core_ticks,_HullGate):
        if any(core_ticks.segment_enters(a,b) for a,b in _segments(bridge)):
            raise ValueError('Bridge footprint enters the forbidden interior core')
        core_point=core_ticks.core_reference
        if core_point is not None and _rational_point_in_filled_tree(core_point,bridge):
            raise ValueError('Bridge footprint contains the forbidden interior core')
        if _topology_counts(_filled_intersection(bridge,source))!=(1,0):
            raise ValueError('Bridge does not have one contractible outer source contact')
        return
    if any(pyclipper.PointInPolygon(center,path)>0 and not hole
           for path,hole in _contours(bridge)):
        raise ValueError('Bridge footprint contains the interior reference center')
    if any(not _point_segment_distance_at_least(center,a,b,core_ticks**2)
           for a,b in _segments(bridge)):
        raise ValueError('Bridge footprint enters the forbidden interior core')
    contact=_filled_intersection(bridge,source)
    if _topology_counts(contact)!=(1,0):
        raise ValueError('Bridge does not have one contractible outer source contact')


def _audit_outer_width_path(points,source,center,launch_ticks):
    # Clipper rounds line/polygon intersection coordinates to native ticks.
    # Two additional ticks give a conservative bound on that rounding error;
    # squared-distance comparisons themselves are exact Python integers.
    clip=pyclipper.Pyclipper()
    clip.AddPath(points,pyclipper.PT_SUBJECT,False)
    clip.AddPaths(_paths(source),pyclipper.PT_CLIP,True)
    outside=pyclipper.OpenPathsFromPolyTree(clip.Execute2(
        pyclipper.CT_DIFFERENCE,pyclipper.PFT_NONZERO,pyclipper.PFT_NONZERO))
    if len(outside)!=1 or len(outside[0])<2:
        raise ValueError('Functional wire leaves and reenters the original support')
    if tuple(points[-1]) not in (tuple(outside[0][0]),tuple(outside[0][-1])):
        raise ValueError('Original-support exterior segment is not the final Pad connection')
    if isinstance(launch_ticks,_HullGate):
        if any(launch_ticks.segment_enters(tuple(a),tuple(b)) for path in outside for a,b in zip(path[:-1],path[1:])):
            raise ValueError('Functional wire leaves original support before the outer exit window')
        return
    threshold2=(launch_ticks+2)**2
    if any(not _point_segment_distance_at_least(center,tuple(a),tuple(b),threshold2)
           for path in outside for a,b in zip(path[:-1],path[1:])):
        raise ValueError('Functional wire leaves original support before the outer exit window')


def _count_outer(tree):
    return sum(not hole for _,hole in _contours(tree))


def _topology_counts(tree):
    contours=_contours(tree)
    return (sum(not hole for _,hole in contours),
            sum(hole for _,hole in contours))


def _independent_polygon_topology(paths):
    """Cross-check source/island topology on raw native-grid polygons.

    Clipper can split a component at coincident slivers in a heavily
    overlapping source.  GEOS union is evaluated separately on integer grid
    coordinates.  This is an independent GDS readback cross-check, not a
    claim of exact topological certification for arbitrary degeneracies.
    """
    polygons=[]
    for path in paths:
        polygon=ShapelyPolygon(path)
        if not polygon.is_valid:
            # The source GDS may contain a self-touching ring.  This
            # independent cross-check uses GEOS filled-region repair; its
            # result is explicitly weaker than exact topology certification.
            polygon=polygon.buffer(0)
        if polygon.is_empty or polygon.geom_type not in ('Polygon','MultiPolygon'):
            raise ValueError('Invalid support polygon in topology audit')
        polygons.append(polygon)
    union=unary_union(polygons)
    components=([union] if union.geom_type=='Polygon' else
                list(union.geoms) if union.geom_type=='MultiPolygon' else None)
    if components is None:raise ValueError('Nonpolygonal support in topology audit')
    return (len(components),sum(len(p.interiors) for p in components))


def _segments(tree):
    return [(a,b) for path,_ in _contours(tree)
            for a,b in zip(path,path[1:]+path[:1]) if a!=b]


def _orientation(a,b,c):
    return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])


def _on_segment(a,b,p):
    return (_orientation(a,b,p)==0 and min(a[0],b[0])<=p[0]<=max(a[0],b[0])
            and min(a[1],b[1])<=p[1]<=max(a[1],b[1]))


def _segments_intersect(a,b,c,d):
    o1,o2=_orientation(a,b,c),_orientation(a,b,d)
    o3,o4=_orientation(c,d,a),_orientation(c,d,b)
    return (((o1==0 and _on_segment(a,b,c)) or
             (o2==0 and _on_segment(a,b,d)) or
             (o3==0 and _on_segment(c,d,a)) or
             (o4==0 and _on_segment(c,d,b))) or
            ((o1>0)!=(o2>0) and (o3>0)!=(o4>0)))


def _point_segment_distance_at_least(p,a,b,threshold2):
    vx,vy=b[0]-a[0],b[1]-a[1]
    wx,wy=p[0]-a[0],p[1]-a[1]
    length2=vx*vx+vy*vy
    projection=wx*vx+wy*vy
    if projection<=0:return wx*wx+wy*wy>=threshold2
    if projection>=length2:
        dx,dy=p[0]-b[0],p[1]-b[1]
        return dx*dx+dy*dy>=threshold2
    cross=wx*vy-wy*vx
    return cross*cross>=threshold2*length2


def _segment_distance_at_least(first,second,threshold2):
    a,b=first;c,d=second
    if _segments_intersect(a,b,c,d):return False
    return (_point_segment_distance_at_least(a,c,d,threshold2) and
            _point_segment_distance_at_least(b,c,d,threshold2) and
            _point_segment_distance_at_least(c,a,b,threshold2) and
            _point_segment_distance_at_least(d,a,b,threshold2))


def _segment_index(segments):
    if not segments:return None
    lines=shapely.linestrings(np.asarray(segments,dtype=float))
    return STRtree(lines)


def _nearby_pairs_clear(left_segments,right_segments,right_tree,required_ticks):
    """Integer exact distance; float STRtree only indexes exact integer boxes."""
    if not left_segments or not right_segments:return False,0
    threshold2=required_ticks*required_ticks
    # A zero-width box can be empty for horizontal/vertical segments. Expand
    # the broad phase by one tick; the exact predicate still uses threshold 0.
    broad_ticks=max(1,required_ticks)
    checked=0
    for start in range(0,len(left_segments),256):
        chunk=left_segments[start:start+256]
        # Native GDS integer boxes are exactly represented by these doubles.
        # Batch only the broad phase; every accepted distance predicate below
        # still uses arbitrary-precision integer orientation and squares.
        bounds=np.asarray([(min(a[0],b[0])-broad_ticks,min(a[1],b[1])-broad_ticks,
                            max(a[0],b[0])+broad_ticks,max(a[1],b[1])+broad_ticks)
                           for a,b in chunk],dtype=float)
        pairs=right_tree.query(shapely.box(*bounds.T))
        for left,index in zip(*pairs):
            checked+=1
            if not _segment_distance_at_least(chunk[int(left)],right_segments[int(index)],threshold2):
                return False,checked
    return True,checked


def _certified_marker_center(marker,diameter_ticks):
    """Return a twice-grid center whose full ideal disk is inside the marker."""
    contours=_contours(marker)
    if len(contours)!=1 or contours[0][1]:
        raise ValueError('Electrode marker is not one simply connected polygon')
    path=contours[0][0]
    xs=[p[0] for p in path];ys=[p[1] for p in path]
    center2=(min(xs)+max(xs),min(ys)+max(ys))
    doubled=[(2*x,2*y) for x,y in path]
    if pyclipper.PointInPolygon(center2,doubled)<=0:
        raise ValueError('Electrode marker bounding-box center is outside metal')
    threshold2=diameter_ticks*diameter_ticks
    if any(not _point_segment_distance_at_least(center2,a,b,threshold2)
           for a,b in zip(doubled,doubled[1:]+doubled[:1])):
        raise ValueError('Electrode marker fails minimum inscribed diameter')
    return center2


def _center_in_original_paths(center2,raw_paths):
    """A half-grid electrode center must lie in an original filled GDS polygon.

    GDS polygons have no holes; source holes are empty regions between their
    filled polygons.  Testing the raw polygons avoids a union round trip and
    accepts a center on an original polygon boundary.
    """
    for path in raw_paths:
        doubled=[(2*x,2*y) for x,y in path]
        if pyclipper.PointInPolygon(center2,doubled)!=0:
            return True
    return False


def _point_in_filled_tree(point,tree):
    crossings=0
    for path,_ in _contours(tree):
        state=pyclipper.PointInPolygon(point,path)
        if state==-1:return True
        crossings+=state==1
    return crossings%2==1


def _verify_width_path(points,metal,electrode,pad,radius_ticks):
    if len(points)<2 or not _point_in_filled_tree(points[0],electrode):
        raise ValueError('Wire width witness does not start in its electrode')
    if not _point_in_filled_tree(points[-1],pad):
        raise ValueError('Wire width witness does not end in its Pad')
    if not _point_in_filled_tree(points[0],metal):
        raise ValueError('Wire width witness starts outside metal')
    boundary=_segments(metal)
    index=_segment_index(boundary)
    segments=[(a,b) for a,b in zip(points,points[1:]) if a!=b]
    clear,checked=_nearby_pairs_clear(segments,boundary,index,radius_ticks)
    if not clear:raise ValueError('Wire width witness has less than minimum radius')
    return checked


def _certified_support_disk(center2,support,radius_ticks):
    """Prove an ideal centered disk stays in the final support polygon."""
    # Work in twice-grid coordinates to accept half-tick electrode centers.
    doubled_center=tuple(center2)
    doubled_paths=[([(2*x,2*y) for x,y in path],hole)
                   for path,hole in _contours(support)]
    crossings=0
    for path,_ in doubled_paths:
        state=pyclipper.PointInPolygon(doubled_center,path)
        if state==-1:
            return False
        crossings += state==1
    if crossings%2!=1:
        return False
    threshold2=(2*radius_ticks)**2
    return all(_point_segment_distance_at_least(
        doubled_center,a,b,threshold2)
        for path,_ in doubled_paths
        for a,b in zip(path,path[1:]+path[:1]))


def _certified_island_envelope(center2,island,radius_ticks):
    """A polygon whose every vertex is in a disk lies in that convex disk."""
    contours=_contours(island)
    if len(contours)!=1 or contours[0][1]:return False
    limit2=(2*radius_ticks)**2
    return all((2*x-center2[0])**2+(2*y-center2[1])**2<=limit2
               for x,y in contours[0][0])


def _certified_pad_bbox(pad,short_ticks,long_ticks):
    contours=_contours(pad)
    if len(contours)!=1 or contours[0][1]:
        return False
    points=contours[0][0]
    width=max(p[0] for p in points)-min(p[0] for p in points)
    height=max(p[1] for p in points)-min(p[1] for p in points)
    # An elongated diagonal or L-shaped contact can have the same bounding
    # box as a real Pad.  Require the exported marker itself to be a complete
    # axis-aligned rectangle, as specified by the four-side Pad model.
    corners={(min(p[0] for p in points),min(p[1] for p in points)),
             (min(p[0] for p in points),max(p[1] for p in points)),
             (max(p[0] for p in points),min(p[1] for p in points)),
             (max(p[0] for p in points),max(p[1] for p in points))}
    return (len(points)==4 and set(points)==corners and
            min(width,height)>=short_ticks and
            max(width,height)>=long_ticks)


def audit_support_polygon_provenance(input_gds,output_gds,layers,*,expected_precision_um=.001):
    """Exact raw-grid polygon multiset proof; no floating overlay tolerance."""
    original,input_hash,input_tick=_read_layer_paths(input_gds,expected_precision_um)
    exported,output_hash,tick=_read_layer_paths(output_gds,expected_precision_um)
    if input_tick!=tick:raise ValueError('Input and output GDS grids differ')
    support_key=tuple(layers['support_layer'])
    original_paths=Counter(map(_canonical_path,original.get(support_key,())))
    output_paths=Counter(map(_canonical_path,exported.get(support_key,())))
    if not original_paths:raise ValueError('Missing original support polygons')
    if any(output_paths[p]<count for p,count in original_paths.items()):
        raise ValueError('Original support is not exactly preserved')
    expected=original_paths.copy()
    for key in ('island_marker_layer','bridge_marker_layer','shell_marker_layer'):
        layer=layers[key]
        if not any(item_layer==layer for item_layer,_ in exported):
            raise ValueError('Missing generated support provenance marker')
        for (item_layer,_),paths in exported.items():
            if item_layer==layer:expected.update(map(_canonical_path,paths))
    if expected!=output_paths:raise ValueError('Exported support contains unmarked or altered geometry')
    return {'passed':True,'input_sha256':input_hash,'output_sha256':output_hash,
            'native_grid_um':tick,'original_support_exactly_preserved':True,
            'all_added_support_polygons_exactly_match_markers':True,
            'predicate':'canonical native-grid polygon multisets including multiplicities',
            'scope':'support source and addition provenance only; shape, spacing and width audit remains mandatory'}


def audit_exported_gds(input_gds,output_gds,layers,*,
                       wire_spacing_um,metal_support_margin_um,
                       expected_nets=None,expected_precision_um=.001,
                       require_exact_source_preservation=True,
                       minimum_electrode_diameter_um=None,
                       minimum_electrode_center_distance_um=None,
                       maximum_electrode_center_radius_um=None,
                       minimum_substrate_disk_radius_um=None,
                       maximum_substrate_disk_radius_um=None,
                       minimum_island_spacing_um=None,
                       minimum_pad_short_side_um=None,
                       minimum_pad_long_side_um=None,
                       minimum_wire_width_um=None,progress=None):
    """Audit exact native-grid polygons against modeled subset/distance rules."""
    for name,value in (('wire_spacing_um',wire_spacing_um),
                       ('metal_support_margin_um',metal_support_margin_um),
                       ('minimum_island_spacing_um',minimum_island_spacing_um),
                       ('minimum_electrode_center_distance_um',minimum_electrode_center_distance_um)):
        if value is not None and (isinstance(value,bool) or not isinstance(value,(int,float)) or
                                  not math.isfinite(value) or value<0):
            raise ValueError(f'{name} must be finite and nonnegative')
    if maximum_electrode_center_radius_um is not None and (
            isinstance(maximum_electrode_center_radius_um,bool) or
            not isinstance(maximum_electrode_center_radius_um,(int,float)) or
            not math.isfinite(maximum_electrode_center_radius_um) or maximum_electrode_center_radius_um<=0):
        raise ValueError('Electrode region radius must be finite and positive')
    if progress:progress('准备独立整数 GDS 审计与原始支撑证明',.89)
    original,input_hash,input_tick=_read_layer_paths(input_gds,expected_precision_um)
    exported,output_hash,tick_um=_read_layer_paths(output_gds,expected_precision_um)
    if input_tick!=tick_um:raise ValueError('Input and output GDS grids differ')
    support_key=tuple(layers['support_layer'])
    raw_source_paths=original.get(support_key,[])
    source_support=_union(raw_source_paths)
    output_support=_union(exported.get(support_key,[]))
    if not source_support or not output_support:
        raise ValueError('Missing source or exported support')
    region_origin_ticks=None
    if maximum_electrode_center_radius_um is not None:
        if minimum_electrode_diameter_um is None:
            raise ValueError('Electrode region audit requires a certified marker center')
        # Recompute from original vertices, not a supplied center or Pad frame.
        hull=shapely.MultiPoint([p for path in raw_source_paths for p in path]).convex_hull
        origin=geometric_origin(hull,1.0)
        region_origin_ticks=tuple(int(v) for v in origin)
    # Exact raw-polygon preservation is stronger and avoids a Clipper union
    # degeneracy at coincident 1 nm slivers in heavily overlapping support.
    original_paths=Counter(map(_canonical_path,raw_source_paths))
    output_paths=Counter(map(_canonical_path,exported.get(support_key,[])))
    source_preserved=all(output_paths[path]>=count
                         for path,count in original_paths.items())
    if require_exact_source_preservation and not source_preserved:
        raise ValueError('Original support is not exactly contained in exported support')
    metal_layer=int(layers['metal_layer'])
    marker_layer=int(layers['electrode_marker_layer'])
    pad_layer=int(layers['pad_contact_layer'])
    bridge_layer=int(layers['bridge_marker_layer'])
    has_support_provenance=('island_marker_layer' in layers and
                            'shell_marker_layer' in layers)
    if ('island_marker_layer' in layers)!=('shell_marker_layer' in layers):
        raise ValueError('Island and shell provenance layers must be supplied together')
    if has_support_provenance:
        island_layer=int(layers['island_marker_layer'])
        shell_layer=int(layers['shell_marker_layer'])
    ids=sorted({datatype for layer,datatype in exported if layer==metal_layer})
    if expected_nets is not None and len(ids)!=expected_nets:
        raise ValueError(f'Expected {expected_nets} metal nets, found {len(ids)}')
    if ids!=sorted({datatype for layer,datatype in exported if layer==marker_layer}) or ids!=sorted({datatype for layer,datatype in exported if layer==pad_layer}) or ids!=sorted({datatype for layer,datatype in exported if layer==bridge_layer}):
        raise ValueError('Electrode, Pad and bridge markers do not match metal ids')
    if has_support_provenance:
        if ids!=sorted({datatype for layer,datatype in exported if layer==island_layer}):
            raise ValueError('Island provenance markers do not match metal ids')
        if sorted(datatype for layer,datatype in exported if layer==shell_layer)!=[0]:
            raise ValueError('Expected exactly one outer shell marker datatype')
        shell=_union(exported[(shell_layer,0)])
        if _count_outer(shell)!=1 or _filled_overlap(shell,source_support):
            raise ValueError('Outer shell is fragmented or overlaps original support')
        expected_support_paths=original_paths.copy()
        for layer in (island_layer,bridge_layer,shell_layer):
            for (item_layer,_),paths in exported.items():
                if item_layer==layer:
                    expected_support_paths.update(map(_canonical_path,paths))
        if output_paths!=expected_support_paths:
            raise ValueError('Exported support contains unmarked or altered geometry')
        raw_islands=[path for network in ids
                     for path in exported[(island_layer,network)]]
        source_topology=_independent_polygon_topology(raw_source_paths)
        islands_topology=_independent_polygon_topology([*raw_source_paths,*raw_islands])
        if islands_topology!=source_topology:
            raise ValueError('Added islands change original support topology')
        protected_holes=[]
        for path,hole in _contours(source_support):
            if not hole:continue
            bounds=(min(p[0] for p in path),min(p[1] for p in path),
                    max(p[0] for p in path),max(p[1] for p in path))
            protected_holes.append((_union([path if _signed_area2(path)>0 else path[::-1]]),bounds))
    required_spacing=math.ceil(wire_spacing_um/tick_um-1e-9)
    required_margin=math.ceil(metal_support_margin_um/tick_um-1e-9)
    support_segments=_segments(output_support)
    support_tree=_segment_index(support_segments)
    nets={}
    islands={}
    checked_support_pairs=0
    checked_width_pairs=0
    electrode_centers_twice_ticks=[]
    diameter_ticks=(math.ceil(minimum_electrode_diameter_um/tick_um-1e-9)
                    if minimum_electrode_diameter_um is not None else None)
    if diameter_ticks is not None and diameter_ticks<=0:
        raise ValueError('Minimum electrode diameter must be positive')
    substrate_radius_ticks=(math.ceil(minimum_substrate_disk_radius_um/tick_um-1e-9)
                            if minimum_substrate_disk_radius_um is not None else None)
    if substrate_radius_ticks is not None and diameter_ticks is None:
        raise ValueError('Substrate disk audit requires electrode center audit')
    island_envelope_ticks=(math.ceil(maximum_substrate_disk_radius_um/tick_um-1e-9)
                           if maximum_substrate_disk_radius_um is not None else None)
    if island_envelope_ticks is not None and not has_support_provenance:
        raise ValueError('Island envelope audit requires support provenance layers')
    if (island_envelope_ticks is not None and substrate_radius_ticks is not None
            and island_envelope_ticks<substrate_radius_ticks):
        raise ValueError('Maximum island radius is below the required disk radius')
    short_ticks=(math.ceil(minimum_pad_short_side_um/tick_um-1e-9)
                 if minimum_pad_short_side_um is not None else None)
    long_ticks=(math.ceil(minimum_pad_long_side_um/tick_um-1e-9)
                if minimum_pad_long_side_um is not None else None)
    if (short_ticks is None)!=(long_ticks is None):
        raise ValueError('Both Pad side minima are required together')
    width_witness=None
    exit_record=layers.get('outer_exit_policy')
    exit_gate=None
    if exit_record is not None:
        if not has_support_provenance:
            raise ValueError('Outer exit audit requires bridge and shell provenance')
        exit_gate=_audit_outer_policy_record(exit_record,raw_source_paths,tick_um,
                                            metal_support_margin_um,minimum_wire_width_um)
    if minimum_wire_width_um is not None or exit_gate is not None:
        witness_name=layers.get('wire_width_witness_path')
        if not witness_name:raise ValueError('Minimum wire width audit requires a path witness')
        witness_bytes=Path(witness_name).read_bytes()
        width_witness=json.loads(witness_bytes)
        if (width_witness.get('format')!='exact_grid_wire_centerline_v1' or
                width_witness.get('output_sha256')!=output_hash or
                abs(width_witness.get('native_grid_um',0)-tick_um)>1e-10):
            raise ValueError('Wire width witness does not match exported GDS')
        witness_routes=width_witness.get('networks',[])
        if exit_gate is not None and width_witness.get('outer_exit_policy')!=exit_record:
            raise ValueError('Wire path witness does not bind the outer exit policy')
        if ([route.get('network') for route in witness_routes]!=ids or
                any(not isinstance(route.get('points_grid_ticks'),list)
                    for route in witness_routes)):
            raise ValueError('Wire width witnesses do not match network ids')
        width_witness_hash=sha256(witness_bytes).hexdigest()
        if minimum_wire_width_um is not None:
            wire_radius_ticks=math.ceil(minimum_wire_width_um/(2*tick_um)-1e-9)
            if wire_radius_ticks<=0:raise ValueError('Minimum wire width must be positive')
    if minimum_island_spacing_um is not None and not has_support_provenance:
        raise ValueError('Island spacing audit requires support provenance layers')
    island_gap_ticks=(math.ceil(minimum_island_spacing_um/tick_um-1e-9)
                      if minimum_island_spacing_um is not None else None)
    for ordinal,network in enumerate(ids,1):
        if progress and (ordinal==1 or ordinal%5==0):
            progress(f'独立整数审计：网络 {ordinal}/{len(ids)}',.90+.06*ordinal/len(ids))
        metal=_union(exported[(metal_layer,network)])
        electrode=_union(exported[(marker_layer,network)])
        pad=_union(exported[(pad_layer,network)])
        bridge=_union(exported[(bridge_layer,network)])
        island=(_union(exported[(island_layer,network)])
                if has_support_provenance else None)
        if (_count_outer(metal)!=1 or _count_outer(electrode)!=1 or
                _count_outer(pad)!=1 or _count_outer(bridge)!=1 or
                (has_support_provenance and _count_outer(island)!=1)):
            raise ValueError(f'Network {network} is disconnected or has fragmented markers')
        if has_support_provenance:
            if not _filled_overlap(island,source_support):
                raise ValueError(f'Network {network} island is not attached to original support')
            if (not _filled_overlap(bridge,source_support) or
                    not _filled_overlap(bridge,shell)):
                raise ValueError(f'Network {network} bridge does not join source and shell')
            bridge_points=[point for path in _paths(bridge) for point in path]
            bridge_bounds=(min(p[0] for p in bridge_points),min(p[1] for p in bridge_points),
                           max(p[0] for p in bridge_points),max(p[1] for p in bridge_points))
            if any(not (bounds[2]<bridge_bounds[0] or bridge_bounds[2]<bounds[0] or
                        bounds[3]<bridge_bounds[1] or bridge_bounds[3]<bounds[1]) and
                   _filled_overlap(bridge,hole) for hole,bounds in protected_holes):
                raise ValueError(f'Network {network} bridge enters a protected source hole')
        if exit_gate is not None:
            try:
                _audit_bridge_outer_window(bridge,source_support,exit_gate[0],exit_gate[2])
            except ValueError as exc:
                raise ValueError(f'Network {network}: {exc}') from exc
        bridge_raw=Counter(map(_canonical_path,exported[(bridge_layer,network)]))
        metal_raw=Counter(map(_canonical_path,exported[(metal_layer,network)]))
        electrode_raw=Counter(map(_canonical_path,exported[(marker_layer,network)]))
        pad_raw=Counter(map(_canonical_path,exported[(pad_layer,network)]))
        electrode_is_metal=all(metal_raw[path]>=count
                               for path,count in electrode_raw.items())
        pad_is_metal=all(metal_raw[path]>=count
                         for path,count in pad_raw.items())
        bridge_is_support=all(output_paths[path]>=count
                              for path,count in bridge_raw.items())
        if (_outside(metal,output_support) or
                (not electrode_is_metal and _outside(electrode,metal)) or
                (not pad_is_metal and _outside(pad,metal)) or
                (not bridge_is_support and _outside(bridge,output_support))):
            raise ValueError(f'Network {network} fails integer polygon containment')
        if diameter_ticks is not None:
            try:
                center2=_certified_marker_center(electrode,diameter_ticks)
            except ValueError as exc:
                raise ValueError(f'Network {network}: {exc}') from exc
            electrode_centers_twice_ticks.append((network,center2))
            if region_origin_ticks is not None and not exported_center_in_disk(
                    center2,region_origin_ticks,maximum_electrode_center_radius_um,tick_um):
                raise ValueError(f'Network {network} electrode center outside permitted disk')
            if not _center_in_original_paths(center2,raw_source_paths):
                raise ValueError(f'Network {network} electrode center is outside original support')
            if (substrate_radius_ticks is not None and
                    not _certified_support_disk(center2,
                                                island if has_support_provenance else output_support,
                                                substrate_radius_ticks)):
                raise ValueError(f'Network {network} lacks the substrate disk')
            if (island_envelope_ticks is not None and
                    not _certified_island_envelope(center2,island,island_envelope_ticks)):
                raise ValueError(f'Network {network} island exceeds its radial envelope')
        if (short_ticks is not None and
                not _certified_pad_bbox(pad,short_ticks,long_ticks)):
            raise ValueError(f'Network {network} Pad is smaller than required')
        if width_witness is not None:
            raw_points=witness_routes[ids.index(network)]['points_grid_ticks']
            if any(not isinstance(point,list) or len(point)!=2 or
                   any(not isinstance(value,int) for value in point)
                   for point in raw_points):
                raise ValueError(f'Network {network} width witness has invalid points')
            try:
                if minimum_wire_width_um is not None:
                    checked_width_pairs+=_verify_width_path(
                        [tuple(point) for point in raw_points],metal,electrode,pad,
                        wire_radius_ticks)
                if exit_gate is not None:
                    _audit_outer_width_path([tuple(point) for point in raw_points],
                                           source_support,exit_gate[0],exit_gate[1])
            except ValueError as exc:
                raise ValueError(f'Network {network}: {exc}') from exc
        metal_segments=_segments(metal)
        # At zero requested margin, containment already proves the rule;
        # touching the support boundary is allowed. Positive margins retain
        # the exact boundary-distance test.
        clear,checked=(_nearby_pairs_clear(metal_segments,support_segments,
                                           support_tree,required_margin)
                       if required_margin>0 else (True,0))
        checked_support_pairs+=checked
        if not clear:
            raise ValueError(f'Network {network} violates metal-to-support margin')
        nets[network]={'metal':metal,'segments':metal_segments,
                       'tree':_segment_index(metal_segments)}
        if has_support_provenance:
            island_segments=_segments(island)
            islands[network]={'shape':island,'segments':island_segments,
                              'tree':_segment_index(island_segments)}
    checked_network_pairs=0
    if minimum_electrode_center_distance_um is not None:
        if diameter_ticks is None:
            raise ValueError('Center distance audit requires electrode diameter audit')
        distance2=math.ceil(2*minimum_electrode_center_distance_um/tick_um-1e-9)
        for index,(network,center) in enumerate(electrode_centers_twice_ticks):
            for other,previous in electrode_centers_twice_ticks[:index]:
                dx,dy=center[0]-previous[0],center[1]-previous[1]
                if dx*dx+dy*dy<distance2*distance2:
                    raise ValueError(f'Electrode centers {other} and {network} violate distance')
    checked_metal_segments=0
    checked_island_segments=0
    if island_gap_ticks is not None:
        for pos,network in enumerate(ids):
            for other in ids[:pos]:
                first=islands[network];second=islands[other]
                if _filled_overlap(first['shape'],second['shape']):
                    raise ValueError(f'Islands {other} and {network} overlap')
                clear,checked=_nearby_pairs_clear(first['segments'],second['segments'],
                                                  second['tree'],island_gap_ticks)
                checked_island_segments+=checked
                if not clear:
                    raise ValueError(f'Islands {other} and {network} violate spacing')
    for pos,network in enumerate(ids):
        for other in ids[:pos]:
            checked_network_pairs+=1
            # Boundary distance alone misses a net contained in another net.
            # With zero gap, overlapping metal must still never certify.
            if required_spacing==0 and _filled_overlap(nets[network]['metal'],nets[other]['metal']):
                raise ValueError(f'Networks {other} and {network} violate spacing: metal overlaps')
            clear,checked=_nearby_pairs_clear(nets[network]['segments'],
                                               nets[other]['segments'],
                                               nets[other]['tree'],
                                               required_spacing)
            checked_metal_segments+=checked
            if not clear:
                raise ValueError(f'Networks {other} and {network} violate spacing')
    return {'passed':True,'scope':'exact polygonal GDS net connectivity, marker contact, bridge containment, metal support margin and inter-net spacing',
            'input_sha256':input_hash,'output_sha256':output_hash,
            'original_support_exactly_preserved':source_preserved,
            'native_grid_um':tick_um,'network_count':len(ids),
            'minimum_inter_net_spacing_um':wire_spacing_um,
            'minimum_metal_support_margin_um':metal_support_margin_um,
            'network_pairs_checked':checked_network_pairs,
            'nearby_metal_segment_pairs_checked':checked_metal_segments,
            'nearby_support_segment_pairs_checked':checked_support_pairs,
            'nearby_island_segment_pairs_checked':checked_island_segments,
            'nearby_wire_width_segment_pairs_checked':checked_width_pairs,
            'minimum_electrode_diameter_um':minimum_electrode_diameter_um,
            'minimum_electrode_center_distance_um':minimum_electrode_center_distance_um,
            'maximum_electrode_center_radius_um':maximum_electrode_center_radius_um,
            'electrode_region_reference_um':([v*tick_um for v in region_origin_ticks] if region_origin_ticks is not None else None),
            'electrode_region_verified':region_origin_ticks is not None,
            'minimum_substrate_disk_radius_um':minimum_substrate_disk_radius_um,
            'maximum_substrate_disk_radius_um':maximum_substrate_disk_radius_um,
            'minimum_island_spacing_um':minimum_island_spacing_um,
            'minimum_pad_short_side_um':minimum_pad_short_side_um,
            'minimum_pad_long_side_um':minimum_pad_long_side_um,
            'electrode_disks_and_center_spacing_verified':diameter_ticks is not None,
            'electrode_centers_in_original_support_verified':diameter_ticks is not None,
            'substrate_disks_verified':substrate_radius_ticks is not None,
            'island_envelopes_verified':island_envelope_ticks is not None,
            'support_provenance_verified':has_support_provenance,
            'source_island_topology_verified':has_support_provenance,
            'bridge_attachment_and_hole_exclusion_verified':has_support_provenance,
            'outer_exit_policy_verified':exit_gate is not None,
            'outer_exit_policy':exit_record,
            'outer_exit_audit_scope':'original-GDS reference and process-sized windows recomputed; whole bridge avoids interior core; one source contact; functional centerline exits once at an outer window' if exit_gate is not None else None,
            'island_island_spacing_verified':island_gap_ticks is not None,
            'pad_dimensions_verified':short_ticks is not None,
            'functional_wire_width_verified':minimum_wire_width_um is not None,
            'minimum_wire_width_um':minimum_wire_width_um,
            'wire_width_witness_sha256':(width_witness_hash
                                         if width_witness is not None else None),
            'predicate':'integer orientation and squared point-to-segment distance; exact integer polygon Boolean',
            'not_covered':', '.join((
                ([] if island_gap_ticks is not None else ['island-island spacing'])+
                ([] if exit_gate is not None else ['outer exit location and bridge shortcut policy'])+
                ['detailed bridge/shell width and other shape policy',
                 'local minimum wire width at unused metal spurs',
                 'complete manufacturing DRC',
                 'exhaustive routing optimality']))}


def rebuild_with_exact_source_and_contacts(input_gds,old_output_gds,
                                           new_output_gds,layers):
    """Re-export an existing routed wrapper against unchanged input cells.

    This is for auditing historical routes without re-solving the MILP. Only
    direct generated wrapper polygons are copied; original input cells come
    from the original GDS. Electrode and Pad marker polygons are duplicated
    into their network's metal layer before the independent integer audit.
    """
    with TemporaryDirectory(prefix='rebuild_exact_gds_') as temporary:
        folder=Path(temporary)
        source=folder/'source.gds'
        routed=folder/'routed.gds'
        output=folder/'output.gds'
        source.write_bytes(Path(input_gds).read_bytes())
        routed.write_bytes(Path(old_output_gds).read_bytes())
        source_lib=gdstk.read_gds(str(source),unit=1e-6)
        routed_lib=gdstk.read_gds(str(routed),unit=1e-6)
        original_tops=source_lib.top_level()
        old_tops=routed_lib.top_level()
        if not original_tops or len(old_tops)!=1:
            raise ValueError('Expected original top cells and one routed wrapper')
        old=old_tops[0]
        original_names={cell.name for cell in original_tops}
        referenced={ref.cell_name for ref in old.references}
        if referenced!=original_names:
            raise ValueError('Historical wrapper references do not match the input cells')
        name=old.name
        existing={cell.name for cell in source_lib.cells}
        while name in existing:name+='X'
        wrapper=source_lib.new_cell(name)
        for cell in original_tops:wrapper.add(gdstk.Reference(cell))
        marker_layers={int(layers['electrode_marker_layer']),
                       int(layers['pad_contact_layer'])}
        metal_layer=int(layers['metal_layer'])
        for polygon in old.polygons:
            wrapper.add(gdstk.Polygon(polygon.points,layer=polygon.layer,
                                      datatype=polygon.datatype))
            if polygon.layer in marker_layers:
                wrapper.add(gdstk.Polygon(polygon.points,layer=metal_layer,
                                          datatype=polygon.datatype))
        source_lib.write_gds(str(output),max_points=4000)
        Path(new_output_gds).parent.mkdir(parents=True,exist_ok=True)
        Path(new_output_gds).write_bytes(output.read_bytes())
    return {'source_sha256':sha256(Path(input_gds).read_bytes()).hexdigest(),
            'historical_output_sha256':sha256(Path(old_output_gds).read_bytes()).hexdigest(),
            'rebuilt_output_sha256':sha256(Path(new_output_gds).read_bytes()).hexdigest(),
            'generated_polygon_count':len(old.polygons),
            'original_top_cell_count':len(original_tops),
            'rebuild_scope':'same generated wrapper polygons; source cells restored byte-equivalent in geometry; electrode and Pad marker polygons duplicated into each metal net'}
