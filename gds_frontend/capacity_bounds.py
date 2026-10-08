"""Candidate-independent center-packing upper bound on a GDS support layer.

The proof uses the original (flattened) GDS polygons, not the navigation
inset, graph, chosen routes, or Pad frame.  When transformed vertices do not
lie on the native GDS grid, the exact-grid certificate is unavailable and the
caller must use another conservative bound.
"""
from __future__ import annotations

from collections import defaultdict
from decimal import Decimal, ROUND_FLOOR
from hashlib import sha256
from pathlib import Path
from tempfile import TemporaryDirectory
import math

import gdstk
import numpy as np


def _orientation(a, b, c):
    return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])


def _on_segment(a, b, p):
    return (_orientation(a, b, p) == 0 and
            min(a[0], b[0]) <= p[0] <= max(a[0], b[0]) and
            min(a[1], b[1]) <= p[1] <= max(a[1], b[1]))


def _segments_intersect(a, b, c, d):
    ab_c = _orientation(a, b, c)
    ab_d = _orientation(a, b, d)
    cd_a = _orientation(c, d, a)
    cd_b = _orientation(c, d, b)
    if (ab_c == 0 and _on_segment(a, b, c)) or (ab_d == 0 and _on_segment(a, b, d)):
        return True
    if (cd_a == 0 and _on_segment(c, d, a)) or (cd_b == 0 and _on_segment(c, d, b)):
        return True
    return (ab_c > 0) != (ab_d > 0) and (cd_a > 0) != (cd_b > 0)


def _simple_integer_polygon(vertices):
    """Exact integer validity screen before using a polygon as a filled set."""
    n=len(vertices)
    if n<3 or len(set(vertices))!=n or _signed_double_area(vertices)==0:
        return False
    for i,b in enumerate(vertices):
        a,c=vertices[i-1],vertices[(i+1)%n]
        if (_orientation(a,b,c)==0 and
                (b[0]-a[0])*(c[0]-b[0])+
                (b[1]-a[1])*(c[1]-b[1])<0):
            return False
    segments=[]
    for i,(a,b) in enumerate(zip(vertices,vertices[1:]+vertices[:1])):
        if a==b:return False
        segments.append((min(a[0],b[0]),max(a[0],b[0]),
                         min(a[1],b[1]),max(a[1],b[1]),i,a,b))
    active=[]
    for current in sorted(segments):
        x0,x1,y0,y1,i,a,b=current
        active=[item for item in active if item[1]>=x0]
        for other in active:
            if other[3]<y0 or y1<other[2]:continue
            j=other[4]
            if (i-j)%n in (1,n-1):continue
            if _segments_intersect(a,b,other[5],other[6]):
                return False
        active.append(current)
    return True


def _point_in_polygon_or_boundary(p, vertices):
    winding = 0
    for a, b in zip(vertices, vertices[1:] + vertices[:1]):
        if _on_segment(a, b, p):
            return True
        turn = _orientation(a, b, p)
        if a[1] <= p[1] < b[1] and turn > 0:
            winding += 1
        elif b[1] <= p[1] < a[1] and turn < 0:
            winding -= 1
    return winding != 0


def _polygon_intersects_closed_cell(vertices, x0, y0, hx, hy=None):
    """Exact integer test for a simple closed polygon and a closed cell."""
    hy = hx if hy is None else hy
    x1, y1 = x0+hx, y0+hy
    if any(x0 <= x <= x1 and y0 <= y <= y1 for x, y in vertices):
        return True
    cell = ((x0, y0), (x1, y0), (x1, y1), (x0, y1))
    for a, b in zip(vertices, vertices[1:] + vertices[:1]):
        if (max(a[0], b[0]) < x0 or min(a[0], b[0]) > x1 or
                max(a[1], b[1]) < y0 or min(a[1], b[1]) > y1):
            continue
        if any(_segments_intersect(a, b, c, d)
               for c, d in zip(cell, cell[1:] + cell[:1])):
            return True
    # If boundaries do not meet and no polygon vertex is inside, the cell
    # can only intersect the polygon by lying completely inside it.
    return _point_in_polygon_or_boundary(cell[0], vertices)


def _grid_vertices(points, grid_um):
    scaled = np.asarray(points, dtype=float) / grid_um
    rounded = np.rint(scaled)
    if (not np.isfinite(scaled).all() or
            np.max(np.abs(scaled)) >= 2**52 or
            np.max(np.abs(scaled-rounded)) > 1e-3):
        raise ValueError('Flattened GDS vertices are not on the native grid')
    vertices = [tuple(map(int, xy)) for xy in rounded]
    if len(vertices) > 1 and vertices[0] == vertices[-1]:
        vertices.pop()
    return vertices


def _floor_grid_ticks(value_um, grid_um):
    # Shrink a positive physical distance by a tiny relative allowance for
    # binary representation of the GDS grid and rule values. This can loosen
    # a bound by a tick, but cannot falsely strengthen one.
    ratio = Decimal(str(value_um))/Decimal(str(grid_um))
    ratio -= Decimal('1e-8')*max(Decimal(1), abs(ratio))
    return int(ratio.to_integral_value(rounding=ROUND_FLOOR))


def _validated_native_grid_um(lib):
    grid_um = lib.precision*1e6
    if not math.isfinite(grid_um) or grid_um <= 0:
        raise ValueError('GDS native grid is invalid')
    return grid_um


def _ceil_isqrt(value):
    root = math.isqrt(value)
    return root if root*root == value else root+1


def _integer_convex_hull(points):
    points = sorted(set(points))
    if len(points) < 3:
        return points
    lower = []
    for point in points:
        while len(lower) >= 2 and _orientation(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    upper = []
    for point in reversed(points):
        while len(upper) >= 2 and _orientation(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return lower[:-1]+upper[:-1]


def _convex_polygon_intersects_box(vertices, x0, y0, x1, y1):
    if any(x0 <= x <= x1 and y0 <= y <= y1 for x, y in vertices):
        return True
    box_vertices = ((x0, y0), (x1, y0), (x1, y1), (x0, y1))
    if any(_point_in_polygon_or_boundary(point, vertices) for point in box_vertices):
        return True
    return any(_segments_intersect(a, b, c, d)
               for a, b in zip(vertices, vertices[1:]+vertices[:1])
               for c, d in zip(box_vertices, box_vertices[1:]+box_vertices[:1]))


def verify_pads_outside_hull_cut(certificate, pad_bounds_um):
    """Certify each rectangular Pad bounding box misses the guarded cut.

    Coordinates are rounded outward by two native grid ticks. A two-tick
    guard around the hull supports matching precision for the input support.
    """
    grid_um = certificate['native_grid_um']
    guard = certificate['square_guard_grid_ticks']
    vertices = [tuple(point) for point in certificate['convex_hull_vertices_grid_ticks']]
    if len(vertices) < 3:
        return False
    for xmin, ymin, xmax, ymax in pad_bounds_um:
        x0 = math.floor(xmin/grid_um)-guard-2
        y0 = math.floor(ymin/grid_um)-guard-2
        x1 = math.ceil(xmax/grid_um)+guard+2
        y1 = math.ceil(ymax/grid_um)+guard+2
        if _convex_polygon_intersects_box(vertices, x0, y0, x1, y1):
            return False
    return True


def verify_pads_outside_radial_cut(certificate, pad_bounds_um):
    """Exact integer bounding-box exclusion from the certified circle."""
    grid_um = certificate['native_grid_um']
    cx2, cy2 = certificate['circle_center_twice_grid_ticks']
    radius2 = certificate['radius_twice_grid_ticks']
    for xmin, ymin, xmax, ymax in pad_bounds_um:
        x0, y0 = math.floor(2*xmin/grid_um)-4, math.floor(2*ymin/grid_um)-4
        x1, y1 = math.ceil(2*xmax/grid_um)+4, math.ceil(2*ymax/grid_um)+4
        dx = max(x0-cx2, cx2-x1, 0)
        dy = max(y0-cy2, cy2-y1, 0)
        if dx*dx+dy*dy <= radius2*radius2:
            return False
    return True


def gds_convex_hull_cut_upper_bound(gds_path, layer, datatype,
                                    wire_width_um, spacing_um):
    """Exact-integer perimeter bound for a guarded convex-hull cut.

    Pad exclusion is a separate applicability check. The support is inside
    its raw-vertex convex hull; adding a two-grid-tick square puts every
    accepted readback vertex strictly inside the cut. The perimeter of a
    convex Minkowski sum is the sum of the two perimeters.
    """
    path = Path(gds_path)
    raw = path.read_bytes()
    with TemporaryDirectory(prefix='capacity_hull_') as temporary:
        local = Path(temporary) / 'source.gds'
        local.write_bytes(raw)
        lib = gdstk.read_gds(str(local), unit=1e-6)
    grid_um = _validated_native_grid_um(lib)
    separation_ticks = _floor_grid_ticks(wire_width_um+spacing_um, grid_um)
    if wire_width_um <= 0 or spacing_um < 0 or separation_ticks < 1:
        raise ValueError('Invalid wire width or spacing')
    vertices = []
    for top in lib.top_level():
        for item in top.get_polygons():
            if (item.layer, item.datatype) == (layer, datatype):
                vertices.extend(_grid_vertices(item.points, grid_um))
    if not vertices:
        raise ValueError(f'GDS has no support polygons on {layer}/{datatype}')
    hull = _integer_convex_hull(vertices)
    if len(hull) < 3:
        raise ValueError('Support convex hull has zero area')
    edge_length_upper = sum(_ceil_isqrt((b[0]-a[0])**2+(b[1]-a[1])**2)
                            for a, b in zip(hull, hull[1:]+hull[:1]))
    guard_ticks = 2
    perimeter_upper_ticks = edge_length_upper+8*guard_ticks
    upper = max(1, perimeter_upper_ticks//separation_ticks)
    xs, ys = [p[0] for p in vertices], [p[1] for p in vertices]
    cx2, cy2 = min(xs)+max(xs), min(ys)+max(ys)
    radius_2ticks = _ceil_isqrt(max((2*x-cx2)**2+(2*y-cy2)**2
                                    for x, y in vertices)) + 8
    return {
        'value': upper,
        'method': 'guarded_integer_convex_hull_perimeter_cut',
        'input_sha256': sha256(raw).hexdigest(),
        'support_layer': [int(layer), int(datatype)],
        'native_grid_um': grid_um,
        'convex_hull_vertex_count': len(hull),
        'convex_hull_vertices_grid_ticks': [list(point) for point in hull],
        'edge_perimeter_upper_grid_ticks': edge_length_upper,
        'square_guard_grid_ticks': guard_ticks,
        'cut_perimeter_upper_grid_ticks': perimeter_upper_ticks,
        'wire_plus_spacing_lower_grid_ticks': separation_ticks,
        'pad_exclusion_circle_center_um': [cx2*grid_um/2, cy2*grid_um/2],
        'pad_exclusion_circle_radius_um': radius_2ticks*grid_um/2,
        'pad_outside_cut_verified': False,
        'applicability_condition': 'every counted electrode center lies in original support, each counted Pad lies outside the guarded convex-hull cut, and each net has an at-least-wire-width centerline path with the specified metal spacing',
        'proof': 'one crossing per net on the guarded convex boundary; consecutive crossing points need at least wire_width+spacing of boundary arc; exact integer upper rounding bounds the hull perimeter',
        'exact_integer_predicates': True,
    }


def gds_radial_cut_upper_bound(gds_path, layer, datatype, wire_width_um,
                               spacing_um):
    """Conditional global net bound for Pads outside the returned circle.

    Each connected, at-least-wire-width route from an interior electrode to
    an exterior Pad has a centerline crossing the circle. Pairwise metal gap
    makes crossing centers at least wire_width+spacing apart. Consecutive
    crossings therefore consume at least that much circle arc length.
    """
    path = Path(gds_path)
    raw = path.read_bytes()
    with TemporaryDirectory(prefix='capacity_cut_') as temporary:
        local = Path(temporary) / 'source.gds'
        local.write_bytes(raw)
        lib = gdstk.read_gds(str(local), unit=1e-6)
    grid_um = _validated_native_grid_um(lib)
    if wire_width_um <= 0 or spacing_um < 0:
        raise ValueError('Wire width must be positive and spacing nonnegative')
    separation_ticks = _floor_grid_ticks(wire_width_um+spacing_um, grid_um)
    if separation_ticks < 1:
        raise ValueError('Wire plus spacing is below the native grid')
    vertices = []
    polygons = 0
    for top in lib.top_level():
        for item in top.get_polygons():
            if (item.layer, item.datatype) == (layer, datatype):
                vertices.extend(_grid_vertices(item.points, grid_um))
                polygons += 1
    if not vertices:
        raise ValueError(f'GDS has no support polygons on {layer}/{datatype}')
    xs, ys = [p[0] for p in vertices], [p[1] for p in vertices]
    cx2, cy2 = min(xs)+max(xs), min(ys)+max(ys)
    radius_2ticks = _ceil_isqrt(max((2*x-cx2)**2+(2*y-cy2)**2
                                    for x, y in vertices)) + 4
    # pi < 22/7. The numerator is an exact upper bound on the perimeter
    # expressed in native-grid ticks, and separation_ticks is a lower bound.
    # A single crossing has no distinct cyclic neighbor; the arc argument
    # applies to N >= 2, while N <= 1 is always conservatively allowed.
    upper = max(1, (22*radius_2ticks)//(7*separation_ticks))
    return {
        'value': upper,
        'method': 'original_support_enclosing_circle_cut',
        'input_sha256': sha256(raw).hexdigest(),
        'support_layer': [int(layer), int(datatype)],
        'native_grid_um': grid_um,
        'circle_center_um': [cx2*grid_um/2, cy2*grid_um/2],
        'circle_center_twice_grid_ticks': [cx2, cy2],
        'circle_radius_um': radius_2ticks*grid_um/2,
        'radius_twice_grid_ticks': radius_2ticks,
        'wire_plus_spacing_lower_grid_ticks': separation_ticks,
        'pi_upper_bound': '22/7',
        'raw_support_polygons': polygons,
        'applicability_condition': 'each counted net has one at-least-wire-width metal path from an electrode center inside this circle to its Pad entirely outside this circle; different nets keep the specified spacing',
        'pad_outside_cut_verified': False,
        'proof': 'choose one centerline crossing per net; radius-wire_width/2 disks at the crossings lie in separate metals, so Euclidean and hence consecutive cyclic arc separation is at least wire_width+spacing; N times this separation does not exceed the circle perimeter',
        'exact_integer_predicates': True,
    }


def gds_enclosing_disk_angular_upper_bound(gds_path, layer, datatype,
                                           center_distance_um):
    """Certify a sharp small-cardinality center bound from one enclosing disk.

    If k points lie in a radius-R disk, two have angular gap at most 2*pi/k.
    Their separation is at most max(R,2*R*sin(pi/k)). Exact integer squared
    comparisons handle k=2,3,4; for k=5,6,7, sin(x)<=x and pi<22/7 give
    a rational upper bound. Returns None when none applies.
    """
    path = Path(gds_path)
    raw = path.read_bytes()
    with TemporaryDirectory(prefix='capacity_angular_') as temporary:
        local = Path(temporary) / 'source.gds'
        local.write_bytes(raw)
        lib = gdstk.read_gds(str(local), unit=1e-6)
    grid_um = _validated_native_grid_um(lib)
    if not math.isfinite(center_distance_um) or center_distance_um <= 0:
        raise ValueError('Center distance must be positive and finite')
    distance_ticks = _floor_grid_ticks(center_distance_um, grid_um)
    if distance_ticks < 1:
        return None
    vertices = []
    for top in lib.top_level():
        for item in top.get_polygons():
            if (item.layer, item.datatype) == (layer, datatype):
                vertices.extend(_grid_vertices(item.points, grid_um))
    if not vertices:
        raise ValueError(f'GDS has no support polygons on {layer}/{datatype}')
    xs, ys = [p[0] for p in vertices], [p[1] for p in vertices]
    cx2, cy2 = min(xs)+max(xs), min(ys)+max(ys)
    radius2 = _ceil_isqrt(max((2*x-cx2)**2+(2*y-cy2)**2
                              for x,y in vertices)) + 4
    # radius2 is twice the radius in native ticks. Hence d > 2R exactly
    # means d_ticks > radius2. The next two tests square 2d versus sqrt(3)R
    # and sqrt(2)R without losing strictness at equality.
    if distance_ticks > radius2:
        value, predicate = 1, 'd_grid > 2R_upper'
    elif 4*distance_ticks*distance_ticks > 3*radius2*radius2:
        value, predicate = 2, '4d_grid^2 > 3(2R_upper)^2'
    elif 4*distance_ticks*distance_ticks > 2*radius2*radius2:
        value, predicate = 3, '4d_grid^2 > 2(2R_upper)^2'
    else:
        value=predicate=None
        for excluded in (5,6,7):
            if (2*distance_ticks>radius2 and
                    14*excluded*distance_ticks>44*radius2):
                value=excluded-1
                predicate=(f'2d_grid > (2R_upper) and '
                           f'14*{excluded}*d_grid > 44*(2R_upper)')
                break
        if value is None:return None
    return {
        'value': value,
        'method': 'native_integer_enclosing_disk_angular_pigeonhole',
        'input_sha256': sha256(raw).hexdigest(),
        'support_layer': [int(layer),int(datatype)],
        'native_grid_um': grid_um,
        'required_center_distance_um': center_distance_um,
        'conservative_distance_grid_ticks': distance_ticks,
        'disk_center_twice_grid_ticks': [cx2,cy2],
        'disk_radius_twice_grid_ticks_upper': radius2,
        'strict_predicate': predicate,
        'pi_upper_bound':'22/7' if value>=4 else None,
        'original_support_vertex_count': len(vertices),
        'exact_integer_predicates': True,
        'scope': 'all electrode centers in the original GDS support, regardless of routes, tracks, bridges and Pads',
        'proof': 'every support polygon lies in the convex enclosing disk; among k points two subtend angle at most 2*pi/k, so their distance is at most max(R,2R sin(pi/k)); for k>=5 use sin(x)<=x and pi<22/7; the strict integer predicate rules out k=value+1'
    }


def verify_gds_angular_certificate(gds_path,certificate):
    """Re-read the original grid and check every premise of the small bound."""
    if certificate is None:return False
    path=Path(gds_path)
    raw=path.read_bytes()
    if sha256(raw).hexdigest()!=certificate.get('input_sha256'):
        return False
    with TemporaryDirectory(prefix='verify_angular_') as temporary:
        local=Path(temporary)/'source.gds'
        local.write_bytes(raw)
        lib=gdstk.read_gds(str(local),unit=1e-6)
    grid=_validated_native_grid_um(lib)
    if grid!=certificate.get('native_grid_um'):
        return False
    layer,datatype=certificate['support_layer']
    vertices=[vertex for top in lib.top_level()
              for polygon in top.get_polygons()
              if (polygon.layer,polygon.datatype)==(layer,datatype)
              for vertex in _grid_vertices(polygon.points,grid)]
    if not vertices or len(vertices)!=certificate.get('original_support_vertex_count'):
        return False
    cx2,cy2=certificate['disk_center_twice_grid_ticks']
    radius2=certificate['disk_radius_twice_grid_ticks_upper']
    if (not all(isinstance(value,int) for value in (cx2,cy2,radius2))
            or radius2<=0):
        return False
    if any((2*x-cx2)**2+(2*y-cy2)**2>radius2*radius2
           for x,y in vertices):
        return False
    d=_floor_grid_ticks(certificate['required_center_distance_um'],grid)
    if d!=certificate.get('conservative_distance_grid_ticks'):
        return False
    value=certificate.get('value')
    return ((value==1 and d>radius2) or
            (value==2 and 4*d*d>3*radius2*radius2) or
            (value==3 and 4*d*d>2*radius2*radius2) or
            (value in (4,5,6) and 2*d>radius2 and
             14*(value+1)*d>44*radius2 and
             certificate.get('pi_upper_bound')=='22/7'))


def _expanded_cell_union_area(occupied, side, shifts, radius, side_y=None):
    """Exact area (grid ticks squared) of the union of expanded cell boxes."""
    side_y = side if side_y is None else side_y
    by_row = defaultdict(list)
    for ix, iy in occupied:
        by_row[iy].append(ix)
    events = []
    for iy, columns in by_row.items():
        columns.sort()
        first = last = columns[0]
        for ix in columns[1:] + [None]:
            if ix is not None and ix == last+1:
                last = ix
                continue
            x0 = shifts[0]+first*side-radius
            x1 = shifts[0]+(last+1)*side+radius
            y0 = shifts[1]+iy*side_y-radius
            y1 = shifts[1]+(iy+1)*side_y+radius
            events.append((x0, 1, y0, y1))
            events.append((x1, -1, y0, y1))
            if ix is not None:
                first = last = ix
    events.sort(key=lambda event: (event[0], -event[1]))
    active = defaultdict(int)
    area = 0
    previous_x = events[0][0]
    index = 0
    while index < len(events):
        x = events[index][0]
        if x > previous_x and active:
            intervals = sorted(active)
            lo, hi = intervals[0]
            covered = 0
            for a, b in intervals[1:]:
                if a > hi:
                    covered += hi-lo
                    lo, hi = a, b
                else:
                    hi = max(hi, b)
            covered += hi-lo
            area += (x-previous_x)*covered
        while index < len(events) and events[index][0] == x:
            _, sign, y0, y1 = events[index]
            key = (y0, y1)
            active[key] += sign
            if active[key] == 0:
                del active[key]
            index += 1
        previous_x = x
    return area


def gds_cell_cover_upper_bound(gds_path, layer, datatype, center_distance_um,
                               *, offset_fraction=(0, 0), cell_aspect=(1, 1),
                               _return_occupied=False):
    """Bound all centers in the selected original support layer.

    Every occupied closed cell has diameter strictly less than the required
    center distance, so each cell contains at most one center.  Closed-cell
    intersection overcounts boundary cells and is therefore safe.  Invalid
    input polygons use their entire bounding-box cell set conservatively.
    """
    path = Path(gds_path)
    raw = path.read_bytes()
    with TemporaryDirectory(prefix='capacity_gds_') as temporary:
        local = Path(temporary) / 'source.gds'
        local.write_bytes(raw)
        lib = gdstk.read_gds(str(local), unit=1e-6)
    grid_um = _validated_native_grid_um(lib)
    if not math.isfinite(center_distance_um) or center_distance_um <= 0:
        raise ValueError('Center distance must be positive and finite')
    distance_ticks = _floor_grid_ticks(center_distance_um, grid_um)
    if distance_ticks < 3:
        raise ValueError('Center distance is too small for a grid-cell bound')
    ax, ay = cell_aspect
    if not isinstance(ax, int) or not isinstance(ay, int) or ax < 1 or ay < 1:
        raise ValueError('Cell aspect must be a pair of positive integers')
    # (ax*k)^2+(ay*k)^2 < d_grid^2, even at the farthest cell corners.
    k = math.isqrt((distance_ticks*distance_ticks-1)//(ax*ax+ay*ay))
    hx, hy = ax*k, ay*k
    if hx < 1 or hy < 1 or hx*hx+hy*hy >= distance_ticks*distance_ticks:
        raise RuntimeError('Cell diameter is not strictly below center distance')
    if len(offset_fraction) != 2 or any(value not in (0, .5) for value in offset_fraction):
        raise ValueError('Only zero and half-cell shifts are supported')
    shifts = tuple(0 if value == 0 else side//2
                   for value, side in zip(offset_fraction, (hx, hy)))
    # A flattened coordinate can differ slightly from its integer GDS grid
    # representation because of floating point readback. Expand every tested
    # cell by two grid ticks, much more than the admitted 0.001 tick error.
    guard_ticks = 2
    occupied = set()
    polygons = 0
    invalid_bbox_polygons = 0
    for top in lib.top_level():
        for item in top.get_polygons():
            if (item.layer, item.datatype) != (layer, datatype):
                continue
            vertices = _grid_vertices(item.points, grid_um)
            if len(set(vertices)) < 3:
                continue
            polygons += 1
            xs = [p[0] for p in vertices]
            ys = [p[1] for p in vertices]
            ix0, ix1 = ((min(xs)-guard_ticks-shifts[0])//hx,
                        (max(xs)+guard_ticks-shifts[0])//hx)
            iy0, iy1 = ((min(ys)-guard_ticks-shifts[1])//hy,
                        (max(ys)+guard_ticks-shifts[1])//hy)
            # The exact fill of a self-intersecting GDS boundary is ambiguous
            # across readers.  Its complete bounding box is a safe cover.
            simple = _simple_integer_polygon(vertices)
            if not simple:
                invalid_bbox_polygons += 1
            for ix in range(ix0, ix1+1):
                for iy in range(iy0, iy1+1):
                    if (ix, iy) in occupied:
                        continue
                    if not simple or _polygon_intersects_closed_cell(
                            vertices, shifts[0]+ix*hx-guard_ticks,
                            shifts[1]+iy*hy-guard_ticks,
                            hx+2*guard_ticks, hy+2*guard_ticks):
                        occupied.add((ix, iy))
    if polygons == 0:
        raise ValueError(f'GDS has no support polygons on {layer}/{datatype}')
    # Disks of radius d_grid/2 have disjoint interiors. Their union is inside
    # the union of occupied cells expanded by ceil(d_grid/2) in L-infinity.
    # The latter's area is obtained exactly from integer rectangles. Use a
    # rational lower bound on pi when dividing, so the quotient rounds up.
    expanded_radius = (distance_ticks+1)//2
    expanded_area = _expanded_cell_union_area(occupied, hx, shifts,
                                               expanded_radius, hy)
    # Archimedes' classical rational enclosure gives 223/71 < pi.
    pi_lower_numerator = 223
    pi_denominator = 71
    area_upper = (expanded_area*4*pi_denominator)//(
        pi_lower_numerator*distance_ticks*distance_ticks)
    upper = min(len(occupied), area_upper)
    digest = sha256()
    for ix, iy in sorted(occupied):
        digest.update(f'{ix},{iy};'.encode('ascii'))
    certificate = {
        'value': upper,
        'method': 'exact_integer_square_cover_and_expanded_area',
        'input_sha256': sha256(raw).hexdigest(),
        'support_layer': [int(layer), int(datatype)],
        'native_grid_um': grid_um,
        'required_center_distance_um': center_distance_um,
        'conservative_distance_grid_ticks': distance_ticks,
        'cell_aspect': list(cell_aspect),
        'cell_side_grid_ticks': [hx, hy],
        'cell_origin_shift_grid_ticks': list(shifts),
        'readback_guard_grid_ticks': guard_ticks,
        'cell_count': len(occupied),
        'expanded_cell_union_area_grid_ticks2': expanded_area,
        'expanded_radius_grid_ticks': expanded_radius,
        'area_packing_upper_bound': area_upper,
        'pi_lower_bound': '223/71',
        'occupied_cell_index_sha256': digest.hexdigest(),
        'raw_support_polygons': polygons,
        'invalid_polygons_covered_by_bbox': invalid_bbox_polygons,
        'exact_integer_predicates': True,
        'scope': 'electrode centers lie in union of selected original GDS support polygons',
        'proof': 'support is covered by intersecting closed cells, each with diameter below d; alternatively disjoint radius-d_grid/2 disks lie in the exact integer rectangle union obtained by expanding those cells',
    }
    return (certificate, occupied) if _return_occupied else certificate


def _signed_double_area(vertices):
    return sum(a[0]*b[1]-b[0]*a[1]
               for a, b in zip(vertices, vertices[1:]+vertices[:1]))


def _integer_rectangle_union_area(rectangles):
    """Exact union area of closed axis-aligned integer rectangles."""
    events = []
    ys = set()
    for x0, y0, x1, y1 in rectangles:
        if x1 <= x0 or y1 <= y0:
            continue
        events.extend(((x0, 1, y0, y1), (x1, -1, y0, y1)))
        ys.update((y0, y1))
    if not events:
        return 0
    ys = sorted(ys)
    y_index = {y: i for i, y in enumerate(ys)}
    covered = [0]*(4*len(ys))
    lengths = [0]*(4*len(ys))
    def update(node, start, stop, lo, hi, change):
        if lo <= start and stop <= hi:
            covered[node] += change
        else:
            midpoint = (start+stop)//2
            if lo < midpoint:
                update(2*node, start, midpoint, lo, hi, change)
            if hi > midpoint:
                update(2*node+1, midpoint, stop, lo, hi, change)
        if covered[node]:
            lengths[node] = ys[stop]-ys[start]
        elif stop-start == 1:
            lengths[node] = 0
        else:
            lengths[node] = lengths[2*node]+lengths[2*node+1]
    events.sort(key=lambda event: (event[0], -event[1]))
    area = 0
    previous_x = events[0][0]
    for x, change, y0, y1 in events:
        area += (x-previous_x)*lengths[1]
        update(1, 0, len(ys)-1, y_index[y0], y_index[y1], change)
        previous_x = x
    return area


def gds_square_packing_upper_bound(gds_path, layer, datatype,
                                   center_distance_um):
    """Candidate-independent safe packing area on the original GDS.

    A center square of half-side r is disjoint from every other center square
    when 2*sqrt(2)*r < the required Euclidean separation. A polygon's square
    dilation is covered by the original polygon plus the union of each edge's
    bounding rectangle enlarged by r. Both areas are bounded using integers;
    overlap between the polygon and edge cover is deliberately counted twice.
    """
    path = Path(gds_path)
    raw = path.read_bytes()
    with TemporaryDirectory(prefix='capacity_square_') as temporary:
        local = Path(temporary) / 'source.gds'
        local.write_bytes(raw)
        lib = gdstk.read_gds(str(local), unit=1e-6)
    grid_um = _validated_native_grid_um(lib)
    if not math.isfinite(center_distance_um) or center_distance_um <= 0:
        raise ValueError('Center distance must be positive and finite')
    distance_ticks = _floor_grid_ticks(center_distance_um, grid_um)
    radius_ticks = math.isqrt((distance_ticks*distance_ticks-1)//8)
    if radius_ticks < 1:
        raise ValueError('Center distance is too small for a grid-square bound')
    # Flattened reference transforms can be within 0.001 native tick of the
    # rounded coordinate accepted by _grid_vertices. Two ticks cover that
    # readback discrepancy on either coordinate and keep this an outer bound.
    readback_guard_ticks = 2
    edge_radius_ticks = radius_ticks+readback_guard_ticks
    edge_rectangles = []
    source_area_twice = 0
    raw_count = 0
    invalid_bbox_count = 0
    for top in lib.top_level():
        for item in top.get_polygons():
            if (item.layer, item.datatype) != (layer, datatype):
                continue
            vertices = _grid_vertices(item.points, grid_um)
            if len(set(vertices)) < 3:
                continue
            raw_count += 1
            # A self-intersecting GDS polygon has reader-dependent fill. Its
            # full bounding box is a deliberately loose but safe superset.
            if not _simple_integer_polygon(vertices):
                invalid_bbox_count += 1
                xs, ys = [p[0] for p in vertices], [p[1] for p in vertices]
                vertices = [(min(xs), min(ys)), (max(xs), min(ys)),
                            (max(xs), max(ys)), (min(xs), max(ys))]
            source_area_twice += abs(_signed_double_area(vertices))
            for a, b in zip(vertices, vertices[1:]+vertices[:1]):
                edge_rectangles.append((min(a[0], b[0])-edge_radius_ticks,
                                        min(a[1], b[1])-edge_radius_ticks,
                                        max(a[0], b[0])+edge_radius_ticks,
                                        max(a[1], b[1])+edge_radius_ticks))
    if raw_count == 0:
        raise ValueError(f'GDS has no support polygons on {layer}/{datatype}')
    rectangle_union_area = _integer_rectangle_union_area(edge_rectangles)
    area_twice = source_area_twice+2*rectangle_union_area
    square_area_twice = 8*radius_ticks*radius_ticks
    if area_twice <= 0:
        raise RuntimeError('Nonempty GDS support yielded zero Minkowski area')
    return {
        'value': area_twice//square_area_twice,
        'method': 'native_integer_edge_rectangle_square_packing',
        'input_sha256': sha256(raw).hexdigest(),
        'support_layer': [int(layer), int(datatype)],
        'native_grid_um': grid_um,
        'required_center_distance_um': center_distance_um,
        'conservative_distance_grid_ticks': distance_ticks,
        'square_half_side_grid_ticks': radius_ticks,
        'readback_guard_grid_ticks': readback_guard_ticks,
        'square_diameter_strictly_below_distance':
            8*radius_ticks*radius_ticks < distance_ticks*distance_ticks,
        'original_polygon_area_sum_twice_grid_ticks2': source_area_twice,
        'edge_rectangle_union_area_grid_ticks2': rectangle_union_area,
        'dilation_area_upper_twice_grid_ticks2': area_twice,
        'one_square_area_twice_grid_ticks2': square_area_twice,
        'raw_support_polygons': raw_count,
        'invalid_polygons_covered_by_bbox': invalid_bbox_count,
        'exact_integer_predicates': True,
        'scope': 'all electrode centers in the original selected GDS support, independent of route candidates, tracks, portals, bridges and Pads',
        'proof': 'strictly separated centers have disjoint-interior axis-aligned squares; every point in support dilated by one square lies either in a rounded original polygon or within square radius plus readback guard of one polygon edge, hence inside that edge bounding rectangle; summing rounded original polygon areas and the exact integer edge-rectangle union area is an upper area bound',
    }


def gds_radial_partition_upper_bound(gds_path, layer, datatype,
                                     center_distance_um, wire_width_um,
                                     spacing_um, pad_bounds_um):
    """Certified circular cut plus packing bound outside the circle.

    Any electrode center inside the selected closed circle must route across
    it to a Pad outside. Centers outside are covered by native-grid cells not
    entirely inside the open circle. The bound includes bridges crossing the
    circle outside the original support; it never assumes support-only arcs.
    """
    if not pad_bounds_um:
        raise ValueError('Actual or proposed Pad bounds are required')
    cut = gds_radial_cut_upper_bound(gds_path, layer, datatype,
                                     wire_width_um, spacing_um)
    cover, occupied = gds_cell_cover_upper_bound(
        gds_path, layer, datatype, center_distance_um,
        _return_occupied=True)
    if cover['input_sha256'] != cut['input_sha256']:
        raise RuntimeError('GDS changed while constructing radial partition')
    grid_um = cut['native_grid_um']
    cx2, cy2 = cut['circle_center_twice_grid_ticks']
    gap_ticks = cut['wire_plus_spacing_lower_grid_ticks']
    hx, hy = cover['cell_side_grid_ticks']
    shift_x, shift_y = cover['cell_origin_shift_grid_ticks']
    # Exactly the outward Pad-box rounding used by the independent verifier.
    max_radius2 = None
    for xmin, ymin, xmax, ymax in pad_bounds_um:
        x0, y0 = math.floor(2*xmin/grid_um)-4, math.floor(2*ymin/grid_um)-4
        x1, y1 = math.ceil(2*xmax/grid_um)+4, math.ceil(2*ymax/grid_um)+4
        dx = max(x0-cx2, cx2-x1, 0)
        dy = max(y0-cy2, cy2-y1, 0)
        distance_squared = dx*dx+dy*dy
        if distance_squared <= 1:
            raise ValueError('A Pad overlaps the origin of every admissible circle')
        allowed = math.isqrt(distance_squared-1)
        max_radius2 = allowed if max_radius2 is None else min(max_radius2, allowed)
    events = defaultdict(int)
    for ix, iy in occupied:
        x0 = 2*(shift_x+ix*hx)-cx2
        x1 = x0+2*hx
        y0 = 2*(shift_y+iy*hy)-cy2
        y1 = y0+2*hy
        far2 = max(x0*x0, x1*x1)+max(y0*y0, y1*y1)
        threshold = math.isqrt(far2)+1
        if threshold <= max_radius2:
            events[threshold] += 1
    outside = len(occupied)
    def circle_capacity(radius2):
        # pi < 22/7; radius2 is in half native-grid ticks.
        return max(1, (22*radius2)//(7*gap_ticks))
    best = (outside+circle_capacity(1), 1, outside)
    for radius2, removed in sorted(events.items()):
        outside -= removed
        candidate = (outside+circle_capacity(radius2), radius2, outside)
        if candidate < best:
            best = candidate
    value, radius2, outside = best
    certificate = {
        'value': value,
        'method': 'native_grid_radial_inside_cut_plus_outside_cell_cover',
        'input_sha256': cut['input_sha256'],
        'support_layer': [int(layer), int(datatype)],
        'native_grid_um': grid_um,
        'circle_center_twice_grid_ticks': [cx2, cy2],
        'circle_center_um': cut['circle_center_um'],
        'radius_twice_grid_ticks': radius2,
        'circle_radius_um': radius2*grid_um/2,
        'maximum_pad_excluding_radius_twice_grid_ticks': max_radius2,
        'wire_plus_spacing_lower_grid_ticks': gap_ticks,
        'required_center_distance_um': center_distance_um,
        'center_cell_diameter_strictly_below_required_distance': True,
        'source_cover_cell_count': len(occupied),
        'source_cover_cell_index_sha256': cover['occupied_cell_index_sha256'],
        'outside_possible_cell_count': outside,
        'inside_circle_routing_cut_upper_bound': circle_capacity(radius2),
        'candidate_radii_checked': len(events)+1,
        'pad_outside_cut_verified': False,
        'exact_integer_predicates': True,
        'proof': 'inside-circle nets cross the full circle and consume at least wire width plus spacing; every outside-circle center lies in an occupied closed cell not wholly inside the open circle, and each cell has diameter below required center separation',
    }
    if not verify_pads_outside_radial_cut(certificate, pad_bounds_um):
        raise RuntimeError('Selected radial partition does not exclude every Pad')
    certificate['pad_outside_cut_verified'] = True
    return certificate
