"""A hard disk constraint on electrode centers, never on wires or exits.

The reference is the original support's minimum enclosing circle center,
rounded to its native grid (the same convention as the outer-port policy).
Polygon clipping proposes sources conservatively; analytic distance and the
exported marker's exact half-grid center independently enforce the rule.
"""
from dataclasses import dataclass
from fractions import Fraction
from functools import cached_property
import math

import numpy as np
import shapely
from shapely.geometry import Point, Polygon, LineString
from shapely.ops import substring


def radius_from_rules(rules):
    # Missing fields in saved reports mean the old run was unconstrained.
    return (rules.get('electrode_region_radius_um') if isinstance(rules, dict)
            else getattr(rules, 'electrode_region_radius_um', None))


def geometric_origin(support, grid_um=.001):
    center = shapely.minimum_bounding_circle(support).centroid.coords[0]
    return tuple(round(v / grid_um) * grid_um for v in center)


@dataclass(frozen=True)
class ElectrodeRegion:
    origin: tuple
    radius_um: float | None
    grid_um: float = .001

    def __post_init__(self):
        if self.radius_um is not None and (
                isinstance(self.radius_um, bool) or not isinstance(self.radius_um, (int, float))
                or not math.isfinite(self.radius_um) or self.radius_um <= 0):
            raise ValueError('Electrode region radius must be finite and positive')

    def contains(self, xy, *, export_safe=False):
        if self.radius_um is None:
            return True
        radius = self.radius_um - (2 * self.grid_um if export_safe else 0)
        return radius >= 0 and math.hypot(xy[0]-self.origin[0], xy[1]-self.origin[1]) <= radius

    @cached_property
    def construction_disk(self):
        if self.radius_um is None:
            raise ValueError('Unrestricted region has no finite disk')
        radius = self.radius_um - 2*self.grid_um
        if radius <= 0:
            return Polygon()
        # Inscribed disk with at most a quarter-grid chord loss. This is a
        # construction restriction only, not an upper-bound domain.
        ratio = min(.25, self.grid_um / (4*radius))
        angle = math.acos(max(-1, 1-ratio))
        quad = max(16, math.ceil(math.pi/(4*angle))) if angle else 4096
        return Point(self.origin).buffer(radius, quad_segs=quad)

    def clip(self, domain):
        if domain.is_empty or self.radius_um is None:
            return domain
        x0,y0,x1,y1 = domain.bounds
        if all(self.contains(p, export_safe=True) for p in ((x0,y0),(x1,y0),(x1,y1),(x0,y1))):
            return domain
        return domain.intersection(self.construction_disk)

    def record(self):
        return {'enabled': self.radius_um is not None, 'radius_um': self.radius_um,
                'diameter_um': 2*self.radius_um if self.radius_um is not None else None,
                'reference_um': list(self.origin), 'native_grid_um': self.grid_um,
                'reference_definition': 'original support minimum enclosing circle center, rounded to native GDS grid',
                'constrained_object': 'electrode_center',
                'wire_and_pad_may_leave_disk': True,
                'construction_rounding_reserve_um': 2*self.grid_um,
                'scope': 'hard center-position constraint; electrode rims and routes are not clipped to this disk'}

    def audit(self, routes):
        radii = [math.dist(route['source_um'], self.origin) for route in routes]
        failed = [i+1 for i,r in enumerate(radii) if self.radius_um is not None and r > self.radius_um]
        return {**self.record(), 'checked_electrodes': len(radii),
                'maximum_selected_radius_um': max(radii, default=None),
                'outside_networks': failed, 'passed': not failed}

    def packing_upper_bound(self, minimum_spacing):
        if self.radius_um is None or minimum_spacing <= 0:
            return None
        radius, spacing = Fraction(str(self.radius_um)), Fraction(str(minimum_spacing))
        # Disjoint open d/2 disks lie in a disk R+d/2; pi cancels exactly.
        quotient = ((2*radius+spacing)/spacing)**2
        value = quotient.numerator // quotient.denominator
        if spacing > 2*radius:
            value = 1
        return {'value': value, 'radius_um': self.radius_um,
                'minimum_center_spacing_um': minimum_spacing,
                'area_ratio_numerator': quotient.numerator, 'area_ratio_denominator': quotient.denominator,
                'uses_exact_rational_arithmetic': True,
                'proof': 'disjoint d/2 center disks inside radius R+d/2; N <= floor((1+2R/d)^2); d>2R implies at most one',
                'scope': 'all centers in the declared analytic disk, independent of sampled candidates and route search'}


def region_for(support, rules, grid_um=.001, *, origin=None):
    return ElectrodeRegion(tuple(origin) if origin is not None else geometric_origin(support,grid_um),
                           radius_from_rules(rules), grid_um)


def region_for_navigation(nav, rules):
    grid = nav.summary.get('gds_native_precision_m',1e-9)*1e6
    origin = nav.summary.get('outer_exit_policy',{}).get('reference_um') if nav.summary.get('outer_exit_policy') else None
    return region_for(nav.support,rules,grid,origin=origin)


def add_region_candidates(graph, candidates, region):
    """Add a source in every in-disk edge portion without clipping graph paths.

    Existing corridor samples may all miss a small placement disk. Splitting
    an edge at an interior source preserves all of its out-of-disk geometry
    and connectivity to distant terminals. This is still a finite search.
    """
    if region.radius_um is None:
        return candidates, 0
    added = []; next_id = max(graph.nodes, default=-1)+1; candidate_set=set(candidates)
    for u,v,data in list(graph.edges(data=True)):
        points = np.asarray(data['points_um'])
        if np.linalg.norm(points[-1]-graph.nodes[u]['xy_um']) < np.linalg.norm(points[0]-graph.nodes[u]['xy_um']):
            points = points[::-1]
        line = LineString(points)
        intersection = region.clip(line)
        pieces = [intersection] if intersection.geom_type=='LineString' else getattr(intersection,'geoms',())
        cuts = []
        for piece in pieces:
            if piece.geom_type!='LineString' or piece.length <= 4*region.grid_um:
                continue
            # Keep existing legal sources at the portion's endpoints.
            if any(n in candidate_set and region.contains(graph.nodes[n]['xy_um'],export_safe=True)
                   and piece.distance(Point(graph.nodes[n]['xy_um'])) < region.grid_um for n in (u,v)):
                continue
            d = line.project(piece.interpolate(.5,normalized=True))
            if 2*region.grid_um < d < line.length-2*region.grid_um:
                cuts.append(d)
        cuts = sorted(set(cuts))
        if not cuts:
            continue
        graph.remove_edge(u,v);chain=[u]
        for d in cuts:
            xy=list(line.interpolate(d).coords[0]);node=next_id;next_id+=1
            graph.add_node(node,xy_um=xy,kind='candidate')
            chain.append(node);added.append(node)
        chain.append(v)
        for a,b,lo,hi in zip(chain,chain[1:],[0,*cuts],[*cuts,line.length]):
            segment=substring(line,lo,hi)
            graph.add_edge(a,b,**{**data,'points_um':np.asarray(segment.coords),'weight':segment.length})
    return [n for n in [*candidates,*added] if region.contains(graph.nodes[n]['xy_um'],export_safe=True)], len(added)


def exported_center_in_disk(center_twice_ticks, origin_ticks, radius_um, grid_um):
    """Exact rational squared-distance test, including half-grid marker centers."""
    dx = int(center_twice_ticks[0])-2*int(origin_ticks[0])
    dy = int(center_twice_ticks[1])-2*int(origin_ticks[1])
    return (dx*dx+dy*dy)*Fraction(str(grid_um))**2 <= 4*Fraction(str(radius_um))**2
