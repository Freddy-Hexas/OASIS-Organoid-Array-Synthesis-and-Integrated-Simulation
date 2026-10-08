"""Process-sized outer ports, derived only from the original GDS geometry.

An exterior face is not an electrical exit. Pad bridges start in a
fabrication-sized band of the original convex envelope, and their whole
footprint avoids the interior core. Legacy circular windows remain readable.
This is an explicit design constraint, not a theorem that every excluded
concave boundary point is physically impossible to connect.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property,lru_cache
import math

import numpy as np
import shapely
from shapely.geometry import LineString, Point, Polygon


LEGACY_EXIT_POLICY_REVISION = 'global_radial_outer_ports_v1'
EXIT_POLICY_REVISION = 'directional_convex_envelope_outer_ports_v2'
OUTER_BRIDGE_CONTACT_POLICY = (
    'added bridges stay outside the original-support interior core and contact one contractible global outer-extremity window; functional wires leave original support only at those windows')


def geometry_min_radius(geometry, origin):
    """Distance to a point includes segment interiors, not just vertices."""
    return float(geometry.distance(Point(origin))) if not geometry.is_empty else math.inf


@dataclass(frozen=True)
class OuterExitPolicy:
    grid_um: float
    center_ticks: tuple[int, int]
    radius_ticks: int
    maximum_vertex_radius_squared_ticks: int
    wire_width_um: float
    margin_um: float
    numeric_guard_um: float
    bridge_width_um: float

    @property
    def origin(self):
        return np.asarray(self.center_ticks, dtype=float) * self.grid_um

    @property
    def center_inset_um(self):
        return self.wire_width_um / 2 + self.margin_um + self.numeric_guard_um

    @property
    def launch_depth_ticks(self):
        # Two native ticks cover port/line rounding at export.  This is a
        # numerical allowance, not a structure-dependent radial percentage.
        return math.ceil(self.center_inset_um / self.grid_um - 1e-9) + 2

    @property
    def contact_depth_ticks(self):
        return self.launch_depth_ticks + math.ceil(
            self.bridge_width_um / (2 * self.grid_um) - 1e-9) + 2

    @property
    def launch_radius_um(self):
        return max(0, self.radius_ticks - self.launch_depth_ticks) * self.grid_um

    @property
    def core_radius_um(self):
        return max(0, self.radius_ticks - self.contact_depth_ticks) * self.grid_um

    @property
    def portal_search_radius_um(self):
        # A center-domain portal is an inset away from the real source edge.
        # This is only a superset for searching; launch and footprint tests
        # below are mandatory before the shortest-path search is performed.
        return max(0.0, self.launch_radius_um - self.center_inset_um - 4*self.grid_um)

    def allows_exit(self, point):
        return float(np.linalg.norm(np.asarray(point) - self.origin)) >= self.launch_radius_um

    def check_bridge(self, source, bridge, *, exit_point=None):
        if exit_point is not None and not self.allows_exit(exit_point):
            return False, 'exit_not_in_outer_extremity_window'
        # Testing the entire bridge also excludes a curved bridge that dives
        # into an empty internal wedge without intersecting the source.
        if geometry_min_radius(bridge, self.origin) < self.core_radius_um + 2*self.grid_um:
            return False, 'bridge_footprint_enters_interior_core'
        contact = bridge.intersection(source)
        if contact.geom_type != 'Polygon' or contact.is_empty or len(contact.interiors):
            return False, 'bridge_source_contact_not_one_contractible_window'
        if contact.area <= self.grid_um**2:
            return False, 'bridge_has_no_source_contact'
        return True, None

    def check_exterior_line(self, source, line):
        outside = line.difference(source)
        if geometry_min_radius(outside, self.origin) < self.launch_radius_um + 2*self.grid_um:
            return False, 'wire_leaves_source_before_outer_window'
        return True, None

    def record(self):
        return {'revision': LEGACY_EXIT_POLICY_REVISION,
                'native_grid_um': self.grid_um,
                'center_grid_ticks': list(self.center_ticks),
                'reference_um': self.origin.tolist(),
                'reference_definition': 'minimum enclosing circle center, rounded to the native GDS grid',
                'enclosing_radius_grid_ticks': self.radius_ticks,
                'maximum_vertex_radius_squared_ticks': self.maximum_vertex_radius_squared_ticks,
                'source_max_radius_um': self.radius_ticks*self.grid_um,
                'wire_width_um': self.wire_width_um, 'margin_um': self.margin_um,
                'numeric_guard_um': self.numeric_guard_um,
                'bridge_width_um': self.bridge_width_um,
                'launch_depth_grid_ticks': self.launch_depth_ticks,
                'contact_depth_grid_ticks': self.contact_depth_ticks,
                'minimum_launch_radius_um': self.launch_radius_um,
                'minimum_bridge_footprint_radius_um': self.core_radius_um,
                'window_derivation': 'launch depth = wire half-width + edge margin + numerical guard + 2 GDS ticks; contact depth adds bridge half-width + 2 GDS ticks',
                'scope': 'mandatory outer-extremity design constraint; no per-component relaxation, filename dispatch, or internal bridge fallback'}

    def edge_intervals(self, a, b, *, search=False, extra_guard_um=0):
        radius = self.portal_search_radius_um if search else self.launch_radius_um + extra_guard_um
        return outside_circle_edge_intervals(a, b, self.origin, radius)


@dataclass(frozen=True)
class ConvexEnvelopeExitPolicy(OuterExitPolicy):
    hull_ticks: tuple[tuple[int, int], ...]

    @cached_property
    def hull(self):
        return Polygon(np.asarray(self.hull_ticks, dtype=float) * self.grid_um)

    @lru_cache(maxsize=32)
    def _core(self, depth):
        return self.hull.buffer(-depth, join_style='mitre')

    @cached_property
    def launch_core(self):
        return self._core(self.launch_depth_ticks * self.grid_um)

    @cached_property
    def bridge_core(self):
        return self._core(self.contact_depth_ticks * self.grid_um)

    def allows_exit(self, point):
        return not self.launch_core.buffer(2 * self.grid_um).covers(Point(point))

    def edge_intervals(self, a, b, *, search=False, extra_guard_um=0):
        depth = self.launch_depth_ticks * self.grid_um - extra_guard_um
        if search:
            depth += self.center_inset_um + 4 * self.grid_um
        clipped = LineString([a, b]).difference(self._core(max(0, depth)))
        pieces = [clipped] if clipped.geom_type == 'LineString' else getattr(clipped, 'geoms', ())
        return [(np.asarray(p.coords[0]), np.asarray(p.coords[-1])) for p in pieces
                if p.geom_type == 'LineString' and p.length > 1e-9]

    def check_bridge(self, source, bridge, *, exit_point=None):
        if exit_point is not None and not self.allows_exit(exit_point):
            return False, 'exit_not_in_outer_extremity_window'
        if bridge.intersects(self.bridge_core.buffer(2 * self.grid_um)):
            return False, 'bridge_footprint_enters_interior_core'
        contact = bridge.intersection(source)
        if contact.geom_type != 'Polygon' or contact.is_empty or len(contact.interiors):
            return False, 'bridge_source_contact_not_one_contractible_window'
        if contact.area <= self.grid_um ** 2:
            return False, 'bridge_has_no_source_contact'
        return True, None

    def check_exterior_line(self, source, line):
        outside = line.difference(source)
        if outside.intersects(self.launch_core.buffer(2 * self.grid_um)):
            return False, 'wire_leaves_source_before_outer_window'
        return True, None

    def record(self):
        return {**super().record(), 'revision': EXIT_POLICY_REVISION,
                'envelope_vertices_grid_ticks': [list(p) for p in self.hull_ticks],
                'minimum_launch_radius_um': None, 'minimum_bridge_footprint_radius_um': None,
                'maximum_launch_depth_from_envelope_um': self.launch_depth_ticks * self.grid_um,
                'maximum_bridge_depth_from_envelope_um': self.contact_depth_ticks * self.grid_um,
                'window_derivation': 'original-GDS convex envelope; one outermost support line per direction; fabrication-sized inward band; enclosed holes and internal concavities excluded',
                'scope': 'uniform directional outer-envelope constraint; bridges cannot enter the eroded hull interior; full source contact and wire checks remain mandatory; no filename dispatch'}


def make_outer_exit_policy(source, *, wire_width_um=5.0, margin_um=4.0,
                           numeric_guard_um=.05, bridge_width_um=40.0,
                           grid_um=.001):
    values = (wire_width_um, margin_um, numeric_guard_um, bridge_width_um, grid_um)
    if (not all(math.isfinite(v) for v in values) or wire_width_um <= 0 or
            margin_um < 0 or numeric_guard_um < 0 or bridge_width_um <= 0 or grid_um <= 0):
        raise ValueError('Invalid outer-port fabrication dimensions')
    if source.is_empty:
        raise ValueError('Empty source cannot define outer ports')
    circle = shapely.minimum_bounding_circle(source)
    center = tuple(int(round(v/grid_um)) for v in circle.centroid.coords[0])
    # Hull vertices are native GDS vertices; internal overlay intersections
    # cannot change the enclosing radius.  Compute squared distances with
    # Python integers and round the enclosing radius outward by at most a tick.
    vertices = np.rint(shapely.get_coordinates(source.convex_hull)/grid_um).astype(np.int64)
    maximum = max((int(x)-center[0])**2 + (int(y)-center[1])**2 for x,y in vertices)
    radius = math.isqrt(maximum)
    if radius*radius < maximum:
        radius += 1
    hull = tuple(map(tuple, vertices[:-1].tolist()))
    return ConvexEnvelopeExitPolicy(grid_um, center, radius, maximum, wire_width_um,
                                   margin_um, numeric_guard_um, bridge_width_um, hull)


def policy_from_record(record):
    if record.get('revision') not in (EXIT_POLICY_REVISION, LEGACY_EXIT_POLICY_REVISION):
        raise ValueError('Unknown outer exit policy revision')
    cls = ConvexEnvelopeExitPolicy if record['revision'] == EXIT_POLICY_REVISION else OuterExitPolicy
    extra = (tuple(map(tuple, record['envelope_vertices_grid_ticks'])),) if cls is ConvexEnvelopeExitPolicy else ()
    return cls(record['native_grid_um'], tuple(record['center_grid_ticks']),
                           record['enclosing_radius_grid_ticks'],
                           record['maximum_vertex_radius_squared_ticks'],
                           record['wire_width_um'], record['margin_um'],
                           record['numeric_guard_um'], record['bridge_width_um'], *extra)


def outside_circle_edge_intervals(a, b, origin, radius):
    """Exact quadratic formula for candidate edge portions outside a disk.

    Float roots construct search proposals only.  Accepted/exported geometry
    is checked separately, with exact native-grid distance predicates.
    """
    a = np.asarray(a, dtype=float); b = np.asarray(b, dtype=float)
    d = b-a; p = a-np.asarray(origin)
    aa = float(d@d)
    if aa <= 1e-24:
        return []
    bb = 2*float(p@d); cc = float(p@p)-radius*radius
    discriminant = bb*bb - 4*aa*cc
    cuts = [0.0, 1.0]
    if discriminant > 0:
        root = math.sqrt(discriminant)
        cuts.extend(t for t in ((-bb-root)/(2*aa), (-bb+root)/(2*aa)) if 0<t<1)
    cuts.sort()
    return [(a+t0*d, a+t1*d) for t0,t1 in zip(cuts[:-1], cuts[1:])
            if np.linalg.norm(p+(t0+t1)/2*d) >= radius]
