"""Certified inner navigation meshes; original polygons remain the DRC oracle.

A simplification is a construction accelerator, never a capacity bound. It is
accepted only if it is contained in the original domain, each original hole
has its own surviving hole, and every component and protected point survives.
"""
from __future__ import annotations

import time
import shapely
from shapely.geometry import Point, Polygon
from shapely.ops import unary_union
from shapely.strtree import STRtree


class SearchDeadlineExceeded(RuntimeError):
    pass


def check_deadline(deadline, stage):
    if deadline is not None and time.perf_counter() >= deadline:
        raise SearchDeadlineExceeded(stage)


def _parts(domain):
    if domain.geom_type == 'Polygon':
        return [domain] if not domain.is_empty else []
    return [p for p in getattr(domain, 'geoms', ()) if p.geom_type == 'Polygon' and not p.is_empty]


def simplify_inside(domain, tolerance, grid_um=.001, *, protected_points=(), deadline=None):
    check_deadline(deadline, 'navigation simplification')
    before = int(shapely.get_num_coordinates(domain))
    records = []
    result = []
    for original in _parts(domain):
        check_deadline(deadline, 'navigation component simplification')
        accepted = original
        selected = 0.0
        trial_tolerance = float(tolerance)
        holes = [Polygon(r) for r in original.interiors]
        while trial_tolerance >= 2 * grid_um:
            # The erosion reserves the simplifier's displacement. Containment
            # is independently checked instead of assuming a library contract.
            simple = original.simplify(trial_tolerance, preserve_topology=True)
            proposal = simple.buffer(-trial_tolerance - 2 * grid_um, join_style='mitre')
            targets = [Point(p) for p in protected_points if original.covers(Point(p))]
            if targets:
                patches = [original.intersection(p.buffer(3 * trial_tolerance + 4 * grid_um)) for p in targets]
                proposal = unary_union([proposal, *patches])
            valid = (proposal.geom_type == 'Polygon' and not proposal.is_empty and proposal.is_valid
                     and len(proposal.interiors) == len(holes) and original.covers(proposal)
                     and all(proposal.covers(p) for p in targets))
            if valid and holes:
                new_holes = [Polygon(r) for r in proposal.interiors]
                tree = STRtree(new_holes)
                matched = set()
                for hole in holes:
                    indices = [int(i) for i in tree.query(hole) if new_holes[int(i)].covers(hole)]
                    if len(indices) != 1 or indices[0] in matched:
                        valid = False
                        break
                    matched.add(indices[0])
            check_deadline(deadline, 'navigation topology and containment certificate')
            if valid:
                # A purported acceleration that adds more vertices is useless.
                if shapely.get_num_coordinates(proposal) < shapely.get_num_coordinates(original):
                    accepted = proposal
                    selected = trial_tolerance
                break
            trial_tolerance /= 2
        result.append(accepted)
        records.append({'tolerance_um': selected, 'holes_bijectively_preserved': True,
                        'coordinates_before': int(shapely.get_num_coordinates(original)),
                        'coordinates_after': int(shapely.get_num_coordinates(accepted))})
    navigation = unary_union(result) if result else domain
    return navigation, {'method': 'verified_inner_simplification_with_hole_correspondence',
                        'coordinates_before': before,
                        'coordinates_after': int(shapely.get_num_coordinates(navigation)),
                        'components': records, 'original_covers_navigation': True,
                        'scope': 'conservative construction mesh; no completeness or capacity claim'}
