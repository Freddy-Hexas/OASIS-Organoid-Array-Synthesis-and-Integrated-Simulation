"""Attachment-aware adaptive search over electrode-center configurations.

Cells cover a guarded construction domain. A failed point is never a proof
that its cell is infeasible: the cell is split, or kept as unresolved at the
requested construction resolution. Only strict whole-cell distance intervals
can exclude a cell because of a second local substrate branch. Full route and
joint-island checks remain the caller's acceptance oracle.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from heapq import heappop, heappush
import math

import numpy as np
from shapely.geometry import Point, box
from shapely.ops import nearest_points

from frontend import _polygonal
from island_router import _attachment, disk
from navigation_simplification import check_deadline


PLACEMENT_SEARCH_REVISION = "attachment_aware_adaptive_cells_v3_fair_coverage_enclosure"


@dataclass
class _Cell:
    geometry: object
    lower_radius: float
    closest: object
    depth: int


class AttachmentAwarePlacementSearch:
    """Fair best-bound subdivision; budget counts attachment-valid proposals.

    The coordinate system only determines the subdivision boxes. Geometry
    predicates, objective and acceptance rules are independent of filenames,
    known generator points, and arm labels. Rotation may change finite-budget
    sampling order, not the physical feasibility checks.
    """
    def __init__(self, domain, support, origin, *, island_radius_um,
                 minimum_substrate_gap_um, grid_um, resolution_um,
                 candidate_budget, maximum_radius_um=math.inf, deadline=None):
        if (grid_um <= 0 or resolution_um <= 0 or candidate_budget < 1 or
                island_radius_um <= 0 or minimum_substrate_gap_um < 0):
            raise ValueError("Invalid attachment placement search limits")
        self.support = support
        self.origin = np.asarray(origin, dtype=float)
        self.reference = Point(self.origin)
        self.radius = float(island_radius_um)
        self.gap = float(minimum_substrate_gap_um)
        self.grid = float(grid_um)
        self.resolution = float(resolution_um)
        self.candidate_budget = int(candidate_budget)
        # Computational effort, not a fabrication or exclusion threshold.
        self.cell_budget = 64 * self.candidate_budget
        self.maximum_radius = float(maximum_radius_um)
        self.deadline = deadline
        self.queue = []
        self.coverage_queue = []
        self.cells = {}
        self.serial = 0
        self.pending = None
        self.seen = set()
        self.visited = 0
        self.coverage_visits = 0
        self.proposed = 0
        self.attachment_checks = 0
        self.attachment_rejections = Counter()
        self.full_rejections = Counter()
        self.excluded = []
        self.excluded_count = 0
        self.excluded_area = 0.
        self.unresolved = []
        self.unresolved_count = 0
        self.unresolved_lower = math.inf
        self.status = "searching"
        self.best_radius = None
        self.radius_inner = float(disk((0, 0), self.radius).boundary.distance(Point(0, 0)))
        self.radius_outer = float(np.linalg.norm(
            np.asarray(disk((0, 0), self.radius).exterior.coords), axis=1).max())
        # The attachment oracle uses circumscribed polygons, whose radial
        # surplus grows with the requested process dimensions. Enclose the
        # actual polygon, rather than assuming a fixed native-grid guard
        # covers that surplus for arbitrarily large electrodes or gaps.
        self.attachment_window_outer = float(np.linalg.norm(np.asarray(
            disk((0, 0), self.radius + self.gap + .1).exterior.coords), axis=1).max())
        scale = max(1., self.radius, *map(abs, domain.bounds)) if not domain.is_empty else 1.
        self.distance_guard = max(2 * self.grid, 64 * math.ulp(scale))
        self.components = 0
        self.domain_area = 0.
        for piece in _polygonal(domain):
            check_deadline(deadline, "attachment position domain construction")
            # Match the old placement search's native-grid interior reserve.
            for interior in _polygonal(piece.buffer(-2 * self.grid, quad_segs=16)):
                self.components += 1
                self.domain_area += interior.area
                self._push(interior, 0)
        self.initial_lower = min((c.lower_radius for c in self.cells.values()), default=None)

    def _push(self, geometry, depth):
        if geometry.is_empty:
            return
        closest = nearest_points(self.reference, geometry)[1]
        # Reserve floating and native-grid uncertainty downward, not upward.
        lower = max(0., closest.distance(self.reference) - self.distance_guard)
        if lower >= self.maximum_radius:
            return
        self.serial += 1
        self.cells[self.serial] = _Cell(geometry, lower, closest, depth)
        heappush(self.queue, (lower, self.serial))
        heappush(self.coverage_queue, (depth, lower, self.serial))

    def _pop_cell(self):
        # One in four visits goes to the least-refined remaining region.
        # This effort policy is independent of process dimensions. It prevents
        # a low-radius rejected patch from consuming the finite search budget
        # before other connected portions of the same domain get inspected.
        coverage = self.visited % 4 == 3
        heap = self.coverage_queue if coverage else self.queue
        while heap:
            entry = heappop(heap)
            key = entry[-1]
            cell = self.cells.pop(key, None)
            if cell is None:
                continue
            if cell.lower_radius >= self.maximum_radius:
                continue
            self.coverage_visits += int(coverage)
            return cell
        return None

    def _remember_unresolved(self, cell, reason):
        self.unresolved_count += 1
        self.unresolved_lower = min(self.unresolved_lower, cell.lower_radius)
        if len(self.unresolved) < 24:
            self.unresolved.append({"bounds_um": list(cell.geometry.bounds),
                "radial_lower_bound_um": cell.lower_radius,
                "reason": reason, "depth": cell.depth})

    def _split(self, cell, reason):
        xmin, ymin, xmax, ymax = cell.geometry.bounds
        width, height = xmax - xmin, ymax - ymin
        if max(width, height) <= self.resolution:
            self._remember_unresolved(cell, "resolution_limit_after_" + reason)
            return
        if width >= height:
            mid = (xmin + xmax) / 2
            cutters = [box(xmin, ymin, mid, ymax), box(mid, ymin, xmax, ymax)]
        else:
            mid = (ymin + ymax) / 2
            cutters = [box(xmin, ymin, xmax, mid), box(xmin, mid, xmax, ymax)]
        children = [p for cutter in cutters for p in _polygonal(cell.geometry.intersection(cutter))]
        if len(children) < 2 or any(p.area >= cell.geometry.area for p in children):
            self._remember_unresolved(cell, "unsplittable_numerical_cell")
            return
        for child in children:
            self._push(child, cell.depth + 1)

    def _substrate_exclusion(self, cell):
        """Strict sufficient interval test; ambiguous contacts stay unresolved.

        For every center x in the cell, ||x-c|| <= h. All candidate disks and
        the attachment oracle's local neighborhoods fit in one larger window.
        A second connected substrate piece in that window cannot join the
        electrode's own piece through the candidate disk. It either stays
        strictly outside at a sub-limit gap for the entire cell, or has a
        positive-area intersection with every candidate disk. Tangency bands
        and topology predicates not covered by these tests are never excluded.
        """
        xmin, ymin, xmax, ymax = cell.geometry.bounds
        center = Point((xmin + xmax) / 2, (ymin + ymax) / 2)
        h = math.hypot(xmax - xmin, ymax - ymin) / 2
        # Large cells cannot satisfy a strict small-gap interval; split first.
        if h >= max(self.gap, self.radius) / 2:
            return None
        query_radius = self.attachment_window_outer + h + self.distance_guard
        window = disk(center.coords[0], query_radius)
        parts = _polygonal(self.support.intersection(window))
        # The entire center cell must lie in one local substrate piece. This
        # establishes the first contact for all its possible centers.
        own = [i for i, p in enumerate(parts) if p.covers(cell.geometry)]
        if len(own) != 1:
            return None
        for i, branch in enumerate(parts):
            if i == own[0]:
                continue
            d = center.distance(branch)
            low = max(0., d - h - self.distance_guard)
            high = d + h + self.distance_guard
            if low > self.radius_outer and high < self.radius_inner + self.gap - 1e-6:
                reason = "substrate_gap_for_all_cell_centers"
            elif high < self.radius_inner - self.distance_guard:
                reason = "two_disconnected_substrate_contacts_for_all_cell_centers"
            else:
                continue
            return {"reason": reason, "bounds_um": list(cell.geometry.bounds),
                "radial_lower_bound_um": cell.lower_radius, "area_um2": cell.geometry.area,
                "window_center_um": list(center.coords[0]), "cell_enclosing_radius_um": h,
                "window_radius_um": query_radius, "local_substrate_components": len(parts),
                "attachment_neighborhood_circumradius_um": self.attachment_window_outer,
                "branch_distance_interval_um": [low, high],
                "island_inradius_um": self.radius_inner, "island_circumradius_um": self.radius_outer,
                "minimum_substrate_gap_um": self.gap, "distance_guard_um": self.distance_guard,
                "predicate": "whole center cell covered by one local support piece; another window component satisfies strict distance intervals",
                "scope": "sufficient exclusion in floating polygon construction with numerical guard; not exact global manufacturing or routing infeasibility"}
        return None

    def next_proposal(self):
        if self.pending is not None:
            raise RuntimeError("Record a proposal result before requesting another")
        while self.cells and self.visited < self.cell_budget and self.proposed < self.candidate_budget:
            check_deadline(self.deadline, "adaptive attachment position search")
            cell = self._pop_cell()
            if cell is None:break
            self.visited += 1
            exclusion = self._substrate_exclusion(cell)
            if exclusion is not None:
                self.excluded_count += 1
                self.excluded_area += cell.geometry.area
                if len(self.excluded) < 24:
                    self.excluded.append(exclusion)
                continue
            xy = np.asarray(cell.closest.coords[0], dtype=float)
            key = tuple(int(round(v / self.grid)) for v in xy)
            if key in self.seen:
                self._split(cell, "repeated_point")
                continue
            self.seen.add(key)
            self.attachment_checks += 1
            island, reason = _attachment(self.support, xy, self.radius, self.gap)
            if island is None:
                self.attachment_rejections[reason] += 1
                self._split(cell, reason)
                continue
            self.proposed += 1
            self.pending = cell
            return xy, island
        self.status = ("cell_budget_reached" if self.cells and self.visited >= self.cell_budget else
            "candidate_budget_reached" if self.cells and self.proposed >= self.candidate_budget else
            "resolution_limited_cells_remain" if self.unresolved_count else
            "covered_construction_cells_processed")
        return None

    def record_result(self, accepted, *, radius_um=None, reason=None):
        cell, self.pending = self.pending, None
        if cell is None:
            raise RuntimeError("No pending attachment-valid proposal")
        if accepted:
            self.best_radius = float(radius_um)
            # This proposal is the cell's radial minimizer (up to the guard),
            # so no further point within this cell can improve it materially.
            self.maximum_radius = min(self.maximum_radius, self.best_radius - self.distance_guard)
        else:
            self.full_rejections[reason or "full_route_check"] += 1
            self._split(cell, reason or "full_route_check")

    def certificate(self, *, deadline_hit=False):
        pending = [c for c in self.cells.values() if c.lower_radius < self.maximum_radius]
        lower = min([cell.lower_radius for cell in pending] +
            ([self.pending.lower_radius] if self.pending is not None else []) +
            ([self.unresolved_lower] if math.isfinite(self.unresolved_lower) else []), default=None)
        if self.best_radius is not None:
            lower = min(lower, self.best_radius) if lower is not None else self.best_radius
        return {"revision": PLACEMENT_SEARCH_REVISION,
            "method": "attachment_oracle_guided_fair_best_bound_cell_subdivision",
            "feasible_components": self.components,
            "minimum_projected_radius_um": self.initial_lower,
            "construction_domain_area_um2": self.domain_area,
            "candidate_budget": self.candidate_budget, "cell_budget": self.cell_budget,
            "position_resolution_um": self.resolution, "native_grid_um": self.grid,
            "visited_cells": self.visited, "pending_cells": len(pending),
            "coverage_priority_visits": self.coverage_visits,
            "coverage_policy": "three radial-bound visits alternate with one least-refined-region visit",
            "attachment_checks": self.attachment_checks, "attachment_valid_proposals": self.proposed,
            "attachment_point_rejections": dict(self.attachment_rejections),
            "full_proposal_rejections": dict(self.full_rejections),
            "whole_cell_exclusion_count": self.excluded_count,
            "whole_cell_excluded_area_um2": self.excluded_area,
            "whole_cell_exclusion_examples": self.excluded,
            "unresolved_cell_count": self.unresolved_count,
            "unresolved_cell_examples": self.unresolved,
            "remaining_radial_lower_bound_um": lower,
            "best_verified_radius_um": self.best_radius,
            "construction_radial_gap_um": (max(0., self.best_radius - lower)
                if self.best_radius is not None and lower is not None else None),
            "deadline_hit": deadline_hit, "status": "deadline_reached" if deadline_hit else self.status,
            "global_center_optimality_proven": False,
            "scope": "guarded simplified residual construction domain with other nets fixed; unresolved cells and finite budgets are not infeasibility certificates; the radial gap is not a bound over the entire original continuous routing problem"}

