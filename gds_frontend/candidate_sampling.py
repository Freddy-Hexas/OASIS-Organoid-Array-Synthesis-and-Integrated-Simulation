"""Finite electrode candidates sampled from every recovered support corridor.

The samples are graph search hints, not a claim that the continuous feasible
electrode region has been exhausted. Every eventual disk and metal route still
needs the vector and exported-GDS checks in demo_router.
"""
from __future__ import annotations

import math

import numpy as np


def _pieces_at_distances(points, distances):
    """Split an ordered polyline at true arc lengths, retaining its vertices."""
    vertices = np.asarray(points, dtype=float)
    lengths = np.linalg.norm(np.diff(vertices, axis=0), axis=1)
    cumulative = np.r_[0.0, np.cumsum(lengths)]
    total = float(cumulative[-1])
    if total <= 0:
        return [], []
    cuts = [0.0, *distances, total]
    pieces = []
    centers = []

    def point_at(d):
        if d <= 0:
            return vertices[0]
        if d >= total:
            return vertices[-1]
        index = int(np.searchsorted(cumulative, d, side='right') - 1)
        index = min(index, len(lengths) - 1)
        alpha = (d - cumulative[index]) / lengths[index] if lengths[index] > 0 else 0.0
        return vertices[index] * (1 - alpha) + vertices[index + 1] * alpha

    centers = [point_at(d).tolist() for d in distances]
    for a, b in zip(cuts[:-1], cuts[1:]):
        middle = vertices[(cumulative > a + 1e-8) & (cumulative < b - 1e-8)]
        part = np.vstack((point_at(a), middle, point_at(b)))
        if np.linalg.norm(np.diff(part, axis=0), axis=1).sum() <= 1e-8:
            raise RuntimeError('Degenerate corridor split')
        pieces.append(part)
    return pieces, centers


def add_corridor_candidates(graph, *, electrode_diameter_um, spacing_um,
                            target_step_um=60.0, max_per_corridor=12):
    """Add one or more source nodes along every edge, including trees and loops.

    The minimum arc-length separation is based on disk diameter plus spacing.
    Sampling density is capped so huge source pools do not obscure the current
    finite-candidate nature of the solver. A candidate on a narrow strip is
    later rejected by the exact disk-plus-margin containment test.
    """
    if target_step_um <= 0 or max_per_corridor < 1:
        raise ValueError('Invalid candidate sampling density')
    out = graph.copy()
    next_node = max(out.nodes, default=-1) + 1
    added = 0
    by_kind = {'tree_or_open': 0, 'cycle_or_loop': 0}
    min_separation = electrode_diameter_um + max(0.0, spacing_um)
    for u, v, key, data in list(out.edges(keys=True, data=True)):
        points = np.asarray(data['points_um'], dtype=float)
        source = np.asarray(out.nodes[u]['xy_um'], dtype=float)
        if np.linalg.norm(points[-1] - source) < np.linalg.norm(points[0] - source):
            points = points[::-1]
        length = float(np.linalg.norm(np.diff(points, axis=0), axis=1).sum())
        if length <= 1e-8:
            continue
        # Keep candidate centers away from abstract graph nodes. These node
        # positions may represent junction windows, not a particular lane.
        end_guard = min(length / 4, max(min_separation / 2, electrode_diameter_um / 2))
        usable = length - 2 * end_guard
        count = min(max_per_corridor, max(1, int(math.floor(usable / target_step_um)) + 1))
        if count > 1 and usable / (count + 1) < min_separation:
            count = max(1, int(math.floor(usable / min_separation)) - 1)
        distances = np.linspace(end_guard, length - end_guard, count + 2)[1:-1]
        pieces, centers = _pieces_at_distances(points, distances.tolist())
        out.remove_edge(u, v, key)
        chain = [u]
        for center in centers:
            mid = next_node; next_node += 1
            out.add_node(mid, xy_um=center, kind='corridor_candidate',
                         collector_ids=[], window_radius_um=data['min_raster_clearance_um'])
            chain.append(mid)
        chain.append(v)
        for a, b, segment in zip(chain[:-1], chain[1:], pieces):
            attrs = dict(data)
            attrs['points_um'] = segment
            attrs['length_um'] = float(np.linalg.norm(np.diff(segment, axis=0), axis=1).sum())
            out.add_edge(a, b, **attrs)
        added += len(centers)
        kind = 'cycle_or_loop' if u == v else 'tree_or_open'
        by_kind[kind] += len(centers)
    return out, {'candidates_added': added, 'by_edge_type': by_kind,
                 'target_step_um': target_step_um,
                 'max_per_corridor': max_per_corridor,
                 'minimum_nominal_separation_um': min_separation}
