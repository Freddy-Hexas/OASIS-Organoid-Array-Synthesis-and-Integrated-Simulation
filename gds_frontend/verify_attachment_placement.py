"""Meaningful search regression: hidden legal bands and honest exclusions."""
import json
import math
import time

import numpy as np
from shapely.affinity import rotate, translate
from shapely.geometry import Point, box

from attachment_placement import AttachmentAwarePlacementSearch
from island_router import _attachment, disk
from navigation_simplification import SearchDeadlineExceeded


def test_hidden_legal_band():
    # A nearby branch makes the closest positions illegal. There is an
    # intermediate legal part of the same long, narrow electrode-center set.
    support = box(-5, -6, 180, 6).union(box(-5, 17, 50, 27))
    domain = box(0, -1, 170, 1)
    reference = np.array([0., 0.])
    records = []
    for angle, offset in [(0, (0, 0)), (37, (321, -172)), (90, (-47, 152))]:
        origin = np.asarray(offset, dtype=float)
        local_support = translate(rotate(support, angle, origin=(0, 0)), *offset)
        local_domain = translate(rotate(domain, angle, origin=(0, 0)), *offset)
        search = AttachmentAwarePlacementSearch(local_domain, local_support, origin,
            island_radius_um=16.05, minimum_substrate_gap_um=4., grid_um=.001,
            resolution_um=1.25, candidate_budget=32)
        best = None
        while True:
            proposal = search.next_proposal()
            if proposal is None:
                break
            xy, island = proposal
            assert _attachment(local_support, xy, 16.05, 4.)[0] is not None
            radius = math.dist(xy, origin)
            best = xy
            search.record_result(True, radius_um=radius)
        certificate = search.certificate()
        assert best is not None and 45 < math.dist(best, origin) < 90, certificate
        assert certificate['attachment_point_rejections'], certificate
        assert not certificate['global_center_optimality_proven']
        # Every recorded whole-cell exclusion really excludes all a fine
        # lattice of center-domain points from the original attachment oracle.
        for exclusion in certificate['whole_cell_exclusion_examples']:
            x0, y0, x1, y1 = exclusion['bounds_um']
            for x in np.linspace(x0, x1, 5):
                for y in np.linspace(y0, y1, 5):
                    if local_domain.covers(Point(x, y)):
                        assert _attachment(local_support, (x, y), 16.05, 4.)[0] is None, exclusion
        records.append({'rotation_deg': angle, 'radius_um': math.dist(best, origin),
            'visited_cells': certificate['visited_cells'],
            'whole_cell_exclusions': certificate['whole_cell_exclusion_count']})
    return records


def test_failed_point_does_not_discard_cell():
    support = box(-5, -8, 180, 8)
    domain = box(0, -2, 170, 2)
    search = AttachmentAwarePlacementSearch(domain, support, (0, 0),
        island_radius_um=16.05, minimum_substrate_gap_um=0, grid_um=.001,
        resolution_um=1.25, candidate_budget=32)
    accepted = None
    # Simulate a later route predicate that only permits the middle band.
    # Failed centers must result in subdivision, never whole-cell rejection.
    while True:
        proposal = search.next_proposal()
        if proposal is None:
            break
        xy, _ = proposal
        if xy[0] < 40:
            search.record_result(False, reason='route_oracle_fixture')
        else:
            accepted = xy
            search.record_result(True, radius_um=float(np.linalg.norm(xy)))
    certificate = search.certificate()
    assert accepted is not None and 40 <= accepted[0] < 90, certificate
    assert certificate['unresolved_cell_count'] > 0
    assert certificate['construction_radial_gap_um'] > 0
    assert certificate['whole_cell_exclusion_count'] == 0
    assert certificate['full_proposal_rejections'].get('route_oracle_fixture')
    # Explicit low budgets cannot masquerade as a no-solution proof.
    limited = AttachmentAwarePlacementSearch(domain, support, (0, 0),
        island_radius_um=16.05, minimum_substrate_gap_um=0, grid_um=.001,
        resolution_um=1.25, candidate_budget=1)
    assert limited.next_proposal() is not None
    limited.record_result(False, reason='route_oracle_fixture')
    assert limited.next_proposal() is None
    cert = limited.certificate()
    assert cert['status'] == 'candidate_budget_reached' and cert['pending_cells'] > 0
    assert not cert['global_center_optimality_proven']
    expired = AttachmentAwarePlacementSearch(domain, support, (0, 0),
        island_radius_um=16.05, minimum_substrate_gap_um=0, grid_um=.001,
        resolution_um=1.25, candidate_budget=32)
    expired.deadline = time.perf_counter() - 1
    try:
        expired.next_proposal()
    except SearchDeadlineExceeded:
        assert expired.certificate(deadline_hit=True)['deadline_hit']
    else:
        raise AssertionError('Expected deadline')


def test_scaled_attachment_windows_are_enclosed():
    # Large user-specified electrodes magnify the circumradius of the oracle's
    # polygonal neighborhood. Each exclusion window must still contain that
    # entire neighborhood at every tested cell corner, including translations.
    for radius in [16.05, 1605., 16050.]:
        for offset in [(0., 0.), (300000., -700000.)]:
            support = translate(box(-20, -6, 20, 6).union(
                box(-20, radius + 2, 20, radius + 12)), *offset)
            domain = translate(box(-.25, -.25, .25, .25), *offset)
            search = AttachmentAwarePlacementSearch(domain, support, offset,
                island_radius_um=radius, minimum_substrate_gap_um=4., grid_um=.001,
                resolution_um=1.25, candidate_budget=8)
            assert search.next_proposal() is None
            exclusions = search.certificate()['whole_cell_exclusion_examples']
            assert exclusions, (radius, offset)
            for item in exclusions:
                window = disk(item['window_center_um'], item['window_radius_um'])
                x0, y0, x1, y1 = item['bounds_um']
                for x in np.linspace(x0, x1, 3):
                    for y in np.linspace(y0, y1, 3):
                        assert window.covers(disk((x, y), radius + 4.1)), (radius, offset)
                        assert _attachment(support, (x, y), radius, 4.)[0] is None


if __name__ == '__main__':
    records = test_hidden_legal_band()
    test_failed_point_does_not_discard_cell()
    test_scaled_attachment_windows_are_enclosed()
    print(json.dumps({'passed': True, 'hidden_legal_band_rigid_transforms': records,
        'whole_cell_exclusions_checked': True, 'full_rejection_subdivision': True,
        'budget_and_resolution_are_not_infeasibility': True,
        'scaled_attachment_window_enclosure': True}, ensure_ascii=False))

