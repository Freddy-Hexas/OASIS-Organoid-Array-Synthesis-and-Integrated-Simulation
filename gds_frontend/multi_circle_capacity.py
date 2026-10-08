"""Candidate-independent multiple-circle routing bound with rational dual.

Every counted electrode center lies in exactly one of the conservative closed
support-cover cells (choose an arbitrary owning cell at boundaries). A cell has
diameter below the required center separation, so its occupancy is at most one.
Every center in a cell wholly inside one of the certified Pad-excluding circles
needs its route to cross that circle. The integer primal relaxation is bounded
by any nonnegative dual multipliers. A floating LP only proposes multipliers;
the returned capacity is recalculated with integer arithmetic.
"""
from __future__ import annotations

from hashlib import sha256
from math import floor, isqrt

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix

from capacity_bounds import (gds_cell_cover_upper_bound,
                             gds_radial_cut_upper_bound,
                             verify_pads_outside_radial_cut)


def _circle_members(ordered_cells, hx, hy, shift_x, shift_y,
                    cx2, cy2, radius2):
    inside = []
    limit2 = radius2*radius2
    for index, (ix, iy) in enumerate(ordered_cells):
        x0 = 2*(shift_x+ix*hx)-cx2
        x1 = x0+2*hx
        y0 = 2*(shift_y+iy*hy)-cy2
        y1 = y0+2*hy
        far2 = max(x0*x0, x1*x1)+max(y0*y0, y1*y1)
        if far2 < limit2:
            inside.append(index)
    return inside


def _member_digest(members):
    digest = sha256()
    for index in members:
        digest.update(f'{index};'.encode('ascii'))
    return digest.hexdigest()


def _circle_cut(cx2, cy2, radius2, gap_ticks, grid_um, pad_bounds_um,
                members):
    cut = {'circle_center_twice_grid_ticks': [cx2, cy2],
           'radius_twice_grid_ticks': radius2,
           'native_grid_um': grid_um}
    if not verify_pads_outside_radial_cut(cut, pad_bounds_um):
        return None
    capacity = max(1, (22*radius2)//(7*gap_ticks))
    if not members or len(members) <= capacity:
        return None
    return {'center_twice_grid_ticks': [cx2, cy2],
            'radius_twice_grid_ticks': radius2,
            'crossing_capacity': capacity,
            'fully_inside_cell_count': len(members),
            'fully_inside_cell_index_sha256': _member_digest(members)}, members


def _proposed_cuts(cover, ordered_cells, gap_ticks, pad_bounds_um,
                   *, centers_per_axis=5, radius_quantiles=(.2,.4,.6,.8,1.0)):
    hx, hy = cover['cell_side_grid_ticks']
    shift_x, shift_y = cover['cell_origin_shift_grid_ticks']
    grid_um = cover['native_grid_um']
    ix_values = [p[0] for p in ordered_cells]
    iy_values = [p[1] for p in ordered_cells]
    def axis_centers(values, side, shift):
        low = min(values)*side+shift
        high = (max(values)+1)*side+shift
        if centers_per_axis == 1:
            return [(low+high)]
        return sorted({2*low + (2*(high-low)*index)//(centers_per_axis-1)
                       for index in range(centers_per_axis)})
    center_xs = axis_centers(ix_values, hx, shift_x)
    center_ys = axis_centers(iy_values, hy, shift_y)
    proposals = {}
    for cx2 in center_xs:
        for cy2 in center_ys:
            far = []
            for ix, iy in ordered_cells:
                x0 = 2*(shift_x+ix*hx)-cx2
                x1 = x0+2*hx
                y0 = 2*(shift_y+iy*hy)-cy2
                y1 = y0+2*hy
                far.append(isqrt(max(x0*x0,x1*x1)+max(y0*y0,y1*y1))+1)
            far.sort()
            radii = {far[min(len(far)-1, floor((len(far)-1)*q))]
                     for q in radius_quantiles}
            for radius2 in radii:
                members = _circle_members(ordered_cells,hx,hy,shift_x,shift_y,
                                          cx2,cy2,radius2)
                candidate = _circle_cut(cx2,cy2,radius2,gap_ticks,grid_um,
                                        pad_bounds_um,members)
                if candidate is not None:
                    proposals[(cx2,cy2,radius2)] = candidate
    return list(proposals.values())


def _exact_dual_value(cuts, members, cell_count, multipliers, denominator):
    cover = [0]*cell_count
    weighted_capacity = 0
    for cut, indices, y in zip(cuts, members, multipliers):
        if y < 0:
            raise ValueError('Dual multipliers must be nonnegative')
        weighted_capacity += y*cut['crossing_capacity']
        for index in indices:
            cover[index] += y
    numerator = weighted_capacity+sum(max(0,denominator-value) for value in cover)
    return numerator//denominator, numerator


def gds_multi_circle_upper_bound(gds_path, layer, datatype,
                                 center_distance_um, wire_width_um,
                                 spacing_um, pad_bounds_um, *,
                                 centers_per_axis=5,
                                 radius_quantiles=(.2,.4,.6,.8,1.0),
                                 denominator=1_000_000):
    """Compute a safe integer upper bound despite a floating LP proposal."""
    if not pad_bounds_um:
        raise ValueError('Actual or universally allowed Pad bounds are required')
    if not isinstance(centers_per_axis,int) or centers_per_axis<1:
        raise ValueError('Circle center count must be a positive integer')
    if not radius_quantiles or any(not 0<q<=1 for q in radius_quantiles):
        raise ValueError('Radius quantiles must lie in (0,1]')
    cover, occupied = gds_cell_cover_upper_bound(
        gds_path,layer,datatype,center_distance_um,_return_occupied=True)
    global_cut = gds_radial_cut_upper_bound(gds_path,layer,datatype,
                                            wire_width_um,spacing_um)
    if cover['input_sha256'] != global_cut['input_sha256']:
        raise RuntimeError('GDS changed during multiple-circle audit')
    ordered = sorted(occupied)
    gap_ticks = global_cut['wire_plus_spacing_lower_grid_ticks']
    proposals = _proposed_cuts(cover,ordered,gap_ticks,pad_bounds_um,
                               centers_per_axis=centers_per_axis,
                               radius_quantiles=radius_quantiles)
    cuts = [p[0] for p in proposals]
    members = [p[1] for p in proposals]
    denominator = int(denominator)
    if denominator < 1:
        raise ValueError('Dual denominator must be positive')
    if cuts:
        row = [i for i,indices in enumerate(members) for _ in indices]
        col = [index for indices in members for index in indices]
        A = csr_matrix((np.ones(len(row)),(row,col)),
                       shape=(len(cuts),len(ordered)))
        result = linprog(-np.ones(len(ordered)),A_ub=A,
                         b_ub=np.array([cut['crossing_capacity'] for cut in cuts]),
                         bounds=(0,1),method='highs')
        proposed = -result.ineqlin.marginals if result.success else np.zeros(len(cuts))
        multipliers = [max(0, floor(float(y)*denominator)) for y in proposed]
    else:
        result = None
        multipliers = []
    value, numerator = _exact_dual_value(cuts,members,len(ordered),
                                          multipliers,denominator)
    return {
        'value': min(value,cover['value']),
        'multi_circle_dual_upper_bound': value,
        'independent_center_cover_upper_bound': cover['value'],
        'method': 'integer_cell_multi_circle_cuts_with_exact_rational_dual',
        'input_sha256': cover['input_sha256'],
        'support_layer': cover['support_layer'],
        'native_grid_um': cover['native_grid_um'],
        'required_center_distance_um': center_distance_um,
        'wire_width_um': wire_width_um,
        'spacing_um': spacing_um,
        'cell_cover': cover,
        'wire_plus_spacing_lower_grid_ticks': gap_ticks,
        'pad_bounds_um': [list(bounds) for bounds in pad_bounds_um],
        'circle_cuts': cuts,
        'dual_multipliers_numerator': multipliers,
        'dual_denominator': denominator,
        'dual_objective_numerator': numerator,
        'proposed_circle_count': len(cuts),
        'floating_lp_success': bool(result is not None and result.success),
        'floating_lp_status': None if result is None else result.message,
        'exact_integer_dual_verified': True,
        'proof': 'Each support-cover cell holds at most one center; centers in a fully-inside cell cross each Pad-excluding full circle; each crossing occupies wire width plus spacing on the circle; the integer dual multipliers and exact residual cell costs give an upper bound for every binary or fractional cell occupancy',
    }


def verify_multi_circle_certificate(gds_path, certificate):
    """Recompute every set and integer dual inequality, ignoring LP metadata."""
    layer, datatype = certificate['support_layer']
    cover, occupied = gds_cell_cover_upper_bound(
        gds_path,layer,datatype,certificate['required_center_distance_um'],
        _return_occupied=True)
    if (cover['input_sha256'] != certificate['input_sha256'] or
            cover['occupied_cell_index_sha256'] !=
            certificate['cell_cover']['occupied_cell_index_sha256']):
        raise ValueError('Input or support cell cover differs from certificate')
    if not certificate['pad_bounds_um']:
        raise ValueError('Pad bounds are missing')
    source_cut = gds_radial_cut_upper_bound(
        gds_path,layer,datatype,certificate['wire_width_um'],
        certificate['spacing_um'])
    if source_cut['wire_plus_spacing_lower_grid_ticks'] != certificate['wire_plus_spacing_lower_grid_ticks']:
        raise ValueError('Wire spacing denominator differs from certificate')
    ordered = sorted(occupied)
    hx, hy = cover['cell_side_grid_ticks']
    shift_x, shift_y = cover['cell_origin_shift_grid_ticks']
    gap_ticks = certificate['wire_plus_spacing_lower_grid_ticks']
    pad_bounds = certificate['pad_bounds_um']
    cuts = certificate['circle_cuts']
    members = []
    for cut in cuts:
        cx2,cy2 = cut['center_twice_grid_ticks']
        radius2 = cut['radius_twice_grid_ticks']
        indices = _circle_members(ordered,hx,hy,shift_x,shift_y,
                                  cx2,cy2,radius2)
        checked = _circle_cut(cx2,cy2,radius2,gap_ticks,
                              cover['native_grid_um'],pad_bounds,indices)
        if (checked is None or checked[0] != cut):
            raise ValueError('Circle member set, capacity or Pad exclusion differs')
        members.append(indices)
    denominator = certificate['dual_denominator']
    if not isinstance(denominator,int) or denominator<1:
        raise ValueError('Dual denominator must be a positive integer')
    multipliers = certificate['dual_multipliers_numerator']
    if len(cuts) != len(multipliers):
        raise ValueError('Dual multiplier count differs')
    value, numerator = _exact_dual_value(cuts,members,len(ordered),
                                          multipliers,denominator)
    if (numerator != certificate['dual_objective_numerator'] or
            min(value,cover['value']) != certificate['value']):
        raise ValueError('Exact dual objective differs')
    return True


def verify_pads_outside_multi_cuts(certificate, pad_bounds_um):
    """Check every actual Pad against every circle used by the dual."""
    grid_um=certificate['native_grid_um']
    for cut in certificate['circle_cuts']:
        item={'native_grid_um':grid_um,
              'circle_center_twice_grid_ticks':cut['center_twice_grid_ticks'],
              'radius_twice_grid_ticks':cut['radius_twice_grid_ticks']}
        if not verify_pads_outside_radial_cut(item,pad_bounds_um):
            return False
    return True
