"""Presentation-only spatial statistics on a fixed analytic target disk.

No routing module is imported. Pixel cells have analytic disk-intersection
areas; Lipschitz distance bounds expose spatial quadrature uncertainty.
"""
from __future__ import annotations

import base64
import math

import numpy as np
from scipy.spatial import cKDTree

REVISION = 'fixed_disk_distribution_v1'
RESOLUTIONS = (192, 384, 768)


def validate_target(radius_um, center_um, coverage_distance_um, resolution=384):
    values = [radius_um, *center_um, coverage_distance_um]
    if len(center_um) != 2 or any(isinstance(x, bool) or not isinstance(x, (int, float))
                                  or not math.isfinite(x) for x in values):
        raise ValueError('评价半径、圆心和覆盖距离须为有限数值')
    if radius_um <= 0 or coverage_distance_um < 0:
        raise ValueError('评价半径须大于 0，覆盖距离可为 0')
    if isinstance(resolution, bool) or resolution not in RESOLUTIONS:
        raise ValueError('评价精度须为 192、384 或 768')
    if not math.isfinite(2*radius_um) or any(not math.isfinite(abs(x)+radius_um) for x in center_um):
        raise ValueError('评价区域超出数值范围')
    return {'shape': 'analytic_disk', 'radius_um': float(radius_um),
            'center_um': list(map(float, center_um)),
            'coverage_distance_um': float(coverage_distance_um), 'resolution': int(resolution)}


def _quadrant_area(x, y):
    """Signed area [0,x] x [0,y] intersected with the unit disk."""
    a = np.clip(np.abs(x), 0, 1); b = np.clip(np.abs(y), 0, 1)
    cutoff = np.sqrt(np.maximum(0, 1-b*b))
    start = np.minimum(a, cutoff)
    def primitive(t):
        return .5*(t*np.sqrt(np.maximum(0, 1-t*t))+np.arcsin(t))
    return np.sign(x)*np.sign(y)*(b*start+primitive(a)-primitive(start))


def _disk_cells(resolution):
    edges = np.linspace(-1, 1, resolution+1)
    integral = _quadrant_area(edges[None, :], edges[:, None])
    weights = np.maximum(0, np.diff(np.diff(integral, axis=0), axis=1))
    centers = (edges[:-1]+edges[1:])/2
    xx, yy = np.meshgrid(centers, centers)
    near_x = np.maximum(np.abs(xx)-1/resolution, 0)
    near_y = np.maximum(np.abs(yy)-1/resolution, 0)
    active = (near_x*near_x+near_y*near_y <= 1) & (weights > 0)
    reps = np.column_stack((xx[active], yy[active]))
    lengths = np.linalg.norm(reps, axis=1)
    reps /= np.maximum(1, lengths)[:, None]
    # Projection onto a closed convex disk is nonexpansive. All points in
    # each intersecting square lie within sqrt(2)/resolution of its rep.
    area_error = abs(float(weights.sum())-math.pi)/math.pi
    if area_error > 1e-10:
        raise ValueError('解析圆与像素交集面积校验失败')
    return reps, weights[active]/math.pi, active, area_error


def evaluate_distribution(points_um, *, radius_um, center_um=(0, 0),
                          coverage_distance_um=100, resolution=384):
    target = validate_target(radius_um, center_um, coverage_distance_um, resolution)
    points = np.asarray(points_um, dtype=float).reshape(-1, 2)
    if not np.isfinite(points).all():
        raise ValueError('电极坐标含非有限数值')
    center = np.asarray(center_um, dtype=float); radius = float(radius_um)
    relative = (points-center)/radius
    if not np.isfinite(relative).all():
        raise ValueError('电极与评价区域的比例超出数值范围')
    count = len(points)
    report = {'revision': REVISION, 'target': target, 'electrode_count': count,
              'inside_target_count': int(np.count_nonzero(np.linalg.norm(relative, axis=1) <= 1)),
              'electrode_centers_um': points.tolist(),
              'distance_definition': 'Euclidean distance to the nearest electrode center; all accepted electrodes contribute, including centers outside the target',
              'scope': '2D geometric coverage statistics; not a calibrated recording range or a routing optimality certificate'}
    nearest = {'cv': None, 'mean_um': None, 'std_um': None,
               'minimum_um': None, 'maximum_um': None, 'status': 'fewer_than_two_electrodes'}
    if count >= 2:
        nn = cKDTree(points).query(points, k=2)[0][:, 1]
        mean = float(nn.mean()); std = float(nn.std(ddof=0))
        nearest = {'cv': std/mean if mean > 0 else None, 'mean_um': mean,
                   'std_um': std, 'minimum_um': float(nn.min()), 'maximum_um': float(nn.max()),
                   'status': 'defined' if mean > 0 else 'zero_mean_spacing',
                   'population_standard_deviation': True}
    report['nearest_neighbor'] = nearest
    if not count:
        report.update({'maximum_uncovered': {'estimate_um': None, 'lower_um': None,
                          'upper_um': None, 'unbounded': True, 'status': 'no_electrodes'},
                       'coverage': {'distance_um': float(coverage_distance_um),
                                    'estimate': 0.0, 'lower': 0.0, 'upper': 0.0},
                       'curve': {'distance_um': [0.0, max(2*radius, coverage_distance_um)],
                                 'estimate': [0.0, 0.0], 'lower': [0.0, 0.0], 'upper': [0.0, 0.0]},
                       'heatmap': None, 'numerics': {'status': 'empty_layout'}})
        return report
    reps, weights, active, area_error = _disk_cells(resolution)
    distances = cKDTree(relative).query(reps, workers=1)[0]*radius
    if not np.isfinite(distances).all():
        raise ValueError('评价距离超出数值范围')
    epsilon = math.sqrt(2)*radius/resolution
    guard = max(1e-9, float(distances.max())*1e-12, radius*1e-12)
    epsilon += guard
    order = np.argsort(distances); ordered = distances[order]
    cumulative = np.r_[0, np.cumsum(weights[order])]
    cumulative /= cumulative[-1]
    def cdf(thresholds):
        return cumulative[np.searchsorted(ordered, thresholds, side='right')]
    def covered(thresholds):
        return cdf(thresholds), cdf(thresholds-epsilon), cdf(thresholds+epsilon)
    worst = int(np.argmax(distances)); maximum = float(distances[worst])
    levels = np.unique(np.r_[np.linspace(0, max(2*radius, coverage_distance_um), 401), coverage_distance_um])
    estimate, lower, upper = covered(levels)
    e, lo, hi = covered(np.array([coverage_distance_um], dtype=float))
    # Finite point sets have zero area at distance 0. Handle it analytically.
    estimate[levels == 0] = lower[levels == 0] = upper[levels == 0] = 0
    if coverage_distance_um == 0: e[:] = lo[:] = hi[:] = 0
    field = np.zeros(active.shape, dtype='<u2')
    scale = max(float(distances.max()), 1e-12)
    field[active] = np.clip(np.rint(distances/scale*65534)+1, 1, 65535).astype('<u2')
    report.update({'maximum_uncovered': {'estimate_um': maximum,
                   'lower_um': max(0, maximum-guard), 'upper_um': maximum+epsilon,
                   'witness_um': (reps[worst]*radius+center).tolist(), 'unbounded': False,
                   'status': 'spatial_bracket'},
                  'coverage': {'distance_um': float(coverage_distance_um), 'estimate': float(e[0]),
                               'lower': float(lo[0]), 'upper': float(hi[0])},
                  'curve': {'distance_um': levels.tolist(), 'estimate': estimate.tolist(),
                            'lower': lower.tolist(), 'upper': upper.tolist()},
                  'heatmap': {'width': resolution, 'height': resolution,
                              'encoding': 'base64_uint16_le; 0=outside; distance=(value-1)/65534*maximum_um',
                              'values_base64': base64.b64encode(field.tobytes()).decode('ascii'),
                              'maximum_um': scale, 'quantization_error_um': scale/65534/2,
                              'rows': 'increasing original GDS y; canvas must flip vertically'},
                  'numerics': {'resolution': resolution, 'cell_step_um': 2*radius/resolution,
                               'distance_uncertainty_um': epsilon, 'sample_count': len(reps),
                               'relative_disk_area_error': area_error,
                               'method': 'analytic circle-cell areas; projected cell representatives; 1-Lipschitz distance envelope',
                               'roundoff_distance_guard_um': guard,
                               'floating_point_bounds': True}})
    return report
