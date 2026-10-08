"""Construct G1 centerlines with conservative certificates in an eroded domain.

Each resolved polyline corner becomes a quadratic Bezier whose tangent matches
the adjacent straight pieces.  The Bezier lies in its control triangle.  Requiring
that entire triangle to lie in the wire-center domain certifies the *continuous*
curve, rather than only a set of sampled vertices.  GDS uses a bounded-error
polygonal approximation of that analytic curve.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
import shapely
from shapely.geometry import LineString, Polygon


@dataclass(frozen=True)
class CurveSettings:
    max_tangent_fraction: float = 0.45
    min_tangent_um: float = 0.25
    min_corner_angle_deg: float = 0.01
    max_chord_error_um: float = 0.01
    max_samples_per_bend: int = 256


class UncertifiableCurve(ValueError):
    """A route corner cannot be rounded inside the specified center domain."""


def _clean_points(points):
    clean = []
    for point in np.asarray(points, dtype=float):
        if not np.isfinite(point).all():
            raise ValueError('Centerline contains non-finite coordinates')
        if not clean or np.linalg.norm(point - clean[-1]) > 1e-7:
            clean.append(point)
    if len(clean) < 2:
        raise ValueError('Centerline needs two distinct endpoints')
    return np.asarray(clean)


def _prune_short_legs(points, domain, threshold_um=2.0):
    """Drop numerical micro-legs only when their replacement chord is legal."""
    work=[point for point in points]
    removed=0
    index=1
    while index<len(work)-1:
        before,vertex,after=work[index-1:index+2]
        if (min(np.linalg.norm(vertex-before),np.linalg.norm(after-vertex))<threshold_um
                and shapely.covers(domain,LineString([before,after]))):
            work.pop(index)
            removed+=1
            index=max(1,index-1)
        else:
            index+=1
    return np.asarray(work),removed


def _quadratic(p0, p1, p2, parameter):
    t = parameter
    return (1-t)**2*p0 + 2*(1-t)*t*p1 + t*t*p2


def _certified_tangent(vertex, incoming, outgoing, cap, domain, settings):
    """Maximise a local bend while its nested control triangle remains legal."""
    def fits(tangent):
        p0 = vertex - tangent*incoming
        p2 = vertex + tangent*outgoing
        triangle = Polygon([p0, vertex, p2])
        return shapely.covers(domain, triangle)

    if fits(cap):
        return cap
    if not fits(settings.min_tangent_um):
        raise UncertifiableCurve('No legal tangent-continuous bend at a route corner')
    low, high = settings.min_tangent_um, cap
    for _ in range(18):
        middle = (low+high)/2
        if fits(middle):
            low = middle
        else:
            high = middle
    return low


def smooth_centerline(points, center_domain, settings=CurveSettings()):
    """Return sampled line and exact Bezier controls for a certified G1 path.

    All significant internal corners are rounded.  Endpoints stay fixed.  The
    maximum deviation of each quadratic from its sampled chords is bounded by
    ``max_chord_error_um`` in the continuous plane.  If a resolved bend cannot
    fit the eroded domain, reject this path column rather than silently retain
    a kink.
    """
    if not (0 < settings.max_tangent_fraction < 0.5):
        raise ValueError('max_tangent_fraction must be in (0, 0.5)')
    if settings.min_tangent_um <= 0 or settings.max_chord_error_um <= 0:
        raise ValueError('Curve lengths and approximation error must be positive')

    shapely.prepare(center_domain)
    input_points = _clean_points(points)
    source,removed_short_legs = _prune_short_legs(input_points,center_domain)
    output = [source[0]]
    controls = []
    tangent_lengths = []
    bend_radii = []
    min_angle = math.radians(settings.min_corner_angle_deg)
    for index in range(1, len(source)-1):
        before, vertex, after = source[index-1:index+2]
        incoming_vector = vertex-before
        outgoing_vector = after-vertex
        length_before = float(np.linalg.norm(incoming_vector))
        length_after = float(np.linalg.norm(outgoing_vector))
        incoming = incoming_vector/length_before
        outgoing = outgoing_vector/length_after
        cross = incoming[0]*outgoing[1]-incoming[1]*outgoing[0]
        turn = math.atan2(float(cross), float(np.dot(incoming, outgoing)))
        angle = abs(turn)
        if angle < min_angle:
            if np.linalg.norm(vertex-output[-1]) > 1e-7:
                output.append(vertex)
            continue
        if angle >= math.pi-1e-5:
            raise UncertifiableCurve('A near-reversal has no regular local quadratic bend')
        cap = settings.max_tangent_fraction*min(length_before,length_after)
        if cap < settings.min_tangent_um:
            raise UncertifiableCurve('Not enough adjacent line length for a certified bend')
        tangent = _certified_tangent(vertex,incoming,outgoing,cap,center_domain,settings)
        p0 = vertex-tangent*incoming
        p1 = vertex
        p2 = vertex+tangent*outgoing
        if np.linalg.norm(p0-output[-1]) > 1e-7:
            output.append(p0)
        # For a quadratic Bezier the chord deviation on 1/N parameter spans is
        # at most ||2 P1 - P0 - P2|| / (4 N^2).
        numerator = float(np.linalg.norm(2*p1-p0-p2))
        samples = max(2,math.ceil(math.sqrt(numerator/(4*settings.max_chord_error_um))))
        if samples > settings.max_samples_per_bend:
            raise UncertifiableCurve('Curve requires more samples than the declared limit')
        for step in range(1,samples+1):
            output.append(_quadratic(p0,p1,p2,step/samples))
        controls.append({'start_um':p0.tolist(),'control_um':p1.tolist(),
                         'end_um':p2.tolist(),'turn_angle_deg':math.degrees(turn),
                         'tangent_length_um':tangent,'samples':samples})
        tangent_lengths.append(tangent)
        # The largest curvature of this symmetric tangent construction occurs
        # at the middle of the quadratic bend.
        half=angle/2
        bend_radii.append(tangent*math.cos(half)**2/math.sin(half))
    if np.linalg.norm(source[-1]-output[-1]) > 1e-7:
        output.append(source[-1])
    line = LineString(np.asarray(output))
    if not shapely.covers(center_domain,line):
        raise UncertifiableCurve('Flattened centerline leaves the legal center domain')
    return line, {'representation':'piecewise_quadratic_bezier_G1',
                  'curve_segments':controls,'curved_corner_count':len(controls),
                  'minimum_tangent_um':min(tangent_lengths,default=None),
                  'minimum_bend_radius_um':min(bend_radii,default=None),
                  'maximum_chord_error_um':settings.max_chord_error_um,
                  'continuous_control_hulls_inside_eroded_domain':True,
                  'removed_short_legs':removed_short_legs,
                  'input_waypoint_count':len(input_points),
                  'original_waypoint_count':len(source),
                  'flattened_waypoint_count':len(output)}
