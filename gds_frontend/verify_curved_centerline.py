"""Small geometric regression checks for v2's continuous bend construction."""

from __future__ import annotations

import math

import numpy as np
from shapely.geometry import box
from shapely.ops import unary_union

from curved_centerline import CurveSettings, UncertifiableCurve, smooth_centerline


def main():
    domain=unary_union([box(-5,-5,35,5),box(25,-5,35,35),box(25,25,65,35)])
    route=[(0,0),(30,0),(30,30),(60,30)]
    line,certificate=smooth_centerline(route,domain)
    assert certificate['curved_corner_count']==2
    assert certificate['continuous_control_hulls_inside_eroded_domain']
    assert line.coords[0]==route[0] and line.coords[-1]==route[-1]
    assert domain.covers(line)
    for bend in certificate['curve_segments']:
        start=np.asarray(bend['start_um'])
        control=np.asarray(bend['control_um'])
        end=np.asarray(bend['end_um'])
        in_tangent=control-start
        out_tangent=end-control
        assert np.linalg.norm(in_tangent)>0 and np.linalg.norm(out_tangent)>0
        assert bend['samples']>=2
        bound=np.linalg.norm(2*control-start-end)/(4*bend['samples']**2)
        assert bound<=certificate['maximum_chord_error_um']+1e-12
    narrow=unary_union([box(-0.01,-0.01,30.01,0.01),
                         box(29.99,-0.01,30.01,30.01),
                         box(29.99,29.99,60.01,30.01)])
    try:
        smooth_centerline(route,narrow,CurveSettings(min_tangent_um=0.25))
    except UncertifiableCurve:
        pass
    else:
        raise AssertionError('Uncertified tight corner was accepted')
    print('G1 tangent bends, convex-hull containment, endpoint and chord bounds passed')


if __name__=='__main__':
    main()
