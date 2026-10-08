"""Analytic reference cases for descriptive spatial coverage, no routing."""
import argparse
import base64
import json
import math
import time
from pathlib import Path

import numpy as np
from distribution_metrics import evaluate_distribution, validate_target


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();checks=[];start=time.perf_counter()
    for radius in (0,-1,math.nan,math.inf,True):
        try:validate_target(radius,[0,0],100)
        except ValueError:pass
        else:raise AssertionError('invalid target accepted')
    for center in ([0,math.inf],[0],[0,0,0]):
        try:validate_target(1,center,0)
        except ValueError:pass
        else:raise AssertionError('invalid center accepted')
    checks.append('invalid targets rejected; zero coverage distance allowed')
    empty=evaluate_distribution([],radius_um=1500)
    assert empty['electrode_count']==0 and empty['nearest_neighbor']['cv'] is None
    assert empty['maximum_uncovered']['unbounded'] and empty['coverage']['upper']==0
    checks.append('empty layout: undefined CV, infinite maximum distance, zero coverage')
    center=[1234.5,-432.5]
    single=evaluate_distribution([center],radius_um=1500,center_um=center,coverage_distance_um=100)
    assert single['nearest_neighbor']['cv'] is None
    assert single['maximum_uncovered']['lower_um']<=1500<=single['maximum_uncovered']['upper_um']
    exact=(100/1500)**2
    assert single['coverage']['lower']<=exact<=single['coverage']['upper']
    for l,lo,hi in zip(single['curve']['distance_um'],single['curve']['lower'],single['curve']['upper']):
        assert lo-1e-12<=min((l/1500)**2,1)<=hi+1e-12,(l,lo,hi)
    checks.append('translated single central electrode: h=R and entire analytic C(l)=(l/R)^2 bracketed')
    outer=evaluate_distribution([[3000,0]],radius_um=1500,coverage_distance_um=100)
    assert outer['inside_target_count']==0 and outer['electrode_count']==1
    assert outer['maximum_uncovered']['lower_um']<=4500<=outer['maximum_uncovered']['upper_um']
    assert outer['coverage']['upper']==0
    checks.append('outside electrodes retained; actual target area is never clipped to support')
    theta=np.arange(32)*2*math.pi/32
    points=np.c_[1200*np.cos(theta),1200*np.sin(theta)]
    ring=evaluate_distribution(points,radius_um=1500)
    assert ring['nearest_neighbor']['cv']<1e-12
    assert abs(ring['nearest_neighbor']['mean_um']-2400*math.sin(math.pi/32))<1e-9
    assert ring['maximum_uncovered']['lower_um']<=1200<=ring['maximum_uncovered']['upper_um']
    checks.append('equal ring: CV=0 yet h=1200 um; no single-score uniformity shortcut')
    angle=.731;rotation=np.array([[math.cos(angle),-math.sin(angle)],[math.sin(angle),math.cos(angle)]])
    transformed=points@rotation.T+center
    shifted=evaluate_distribution(transformed,radius_um=1500,center_um=center)
    assert abs(shifted['nearest_neighbor']['mean_um']-ring['nearest_neighbor']['mean_um'])<1e-9
    assert shifted['maximum_uncovered']['lower_um']<=1200<=shifted['maximum_uncovered']['upper_um']
    assert abs(shifted['coverage']['estimate']-ring['coverage']['estimate'])<.005
    checks.append('translation/rotation invariance within declared grid uncertainty')
    detailed=evaluate_distribution([center],radius_um=1500,center_um=center,resolution=768)
    assert detailed['numerics']['distance_uncertainty_um']<single['numerics']['distance_uncertainty_um']
    assert detailed['numerics']['relative_disk_area_error']<1e-10
    for array in (ring['curve']['estimate'],ring['curve']['lower'],ring['curve']['upper']):
        assert min(array)>=0 and max(array)<=1 and np.all(np.diff(array)>=-1e-12)
    raster=np.frombuffer(base64.b64decode(ring['heatmap']['values_base64']),dtype='<u2')
    assert len(raster)==384**2 and np.count_nonzero(raster)==ring['numerics']['sample_count']
    assert ring['heatmap']['quantization_error_um']<.01
    checks.append('analytic disk area, monotone curves, heatmap encoding, finer-grid uncertainty reduction')
    result={'status':'passed','checks':checks,'seconds':time.perf_counter()-start,
            'single_electrode':{key:single[key] for key in ('nearest_neighbor','maximum_uncovered','coverage','numerics')},
            'ring_example':{key:ring[key] for key in ('electrode_count','nearest_neighbor','maximum_uncovered','coverage')}}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({'status':'passed','checks':checks,'seconds':result['seconds']},ensure_ascii=False))


if __name__=='__main__':main()
