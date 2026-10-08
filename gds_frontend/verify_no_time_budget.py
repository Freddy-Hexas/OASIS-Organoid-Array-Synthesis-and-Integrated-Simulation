"""Check unlimited defaults, solver options and real geometry under clock jumps.

No workbench listener is started. Synthetic geometry uses the same production
flow, exported GDS, independent integer audit and incumbent recertification.
"""
from inspect import signature
from itertools import count
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
import argparse
import json
import math
import time

import gdstk
import numpy as np
from shapely.geometry import box
from shapely.ops import unary_union

import island_router
import joint_path_flow
from center_compaction import CompactionSettings, compact_centers
from exact_gds_audit import audit_exported_gds
from island_router import IslandSettings, build_navigation
from joint_port_augment import AugmentSettings
from outer_exit_policy import make_outer_exit_policy
from pad_router import (PadSettings, audit_pad_gds, make_pad_frame,
                        solve_with_pads, write_pad_gds)
from process_geometry import ProcessRules
from route_incumbent import recertify_incumbent
from solver_time_policy import deadline_after, deadline_expired, solver_options
from verify_center_compaction import simple_net


def check_defaults_and_milp_options():
    assert CompactionSettings().time_limit_s is None
    assert AugmentSettings().time_limit_s is None
    assert IslandSettings().milp_time_limit_s is None
    assert signature(joint_path_flow.joint_path_routes).parameters['time_limit_s'].default is None
    assert signature(recertify_incumbent).parameters['time_limit_s'].default is None
    assert deadline_after(None) is None and not deadline_expired(None)
    assert solver_options(None, mip_rel_gap=0) == {'mip_rel_gap': 0}
    assert solver_options(5)['time_limit'] == 5
    for bad in [0, -1, math.inf, math.nan, True]:
        try:
            solver_options(bad)
        except ValueError:
            pass
        else:
            raise AssertionError(bad)
    columns=[]
    for i,x in enumerate([100,115,150,300]):
        columns.append({'source_node':i//2,'source_um':(x,0),'depth_um':1,
            'metal':box(x-2,-2,x+2,2),'island':box(x-5,-5,x+5,5),
            'pad_id':f'top-{i//2+1:02d}','pad_side':'top','pad_index':i//2+1})
    observed=[]
    original=island_router.milp
    def spy(*args,**kwargs):
        observed.append(dict(kwargs.get('options',{})))
        assert 'time_limit' not in observed[-1]
        return original(*args,**kwargs)
    with patch.object(island_router,'milp',side_effect=spy):
        chosen,record=island_router._select(columns,ProcessRules(),IslandSettings(),lambda *_:None,
            pad_order={'top':['top-01','top-02']},center_origin=(0,0),baseline_candidate_ids=[1,3])
    assert len(observed)==3,observed
    assert [c['source_um'][0] for c in chosen]==[100,300]
    assert record['finite_model_proven_optimal'] and record['secondary_proven_optimal']
    return {'baseline_and_primary_and_secondary_milp_options':observed}


def check_all_centers_survive_clock_jumps():
    jumps=count(10000,10000)
    with patch('time.perf_counter',side_effect=lambda:next(jumps)):
        result,record=compact_centers(box(-500,-90,500,90),
            [simple_net(300,35,1),simple_net(300,-35,2)],ProcessRules(),IslandSettings(),(0,0))
    assert len(result)==2 and not record['time_budget_hit']
    assert record['search_limits']['time_limit_s'] is None
    assert record['wall_clock_limit_enabled'] is False
    assert {a['network_index'] for a in record['attempts']}=={1,2}
    assert all(not a.get('search_deadline_hit') for a in record['attempts'])
    assert record['after']['sum_squared_radius_um2']<record['before']['sum_squared_radius_um2']
    return {'every_net_inspected':True,'moves_retained':len(record['accepted_updates']),
            'clock_step_seconds':10000,'time_budget_hit':False}


def check_real_pipeline(output_dir):
    rules=ProcessRules();settings=IslandSettings();pads=PadSettings()
    source=output_dir/'unlabeled_cross.gds';output=output_dir/'routing.gds'
    support=unary_union([box(-700,-20,700,20),box(-20,-700,20,700)])
    lib=gdstk.Library(unit=1e-6,precision=1e-9);cell=lib.new_cell('NO_METADATA')
    cell.add(gdstk.Polygon(np.asarray(support.exterior.coords),layer=42))
    with TemporaryDirectory(prefix='no_time_limit_source_') as temp:
        ascii_path=Path(temp)/'source.gds'
        lib.write_gds(str(ascii_path));source.write_bytes(ascii_path.read_bytes())
    policy=make_outer_exit_policy(support)
    nav=build_navigation(source,rules,layer=42,outlet_mode='external_boundary',
        settings=settings,outer_exit_policy=policy)
    milp_calls=[];lp_calls=[]
    original_milp=island_router.milp;original_lp=joint_path_flow.linprog
    def milp_spy(*args,**kwargs):
        milp_calls.append(kwargs.get('options',{}))
        assert 'time_limit' not in milp_calls[-1]
        return original_milp(*args,**kwargs)
    def lp_spy(*args,**kwargs):
        lp_calls.append(kwargs.get('options',{}))
        assert 'time_limit' not in lp_calls[-1]
        return original_lp(*args,**kwargs)
    jumps=count(10000,10000)
    with patch('time.perf_counter',side_effect=lambda:next(jumps)), \
         patch.object(island_router,'milp',side_effect=milp_spy), \
         patch.object(joint_path_flow,'linprog',side_effect=lp_spy):
        routing=solve_with_pads(nav,rules,settings=settings,pad_settings=pads)
    assert routing['retained_routes']>=4,routing['retained_routes']
    assert lp_calls,'The joint path LP must run even after the simulated former source deadline'
    for name in ['center_compaction','port_augmentation','joint_path_flow']:
        phase=routing[name]
        assert phase['wall_clock_limit_enabled'] is False,(name,phase)
        assert phase.get('time_budget_hit',False) is False,name
    layers=write_pad_gds(source,output,routing,(42,0))
    readback=audit_pad_gds(output,nav,routing,layers)
    assert readback['passed']
    audit=audit_exported_gds(source,output,layers,wire_spacing_um=rules.spacing_um,
        metal_support_margin_um=rules.margin_um,expected_nets=routing['retained_routes'],
        minimum_electrode_diameter_um=rules.electrode_diameter_um,
        minimum_electrode_center_distance_um=rules.minimum_center_spacing_um,
        maximum_electrode_center_radius_um=rules.electrode_region_radius_um,
        minimum_substrate_disk_radius_um=rules.electrode_diameter_um/2+rules.margin_um+settings.numeric_guard_um,
        maximum_substrate_disk_radius_um=rules.electrode_diameter_um/2+rules.margin_um+settings.numeric_guard_um+.01,
        minimum_island_spacing_um=settings.pad_gap_um,
        minimum_pad_short_side_um=pads.pad_width_um,minimum_pad_long_side_um=pads.pad_length_um,
        minimum_wire_width_um=rules.wire_width_um)
    assert audit['passed']
    (output_dir/'integer_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
    records=[{'job_id':str(i),'report':{'geometry':{'sha256':nav.summary['sha256']},
        'selected_support_layer':nav.summary['support_layer'],'routing':routing}} for i in range(2)]
    frame=make_pad_frame(nav.support,pads,exit_policy=policy)
    jumps=count(10000,10000)
    with patch('time.perf_counter',side_effect=lambda:next(jumps)):
        chosen,_,record=recertify_incumbent(records,nav,rules,settings,pads,frame,policy)
    assert len(chosen)==routing['retained_routes'],record
    assert len(record['records_checked'])==2 and not record['wall_clock_limit_enabled'],record
    return {'connected_count':len(chosen),'export_roundtrip_passed':True,'independent_integer_audit_passed':True,
        'lp_calls_without_time_limit':len(lp_calls),'milp_calls_without_time_limit':len(milp_calls),
        'all_default_phases_survive_clock_jumps':True,'both_historical_constructions_rechecked':True,
        'clock_step_seconds':10000}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path)
    args=parser.parse_args();started=time.perf_counter()
    defaults=check_defaults_and_milp_options();centers=check_all_centers_survive_clock_jumps()
    if args.output_dir is None:
        with TemporaryDirectory(prefix='no_layout_time_limit_') as temp:
            pipeline=check_real_pipeline(Path(temp))
    else:
        args.output_dir.mkdir(parents=True,exist_ok=True)
        pipeline=check_real_pipeline(args.output_dir)
    result={'passed':True,'defaults_and_milp':defaults,'center_clock_jumps':centers,
        'synthetic_unlabeled_gds':pipeline,'real_elapsed_seconds':time.perf_counter()-started,
        'starts_or_restarts_web_server':False}
    if args.output_dir:
        (args.output_dir/'verification.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
