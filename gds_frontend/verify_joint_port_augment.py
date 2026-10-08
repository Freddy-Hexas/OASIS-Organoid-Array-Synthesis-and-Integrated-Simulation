"""Ordered exterior routing and complete multi-track residual additions."""
from dataclasses import asdict
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import patch
import json
import math

import gdstk
import networkx as nx
import numpy as np
from shapely.geometry import LineString,Point,box
from shapely.ops import unary_union

from exact_gds_audit import audit_exported_gds
from island_router import IslandSettings,build_navigation,betti
from joint_port_augment import AugmentSettings,_interior_seed,augment_joint_ports
from ordered_pad_fanout import ordered_pad_assignment,polar_fanout,rebuild_ordered_outer
from outer_exit_policy import make_outer_exit_policy
from pad_router import PadSettings,make_pad_frame,_pad_occupancy,write_pad_gds,audit_pad_gds
from pad_router import solve_with_pads
import pad_router
from process_geometry import ProcessRules


def assignment_checks():
    source=Point(0,0).buffer(3500)
    frame=make_pad_frame(source)
    angles=np.arange(0,360,15,dtype=float)+.1
    exits=np.column_stack((3500*np.sin(np.radians(angles)),3500*np.cos(np.radians(angles))))
    assignment,record=ordered_pad_assignment(exits,frame)
    assert assignment is not None and record['minimum_endpoint_angular_gap_rad']>0
    for bank in record['banks'].values():
        slots=bank['slot_indices']
        assert not slots or slots[-1]-slots[0]+1==len(slots)
    translated={**frame,'origin_um':[123.0,-456.0],
                'slots':[{**s,'target_um':np.asarray(s['target_um'])+[123,-456]} for s in frame['slots']]}
    shifted,_=ordered_pad_assignment(exits+[123,-456],translated)
    assert [assignment[i][0]['pad_id'] for i in range(len(exits))]==[
        shifted[i][0]['pad_id'] for i in range(len(exits))]
    rotated,_=ordered_pad_assignment(np.column_stack((-exits[:,1],exits[:,0])),frame)
    assert rotated is not None
    a,b=assignment[0][1:]
    points,certificate=polar_fanout(a,b,4100,12800,np.array([0.,0.]),.01)
    assert certificate['maximum_chord_error_um']<=.01+1e-12
    assert np.isclose(np.linalg.norm(points[0]),4100) and np.isclose(np.linalg.norm(points[-1]),12800)
    crowded=exits[[0]*30]
    rejected,reason=ordered_pad_assignment(crowded,frame)
    assert rejected is None and reason['reason']=='ordered_bank_has_too_few_slots'
    duplicated,reason=ordered_pad_assignment(exits[[0,0]],frame)
    assert duplicated is None and reason['reason']=='source_or_target_cyclic_order_is_degenerate'
    # Strict angular order proves only disjoint zero-width curves. These
    # different close outlets still fail the actual finite-width check.
    narrow=box(-1000,-40,1000,40);policy=make_outer_exit_policy(narrow)
    graph=nx.Graph()
    for i,y in enumerate((-.5,.5)):
        graph.add_node(i,xy_um=[993.4,y],terminal_kind='external_boundary')
    nav=SimpleNamespace(support=narrow,graph=graph,summary={'outer_exit_policy':policy.record()})
    rules=ProcessRules();pads=PadSettings()
    seeds=[_interior_seed(nav,i,rules) for i in range(2)]
    rebuilt,reason=rebuild_ordered_outer(nav,seeds,rules,IslandSettings(),pads,make_pad_frame(narrow),policy)
    assert rebuilt is None and reason['reason']=='ordered_fanout_full_metal_conflict',reason


def geometry_checks():
    rules=ProcessRules();pads=PadSettings();settings=IslandSettings()
    with TemporaryDirectory(prefix='joint_ports_') as tmp:
        source=Path(tmp)/'unlabeled_geometry.gds';output=Path(tmp)/'connected.gds'
        # Support only: no generator centerlines, topology labels or metadata.
        support=unary_union([box(-700,-20,700,20),box(-20,-700,20,700)])
        lib=gdstk.Library(unit=1e-6,precision=1e-9);cell=lib.new_cell('UNLABELED')
        for part in (support.geoms if hasattr(support,'geoms') else [support]):
            cell.add(gdstk.Polygon(np.asarray(part.exterior.coords),layer=42))
        lib.write_gds(str(source))
        policy=make_outer_exit_policy(support)
        nav=build_navigation(source,rules,layer=42,outlet_mode='external_boundary',
            settings=settings,outer_exit_policy=policy)
        frame=make_pad_frame(nav.support,pads,exit_policy=policy)
        # An empty initial library is not geometric infeasibility. The public
        # solver must still construct networks, with internally consistent
        # count, status, Pad occupancy and optimization-bound scope.
        with patch.object(pad_router,'_route_columns',return_value=([],{'path_columns':0})), \
             patch.object(pad_router,'joint_path_routes',return_value=([],{'status':'empty_test_proposals'})):
            from_empty=solve_with_pads(nav,rules,settings=settings,pad_settings=pads)
        assert from_empty['retained_routes']>=4
        assert from_empty['status']=='awaiting_export_audit'
        assert from_empty['selected_pad_count']==from_empty['retained_routes']
        assert from_empty['optimization']['finite_upper_bound'] is None
        assert from_empty['initial_finite_library_optimization']['finite_upper_bound']==0
        initial,first=augment_joint_ports(nav,[],rules,settings,pads,frame,policy,
            capacity_upper=2,search=AugmentSettings(time_limit_s=30,max_sweeps=1))
        assert len(initial)==2,first
        frame=first.pop('_frame')
        updated,record=augment_joint_ports(nav,initial,rules,settings,pads,frame,policy,
            search=AugmentSettings(time_limit_s=35,max_sweeps=1))
        frame=record.pop('_frame')
        assert len(updated)>len(initial),record
        assert record['after']['used_window_count']>record['before']['used_window_count']
        assert record['after']['used_window_count']==4,record['after']
        assert not record['one_track_per_window_imposed']
        assert len(updated)>4,'The 40 um support can carry multiple real tracks per window'
        assert record['count_nondecreasing'] and not record['global_maximum_proven']
        assert len({c['pad_id'] for c in updated})==len(updated)
        assert all(bank['internal_empty_slots']==0 for bank in _pad_occupancy(updated,pads).values())
        # The outside may be reassigned, while existing interior geometry and
        # electrodes remain unchanged during an addition.
        for old,new in zip(initial,updated):
            assert old['central_wire'].equals(new['central_wire'])
            assert old['source_um']==new['source_um']
        final=unary_union([nav.support,frame['shell'],*[c['bridge'] for c in updated],
                          *[c['island'] for c in updated]])
        routing={'_chosen':updated,'_final_support':final,'_shell':frame['shell'],
                 'rules':asdict(rules),'pad_settings':asdict(pads),
                 'outer_exit_policy':policy.record()}
        layers=write_pad_gds(source,output,routing,(42,0))
        roundtrip=audit_pad_gds(output,nav,routing,layers)
        exact=audit_exported_gds(source,output,layers,wire_spacing_um=4,
            metal_support_margin_um=4,expected_nets=len(updated),
            minimum_electrode_diameter_um=30,minimum_electrode_center_distance_um=70,
            minimum_substrate_disk_radius_um=19.05,maximum_substrate_disk_radius_um=19.06,
            minimum_island_spacing_um=4,minimum_pad_short_side_um=500,
            minimum_pad_long_side_um=3000,minimum_wire_width_um=5)
        assert exact['passed'] and exact['outer_exit_policy_verified']
        assert exact['original_support_exactly_preserved']
        return {'initial_count':len(initial),'final_count':len(updated),
            'windows_before':record['before']['used_window_count'],
            'windows_after':record['after']['used_window_count'],
            'integer_audit_passed':True,
            'minimum_metal_spacing_um':roundtrip['minimum_inter_net_gap_um']}


def main():
    assignment_checks()
    print(json.dumps({'passed':True,'synthetic_unknown_gds':geometry_checks()},indent=2),flush=True)


if __name__=='__main__':main()
