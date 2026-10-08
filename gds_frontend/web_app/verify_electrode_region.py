"""Targeted hard-radius verification; never reruns the production case set."""
from pathlib import Path
from tempfile import TemporaryDirectory
from dataclasses import asdict
import sys
import json
import math
import argparse

import gdstk
import numpy as np
import networkx as nx
from shapely.geometry import box, Point

import server
from electrode_region import (ElectrodeRegion,region_for,add_region_candidates,
                              exported_center_in_disk)
from process_geometry import ProcessRules
from exact_gds_audit import audit_exported_gds
from center_compaction import _try_replacement
from joint_port_augment import _interior_seed
from island_router import IslandSettings
from types import SimpleNamespace


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.parent.mkdir(parents=True,exist_ok=True)
    checks=[]
    assert server._validated_settings({'method':'four_side_pads'})[1].electrode_region_radius_um==3000
    for value in (0,-1,True,None,math.nan,math.inf,'3000'):
        try:server._validated_settings({'rules':{'electrode_region_radius_um':value}})
        except ValueError:pass
        else:raise AssertionError(('invalid radius accepted',value))
    checks.append('default 3000 um and invalid radius API rejection')
    region=region_for(box(900,-100,1100,100),ProcessRules(electrode_region_radius_um=30))
    assert region.origin==(1000.,0.)
    assert region.contains([1030,0]) and not region.contains([1030.000001,0])
    assert exported_center_in_disk((2060000,0),(1000000,0),30,.001)
    assert not exported_center_in_disk((2060001,0),(1000000,0),30,.001)
    assert exported_center_in_disk((2000001,0),(1000000,0),.0005,.001)
    assert not exported_center_in_disk((2000001,1),(1000000,0),.0005,.001)
    assert region.packing_upper_bound(70)['value']==1
    assert region.packing_upper_bound(60)['value']==4 # bound remains conservative at equality
    assert region.clip(box(0,0,10,10)).is_empty
    checks.append('translated center, inclusive boundary, half-grid exact rejection, disjoint-domain and rational packing')
    g=nx.Graph()
    for n,xy in ((0,[500.,0.]),(1,[1500.,0.])):g.add_node(n,xy_um=xy,kind='leaf')
    g.add_edge(0,1,points_um=np.array([[500.,0.],[1500.,0.]]),weight=1000.,corridor_id=0)
    candidates,added=add_region_candidates(g,[0,1],region)
    assert added==1 and len(candidates)==1 and region.contains(g.nodes[candidates[0]]['xy_um'])
    assert nx.has_path(g,0,1) and math.isclose(nx.shortest_path_length(g,0,1,weight='weight'),1000.)
    assert g.nodes[0]['xy_um']==[500.,0.] and g.nodes[1]['xy_um']==[1500.,0.]
    checks.append('small disk adds missing interior candidate without cutting out-of-disk wire graph')
    residual_graph=nx.Graph();residual_graph.add_node(1,xy_um=[1690.,-432.])
    residual_nav=SimpleNamespace(graph=residual_graph)
    residual_rules=ProcessRules(electrode_region_radius_um=40,margin_um=1)
    residual_support=box(534,-452,1934,-412)
    seed=_interior_seed(residual_nav,1,residual_rules)
    proposal,detail=_try_replacement(residual_support,seed,[],residual_rules,
        IslandSettings(),np.array([1234.,-432.]),.001,32,inward_only=False)
    assert proposal is not None,detail
    assert math.dist(proposal['source_um'],[1234.,-432.])<=40
    assert math.dist(proposal['points_um'][-1],[1234.,-432.])>400
    checks.append('residual augmentation moves exterior seed into allowed disk while retaining distant exit')

    original_runs=server.RUNS_DIR;tag='radius-targeted';record=None
    try:
        with TemporaryDirectory(prefix='gds_radius_verify_') as temporary:
            base=Path(temporary);source=base/'unfamiliar_translated_strip.gds'
            lib=gdstk.Library(unit=1e-6,precision=1e-9);cell=lib.new_cell('UNKNOWN_TRANSLATED_SUPPORT')
            center=np.array([1234.,-432.])
            cell.add(gdstk.rectangle(center+[-700,-20],center+[-20,20],layer=42),
                     gdstk.rectangle(center+[20,-20],center+[700,20],layer=42));lib.write_gds(str(source))
            server.RUNS_DIR=base/'runs';server.RUNS_DIR.mkdir();folder=server.RUNS_DIR/tag;folder.mkdir()
            values,rules,pads,_,_=server._validated_settings({'method':'four_side_pads','rules':{
                'electrode_region_radius_um':40.,'minimum_center_spacing_um':120.,
                'margin_um':1.,'electrode_diameter_um':30.,'wire_width_um':5.,'spacing_um':4.}})
            server.JOBS[tag]={'id':tag,'file_id':tag,'created_at':server.utc_now(),'input_name':source.name,
                              'method':'four_side_pads','rules':values,'status':'queued','artifacts':[]}
            print('Running one translated support through placement, Pad routing, export and integer audit',flush=True)
            server._run_island_job(tag,source,rules,'external_boundary',(42,0),pad_mode=True,pad_settings=pads)
            status=json.loads((folder/'status.json').read_text(encoding='utf-8'))
            assert status['status']=='complete',(status.get('message'),(folder/'error.log').read_text() if (folder/'error.log').exists() else '')
            result=json.loads((folder/'summary.json').read_text(encoding='utf-8'));routing=result['routing']
            assert routing['retained_routes']==1,routing['candidate_library']
            assert routing['electrode_region']['passed'] and routing['electrode_region']['reference_um']==center.tolist()
            assert all(math.dist(r['source_um'],center)<=40 for r in routing['routes'])
            witness=json.loads((folder/'routing.width_witness.json').read_text(encoding='utf-8'))
            assert any(math.dist(np.array(p)*witness['native_grid_um'],center)>600
                       for r in witness['networks'] for p in r['points_grid_ticks'])
            assert routing['gds_roundtrip_audit']['passed']
            audit=routing['integer_polygon_audit']
            assert audit['passed'] and audit['electrode_region_verified'],audit
            assert audit['electrode_region_reference_um']==center.tolist()
            assert audit['maximum_electrode_center_radius_um']==40
            assert result['capacity_interval']['integer_polygon_lower_verified']
            assert routing['center_distance_upper_bound']['electrode_region_packing_certificate']['value']==1
            # Audit the same actual marker against an impossibly small disk;
            # perturbing report source coordinates cannot bypass this check.
            marker_radius=math.dist(routing['routes'][0]['source_um'],center)
            negative_test='exact half-grid predicate separately tested'
            if marker_radius>.01:
                try:
                    audit_exported_gds(source,folder/'routing.gds',routing['gds_roundtrip_audit']['output_layers'],
                        wire_spacing_um=4,metal_support_margin_um=1,expected_nets=1,
                        minimum_electrode_diameter_um=30,maximum_electrode_center_radius_um=marker_radius/2)
                except ValueError as exc:
                    assert 'outside permitted disk' in str(exc),str(exc)
                    negative_test=str(exc)
                else:raise AssertionError('Exported out-of-circle marker passed radius audit')
            for name in ('routing.gds','routing.width_witness.json','summary.json','integer_polygon_audit.json'):
                (args.output.parent/('synthetic_'+name)).write_bytes((folder/name).read_bytes())
            record={'rules':values,'connected_count':routing['retained_routes'],
                    'electrode_region':routing['electrode_region'],'integer_audit':audit,
                    'wires_leave_placement_disk':True,'export_negative_test':negative_test,
                    'capacity_interval':result['capacity_interval']}
            checks.append('end-to-end translated-source electrode inside disk, distant physical Pad, and exact GDS radius audit')
            # The entire requested disk lies in the gap between the two strips.
            # Empty placement must remain empty, rather than placing outside.
            empty_tag='radius-empty';empty_folder=server.RUNS_DIR/empty_tag;empty_folder.mkdir()
            empty_rules=ProcessRules(**{**asdict(rules),'electrode_region_radius_um':10.})
            server.JOBS[empty_tag]={'id':empty_tag,'file_id':empty_tag,'created_at':server.utc_now(),
                'input_name':source.name,'method':'four_side_pads','rules':asdict(empty_rules),
                'status':'queued','artifacts':[]}
            server._run_island_job(empty_tag,source,empty_rules,'external_boundary',(42,0),pad_mode=True,pad_settings=pads)
            empty=json.loads((empty_folder/'status.json').read_text(encoding='utf-8'))
            assert empty['status']=='complete',empty.get('message')
            empty_result=json.loads((empty_folder/'summary.json').read_text(encoding='utf-8'))
            assert empty_result['routing']['retained_routes']==0
            assert not (empty_folder/'routing.gds').exists()
            assert not empty_result['capacity_interval']['global_optimality_proven']
            record['empty_region_result']={'retained_routes':0,'exported_gds':False,'status':'complete'}
            checks.append('disk disjoint from support yields zero routes and no out-of-region fallback')
    finally:
        server.RUNS_DIR=original_runs
        server.JOBS.pop(tag,None);server.RECENT_JOBS.pop(tag,None)
        server.JOBS.pop('radius-empty',None);server.RECENT_JOBS.pop('radius-empty',None)
    args.output.write_text(json.dumps({'status':'passed','checks':checks,'real_pipeline':record},ensure_ascii=False,indent=2),encoding='utf-8')
    print('PASS',*checks,sep='\n',flush=True)


if __name__=='__main__':main()
