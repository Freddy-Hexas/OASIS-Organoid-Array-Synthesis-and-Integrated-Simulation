"""Custom process rules through validation, geometry, export and strict audit."""
from dataclasses import asdict
from pathlib import Path
from tempfile import TemporaryDirectory
import argparse
import json
import math

import gdstk
import numpy as np
from shapely.geometry import box, GeometryCollection

import server
from island_router import _physical_conflicts


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    checks=[]
    for changes in ({'wire_width_um':0},{'electrode_diameter_um':0},
                    {'spacing_um':-1},{'margin_um':-1},
                    {'minimum_center_spacing_um':-1},{'margin_um':math.nan},
                    {'wire_width_um':math.inf},{'wire_width_um':True}):
        try:
            server._validated_settings({'method':'four_side_pads','rules':changes})
        except ValueError:
            checks.append({'invalid_rule':str(changes),'rejected':True})
        else:
            raise AssertionError(f'Invalid rule accepted: {changes}')
    _,wide,pads,_,_=server._validated_settings({'method':'four_side_pads',
        'rules':{'electrode_diameter_um':45,'wire_width_um':80,'margin_um':130,
                 'spacing_um':0,'minimum_center_spacing_um':150000}})
    assert pads.bridge_width_um>=wide.wire_width_um+2*wide.margin_um
    assert pads.pad_outer_setback_um>=wide.margin_um
    _,_,tiny,_,_=server._validated_settings({'method':'four_side_pads',
        'rules':{'wire_width_um':.5,'spacing_um':0,'margin_um':0},
        'pad_settings':{'square_side_um':.1,'pad_width_um':3,'pad_length_um':4,
                        'pad_pitch_um':3,'minimum_pad_width_um':2,'minimum_pad_length_um':3}})
    assert tiny.pad_width_um==3 and tiny.square_side_um==.1
    # Arbitrarily small positive gaps, and zero gaps, must detect shorts.
    codes=np.array([0,1],dtype=np.int32)
    pad_codes=np.array([-1,-1],dtype=np.int32)
    for collection in (False,True):
        for gap in (0.,1e-9):
            shapes=[box(0,0,1,1),box(1,0,2,1)]
            if collection:shapes=[GeometryCollection([shape]) for shape in shapes]
            assert _physical_conflicts(shapes,gap,codes,pad_codes)=={(0,1)}
    original_runs=server.RUNS_DIR
    records=[]
    try:
        with TemporaryDirectory(prefix='gds_custom_rules_') as temporary:
            base=Path(temporary)
            source=base/'unfamiliar_support.gds'
            lib=gdstk.Library(unit=1e-6,precision=1e-9)
            cell=lib.new_cell('CUSTOM_RULES_UNKNOWN_SOURCE')
            cell.add(gdstk.rectangle((-700,-20),(700,20),layer=42),
                     gdstk.rectangle((600,-150),(900,150),layer=42))
            lib.write_gds(str(source))
            server.RUNS_DIR=base/'runs';server.RUNS_DIR.mkdir()
            for tag,gap,margin in (('sub_default',1.,.5),('zero_gap_margin',0.,0.)):
                settings=server._validated_settings({'method':'four_side_pads',
                    'rules':{'electrode_diameter_um':20,'wire_width_um':2,
                             'spacing_um':gap,'margin_um':margin,
                             'minimum_center_spacing_um':5000}})
                defaults,rules,pad_settings,_,_=settings
                output=server.RUNS_DIR/tag;output.mkdir()
                server.JOBS[tag]={'id':tag,'file_id':tag,'created_at':server.utc_now(),
                                  'input_name':source.name,'method':'four_side_pads',
                                  'rules':defaults,'status':'queued','artifacts':[]}
                print(f'Running real geometry/export/audit: {tag}',flush=True)
                server._run_island_job(tag,source,rules,'external_boundary',(42,0),
                                       pad_mode=True,pad_settings=pad_settings)
                job=json.loads((output/'status.json').read_text(encoding='utf-8'))
                if job['status']!='complete':
                    raise AssertionError((job,(output/'error.log').read_text(encoding='utf-8')))
                result=job['result'];routing=result['routing']
                if args.output:
                    args.output.parent.mkdir(parents=True,exist_ok=True)
                    (args.output.parent/f'{tag}.summary.json').write_text(
                        json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
                assert routing['retained_routes']>=1,routing['candidate_library']
                assert routing['rules']==asdict(rules)
                assert routing['settings']['pad_gap_um']==gap
                assert routing['gds_roundtrip_audit']['passed']
                audit=routing['integer_polygon_audit']
                assert audit.get('passed') is True and audit['status']=='passed',audit
                assert audit['minimum_inter_net_spacing_um']==gap
                assert audit['minimum_metal_support_margin_um']==margin
                assert audit['minimum_island_spacing_um']==gap
                assert audit['minimum_electrode_diameter_um']==20
                assert audit['minimum_wire_width_um']==2
                assert result['capacity_interval']['integer_polygon_lower_verified']
                assert result['capacity_interval']['lower_bound']==routing['retained_routes']
                records.append({'case':tag,'rules':asdict(rules),
                    'island_gap_um':routing['settings']['pad_gap_um'],
                    'connected_electrodes':routing['retained_routes'],
                    'roundtrip_passed':True,'integer_audit':audit,
                    'capacity_interval':result['capacity_interval']})
                print(f'Passed {tag}: {routing["retained_routes"]} complete electrode-to-Pad nets',flush=True)
                server.JOBS.pop(tag,None)
                server.RECENT_JOBS.pop(tag,None)
    finally:
        server.RUNS_DIR=original_runs
        for tag in ('sub_default','zero_gap_margin'):
            server.JOBS.pop(tag,None)
            server.RECENT_JOBS.pop(tag,None)
    report={'status':'passed','invalid_inputs':checks,'adaptive_ancillary_geometry':asdict(pads),
            'small_pad_reference_accepted':asdict(tiny),'zero_gap_short_conflicts_verified':True,
            'geometry_export_audit_cases':records}
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print('Editable rule validation, zero-gap conflicts and real exported-GDS audits passed',flush=True)


if __name__=='__main__':main()
