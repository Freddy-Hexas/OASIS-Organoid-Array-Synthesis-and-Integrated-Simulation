"""Keep the best independently audited construction across search settings.

Every input result must refer to the same source GDS, support layer and
physical rules. The chosen GDS is re-audited by run_capacity_audit.py, which
also recomputes its candidate-independent continuous upper bounds. Search
resolution and solver limits never enter the upper-bound proof.
"""
from __future__ import annotations

from argparse import ArgumentParser
from dataclasses import asdict
from hashlib import sha256
from pathlib import Path
import json
import subprocess
import sys

from pad_router import PadSettings
from process_geometry import ProcessRules


def _digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _physical_pad_policy(settings):
    search_only={'square_side_um','pad_pitch_um','columns_per_pad',
                 'deep_extra_columns_per_pad','center_extra_columns_per_pad',
                 'pad_choices_per_portal'}
    return {key:value for key,value in settings.items()
            if key not in search_only}


def _load_verified_claim(result_path,source,source_hash):
    result=json.loads(result_path.read_text(encoding='utf-8'))
    report_path=Path(result['report']).resolve(strict=True)
    report=json.loads(report_path.read_text(encoding='utf-8'))
    if (Path(result['input']).resolve(strict=True)!=source or
            Path(report['input']).resolve(strict=True)!=source or
            report['input_sha256']!=source_hash or
            result['rules']!=report['rules']):
        raise ValueError(f'Source or process rules differ: {result_path}')
    output=Path(report['output']).resolve(strict=True)
    if _digest(output)!=report['output_sha256']:
        raise ValueError(f'Exported GDS hash differs: {result_path}')
    if (not isinstance(result.get('certified_lower_bound'),int) or
            result['certified_lower_bound']<0 or
            result['certified_lower_bound']!=report['connected_pad_count'] or
            report.get('integer_polygon_audit',{}).get('passed') is not True):
        raise ValueError(f'Result has no audited constructed lower bound: {result_path}')
    if report['integer_polygon_audit'].get('outer_exit_policy_verified') is not True:
        raise ValueError(f'Result belongs to an older exit model; rerun base routing: {result_path}')
    return {'result_file':str(result_path),'result_sha256':_digest(result_path),
            'report':str(report_path),'output':str(output),
            'output_sha256':report['output_sha256'],
            'count':result['certified_lower_bound'],
            'rules':result['rules'],
            'support_layer':report['support_layer'],
            'search_settings':report.get('island_settings'),
            'pad_settings':report.get('pad_settings'),
            'physical_pad_policy':_physical_pad_policy(report['pad_settings'])}


def main():
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--result',type=Path,action='append',required=True,
                        help='Result JSON from prove_maximum.py; repeat to form a portfolio')
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args()
    source=args.input.resolve(strict=True)
    source_hash=_digest(source)
    claims=[_load_verified_claim(p.resolve(strict=True),source,source_hash)
            for p in args.result]
    physical_rules=claims[0]['rules']
    layer=claims[0]['support_layer']
    pad_policy=claims[0]['physical_pad_policy']
    if any(c['rules']!=physical_rules or c['support_layer']!=layer or
           c['physical_pad_policy']!=pad_policy for c in claims):
        parser.error('Portfolio claims use different physical rules, support layers or Pad policies')
    expected_rules=asdict(ProcessRules(
        minimum_center_spacing_um=physical_rules['minimum_center_spacing_um']))
    expected_pad_policy=_physical_pad_policy(asdict(PadSettings()))
    if physical_rules!=expected_rules or pad_policy!=expected_pad_policy:
        parser.error('The independent capacity auditor currently supports default physical and Pad rules with user-selected center spacing only')
    best=max(claims,key=lambda c:c['count'])
    destination=args.output_dir.resolve()
    destination.mkdir(parents=True,exist_ok=True)
    capacity_path=destination/'best_capacity_audit.json'
    command=[sys.executable,'-B',
             str(Path(__file__).with_name('run_capacity_audit.py')),
             '--input-dir',str(source.parent),'--input-names',source.name,
             '--output',str(capacity_path),'--layer',str(layer[0]),
             '--datatype',str(layer[1]),
             '--center-spacing-um',
             str(physical_rules['minimum_center_spacing_um']),
             '--reference-reports',best['report']]
    subprocess.run(command,check=True)
    record=json.loads(capacity_path.read_text(encoding='utf-8'))['records'][0]
    if (record['input_sha256']!=source_hash or
            record['strict_integer_electrode_lower_bound']!=best['count'] or
            record['strict_integer_electrode_audit']['output_sha256']!=
            best['output_sha256']):
        raise RuntimeError('Independent best-GDS audit differs from portfolio')
    summary={'input':str(source),'input_sha256':source_hash,
             'support_layer':layer,'rules':physical_rules,
             'physical_pad_policy':pad_policy,
             'searches':[{'result_file':c['result_file'],
                          'result_sha256':c['result_sha256'],
                          'certified_lower_bound':c['count'],
                          'search_settings':c['search_settings'],
                          'pad_search_settings':{
                              key:value for key,value in c['pad_settings'].items()
                              if key not in c['physical_pad_policy']}}
                         for c in claims],
             'selected_report':best['report'],
             'selected_routing_gds':best['output'],
             'selected_output_sha256':best['output_sha256'],
             'certified_lower_bound':best['count'],
             'continuous_upper_bound':
                 record['continuous_geometric_upper_bound'],
             'declared_geometric_model_optimality_proven':
                 record['declared_geometric_model_optimality_proven'],
             'full_manufacturing_drc_proven':False,
             'capacity_certificate':str(capacity_path),
             'capacity_certificate_sha256':_digest(capacity_path),
             'selection_rule':'maximum re-audited connected-Pad lower bound; finite-library optima are never continuous upper bounds'}
    output=destination/'portfolio_result.json'
    output.write_text(json.dumps(summary,ensure_ascii=False,indent=2),
                      encoding='utf-8')
    print(json.dumps({'certified_lower_bound':best['count'],
                      'continuous_upper_bound':summary['continuous_upper_bound'],
                      'declared_geometric_model_optimality_proven':
                          summary['declared_geometric_model_optimality_proven'],
                      'result':str(output)},ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
