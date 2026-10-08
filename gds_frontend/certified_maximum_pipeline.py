"""End-to-end GDS construction with an honest, independently checked L/U.

The command never calls a candidate-library optimum a continuous maximum.
It starts from one audited report (or builds one), optionally adds audited
residual routes, then recomputes the candidate-independent capacity bound.
Only equality of the readback lower bound and continuous upper bound sets
`declared_geometric_model_optimality_proven` to true.
"""
from __future__ import annotations

from argparse import ArgumentParser
from dataclasses import asdict
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import subprocess
import sys

from pad_router import PadSettings
from process_geometry import ProcessRules
from outer_exit_policy import OUTER_BRIDGE_CONTACT_POLICY
from audit_per_net_substrate import audit as audit_per_net_substrate
from audit_bridge_contacts import audit_bridge_contacts
from audit_bridge_spine import (generate_bridge_spine_certificate,
                                verify_bridge_spine_certificate)


def _hash(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _run(*args):
    subprocess.run([sys.executable,'-B',str(Path(__file__).with_name(args[0])),
                    *map(str,args[1:])],check=True)


def _validated_report(path,source,source_hash):
    report=json.loads(Path(path).read_text(encoding='utf-8'))
    output=Path(report['output']).resolve(strict=True)
    if (Path(report['input']).resolve(strict=True)!=source or
            report['input_sha256']!=source_hash or
            _hash(output)!=report['output_sha256'] or
            report.get('integer_polygon_audit',{}).get('passed') is not True):
        raise ValueError(f'Report GDS/hash/audit mismatch: {path}')
    if report.get('integer_polygon_audit',{}).get('outer_exit_policy_verified') is not True:
        raise ValueError('Initial report lacks the current outer exit certificate; rerun base Pad routing')
    rules=ProcessRules(minimum_center_spacing_um=
                       report['rules']['minimum_center_spacing_um'])
    if report['rules']!=asdict(rules):
        raise ValueError('The independent capacity auditor supports default rules with user-selected center spacing')
    pads=PadSettings()
    if any(report['pad_settings'][name]!=getattr(pads,name)
           for name in ('pad_width_um','pad_length_um',
                        'minimum_inner_radius_um','source_to_shell_gap_um',
                        'pad_outer_setback_um','bridge_width_um',
                        'bridge_shell_overlap_um')):
        raise ValueError('The independent Pad audit requires the default physical Pad dimensions')
    if (report['pad_settings']['square_side_um']<pads.square_side_um or
            report['pad_settings']['pad_pitch_um']<
                pads.pad_width_um+rules.spacing_um):
        raise ValueError('Pad frame or Pad pitch is outside the declared policy')
    return report


def _normalize_legacy_report(path,destination):
    """Add explicit current defaults to older, already integer-audited reports."""
    report=json.loads(Path(path).read_text(encoding='utf-8'))
    if all(key in report for key in ('rules','pad_settings','support_layer')):
        return Path(path)
    audit=report.get('integer_polygon_audit',{})
    layers=report.get('audit',{}).get('output_layers',{})
    frame=report.get('pad_frame',{})
    if (audit.get('passed') is not True or
            audit.get('minimum_electrode_diameter_um')!=30.0 or
            audit.get('minimum_wire_width_um')!=5.0 or
            audit.get('minimum_inter_net_spacing_um')!=4.0 or
            audit.get('minimum_metal_support_margin_um')!=4.0 or
            audit.get('minimum_pad_short_side_um')!=500.0 or
            audit.get('minimum_pad_long_side_um')!=3000.0 or
            audit.get('minimum_substrate_disk_radius_um')!=19.05 or
            frame.get('square_side_um',0)<32000 or
            frame.get('pad_pitch_um',0)<504 or
            'support_layer' not in layers):
        raise ValueError('Historical report cannot be normalized to current physical rules')
    center_distance=float(audit['minimum_electrode_center_distance_um'])
    # The independent audit records the effective separation. For a legacy
    # report, use that value as the declared center rule. This is conservative
    # for the same output and binds the subsequent upper-bound computation.
    report['rules']=asdict(ProcessRules(
        minimum_center_spacing_um=center_distance))
    report['pad_settings']=asdict(PadSettings(
        square_side_um=frame['square_side_um'],
        pad_pitch_um=frame['pad_pitch_um']))
    report['support_layer']=layers['support_layer']
    report['normalization_provenance']={
        'legacy_report':str(Path(path).resolve()),
        'legacy_report_sha256':_hash(path),
        'inferred_center_rule_from_integer_audit_um':center_distance}
    normalized=destination/'normalized_initial_report.json'
    normalized.write_text(json.dumps(report,ensure_ascii=False,indent=2),
                          encoding='utf-8')
    return normalized


def main():
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--initial-report',type=Path,
                        help='Reuse an audited layout instead of solving a fresh one')
    parser.add_argument('--center-spacing-um',type=float,default=70)
    parser.add_argument('--support-layer',type=int,default=10)
    parser.add_argument('--support-datatype',type=int,default=0)
    parser.add_argument('--candidate-step-um',type=float,default=600)
    parser.add_argument('--external-portal-pitch-um',type=float,default=40)
    parser.add_argument('--external-portal-mode',
                        choices=('radial_envelope','full_boundary'),
                        default='radial_envelope')
    parser.add_argument('--pad-columns-per-pad',type=int,default=12)
    parser.add_argument('--milp-time-limit-s',type=float,default=None,
                        help='Optional explicit offline limit; omitted means no MILP wall-clock limit')
    parser.add_argument('--residual-iterations',type=int,default=0)
    parser.add_argument('--residual-max-sources',type=int,default=200)
    parser.add_argument('--residual-target-trials',type=int,default=1)
    parser.add_argument('--residual-external-portal-pitch-um',type=float,
                        default=9)
    parser.add_argument('--residual-external-portal-mode',
                        choices=('radial_envelope','full_boundary'),
                        default='radial_envelope')
    parser.add_argument('--residual-anchor-order',choices=(
                        'farthest_existing','closest_eligible_exit'),
                        default='closest_eligible_exit')
    args=parser.parse_args()
    if (args.residual_iterations<0 or args.residual_max_sources<1 or
            args.residual_target_trials<1):
        parser.error('Residual iterations must be nonnegative; source and target trial limits positive')
    source=args.input.resolve(strict=True)
    source_hash=_hash(source)
    destination=args.output_dir.resolve()
    destination.mkdir(parents=True,exist_ok=True)
    if args.initial_report is None:
        solver_limit_args=(['--milp-time-limit-s',args.milp_time_limit_s]
                           if args.milp_time_limit_s is not None else [])
        _run('prove_maximum.py','--input',source,'--output-dir',
             destination/'initial','--center-spacing-um',
             args.center_spacing_um,'--support-layer',args.support_layer,
             '--support-datatype',args.support_datatype,
             '--candidate-step-um',args.candidate_step_um,
             '--external-portal-pitch-um',args.external_portal_pitch_um,
             '--external-portal-mode',args.external_portal_mode,
             '--pad-columns-per-pad',args.pad_columns_per_pad,
             *solver_limit_args,
             '--allow-short-exit-connections',
             '--include-guarded-outlet-alternative')
        report_path=destination/'initial'/source.stem/'report.json'
    else:
        report_path=args.initial_report.resolve(strict=True)
    report_path=_normalize_legacy_report(report_path,destination)
    report=_validated_report(report_path,source,source_hash)
    initial_count=report['connected_pad_count']
    search_files=[]
    if args.residual_iterations:
        residual_dir=destination/'residual'
        _run('residual_pad_augment.py','--report',report_path,
             '--output-dir',residual_dir,'--iterations',
             args.residual_iterations,'--max-sources',
             args.residual_max_sources,'--external-portal-pitch-um',
             args.residual_external_portal_pitch_um,
             '--external-portal-mode',args.residual_external_portal_mode,
             '--anchor-order',args.residual_anchor_order,
             '--target-trials',args.residual_target_trials)
        for result_path in sorted(residual_dir.glob('iteration_*/residual_search_result.json')):
            result=json.loads(result_path.read_text(encoding='utf-8'))
            search_files.append(str(result_path))
            if result['status']!='audited_additional_network':
                break
            next_path=Path(result['report']).resolve(strict=True)
            next_report=_validated_report(next_path,source,source_hash)
            if next_report['connected_pad_count']<=report['connected_pad_count']:
                raise RuntimeError('Residual construction did not improve the audited lower bound')
            report_path,report=next_path,next_report
    capacity=destination/'capacity_certificate.json'
    _run('run_capacity_audit.py','--input-dir',source.parent,
         '--input-names',source.name,'--output',capacity,
         '--layer',report['support_layer'][0],
         '--datatype',report['support_layer'][1],
         '--center-spacing-um',
         report['rules']['minimum_center_spacing_um'],
         '--reference-reports',report_path)
    record=json.loads(capacity.read_text(encoding='utf-8'))['records'][0]
    lower=record['strict_integer_electrode_lower_bound']
    upper=record['continuous_geometric_upper_bound']
    if (record['input_sha256']!=source_hash or
            lower!=report['connected_pad_count'] or
            record['strict_integer_electrode_audit']['output_sha256']!=
                report['output_sha256'] or
            upper<lower):
        raise RuntimeError('Independent GDS audit and capacity certificate disagree')
    proven=bool(record['declared_geometric_model_optimality_proven'])
    if proven!=(lower==upper):
        raise RuntimeError('Optimality flag is inconsistent with the proven interval')
    own_support=destination/'per_net_substrate_audit.json'
    own_audit=audit_per_net_substrate(report_path)
    if (own_audit['passed'] is not True or
            own_audit['network_count']!=lower or
            own_audit['output_sha256']!=report['output_sha256']):
        raise RuntimeError('Per-net substrate provenance audit disagrees with the GDS')
    own_support.write_text(json.dumps(own_audit,ensure_ascii=False,indent=2),
                           encoding='utf-8')
    bridge_contact_path=destination/'bridge_contact_diagnostic.json'
    bridge_contacts=audit_bridge_contacts(report_path)
    if (bridge_contacts['network_count']!=lower or
            bridge_contacts['input_sha256']!=source_hash or
            bridge_contacts['output_sha256']!=report['output_sha256']):
        raise RuntimeError('Bridge-contact audit disagrees with selected GDS')
    bridge_contact_path.write_text(json.dumps(
        bridge_contacts,ensure_ascii=False,indent=2),encoding='utf-8')
    bridge_spine_path=destination/'bridge_spine_certificate.json'
    bridge_spine=generate_bridge_spine_certificate(report_path)
    if (bridge_spine['network_count']!=lower or
            not verify_bridge_spine_certificate(report_path,bridge_spine)):
        raise RuntimeError('Bridge spine width witness failed independent replay')
    bridge_spine_path.write_text(json.dumps(
        bridge_spine,ensure_ascii=False,indent=2),encoding='utf-8')
    summary={
        'created_at':datetime.now(timezone.utc).isoformat(),
        'input':str(source),'input_sha256':source_hash,
        'rules':report['rules'],'physical_pad_policy':{
            key:value for key,value in report['pad_settings'].items()
            if key not in ('square_side_um','pad_pitch_um','columns_per_pad',
                           'deep_extra_columns_per_pad','center_extra_columns_per_pad',
                           'pad_choices_per_portal')},
        'initial_connected_pad_count':initial_count,
        'certified_lower_bound':lower,
        'continuous_upper_bound':upper,
        'declared_geometric_model_optimality_proven':proven,
        'full_manufacturing_drc_proven':False,
        'selected_report':str(report_path),
        'selected_routing_gds':report['output'],
        'selected_output_sha256':report['output_sha256'],
        'capacity_certificate':str(capacity),
        'capacity_certificate_sha256':_hash(capacity),
        'per_net_substrate_audit':str(own_support),
        'per_net_substrate_audit_sha256':_hash(own_support),
        'per_net_substrate_provenance_verified':True,
        'bridge_contact_diagnostic':str(bridge_contact_path),
        'bridge_contact_diagnostic_sha256':_hash(bridge_contact_path),
        'multiple_disconnected_source_contact_count':bridge_contacts[
            'multiple_disconnected_source_contact_count'],
        'single_outlet_contact_proven':False,
        'declared_bridge_contact_policy':OUTER_BRIDGE_CONTACT_POLICY,
        'outer_exit_policy':report['audit']['output_layers']['outer_exit_policy'],
        'outer_exit_policy_verified':True,
        'bridge_spine_certificate':str(bridge_spine_path),
        'bridge_spine_certificate_sha256':_hash(bridge_spine_path),
        'certified_continuous_bridge_width_um':bridge_spine[
            'certified_continuous_bridge_width_um'],
        'bridge_spine_width_verified':True,
        'residual_search_results':search_files,
        'proof_rule':'maximum is claimable only when independent exported-GDS lower bound equals candidate-independent continuous upper bound'}
    summary_path=destination/'proof_result.json'
    summary_path.write_text(json.dumps(summary,ensure_ascii=False,indent=2),
                            encoding='utf-8')
    print(json.dumps({'lower':lower,'upper':upper,'maximum_proven':proven,
                      'proof_result':str(summary_path)},
                     ensure_ascii=False,indent=2),flush=True)


if __name__=='__main__':
    main()
