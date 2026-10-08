"""Replay a proof bundle from source GDS and exported metal GDS.

This verifier does not trust the saved lower/upper values. It re-reads the
selected report and both GDS files, runs the integer network audit plus all
continuous bounds again, and compares the saved claim with fresh evidence.
"""
from __future__ import annotations

from argparse import ArgumentParser
from hashlib import sha256
from pathlib import Path
from tempfile import TemporaryDirectory
import json
import subprocess
import sys

from frontend import read_support
from pad_router import PadSettings, make_pad_frame
from process_geometry import ProcessRules
from outer_exit_policy import OUTER_BRIDGE_CONTACT_POLICY
from audit_per_net_substrate import audit as audit_per_net_substrate
from audit_bridge_contacts import audit_bridge_contacts
from audit_bridge_spine import verify_bridge_spine_certificate
from residual_pad_augment import (_grouped_polygons,
                                  _occupancy_record,
                                  _readback_pad_occupancy)
from dataclasses import asdict


def _digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def verify_proof_result(path):
    summary=json.loads(Path(path).read_text(encoding='utf-8'))
    bridge_policy=OUTER_BRIDGE_CONTACT_POLICY
    if summary.get('outer_exit_policy_verified') is not True:
        raise ValueError('Proof bundle predates the mandatory outer exit policy; recompute it under the current model')
    if (('bridge_contact_diagnostic' in summary and
         summary.get('declared_bridge_contact_policy')!=bridge_policy) or
            ('declared_bridge_contact_policy' in summary and
             summary['declared_bridge_contact_policy']!=bridge_policy)):
        raise ValueError('Bridge contact policy differs from verified model')
    if summary.get('full_manufacturing_drc_proven') is not False:
        raise ValueError('Full manufacturing DRC is not certified by this proof bundle')
    source=Path(summary['input']).resolve(strict=True)
    report_path=Path(summary['selected_report']).resolve(strict=True)
    report=json.loads(report_path.read_text(encoding='utf-8'))
    output=Path(summary['selected_routing_gds']).resolve(strict=True)
    saved_certificate=Path(summary['capacity_certificate']).resolve(strict=True)
    if (_digest(source)!=summary['input_sha256'] or
            _digest(output)!=summary['selected_output_sha256'] or
            _digest(saved_certificate)!=summary['capacity_certificate_sha256'] or
            Path(report['input']).resolve(strict=True)!=source or
            Path(report['output']).resolve(strict=True)!=output or
            report['input_sha256']!=summary['input_sha256'] or
            report['output_sha256']!=summary['selected_output_sha256'] or
            report['integer_polygon_audit']['passed'] is not True):
        raise ValueError('Proof bundle input, output, report or certificate hash differs')
    rules=ProcessRules(minimum_center_spacing_um=
                       report['rules']['minimum_center_spacing_um'])
    if report['rules']!=asdict(rules) or summary['rules']!=report['rules']:
        raise ValueError('Proof bundle rule set differs from supported process model')
    default_pads=PadSettings()
    physical=summary['physical_pad_policy']
    if any(physical.get(name)!=getattr(default_pads,name)
           for name in ('minimum_inner_radius_um',
                        'source_to_shell_gap_um','pad_width_um',
                        'pad_length_um','pad_outer_setback_um',
                        'bridge_width_um','bridge_shell_overlap_um')):
        raise ValueError('Proof bundle has an unverified Pad policy')
    if any(report['pad_settings'][name]!=physical[name] for name in physical):
        raise ValueError('Report and summary Pad rules disagree')
    settings=PadSettings(**report['pad_settings'])
    layers=report['audit']['output_layers']
    if (summary.get('outer_exit_policy')!=layers.get('outer_exit_policy') or
            report['integer_polygon_audit'].get('outer_exit_policy_verified') is not True):
        raise ValueError('Proof outer exit policy differs from the audited GDS report')
    groups=_grouped_polygons(output)
    shell=groups.get((layers['shell_marker_layer'],0))
    if shell is None or shell.is_empty:
        raise ValueError('Exported GDS has no Pad shell marker')
    support=read_support(source,report['support_layer'][0],
                         report['support_layer'][1])[0]
    side=float(shell.bounds[2]-shell.bounds[0])
    frame=make_pad_frame(support,settings,minimum_side_um=side,
                         spacing_um=rules.spacing_um,
                         margin_um=rules.margin_um)
    if abs(frame['square_side_um']-side)>2*report['integer_polygon_audit']['native_grid_um']:
        raise ValueError('Exported Pad frame has the wrong square size')
    bank=_readback_pad_occupancy(
        groups,frame,layers,report['connected_pad_count'],
        report['integer_polygon_audit']['native_grid_um'])
    record=_occupancy_record(bank,settings.pad_pitch_um,
                             settings.pad_width_um)
    if record!=report.get('pad_bank_occupancy'):
        raise ValueError('Exported GDS Pad bank differs from report')
    saved=json.loads(saved_certificate.read_text(encoding='utf-8'))
    if len(saved.get('records',[]))!=1:
        raise ValueError('Expected one saved capacity record')
    with TemporaryDirectory(prefix='proof_replay_') as temporary:
        fresh_path=Path(temporary)/'recomputed_capacity.json'
        command=[sys.executable,'-B',
                 str(Path(__file__).with_name('run_capacity_audit.py')),
                 '--input-dir',str(source.parent),'--input-names',source.name,
                 '--output',str(fresh_path),
                 '--layer',str(report['support_layer'][0]),
                 '--datatype',str(report['support_layer'][1]),
                 '--center-spacing-um',str(rules.minimum_center_spacing_um),
                 '--reference-reports',str(report_path)]
        subprocess.run(command,check=True,capture_output=True,text=True)
        fresh=json.loads(fresh_path.read_text(encoding='utf-8'))
    if len(fresh.get('records',[]))!=1:
        raise RuntimeError('Recomputed capacity has the wrong number of records')
    old_record=saved['records'][0]
    record=fresh['records'][0]
    for name in ('input_sha256','strict_integer_electrode_lower_bound',
                 'continuous_geometric_upper_bound',
                 'declared_geometric_model_optimality_proven'):
        if record[name]!=old_record[name]:
            raise ValueError(f'Recomputed {name} differs from the saved certificate')
    lower=record['strict_integer_electrode_lower_bound']
    upper=record['continuous_geometric_upper_bound']
    if (lower is None or lower>upper or
            record['strict_integer_electrode_audit']['output_sha256']!=
                summary['selected_output_sha256'] or
            lower!=summary['certified_lower_bound'] or
            upper!=summary['continuous_upper_bound'] or
            summary['declared_geometric_model_optimality_proven'] is not
                (lower==upper) or
            record['declared_geometric_model_optimality_proven'] is not
                (lower==upper)):
        raise ValueError('Recomputed optimality interval differs from the saved claim')
    own_verified=False
    own_path=summary.get('per_net_substrate_audit')
    if own_path is not None:
        own_path=Path(own_path).resolve(strict=True)
        if (_digest(own_path)!=summary.get('per_net_substrate_audit_sha256') or
                summary.get('per_net_substrate_provenance_verified') is not True):
            raise ValueError('Per-net substrate audit file or claim differs')
        saved_own=json.loads(own_path.read_text(encoding='utf-8'))
        fresh_own=audit_per_net_substrate(report_path)
        for name in ('input_sha256','output_sha256','native_grid_um',
                     'network_count','minimum_own_support_margin_um',
                     'nearby_integer_segment_pairs_checked','passed'):
            if saved_own.get(name)!=fresh_own.get(name):
                raise ValueError(f'Per-net substrate audit {name} differs')
        own_verified=True
    elif summary.get('per_net_substrate_provenance_verified') is not None:
        raise ValueError('Per-net substrate claim has no audit file')
    contacts_path=summary.get('bridge_contact_diagnostic')
    contacts_verified=False
    if contacts_path is not None:
        contacts_path=Path(contacts_path).resolve(strict=True)
        if (_digest(contacts_path)!=summary.get(
                'bridge_contact_diagnostic_sha256') or
                summary.get('single_outlet_contact_proven') is not False):
            raise ValueError('Bridge contact diagnostic file or scope differs')
        saved_contacts=json.loads(contacts_path.read_text(encoding='utf-8'))
        fresh_contacts=audit_bridge_contacts(report_path)
        for name in ('input_sha256','output_sha256','native_grid_um',
                     'network_count',
                     'multiple_disconnected_source_contact_networks',
                     'multiple_disconnected_source_contact_count',
                     'all_bridges_have_one_connected_contact_footprint',
                     'single_outlet_contact_proven','per_network'):
            if saved_contacts.get(name)!=fresh_contacts.get(name):
                raise ValueError(f'Bridge contact diagnostic {name} differs')
        if (summary.get('multiple_disconnected_source_contact_count')!=
                fresh_contacts['multiple_disconnected_source_contact_count']):
            raise ValueError('Bridge contact count differs from proof summary')
        contacts_verified=True
    elif (summary.get('bridge_contact_diagnostic_sha256') is not None or
          summary.get('multiple_disconnected_source_contact_count') is not None or
          summary.get('single_outlet_contact_proven') is not None):
        raise ValueError('Bridge contact claim has no audit file')
    spine_verified=False
    spine_path=summary.get('bridge_spine_certificate')
    if spine_path is not None:
        spine_path=Path(spine_path).resolve(strict=True)
        if (_digest(spine_path)!=summary.get('bridge_spine_certificate_sha256') or
                summary.get('bridge_spine_width_verified') is not True):
            raise ValueError('Bridge spine witness hash or claim differs')
        spine=json.loads(spine_path.read_text(encoding='utf-8'))
        if (not verify_bridge_spine_certificate(report_path,spine) or
                spine['network_count']!=lower or
                spine['certified_continuous_bridge_width_um']!=
                summary.get('certified_continuous_bridge_width_um')):
            raise ValueError('Bridge spine width witness failed replay')
        spine_verified=True
    elif (summary.get('bridge_spine_certificate_sha256') is not None or
          summary.get('bridge_spine_width_verified') is not None or
          summary.get('certified_continuous_bridge_width_um') is not None):
        raise ValueError('Bridge spine width claim has no certificate')
    return {'input_sha256':summary['input_sha256'],
            'output_sha256':summary['selected_output_sha256'],
            'certified_lower_bound':lower,
            'continuous_upper_bound':upper,
            'declared_geometric_model_optimality_proven':lower==upper,
            'full_manufacturing_drc_proven':False,
            'four_side_contiguous_pad_bank_verified':True,
            'outer_exit_policy_verified':True,
            'outer_exit_policy':summary['outer_exit_policy'],
            'per_net_substrate_provenance_verified':own_verified,
            'bridge_contact_diagnostic_replayed':contacts_verified,
            'multiple_disconnected_source_contact_count':(
                summary.get('multiple_disconnected_source_contact_count')
                if contacts_verified else None),
            'bridge_spine_width_verified':spine_verified,
            'certified_continuous_bridge_width_um':(
                summary.get('certified_continuous_bridge_width_um')
                if spine_verified else None),
            'proof_bundle_replayed':True}


def main():
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--proof-result',type=Path,required=True)
    args=parser.parse_args()
    print(json.dumps(verify_proof_result(args.proof_result),
                     ensure_ascii=False,indent=2),flush=True)


if __name__=='__main__':
    main()
