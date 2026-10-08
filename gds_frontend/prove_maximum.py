"""Construct and certify a maximum-electrode claim for an arbitrary GDS.

The command never substitutes a finite candidate optimum for a continuous
upper bound.  It emits either a closed, independently audited interval or an
honest open interval.  All added support and routing are exported as GDS.
"""
from __future__ import annotations

from argparse import ArgumentParser
from dataclasses import asdict
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import math
import subprocess
import sys
import time

from exact_gds_audit import audit_exported_gds
from island_router import IslandSettings, build_navigation, preview_navigation
from pad_router import PadSettings, audit_pad_gds, solve_with_pads, write_pad_gds
from process_geometry import ProcessRules
from frontend import read_support
from outer_exit_policy import make_outer_exit_policy


def main():
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--support-layer',type=int,default=10)
    parser.add_argument('--support-datatype',type=int,default=0)
    parser.add_argument('--center-spacing-um',type=float,default=70)
    parser.add_argument('--candidate-step-um',type=float,default=600)
    parser.add_argument('--routes-per-candidate',type=int,default=1)
    parser.add_argument('--lane-half-steps',type=int,default=1)
    parser.add_argument('--external-portal-pitch-um',type=float)
    parser.add_argument('--external-portal-mode',
                        choices=('radial_envelope','full_boundary'),
                        default='radial_envelope')
    parser.add_argument('--pad-columns-per-pad',type=int,default=12)
    parser.add_argument('--minimum-pad-square-side-um',type=float,default=32000,
                        help='Starting Pad-frame side length; a search choice, not a fixed Pad-count limit')
    parser.add_argument('--pad-pitch-um',type=float,default=1000,
                        help='Center pitch of contiguous Pads; must preserve the required gap')
    parser.add_argument('--deep-pad-columns-per-pad',type=int,default=0)
    parser.add_argument('--milp-time-limit-s',type=float,default=None,
                        help='Optional explicit offline limit; omitted means no MILP wall-clock limit')
    parser.add_argument('--skip-auxiliary-view',action='store_true',
                        help='Compatibility option; Pad routing now uses one full-source outer-port view')
    parser.add_argument('--allow-short-exit-connections',action='store_true',
                        help='Retain short source-to-exit paths supported by an electrode island')
    parser.add_argument('--include-guarded-outlet-alternative',action='store_true',
                        help='Search both the nearest outlet and a farther outlet beyond the island clearance')
    args=parser.parse_args()
    source=args.input.resolve(strict=True)
    output_root=args.output_dir.resolve()
    output=output_root/source.stem
    output.mkdir(parents=True,exist_ok=True)
    rules=ProcessRules(minimum_center_spacing_um=args.center_spacing_um)
    rules.validate()
    settings=IslandSettings(candidate_step_um=args.candidate_step_um,
                            routes_per_candidate=args.routes_per_candidate,
                            lane_half_steps=args.lane_half_steps,
                            external_portal_pitch_um=args.external_portal_pitch_um,
                            external_portal_mode=args.external_portal_mode,
                            navigation_component_exclusion_distance_um=(
                                0.0 if args.allow_short_exit_connections else None),
                            navigation_prune_leaf_length_um=(
                                0.0 if args.allow_short_exit_connections else None),
                            anchor_terminal_search_gap_um=(
                                0.0 if args.allow_short_exit_connections else None),
                            include_guarded_outlet_alternative=(
                                args.include_guarded_outlet_alternative),
                            milp_time_limit_s=args.milp_time_limit_s)
    if (settings.candidate_step_um<=0 or settings.routes_per_candidate<1 or
            settings.lane_half_steps<0 or
            (settings.milp_time_limit_s is not None and
             (not math.isfinite(settings.milp_time_limit_s) or settings.milp_time_limit_s<=0)) or
            (settings.external_portal_pitch_um is not None and
             settings.external_portal_pitch_um<=0)):
        parser.error('Search resolution, route count and time limit must be valid')
    if (args.pad_columns_per_pad<1 or args.deep_pad_columns_per_pad<0 or
            args.minimum_pad_square_side_um<PadSettings().square_side_um or
            args.pad_pitch_um<PadSettings().pad_width_um+rules.spacing_um):
        parser.error('Pad columns must be valid, frame at least the reference minimum, and pitch must meet Pad spacing')
    pads=PadSettings(square_side_um=args.minimum_pad_square_side_um,
                     pad_pitch_um=args.pad_pitch_um,
                     columns_per_pad=args.pad_columns_per_pad,
                     deep_extra_columns_per_pad=args.deep_pad_columns_per_pad)
    started=time.perf_counter()
    def progress(stage,value):
        print(f'{time.perf_counter()-started:.1f}s {value:.0%} {stage}',flush=True)
    support,_,metadata=read_support(source,args.support_layer,args.support_datatype)
    policy=make_outer_exit_policy(support,wire_width_um=rules.wire_width_um,
        margin_um=rules.margin_um,numeric_guard_um=settings.numeric_guard_um,
        bridge_width_um=pads.bridge_width_um,
        grid_um=metadata['gds_native_precision_m']*1e6)
    nav=build_navigation(source,rules,layer=args.support_layer,
                         datatype=args.support_datatype,
                         outlet_mode='external_boundary',settings=settings,
                         progress=progress,outer_exit_policy=policy)
    routing=solve_with_pads(nav,rules,settings=settings,pad_settings=pads,
                            progress=progress)
    report={'created_at':datetime.now(timezone.utc).isoformat(),
            'input':str(source),'input_sha256':sha256(source.read_bytes()).hexdigest(),
            'support_layer':[args.support_layer,args.support_datatype],
            'rules':asdict(rules),'island_settings':asdict(settings),
            'pad_settings':asdict(pads),
            'outer_exit_policy':routing['outer_exit_policy'],
            'connected_pad_count':routing['retained_routes'],
            'continuous_upper_bound':routing['continuous_upper_bound'],
            'optimization':routing['optimization'],
            'geometry_summary':nav.summary,
            'candidate_library':routing['candidate_library']}
    if routing['retained_routes']:
        gds=output/'routing.gds'
        layers=write_pad_gds(source,gds,routing,
                             (args.support_layer,args.support_datatype))
        roundtrip=audit_pad_gds(gds,nav,routing,layers)
        island_radius=(rules.electrode_diameter_um/2+rules.margin_um+
                       settings.numeric_guard_um)
        gap=max(rules.minimum_center_spacing_um,
                rules.electrode_diameter_um+rules.spacing_um,
                2*island_radius+settings.pad_gap_um)
        exact=audit_exported_gds(
            source,gds,layers,wire_spacing_um=rules.spacing_um,
            metal_support_margin_um=rules.margin_um,
            expected_nets=routing['retained_routes'],
            minimum_electrode_diameter_um=rules.electrode_diameter_um,
            minimum_electrode_center_distance_um=gap,
            minimum_substrate_disk_radius_um=island_radius,
            maximum_substrate_disk_radius_um=island_radius+.01,
            minimum_island_spacing_um=settings.pad_gap_um,
            minimum_pad_short_side_um=pads.pad_width_um,
            minimum_pad_long_side_um=pads.pad_length_um,
            minimum_wire_width_um=rules.wire_width_um)
        preview_navigation(output/'routing.png',nav,routing)
        report.update({'output':str(gds.resolve()),
                       'output_sha256':sha256(gds.read_bytes()).hexdigest(),
                       'audit':roundtrip,'integer_polygon_audit':exact,
                       'pad_bank_occupancy':routing['pad_bank_occupancy']})
    report['elapsed_seconds']=time.perf_counter()-started
    report_path=output/'report.json'
    report_path.write_text(json.dumps(report,ensure_ascii=False,indent=2),
                           encoding='utf-8')
    audit_path=output_root/(source.stem+'_capacity_audit.json')
    command=[sys.executable,'-B',str(Path(__file__).with_name('run_capacity_audit.py')),
             '--input-dir',str(source.parent),'--input-names',source.name,
             '--output',str(audit_path),'--layer',str(args.support_layer),
             '--datatype',str(args.support_datatype),
             '--center-spacing-um',str(args.center_spacing_um)]
    if routing['retained_routes']:
        command.extend(('--reference-reports',str(report_path)))
    subprocess.run(command,check=True)
    certificate=json.loads(audit_path.read_text(encoding='utf-8'))['records'][0]
    result={'input':str(source),
            'rules':asdict(rules),
            'constructed_connected_pad_count':routing['retained_routes'],
            'certified_lower_bound':certificate['strict_integer_electrode_lower_bound'],
            'continuous_upper_bound':certificate['continuous_geometric_upper_bound'],
            'declared_geometric_model_optimality_proven':certificate.get(
                'declared_geometric_model_optimality_proven',False),
            'full_manufacturing_drc_proven':False,
            'report':str(report_path),
            'capacity_certificate':str(audit_path),
            'routing_gds':report.get('output'),
            'wire_width_witness':(layers.get('wire_width_witness_path')
                                  if routing['retained_routes'] else None)}
    (output_root/(source.stem+'_result.json')).write_text(
        json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2),flush=True)


if __name__=='__main__':main()
