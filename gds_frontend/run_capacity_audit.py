"""Recompute candidate-independent capacity certificates for any GDS files.

Example:
    python -B gds_frontend/run_capacity_audit.py --input-dir data --output outputs/capacity_audit.json
"""
from __future__ import annotations

from argparse import ArgumentParser
from hashlib import sha256
from pathlib import Path
from tempfile import TemporaryDirectory
import json
import math

import gdstk
from capacity_bounds import (gds_cell_cover_upper_bound,
                             gds_convex_hull_cut_upper_bound,
                             gds_enclosing_disk_angular_upper_bound,
                             verify_gds_angular_certificate,
                             gds_radial_partition_upper_bound,
                             gds_radial_cut_upper_bound,
                             gds_square_packing_upper_bound,
                             verify_pads_outside_hull_cut,
                             verify_pads_outside_radial_cut)
from multi_circle_capacity import (gds_multi_circle_upper_bound,
                                   verify_multi_circle_certificate,
                                   verify_pads_outside_multi_cuts)
from small_cell_conflicts import (gds_small_cell_conflict_upper_bound,
                                  verify_small_cell_conflict_certificate)
from frontend import read_support
from exact_gds_audit import audit_exported_gds
from island_router import IslandSettings
from pad_router import PadSettings, make_pad_frame
from process_geometry import ProcessRules
from outer_exit_policy import make_outer_exit_policy
from shapely.geometry import Point


def _universal_pad_exclusion(frame,settings,grid_um,circles):
    """Certify cuts for every Pad contained in the permitted outer shell.

    Any Pad in a valid frame is outside the shell's circular opening.  The
    minimum opening radius is a process rule, independent of Pad pitch,
    exact Pad dimensions and square expansion.  This is deliberately a
    broader claim than merely checking the current exported Pad markers.
    """
    opening_um=float(settings.minimum_inner_radius_um)
    if opening_um<=0 or frame['inner_radius_um']<opening_um:
        raise ValueError('Pad shell has no certified minimum inner opening')
    origin2=[round(2*coordinate/grid_um) for coordinate in frame['origin_um']]
    # Eight half-grid ticks cover input/output rounding and the frame-origin
    # conversion; a failed comparison merely discards the affected cut.
    opening2=math.floor(2*opening_um/grid_um)-8
    entries=[]
    for name,center2,radius2 in circles:
        dx=center2[0]-origin2[0]
        dy=center2[1]-origin2[1]
        offset2=math.isqrt(dx*dx+dy*dy)
        if offset2*offset2<dx*dx+dy*dy:
            offset2+=1
        outside=opening2>offset2+radius2+8
        entries.append({'cut':name,'center_twice_grid_ticks':list(center2),
                        'radius_twice_grid_ticks':int(radius2),
                        'center_offset_upper_twice_grid_ticks':offset2+8,
                        'all_allowed_pad_rectangles_outside':outside})
    return {'method':'minimum_shell_opening_integer_triangle_inequality',
            'minimum_frame_side_um':frame['square_side_um'],
            'minimum_shell_opening_radius_um':opening_um,
            'actual_reference_shell_opening_radius_um':frame['inner_radius_um'],
            'minimum_pad_radius_lower_twice_grid_ticks':opening2,
            'frame_origin_twice_grid_ticks':origin2,
            'scope':'all Pad rectangles wholly in a shell with the declared minimum circular opening, independent of Pad size, pitch, slot choice and square expansion',
            'cut_checks':entries}


def _verify_universal_pad_exclusion(certificate,frame,settings,grid_um,
                                    circles):
    """Recompute premises instead of trusting serialized cut booleans."""
    fresh=_universal_pad_exclusion(frame,settings,grid_um,circles)
    return fresh==certificate


def main():
    parser = ArgumentParser()
    parser.add_argument('--input-dir', type=Path, required=True)
    parser.add_argument('--input-glob', default='*.gds')
    parser.add_argument('--input-names', nargs='*', help='Optional exact GDS filenames')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--layer', type=int, default=10)
    parser.add_argument('--datatype', type=int, default=0)
    parser.add_argument('--center-spacing-um', type=float, default=70.0)
    parser.add_argument('--reference-report-dir', type=Path)
    parser.add_argument('--reference-reports', type=Path, nargs='*',
                        help='Explicit report.json files, one per input GDS stem')
    args = parser.parse_args()
    rules = ProcessRules(minimum_center_spacing_um=args.center_spacing_um)
    settings = IslandSettings()
    island_radius = (rules.electrode_diameter_um/2 + rules.margin_um +
                     settings.numeric_guard_um)
    gap = max(rules.minimum_center_spacing_um,
              rules.electrode_diameter_um+rules.spacing_um,
              2*island_radius+settings.pad_gap_um)
    records = []
    report_map = {}
    for report_path in args.reference_reports or ():
        declared = json.loads(report_path.read_text(encoding='utf-8'))
        if 'input' not in declared:
            raise ValueError(f'Reference report has no input GDS: {report_path}')
        key = Path(declared['input']).stem
        if key in report_map:
            raise ValueError(f'Duplicate report for {key}')
        report_map[key] = report_path
    paths = ([args.input_dir/name for name in args.input_names]
             if args.input_names else sorted(args.input_dir.glob(args.input_glob)))
    for path in paths:
        candidates = [gds_cell_cover_upper_bound(path, args.layer, args.datatype,
                                                  gap, offset_fraction=offset)
                      for offset in ((0, 0), (0, .5), (.5, 0), (.5, .5))]
        chosen = min(candidates, key=lambda result: result['value'])
        small_cells=(gds_small_cell_conflict_upper_bound(
            path,args.layer,args.datatype,gap)
            if chosen['cell_count']<=18 else None)
        if small_cells is not None and not verify_small_cell_conflict_certificate(path,small_cells):
            raise RuntimeError('Small-cell conflict certificate failed verification')
        angular = gds_enclosing_disk_angular_upper_bound(
            path, args.layer, args.datatype, gap)
        if angular is not None and angular['input_sha256'] != chosen['input_sha256']:
            raise RuntimeError('GDS changed while calculating angular bound')
        if angular is not None and not verify_gds_angular_certificate(path,angular):
            raise RuntimeError('Angular capacity certificate failed independent verification')
        square_envelope = gds_square_packing_upper_bound(
            path, args.layer, args.datatype, gap)
        if square_envelope['input_sha256'] != chosen['input_sha256']:
            raise RuntimeError('GDS changed while calculating packing bounds')
        cut = gds_radial_cut_upper_bound(path, args.layer, args.datatype,
                                         rules.wire_width_um, rules.spacing_um)
        support, _, support_meta = read_support(path, args.layer, args.datatype)
        pad_policy=PadSettings()
        frame = make_pad_frame(support, pad_policy, spacing_um=rules.spacing_um,
                               margin_um=rules.margin_um)
        partition = gds_radial_partition_upper_bound(
            path, args.layer, args.datatype, gap, rules.wire_width_um,
            rules.spacing_um, [slot['polygon'].bounds for slot in frame['slots']])
        multi_circle = gds_multi_circle_upper_bound(
            path, args.layer, args.datatype, gap, rules.wire_width_um,
            rules.spacing_um, [slot['polygon'].bounds for slot in frame['slots']])
        if multi_circle['input_sha256'] != chosen['input_sha256']:
            raise RuntimeError('GDS changed while calculating circle cuts')
        if not verify_multi_circle_certificate(path,multi_circle):
            raise RuntimeError('Multiple-circle dual certificate failed verification')
        origin = Point(cut['circle_center_um'])
        clearance = min(slot['polygon'].distance(origin)-cut['circle_radius_um']
                        for slot in frame['slots'])
        if not verify_pads_outside_radial_cut(
                cut,[slot['polygon'].bounds for slot in frame['slots']]):
            raise RuntimeError(f'Default Pad frame is not outside radial cut: {path}')
        cut['pad_outside_cut_verified'] = True
        cut['minimum_pad_to_cut_clearance_um'] = clearance
        cut['pad_exclusion_method'] = 'exact_integer_circle_vs_outward_rounded_pad_boxes'
        hull_cut = gds_convex_hull_cut_upper_bound(path, args.layer, args.datatype,
                                                   rules.wire_width_um,rules.spacing_um)
        hull_origin = Point(hull_cut['pad_exclusion_circle_center_um'])
        hull_clearance = min(slot['polygon'].distance(hull_origin)-
                             hull_cut['pad_exclusion_circle_radius_um']
                             for slot in frame['slots'])
        if not verify_pads_outside_hull_cut(
                hull_cut,[slot['polygon'].bounds for slot in frame['slots']]):
            raise RuntimeError(f'Default Pad frame is not outside hull cut: {path}')
        hull_cut['pad_outside_cut_verified'] = True
        hull_cut['minimum_pad_to_exclusion_circle_um'] = hull_clearance
        hull_cut['pad_exclusion_method'] = 'exact_integer_hull_vs_outward_rounded_pad_boxes'
        grid_um=cut['native_grid_um']
        hull_center2=[round(2*v/grid_um) for v in
                      hull_cut['pad_exclusion_circle_center_um']]
        hull_radius2=math.ceil(2*hull_cut[
            'pad_exclusion_circle_radius_um']/grid_um)+2
        circles=[('radial',cut['circle_center_twice_grid_ticks'],
                  cut['radius_twice_grid_ticks']),
                 ('radial_partition',partition['circle_center_twice_grid_ticks'],
                  partition['radius_twice_grid_ticks']),
                 ('convex_hull_enclosing_circle',hull_center2,hull_radius2)]
        circles.extend((f'multi_circle_{index}',item['center_twice_grid_ticks'],
                        item['radius_twice_grid_ticks'])
                       for index,item in enumerate(multi_circle['circle_cuts']))
        universal_pads=_universal_pad_exclusion(
            frame,pad_policy,grid_um,circles)
        if not _verify_universal_pad_exclusion(
                universal_pads,frame,pad_policy,grid_um,circles):
            raise RuntimeError('Universal Pad exclusion replay failed')
        valid_cuts={entry['cut']:entry['all_allowed_pad_rectangles_outside']
                    for entry in universal_pads['cut_checks']}
        eligible_bounds=[chosen['value'],square_envelope['value']]
        if angular:eligible_bounds.append(angular['value'])
        if small_cells:eligible_bounds.append(small_cells['value'])
        if valid_cuts['radial']:eligible_bounds.append(cut['value'])
        if valid_cuts['convex_hull_enclosing_circle']:
            eligible_bounds.append(hull_cut['value'])
        if valid_cuts['radial_partition']:
            eligible_bounds.append(partition['value'])
        if all(valid_cuts[f'multi_circle_{index}']
               for index in range(len(multi_circle['circle_cuts']))):
            eligible_bounds.append(multi_circle['value'])
        upper=min(eligible_bounds)
        lower = None
        strict_lower = None
        strict_lower_audit = None
        strict_lower_error = None
        finite_upper = None
        evidence = None
        if args.reference_report_dir or report_map:
            report = report_map.get(path.stem)
            if report is None and args.reference_report_dir:
                report = args.reference_report_dir/path.stem/'report.json'
            if report is not None and report.exists():
                previous = json.loads(report.read_text(encoding='utf-8'))
                output_path = Path(previous.get('output', ''))
                output_verified = (output_path.is_file() and
                                   sha256(output_path.read_bytes()).hexdigest() ==
                                   previous.get('output_sha256'))
                if (previous.get('input_sha256') == chosen['input_sha256'] and
                    previous.get('audit', {}).get('passed') is True and
                    previous.get('audit', {}).get('pad_connection_verified') is True and
                    output_verified):
                    with TemporaryDirectory(prefix='capacity_output_') as temporary:
                        local = Path(temporary)/'output.gds'
                        local.write_bytes(output_path.read_bytes())
                        output_lib = gdstk.read_gds(str(local),unit=1e-6)
                    pad_layer = previous['audit']['output_layers']['pad_contact_layer']
                    actual_pad_bounds = []
                    for top in output_lib.top_level():
                        for polygon in top.get_polygons():
                            if polygon.layer != pad_layer:
                                continue
                            points = polygon.points
                            actual_pad_bounds.append((float(points[:,0].min()),
                                                      float(points[:,1].min()),
                                                      float(points[:,0].max()),
                                                      float(points[:,1].max())))
                    if (len(actual_pad_bounds) < previous['connected_pad_count'] or
                            not verify_pads_outside_radial_cut(cut,actual_pad_bounds) or
                            not verify_pads_outside_hull_cut(hull_cut,actual_pad_bounds) or
                            not verify_pads_outside_radial_cut(partition,actual_pad_bounds) or
                            not verify_pads_outside_multi_cuts(multi_circle,actual_pad_bounds)):
                        raise RuntimeError(f'Output Pads are not certified outside both cuts: {path}')
                    lower = previous['connected_pad_count']
                    try:
                        expected_exit=make_outer_exit_policy(support,
                            wire_width_um=rules.wire_width_um,margin_um=rules.margin_um,
                            numeric_guard_um=settings.numeric_guard_um,
                            bridge_width_um=pad_policy.bridge_width_um,
                            grid_um=support_meta['gds_native_precision_m']*1e6)
                        if previous['audit']['output_layers'].get('outer_exit_policy')!=expected_exit.record():
                            raise ValueError('Lower-bound GDS has a missing or different outer exit rule set')
                        strict_lower_audit = audit_exported_gds(
                            path, output_path, previous['audit']['output_layers'],
                            wire_spacing_um=rules.spacing_um,
                            metal_support_margin_um=rules.margin_um,
                            expected_nets=lower,
                            minimum_electrode_diameter_um=rules.electrode_diameter_um,
                            minimum_electrode_center_distance_um=gap,
                            minimum_substrate_disk_radius_um=island_radius,
                            maximum_substrate_disk_radius_um=(island_radius+.01
                                                              if 'island_marker_layer' in previous['audit']['output_layers']
                                                              else None),
                            minimum_island_spacing_um=(settings.pad_gap_um
                                                       if 'island_marker_layer' in previous['audit']['output_layers']
                                                       else None),
                            minimum_pad_short_side_um=PadSettings().pad_width_um,
                            minimum_pad_long_side_um=PadSettings().pad_length_um,
                            minimum_wire_width_um=(rules.wire_width_um
                                                   if 'wire_width_witness_path' in previous['audit']['output_layers']
                                                   else None))
                        if (strict_lower_audit['input_sha256'] !=
                                previous['input_sha256'] or
                                strict_lower_audit['output_sha256'] !=
                                previous['output_sha256']):
                            raise ValueError('Independent audit input or output hash differs')
                        if strict_lower_audit.get('outer_exit_policy_verified') is not True:
                            raise ValueError('Historical layout is not certified under the mandatory outer-extremity exit policy; rerun Pad routing')
                        strict_lower = lower
                    except (ValueError,RuntimeError) as exc:
                        strict_lower_audit = None
                        strict_lower_error = f'{type(exc).__name__}: {exc}'
                    finite_upper = previous.get('optimization', {}).get('finite_upper_bound')
                    evidence = str(report.resolve())
        if lower is not None and lower > upper:
            raise RuntimeError(f'Constructed layout exceeds geometric bound: {path}')
        model_proven=bool(
            strict_lower is not None and strict_lower==upper and evidence and
            strict_lower_audit and
            all(strict_lower_audit.get(field) is True for field in (
                'electrode_centers_in_original_support_verified',
                'electrode_disks_and_center_spacing_verified',
                'substrate_disks_verified','island_envelopes_verified',
                'support_provenance_verified','source_island_topology_verified',
                'bridge_attachment_and_hole_exclusion_verified',
                'outer_exit_policy_verified',
                'island_island_spacing_verified',
                'pad_dimensions_verified','functional_wire_width_verified')))
        records.append({'input': str(path.resolve()),
                        'input_sha256': chosen['input_sha256'],
                        'connected_pad_lower_bound': lower,
                        'strict_integer_electrode_lower_bound': strict_lower,
                        'strict_integer_electrode_audit': strict_lower_audit,
                        'strict_integer_electrode_error': strict_lower_error,
                        'continuous_geometric_upper_bound': upper,
                        'finite_candidate_upper_bound': finite_upper,
                        'bound_gap_closed': (strict_lower == upper if strict_lower is not None else False),
                        'declared_geometric_model_optimality_proven':model_proven,
                        'proven_globally_optimal': False,
                        'optimality_note': ('The declared geometric model is globally optimal for this input and spacing; complete manufacturing DRC and extra bridge/shell process restrictions are outside the declared model' if model_proven else
                                            'The exact integer audit verifies original-support anchors, 30 um electrode disks, modeled center spacing, substrate disks, and rectangular Pad dimensions; the continuous upper bound has not closed or the full geometric witness is incomplete' if strict_lower is not None else
                                            'The current lower-bound GDS has not passed the independent exact electrode diameter and center spacing audit'),
                        'lower_bound_evidence': evidence,
                        'actual_output_pads_certified_outside_cuts':bool(evidence),
                        'universal_pad_policy_certificate':universal_pads,
                        'upper_bound_certificate': chosen,
                        'square_envelope_packing_certificate': square_envelope,
                        'enclosing_disk_angular_certificate': angular,
                        'small_cell_conflict_certificate':small_cells,
                        'radial_routing_cut_certificate': cut,
                        'radial_partition_certificate': partition,
                        'multi_circle_dual_certificate': multi_circle,
                        'convex_hull_routing_cut_certificate': hull_cut,
                        'four_grid_shift_upper_bounds': [c['value'] for c in candidates]})
        print(f'{path.name}: {strict_lower if strict_lower is not None else "?"} <= N <= {upper}'
              + (' [strict integer electrode audit]' if strict_lower is not None else
                 ' [strict integer electrode audit unavailable]'),
              flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps({'model_scope': 'four-side complete Pad topology; process-sized global outer-extremity exits and bridges outside the interior core; each Pad outside the original-support enclosing circle and guarded convex hull; electrodes anchored in original GDS support',
                                       'rules': {'electrode_diameter_um': rules.electrode_diameter_um,
                                                 'wire_width_um': rules.wire_width_um,
                                                 'spacing_um': rules.spacing_um,
                                                 'margin_um': rules.margin_um,
                                                 'minimum_center_spacing_um': rules.minimum_center_spacing_um,
                                                 'island_gap_um': settings.pad_gap_um,
                                                 'numeric_guard_um': settings.numeric_guard_um,
                                                 'required_center_separation_um': gap},
                                       'records': records}, ensure_ascii=False,
                                      indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
