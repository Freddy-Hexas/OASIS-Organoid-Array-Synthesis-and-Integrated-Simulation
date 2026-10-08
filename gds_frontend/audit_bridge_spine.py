"""Independent integer witness for a wide source-to-shell bridge spine.

The constructor proposes an inner polygon with Clipper's negative offset.
The verifier does *not* trust that offset: it checks the proposed polygon is
inside the exported bridge, is connected, touches original source and outer
shell, and that its entire boundary clears the true bridge boundary by the
declared radius.  The last condition extends from the candidate boundary to
all its interior points: a segment from an interior point to a closer bridge
boundary point would first exit the candidate at an even closer point.

Thus a disk of the certified radius can move continuously from the original
support into the Pad shell while staying in that network's bridge polygon.
This proves an existence width, not exact equality with the nominal 40 um
bridge target and not full substrate manufacturing DRC.
"""
from __future__ import annotations

from argparse import ArgumentParser
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import math

import pyclipper

from capacity_bounds import _simple_integer_polygon
from exact_gds_audit import (_count_outer, _filled_overlap, _nearby_pairs_clear,
                             _outside, _paths, _read_layer_paths, _segment_index,
                             _segments, _union)


def _read_report(report_path):
    path = Path(report_path).resolve(strict=True)
    report = json.loads(path.read_text(encoding='utf-8'))
    source = Path(report['input']).resolve(strict=True)
    output = Path(report['output']).resolve(strict=True)
    source_hash = sha256(source.read_bytes()).hexdigest()
    output_hash = sha256(output.read_bytes()).hexdigest()
    if (source_hash != report['input_sha256'] or
            output_hash != report['output_sha256'] or
            report.get('integer_polygon_audit', {}).get('passed') is not True):
        raise ValueError('Bridge spine source, output or prior audit differs')
    grid_um = report['integer_polygon_audit']['native_grid_um']
    original, read_source_hash, source_grid = _read_layer_paths(source, grid_um)
    exported, read_output_hash, output_grid = _read_layer_paths(output, grid_um)
    if (source_hash != read_source_hash or output_hash != read_output_hash or
            source_grid != output_grid):
        raise ValueError('Bridge spine GDS readback differs')
    layers = report['audit']['output_layers']
    source_shape = _union(original[tuple(layers['support_layer'])])
    shell_shape = _union(exported[(layers['shell_marker_layer'], 0)])
    bridge_layer = layers['bridge_marker_layer']
    network_count = report['connected_pad_count']
    bridge_ids = sorted(datatype for layer, datatype in exported
                        if layer == bridge_layer)
    if (not source_shape or not shell_shape or
            bridge_ids != list(range(1, network_count + 1))):
        raise ValueError('Bridge spine layer markers are incomplete')
    bridges = {network: _union(exported[(bridge_layer, network)])
               for network in bridge_ids}
    return report, source_shape, shell_shape, bridges, grid_um


def _core_passes(raw_paths, bridge, source, shell, radius_ticks):
    if (not isinstance(raw_paths, list) or not raw_paths or
            not all(isinstance(path, list) and
                    all(isinstance(point, list) and len(point)==2 and
                        all(type(value) is int for value in point)
                        for point in path) and
                    _simple_integer_polygon([tuple(point) for point in path])
                    for path in raw_paths)):
        return False, 0
    core = _union(raw_paths)
    if (_count_outer(core) != 1 or _outside(core, bridge) or
            not _filled_overlap(core, source) or
            not _filled_overlap(core, shell)):
        return False, 0
    bridge_segments = _segments(bridge)
    if not bridge_segments:
        return False, 0
    clear, checked = _nearby_pairs_clear(
        _segments(core), bridge_segments,
        _segment_index(bridge_segments), radius_ticks)
    return clear, checked


def generate_bridge_spine_certificate(report_path, *,
                                      nominal_width_um=40.0,
                                      certified_width_um=39.78,
                                      proposed_inset_um=19.9):
    report, source, shell, bridges, grid_um = _read_report(report_path)
    if (not all(math.isfinite(v) and v > 0 for v in
                (nominal_width_um, certified_width_um,
                 proposed_inset_um)) or
            certified_width_um > nominal_width_um or
            proposed_inset_um <= certified_width_um / 2):
        raise ValueError('Invalid nominal or certified bridge width')
    if report['pad_settings']['bridge_width_um'] != nominal_width_um:
        raise ValueError('Nominal bridge width differs from report')
    radius_ticks = math.ceil(certified_width_um / (2 * grid_um) - 1e-9)
    proposal_ticks = math.ceil(proposed_inset_um / grid_um - 1e-9)
    entries = []
    for network, bridge in bridges.items():
        offset = pyclipper.PyclipperOffset()
        offset.AddPaths(_paths(bridge), pyclipper.JT_ROUND,
                        pyclipper.ET_CLOSEDPOLYGON)
        proposal = _union(offset.Execute(-proposal_ticks))
        paths = [[list(point) for point in polygon]
                 for polygon in _paths(proposal)]
        passed, checked = _core_passes(paths, bridge, source, shell,
                                      radius_ticks)
        if not passed:
            raise ValueError(f'Network {network} has no certified wide spine')
        entries.append({'network': network,
                        'core_polygons_grid_ticks': paths,
                        'nearby_boundary_segment_pairs_checked': checked})
    return {'created_at': datetime.now(timezone.utc).isoformat(),
            'input_sha256': report['input_sha256'],
            'output_sha256': report['output_sha256'],
            'native_grid_um': grid_um,
            'nominal_bridge_width_um': nominal_width_um,
            'certified_continuous_bridge_width_um':
                2 * radius_ticks * grid_um,
            'certified_radius_grid_ticks': radius_ticks,
            'candidate_offset_grid_ticks': proposal_ticks,
            'network_count': len(entries),
            'per_network': entries,
            'passed': True,
            'scope': 'existence of a disk-clear connected path inside each exported bridge from original source to outer Pad shell',
            'not_covered': 'exact nominal 40 um width, bridge shape tolerance, detailed substrate manufacturing DRC, continuous global optimality'}


def verify_bridge_spine_certificate(report_path, certificate):
    try:
        report, source, shell, bridges, grid_um = _read_report(report_path)
        radius = certificate['certified_radius_grid_ticks']
        if (type(radius) is not int or radius < 1 or
                certificate['input_sha256'] != report['input_sha256'] or
                certificate['output_sha256'] != report['output_sha256'] or
                certificate['native_grid_um'] != grid_um or
                certificate['nominal_bridge_width_um'] !=
                    report['pad_settings']['bridge_width_um'] or
                certificate['certified_continuous_bridge_width_um'] !=
                    2 * radius * grid_um or
                certificate['certified_continuous_bridge_width_um'] >
                    certificate['nominal_bridge_width_um'] or
                certificate['network_count'] != len(bridges) or
                certificate.get('passed') is not True):
            return False
        entries = certificate['per_network']
        if [entry['network'] for entry in entries] != sorted(bridges):
            return False
        for entry in entries:
            passed, checked = _core_passes(
                entry['core_polygons_grid_ticks'],
                bridges[entry['network']], source, shell, radius)
            if (not passed or checked !=
                    entry['nearby_boundary_segment_pairs_checked']):
                return False
        return True
    except (KeyError, TypeError, ValueError, IndexError,
            pyclipper.ClipperException):
        return False


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--report', type=Path, required=True)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--verify-certificate', type=Path)
    args = parser.parse_args()
    if args.verify_certificate:
        certificate = json.loads(
            args.verify_certificate.read_text(encoding='utf-8'))
        if not verify_bridge_spine_certificate(args.report, certificate):
            raise RuntimeError('Saved bridge spine certificate failed replay')
    else:
        if not args.output:
            parser.error('--output is required when generating')
        certificate = generate_bridge_spine_certificate(args.report)
        if not verify_bridge_spine_certificate(args.report, certificate):
            raise RuntimeError('Generated bridge spine certificate failed replay')
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(certificate, ensure_ascii=False,
                                          indent=2), encoding='utf-8')
    print(json.dumps({key: certificate[key] for key in
                      ('network_count',
                       'certified_continuous_bridge_width_um',
                       'passed')}), flush=True)


if __name__ == '__main__':
    main()
