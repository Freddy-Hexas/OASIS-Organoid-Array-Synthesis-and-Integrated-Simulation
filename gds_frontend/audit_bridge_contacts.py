"""Replay the spatial contacts made by each added substrate bridge.

Multiple disconnected contact footprints do not by themselves prove that a
bridge touches a different source *component*: the source may already be one
connected polygon.  They do prove that a single-outlet-contact condition has
not been established by the existing GDS audit.  This tool records the facts
without silently changing the permitted layout model.
"""
from __future__ import annotations

from argparse import ArgumentParser
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json

import pyclipper

from exact_gds_audit import (_contours, _count_outer, _paths,
                             _read_layer_paths, _signed_area2, _union)


def _boolean(subject, clip, operation):
    solver = pyclipper.Pyclipper()
    solver.AddPaths(_paths(subject), pyclipper.PT_SUBJECT, True)
    solver.AddPaths(_paths(clip), pyclipper.PT_CLIP, True)
    return solver.Execute2(operation, pyclipper.PFT_NONZERO,
                           pyclipper.PFT_NONZERO)


def _filled_area_twice(tree):
    return sum((-1 if hole else 1) * abs(_signed_area2(path))
               for path, hole in _contours(tree))


def audit_bridge_contacts(report_path):
    path = Path(report_path).resolve(strict=True)
    report = json.loads(path.read_text(encoding='utf-8'))
    source = Path(report['input']).resolve(strict=True)
    output = Path(report['output']).resolve(strict=True)
    source_hash = sha256(source.read_bytes()).hexdigest()
    output_hash = sha256(output.read_bytes()).hexdigest()
    if (source_hash != report['input_sha256'] or
            output_hash != report['output_sha256'] or
            report.get('integer_polygon_audit', {}).get('passed') is not True):
        raise ValueError('Bridge-contact report or GDS hash differs')
    grid_um = report['integer_polygon_audit']['native_grid_um']
    original, read_source_hash, original_grid = _read_layer_paths(
        source, grid_um)
    exported, read_output_hash, output_grid = _read_layer_paths(
        output, grid_um)
    if (read_source_hash != source_hash or read_output_hash != output_hash or
            original_grid != output_grid):
        raise ValueError('Bridge-contact GDS readback differs')
    layers = report['audit']['output_layers']
    support_key = tuple(layers['support_layer'])
    bridge_layer = layers['bridge_marker_layer']
    source_paths = original.get(support_key)
    if not source_paths:
        raise ValueError('Original support layer is empty')
    source_support = _union(source_paths)
    bridge_ids = sorted(datatype for layer, datatype in exported
                        if layer == bridge_layer)
    expected_ids = list(range(1, report['connected_pad_count'] + 1))
    if bridge_ids != expected_ids:
        raise ValueError('Bridge markers do not match expected network ids')
    records = []
    for network in bridge_ids:
        bridge = _union(exported[(bridge_layer, network)])
        contact = _boolean(bridge, source_support,
                           pyclipper.CT_INTERSECTION)
        external = _boolean(bridge, source_support,
                            pyclipper.CT_DIFFERENCE)
        positive_contact_areas = sorted((abs(_signed_area2(path))
            for path, hole in _contours(contact)
            if not hole and _signed_area2(path) != 0), reverse=True)
        count = len(positive_contact_areas)
        if count < 1:
            raise ValueError(f'Network {network} bridge misses original support')
        records.append({
            'network': network,
            'source_contact_components': count,
            'source_contact_component_areas_twice_native_grid_ticks2':
                positive_contact_areas,
            'source_contact_area_twice_native_grid_ticks2':
                _filled_area_twice(contact),
            'outside_source_bridge_components': _count_outer(external),
            'outside_source_bridge_area_twice_native_grid_ticks2':
                _filled_area_twice(external),
        })
    multiple = [item['network'] for item in records
                if item['source_contact_components'] > 1]
    return {
        'created_at': datetime.now(timezone.utc).isoformat(),
        'input_sha256': source_hash,
        'output_sha256': output_hash,
        'native_grid_um': grid_um,
        'network_count': len(records),
        'multiple_disconnected_source_contact_networks': multiple,
        'multiple_disconnected_source_contact_count': len(multiple),
        'all_bridges_have_one_connected_contact_footprint':
            len(multiple) == 0,
        'single_outlet_contact_proven': False,
        'per_network': records,
        'audit_scope': 'integer GDS bridge/source intersection components',
        'interpretation_limit':
            'separate contact footprints need not be distinct connected source branches; one footprint also does not prove that a bridge has no second contact with the same source component',
    }


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--report', type=Path, required=True)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = audit_bridge_contacts(args.report)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, ensure_ascii=False,
                                          indent=2), encoding='utf-8')
    print(json.dumps({key: result[key] for key in
                      ('network_count',
                       'multiple_disconnected_source_contact_count',
                       'single_outlet_contact_proven')},
                     ensure_ascii=False), flush=True)


if __name__ == '__main__':
    main()
