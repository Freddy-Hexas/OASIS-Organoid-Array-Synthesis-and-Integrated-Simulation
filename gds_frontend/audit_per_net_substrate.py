"""Integer audit that each network uses only its own added substrate.

The standard exported-GDS audit checks metal against the union of all support.
This optional stronger check uses original GDS support, the common outer shell,
and only the selected network's own island and bridge. It prevents another
network's substrate additions from silently carrying a wire. It is a lower-
bound witness check, not an optimality proof or complete manufacturing DRC.
"""
from __future__ import annotations

from argparse import ArgumentParser
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import math

from exact_gds_audit import (_read_layer_paths, _union, _outside,
                             _segments, _segment_index, _nearby_pairs_clear)


def audit(report_path):
    path = Path(report_path).resolve(strict=True)
    report = json.loads(path.read_text(encoding='utf-8'))
    source = Path(report['input']).resolve(strict=True)
    output = Path(report['output']).resolve(strict=True)
    if (sha256(source.read_bytes()).hexdigest()!=report['input_sha256'] or
            sha256(output.read_bytes()).hexdigest()!=report['output_sha256'] or
            report.get('integer_polygon_audit',{}).get('passed') is not True):
        raise ValueError('Input, output or prior integer audit does not match report')
    layers = report['audit']['output_layers']
    native_grid = report['integer_polygon_audit']['native_grid_um']
    original, original_hash, source_grid = _read_layer_paths(source,native_grid)
    exported, output_hash, output_grid = _read_layer_paths(output,native_grid)
    if (original_hash!=report['input_sha256'] or
            output_hash!=report['output_sha256'] or
            source_grid!=output_grid):
        raise ValueError('Independent GDS readback differs from report')
    support_key = tuple(layers['support_layer'])
    original_paths = original.get(support_key,[])
    shell_paths = exported.get((layers['shell_marker_layer'],0),[])
    if not original_paths or not shell_paths:
        raise ValueError('Missing original support or common shell')
    required = math.ceil(report['rules']['margin_um']/native_grid-1e-9)
    if required<=0:
        raise ValueError('Positive support-edge margin is required')
    checks = 0
    for network in range(1,report['connected_pad_count']+1):
        metal_paths = exported.get((layers['metal_layer'],network),[])
        island_paths = exported.get((layers['island_marker_layer'],network),[])
        bridge_paths = exported.get((layers['bridge_marker_layer'],network),[])
        if not metal_paths or not island_paths or not bridge_paths:
            raise ValueError(f'Network {network} lacks metal, island or bridge')
        metal = _union(metal_paths)
        own_support = _union([*original_paths,*shell_paths,
                              *island_paths,*bridge_paths])
        if _outside(metal,own_support):
            raise ValueError(f'Network {network} uses another net\'s substrate')
        metal_segments = _segments(metal)
        support_segments = _segments(own_support)
        clear,count = _nearby_pairs_clear(
            metal_segments,support_segments,
            _segment_index(support_segments),required)
        checks += count
        if not clear:
            raise ValueError(f'Network {network} violates own-support margin')
    return {
        'created_at':datetime.now(timezone.utc).isoformat(),
        'input_sha256':original_hash,
        'output_sha256':output_hash,
        'native_grid_um':native_grid,
        'network_count':report['connected_pad_count'],
        'minimum_own_support_margin_um':report['rules']['margin_um'],
        'nearby_integer_segment_pairs_checked':checks,
        'passed':True,
        'scope':'each metal polygon clears the original support plus common shell and only its own island and bridge',
        'not_covered':'continuous routing optimality and complete manufacturing DRC',
    }


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--report',type=Path,required=True)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    result = audit(args.report)
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2),
                               encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2),flush=True)


if __name__=='__main__':
    main()
