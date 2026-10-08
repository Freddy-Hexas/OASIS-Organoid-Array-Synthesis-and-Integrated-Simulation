"""Read geometry-only JSON/WKB navigation snapshots and recheck their domain.

No pickle, generator metadata, or capacity certificate is loaded. Each reused
edge is checked in the current original-GDS wire-center domain. Cache failure
just invokes the same geometry frontend again.
"""
from dataclasses import asdict
from hashlib import sha256
import gzip
import json
from pathlib import Path
import re
import time

import networkx as nx
import numpy as np
import shapely
from shapely.geometry import LineString,Polygon

from frontend import read_support
from island_router import VectorNavigation,betti


def load_navigation_snapshot(runs,path,rules,support_spec,settings,policy,*,progress=None):
    start=time.perf_counter();input_sha=sha256(Path(path).read_bytes()).hexdigest()
    folders=sorted(Path(runs).glob('*/regions.json.gz'),key=lambda p:p.stat().st_mtime,reverse=True)
    for region_path in folders:
        try:
            status_path=region_path.parent/'status.json'
            with status_path.open(encoding='utf-8') as stream:header=stream.read(4096)
            if not re.search(r'"input_sha256"\s*:\s*"'+input_sha+'"',header):continue
            status=json.loads(status_path.read_text(encoding='utf-8'))
            if (status.get('method')!='four_side_pads' or status['rules']!=asdict(rules) or
                    status['support_layer']!='/'.join(map(str,support_spec))):continue
            with gzip.open(region_path,'rt',encoding='utf-8') as stream:region=json.load(stream)
            cached_settings=region.get('navigation_settings')
            if cached_settings is None:
                cached_settings=status.get('result',{}).get('routing',{}).get('settings')
            if cached_settings!=asdict(settings):continue
            summary=region['rules']
            if (summary.get('outer_exit_policy')!=policy.record() or
                    summary.get('outlet_mode')!='external_boundary' or summary.get('sha256')!=input_sha):continue
            with gzip.open(region_path.parent/'graph.json.gz','rt',encoding='utf-8') as stream:record=json.load(stream)
            if record['input_sha256']!=input_sha:continue
            if progress:progress('重检已有矢量导航快照与原 GDS 的包含关系',.18)
            support,_,meta=read_support(path,*support_spec)
            cached_support=shapely.from_wkb(bytes.fromhex(region['support_wkb_hex']))
            if support.wkb!=cached_support.wkb:continue
            domain=shapely.from_wkb(bytes.fromhex(region['attachment_anchor_region_wkb_hex']))
            center=support.buffer(-(rules.wire_width_um/2+rules.margin_um+settings.numeric_guard_um),quad_segs=32)
            if (not center.covers(domain) or betti(domain)!=summary['navigation_domain_topology'] or
                    summary.get('extractor_uses_generator_source') is not False):continue
            shapely.prepare(center)
            graph=nx.Graph()
            for node in record['nodes']:
                attrs={k:v for k,v in node.items() if k!='id'}
                graph.add_node(node['id'],**attrs)
            lines=[]
            for edge in record['corridors']:
                points=np.asarray(edge['points_um'],dtype=float);line=LineString(points)
                if (not np.isfinite(points).all() or
                    min(np.linalg.norm(points[0]-graph.nodes[edge['start']]['xy_um'])+
                        np.linalg.norm(points[-1]-graph.nodes[edge['end']]['xy_um']),
                        np.linalg.norm(points[-1]-graph.nodes[edge['start']]['xy_um'])+
                        np.linalg.norm(points[0]-graph.nodes[edge['end']]['xy_um']))>meta['gds_native_precision_m']*1e6):
                    raise ValueError('Cache corridor has invalid endpoints')
                graph.add_edge(edge['start'],edge['end'],points_um=points,weight=line.length,
                               corridor_id=edge['corridor_id']);lines.append(line)
            if not bool(np.all(shapely.covers(center,lines))):continue
            components=nx.number_connected_components(graph)
            if {'components':components,'holes':graph.number_of_edges()-graph.number_of_nodes()+components}!=record['topology_certificate']:
                continue
            if any(n not in graph for n in [*record['outlets'],*record['candidate_anchors']]):continue
            if not bool(np.all(shapely.covers(center,
                    shapely.points([graph.nodes[n]['xy_um'] for n in graph])))):continue
            summary={**summary,'navigation_cache':{'source_job_id':region_path.parent.name,
                'seconds':time.perf_counter()-start,'current_original_gds_rechecked':True,
                'all_corridor_centerlines_in_current_physical_domain':True,
                'scope':'geometry proposal reuse; original-GDS metal and independent export audit remain mandatory'}}
            return VectorNavigation(support,center,domain,Polygon(),graph,summary,
                                    record['candidate_anchors'],record['outlets'])
        except (OSError,ValueError,KeyError,TypeError,IndexError):continue
    return None
