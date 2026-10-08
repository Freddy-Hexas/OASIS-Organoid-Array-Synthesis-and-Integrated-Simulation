"""GDS geometry workbench with a bounded parallel job queue.

Run: python start_workbench.py from the package root. Inputs default to data/.
No path submitted by the browser is ever opened as an input GDS.
"""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict
from datetime import datetime, timezone
from hashlib import sha256
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from tempfile import TemporaryDirectory
from threading import Lock
from urllib.parse import unquote, urlsplit, parse_qs
from xml.sax.saxutils import quoteattr
import math
import gzip
import json
import os
import re
import sys
import traceback
import uuid

HERE=Path(__file__).resolve().parent
FRONTEND_DIR=HERE.parent
PACKAGE_DIR=FRONTEND_DIR.parent
sys.path[:0]=[str(HERE),str(FRONTEND_DIR)]
from workspace_paths import DATA_DIR, RUNS_DIR
INPUT_DIR=DATA_DIR
RUNS_DIR.mkdir(parents=True,exist_ok=True)

import gdstk
from frontend import extract,read_support,_polygonal
from process_geometry import ProcessRules,characterize
from demo_router import route,write_demo_gds,audit_demo_gds
from island_router import (IslandSettings,build_navigation,solve as solve_islands,
                           write_gds as write_island_gds,audit_gds as audit_island_gds,
                           save_navigation_graph,preview_navigation,betti)
from pad_router import (PAD_LAYOUT_REVISION,PadSettings,solve_with_pads,
                        write_pad_gds,audit_pad_gds)
from exact_gds_audit import audit_exported_gds
from outer_exit_policy import make_outer_exit_policy
from route_incumbent import discover_incumbents
from navigation_cache import load_navigation_snapshot
from worker_identity import is_verification_worker_running
from electrode_region import region_for_navigation, region_for
from distribution_metrics import REVISION as DISTRIBUTION_REVISION, validate_target, evaluate_distribution
from run_frontend import save_graph,preview,routing_preview,_plain
from preview_cache import DISPLAY_CACHE,file_stamp,svg_scene,preview_job,compact_status

LOCK=Lock()
SAVE_LOCK=Lock()
PLOT_LOCK=Lock()
MAX_WORKERS=max(1,min(4,int(os.environ.get('GDS_WORKBENCH_WORKERS','2'))))
MAX_PENDING_JOBS=32
POOL=ThreadPoolExecutor(max_workers=MAX_WORKERS,thread_name_prefix='gds-workbench')
JOBS={}
RECENT_JOBS={}
INPUT_LIMIT=32*1024
ARTIFACT_TYPES={'.png':'image/png','.gds':'application/octet-stream',
                '.gz':'application/gzip','.json':'application/json'}
PAD_INPUT_FIELDS=('square_side_um','pad_width_um','pad_length_um','pad_pitch_um')
PAD_SIZING_FIELDS=('sizing_mode','minimum_pad_width_um','minimum_pad_length_um')


def _pad_settings_from_request(raw,rules):
    """Read reference dimensions; the solver derives slot counts from geometry."""
    if raw is None:
        raw={}
    if not isinstance(raw,dict) or set(raw)-set(PAD_INPUT_FIELDS+PAD_SIZING_FIELDS):
        raise ValueError('Unknown Pad geometry or sizing parameter')
    defaults=PadSettings()
    values={}
    for key in PAD_INPUT_FIELDS:
        value=raw.get(key,getattr(defaults,key))
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value):
            raise ValueError(f'{key} 必须是有限的数值，单位 μm')
        values[key]=float(value)
    minimum_contact=rules.wire_width_um+2*rules.margin_um
    if (values['square_side_um']<=0 or values['pad_width_um']<minimum_contact or
            values['pad_length_um']<minimum_contact):
        raise ValueError(f'Pad 宽度和长度至少需容纳导线及两侧余量（本次 {minimum_contact:g} μm）')
    if values['pad_pitch_um']<values['pad_width_um']+rules.spacing_um:
        raise ValueError('Pad 节距必须不小于 Pad 宽度加网络间距')
    values['sizing_mode']=raw.get('sizing_mode','expand_frame')
    for key in PAD_SIZING_FIELDS[1:]:
        value=raw.get(key)
        if value is not None and (isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value)):
            raise ValueError(f'{key} must be a finite dimension or null')
        values[key]=value
    from pad_sizing import validate_pad_sizing
    validate_pad_sizing(PadSettings(**values),wire_width_um=rules.wire_width_um,margin_um=rules.margin_um)
    # Ancillary support dimensions adapt to the supplied process rules.
    # The historical 40/100 um values are defaults, never process ceilings.
    guard=IslandSettings().numeric_guard_um
    return PadSettings(**values,
        bridge_width_um=max(defaults.bridge_width_um,minimum_contact+2*guard),
        pad_outer_setback_um=max(defaults.pad_outer_setback_um,rules.margin_um+guard))


def utc_now():
    return datetime.now(timezone.utc).isoformat(timespec='milliseconds')


def discover_files():
    records=[]
    for path in INPUT_DIR.iterdir():
        if not path.is_file() or path.suffix.lower()!='.gds':
            continue
        stat=path.stat()
        name=(Path('data')/path.relative_to(INPUT_DIR)).as_posix()
        digest=DISPLAY_CACHE.get(('input_hash',file_stamp(path)),
                                 lambda:sha256(path.read_bytes()).hexdigest())
        records.append({'id':sha256(name.encode('utf-8')).hexdigest()[:16],
                        'name':path.name,'relative_path':name,'absolute_path':str(path),
                        'category':'结构输入','size_bytes':stat.st_size,
                        'sha256':digest,
                        'modified_at':datetime.fromtimestamp(stat.st_mtime,timezone.utc).isoformat(timespec='seconds'),
                        'routable_input':True})
    records.sort(key=lambda r:(r['name'].casefold(),r['relative_path']))
    return records


def latest_tasks():
    """One most recent run per input GDS, including runs before a server restart."""
    tasks={}
    ordering={}
    for path in RUNS_DIR.glob('*/status.json'):
        try:
            cached=DISPLAY_CACHE.get(('status_index',file_stamp(path)),lambda:compact_status(path))
            job=dict(cached)
            file_id=job['file_id']
            if (job.get('status') in ('queued','running') and job['id'] not in JOBS and
                    not is_verification_worker_running(job)):
                job={**job,'status':'interrupted','message':'服务器重启后，该运行没有继续。可以重新启动。'}
            compact={k:job.get(k) for k in ('id','file_id','input_name','status','stage','progress',
                     'message','created_at','finished_at','support_layer','outlet_mode','method',
                     'artifacts','rules','pad_settings','solver_revision','batch_id')}
            compact['lower_bound']=job.get('lower_bound')
            compact['input_sha256']=(job.get('input_sha256') or
                                     job.get('result',{}).get('geometry',{}).get('sha256'))
            created=datetime.fromisoformat(compact['created_at']).astimezone(timezone.utc)
            rank=(created,path.stat().st_mtime_ns)
            if file_id not in ordering or rank>ordering[file_id]:
                tasks[file_id]=compact
                ordering[file_id]=rank
        except (OSError,ValueError,KeyError,TypeError):
            continue
    # Compatible historical layouts stay viewable with their recorded revision.
    # Earlier sparse-Pad runs still require a rerun; input hashes alone cannot
    # distinguish those from layouts produced by the contiguous-bank solver.
    for task in tasks.values():
        if task.get('method')=='four_side_pads' and task.get('solver_revision')!=PAD_LAYOUT_REVISION:
            if task.get('solver_revision')=='center_preferred_pad_routing_v11_directional_envelope_joint_paths':
                task['parameter_note']='历史任务未限制电极区域半径；新任务默认 3 mm'
            elif task.get('solver_revision')=='center_preferred_pad_routing_v12_electrode_region':
                task['parameter_note']='历史布局采用旧位置搜索；新任务使用圆岛附着自适应搜索。原结果保留供查看，尚未按新策略重跑。'
            elif task.get('solver_revision')=='center_preferred_pad_routing_v13_attachment_cells':
                task['parameter_note']='历史布局采用限时搜索；新任务取消默认求解时间限制。原结果保留供查看，尚未按不限时策略重跑。'
            else:
                task['stale_reason']='pad_layout_algorithm_changed'
    return tasks


def file_by_id(file_id):
    for record in discover_files():
        if record['id']==file_id:
            return record
    raise FileNotFoundError('File no longer exists in the fixed input directory')


def _read_gds_details(path:Path):
    raw=path.read_bytes()
    with TemporaryDirectory(prefix='gds_web_read_') as temp:
        ascii_path=Path(temp)/'input.gds'
        ascii_path.write_bytes(raw)
        lib=gdstk.read_gds(str(ascii_path),unit=1e-6)
    layers={}
    for top in lib.top_level():
        for p in top.get_polygons():
            key=f'{p.layer}/{p.datatype}'
            layers[key]=layers.get(key,0)+1
    support_candidates={k:v for k,v in layers.items() if v>0}
    if not support_candidates:
        raise ValueError('GDS contains no polygonal layers')
    if '10/0' in support_candidates:
        recommendation='10/0';reason='优先候选层；未知来源必须人工确认材料语义'
    else:
        recommendation=max(support_candidates,key=support_candidates.get)
        reason='按多边形数量选择的候选层；未知 GDS 请人工确认材料语义'
    layer,datatype=(int(x) for x in recommendation.split('/'))
    support,_,geometry=read_support(path,layer,datatype)
    geometry['electrode_region_reference_um']=region_for(support,{},geometry['gds_native_precision_m']*1e6).record()['reference_um']
    return {'sha256':sha256(raw).hexdigest(),'gds_native_unit_m':lib.unit,
            'gds_native_precision_m':lib.precision,
            'top_cells':[c.name for c in lib.top_level()],
            'layers':dict(sorted(layers.items())),
            'suggested_support_layer':recommendation,
            'layer_suggestion_reason':reason,
            'geometry':geometry}


def _cached_gds_details(path):
    # The endpoint adds selected-layer facts; never mutate the cached object.
    key=('gds_details',file_stamp(path))
    def load():
        packed=DISPLAY_CACHE.disk_bytes(RUNS_DIR,key,
            lambda:json.dumps(_read_gds_details(path),ensure_ascii=False,default=_plain).encode('utf-8'))
        return json.loads(gzip.decompress(packed))
    return dict(DISPLAY_CACHE.get(key,load))


def _gds_svg(path:Path,support_spec,metal_layer=None,marker_layer=None,*,highlight_metal=False,focus_bounds=None):
    key=('gds_svg',file_stamp(path),support_spec,metal_layer,marker_layer,highlight_metal,focus_bounds)
    packed=DISPLAY_CACHE.disk_bytes(RUNS_DIR,key,
        lambda:_build_gds_svg(path,support_spec,metal_layer,marker_layer,
                              highlight_metal=highlight_metal,focus_bounds=focus_bounds))
    return gzip.decompress(packed)


def _build_gds_svg(path:Path, support_spec:tuple[int,int], metal_layer=None, marker_layer=None, *, highlight_metal=False, focus_bounds=None):
    """Render the actual GDS polygon vertices as SVG paths in micrometer coordinates."""
    with TemporaryDirectory(prefix='gds_web_vector_') as temp:
        ascii_path=Path(temp)/'input.gds'
        ascii_path.write_bytes(path.read_bytes())
        lib=gdstk.read_gds(str(ascii_path),unit=1e-6)
    polys=[p for cell in lib.top_level() for p in cell.get_polygons()]
    selected=[p for p in polys if (p.layer,p.datatype)==support_spec or
              (metal_layer is not None and p.layer==metal_layer) or
              (marker_layer is not None and p.layer==marker_layer)]
    if not selected:
        raise ValueError(f'No displayable polygons on support layer {support_spec}')
    import numpy as np
    bounds=np.vstack([p.points for p in selected])
    minx,miny=bounds.min(axis=0);maxx,maxy=bounds.max(axis=0)
    pad=max(float(maxx-minx),float(maxy-miny))*.025+1
    x=float(minx-pad);y=float(miny-pad);w=float(maxx-minx+2*pad);h=float(maxy-miny+2*pad)
    if focus_bounds is not None:
        fx1,fy1,fx2,fy2=focus_bounds
        focus_pad=max(fx2-fx1,fy2-fy1)*.14+35
        x=fx1-focus_pad;y=float(miny+maxy-fy2-focus_pad)
        w=fx2-fx1+2*focus_pad;h=fy2-fy1+2*focus_pad
    lines=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{x:.6f} {y:.6f} {w:.6f} {h:.6f}" '
           f'data-original-viewbox="{x:.6f} {y:.6f} {w:.6f} {h:.6f}" role="img" '
           f'aria-label={quoteattr(path.name + " 的 GDS 矢量多边形")}>',
           f'<rect x="{x:.6f}" y="{y:.6f}" width="{w:.6f}" height="{h:.6f}" fill="#f7f8f5"/>',
           f'<g transform="translate(0 {float(miny+maxy):.6f}) scale(1 -1)">']
    if highlight_metal:
        selected.sort(key=lambda p: ((p.layer,p.datatype)!=support_spec,
                                      p.layer==metal_layer,p.layer==marker_layer))
    for p in selected:
        if (p.layer,p.datatype)==support_spec:
            fill='#3b807d';opacity='.24' if highlight_metal else '.76'
        elif p.layer==metal_layer:
            fill='#da683e';opacity='1'
        else:
            fill='#142e51';opacity='.95'
        points=p.points
        if len(points)<3:
            continue
        d='M'+' '.join(f'{float(a):.6f},{float(b):.6f}' for a,b in points)+'Z'
        outline=(f' stroke="{fill}" stroke-width="55" stroke-linejoin="round"'
                 if highlight_metal and p.layer in (metal_layer,marker_layer) else '')
        lines.append(f'<path d="{d}" fill="{fill}" fill-opacity="{opacity}"'
                     f' fill-rule="evenodd"{outline}/>')
    lines.extend(['</g>','</svg>'])
    return ''.join(lines).encode('utf-8')


def _feasible_svg(job):
    """Show the exact vector center region used by the job, over its input GDS."""
    import shapely.wkb
    record=file_by_id(job['file_id'])
    spec=tuple(int(x) for x in job['support_layer'].split('/'))
    base=_gds_svg(Path(record['absolute_path']),spec).decode('utf-8')
    with gzip.open(RUNS_DIR/job['id']/'regions.json.gz','rt',encoding='utf-8') as file:
        regions=json.load(file)
    key=('attachment_anchor_region_wkb_hex' if job.get('method') in ('attached_islands','four_side_pads')
         else 'electrode_center_region_wkb_hex')
    feasible=shapely.wkb.loads(bytes.fromhex(regions[key]))
    if job.get('rules',{}).get('electrode_region_radius_um') is not None:
        source=shapely.wkb.loads(bytes.fromhex(regions['support_wkb_hex']))
        grid=regions.get('rules',{}).get('gds_native_precision_m',1e-9)*1e6
        feasible=region_for(source,job['rules'],grid).clip(feasible)
    paths=[]
    for polygon in _polygonal(feasible):
        rings=[polygon.exterior,*polygon.interiors]
        d=' '.join('M'+' '.join(f'{float(x):.6f},{float(y):.6f}' for x,y in ring.coords)+'Z'
                   for ring in rings)
        paths.append(f'<path d="{d}" fill="#9b67d0" fill-opacity=".72" fill-rule="evenodd"/>')
    return base.replace('</g></svg>',''.join(paths)+'</g></svg>').encode('utf-8')


def _curve_svg(job, *, focus_first_route=False):
    """Overlay exact quadratic centerline commands on the exported GDS view."""
    routes=job.get('result',{}).get('routing',{}).get('routes',[])
    if not routes or not any(route.get('curve',{}).get('curve_segments') for route in routes):
        raise FileNotFoundError('This job has no analytic curved centerlines')
    layers=job['result']['routing']['gds_roundtrip_audit']['output_layers']
    support_spec=tuple(int(x) for x in job.get('support_layer','10/0').split('/'))
    focus=None
    if focus_first_route:
        points=routes[0]['points_um']
        focus=(min(point[0] for point in points),min(point[1] for point in points),
               max(point[0] for point in points),max(point[1] for point in points))
    base=_gds_svg(RUNS_DIR/job['id']/'routing.gds',support_spec,
                  layers.get('metal_layer'),layers.get('electrode_marker_layer'),
                  focus_bounds=focus).decode('utf-8')
    paths=[]
    def xy(point):
        return f'{float(point[0]):.6f},{float(point[1]):.6f}'
    for index,route in enumerate(routes,1):
        points=route['points_um']
        commands=['M'+xy(points[0])]
        cursor=0
        for segment in route.get('curve',{}).get('curve_segments',[]):
            start=segment['start_um']
            hit=next((i for i in range(cursor,len(points))
                      if math.dist(points[i],start)<=1e-6),None)
            if hit is None:
                raise ValueError(f'Curve start not found on sampled route {index}')
            commands.extend('L'+xy(point) for point in points[cursor+1:hit+1])
            commands.append('Q'+xy(segment['control_um'])+' '+xy(segment['end_um']))
            cursor=hit+segment['samples']
            if cursor>=len(points) or math.dist(points[cursor],segment['end_um'])>1e-6:
                raise ValueError(f'Curve end not found on sampled route {index}')
        commands.extend('L'+xy(point) for point in points[cursor+1:])
        path=' '.join(commands)
        paths.append(f'<path d={quoteattr(path)} fill="none" stroke="#142e51" '
                     f'stroke-width="{10 if focus_first_route else 1.2}" stroke-linecap="round" stroke-linejoin="round" '
                     f'aria-label={quoteattr(f"网络 {index} 的解析二次贝塞尔中心线")}/>')
    return base.replace('</g></svg>',''.join(paths)+'</g></svg>').encode('utf-8')


def _save_job(job_id):
    # Serialize snapshots so a slow earlier write cannot replace newer progress.
    with SAVE_LOCK:
        with LOCK:
            data=dict(JOBS[job_id])
        path=RUNS_DIR/job_id/'status.json'
        temp=path.with_name(f'status.{uuid.uuid4().hex}.tmp')
        try:
            temp.write_text(json.dumps(data,ensure_ascii=False,indent=2,default=_plain),encoding='utf-8')
            os.replace(temp,path)
        finally:
            temp.unlink(missing_ok=True)


def _update_job(job_id,**changes):
    with LOCK:
        JOBS[job_id].update(changes)
    _save_job(job_id)
    if changes.get('status') in ('complete','error','interrupted'):
        # Full results live on disk. Keep only a small latest-per-file summary
        # in RAM so repeated batches do not retain every route geometry.
        with LOCK:
            job=JOBS.pop(job_id,None)
            if job:
                compact=_compact_job(job)
                previous=RECENT_JOBS.get(job['file_id'])
                if previous is None or compact['created_at']>=previous['created_at']:
                    RECENT_JOBS[job['file_id']]=compact


def _job_public(job_id):
    if not re.fullmatch(r'[0-9a-f]{16}',job_id):
        raise FileNotFoundError('Unknown job')
    with LOCK:
        job=JOBS.get(job_id)
        if job:
            return dict(job)
    path=RUNS_DIR/job_id/'status.json'
    if path.exists():
        data=json.loads(path.read_text(encoding='utf-8'))
        if data.get('status') in ('queued','running') and not is_verification_worker_running(data):
            data['status']='interrupted'
            data['message']='服务器重启后，该运行没有继续。可以重新启动。'
        return data
    raise FileNotFoundError('Unknown job')


def _preview_job_bytes(job_id):
    if not re.fullmatch(r'[0-9a-f]{16}',job_id):raise FileNotFoundError('Unknown job')
    with LOCK:
        active=JOBS.get(job_id)
        if active:return gzip.compress(json.dumps(preview_job(dict(active)),ensure_ascii=False,default=_plain).encode('utf-8'),mtime=0)
    path=RUNS_DIR/job_id/'status.json'
    key=('preview_job',file_stamp(path))
    return DISPLAY_CACHE.disk_bytes(RUNS_DIR,key,
        lambda:json.dumps(preview_job(_job_public(job_id)),ensure_ascii=False,default=_plain).encode('utf-8'))


def _display_job(job_id):
    return json.loads(gzip.decompress(_preview_job_bytes(job_id)))


def _distribution_source(job_id):
    if not re.fullmatch(r'[0-9a-f]{16}',job_id):
        raise FileNotFoundError('Unknown job')
    stamp=file_stamp(RUNS_DIR/job_id/'status.json')
    def read():
        job=_job_public(job_id)
        if job.get('status')!='complete':
            raise ValueError('电极分布评价须等待任务完成')
        result=job.get('result',{});routing=result.get('routing',{})
        count=routing.get('retained_routes')
        if not isinstance(count,int) or count<0:
            raise ValueError('该历史报告没有可用的最终电极数量')
        audit=routing.get('gds_roundtrip_audit',{})
        if count and audit.get('passed') is not True:
            raise ValueError('该布局尚未通过 GDS 回读，不能作为已布置电极评价')
        routes=routing.get('routes',[])
        if len(routes)!=count or any('source_um' not in route for route in routes):
            raise ValueError('最终报告缺少完整电极坐标，请查看原报告；不会用候选点代替')
        return json.dumps({'job_id':job_id,'input_name':job.get('input_name'),
                           'input_sha256':result.get('geometry',{}).get('sha256'),
                           'method':job.get('method',result.get('method')),
                           'points_um':[route['source_um'] for route in routes],
                           'coordinate_source':'final report source_um after successful GDS round-trip audit; not sampled candidates',
                           'gds_roundtrip_verified':audit.get('passed') is True,
                           'integer_polygon_verified':routing.get('integer_polygon_audit',{}).get('passed') is True},
                          ensure_ascii=False,separators=(',',':')).encode('utf-8')
    return stamp,json.loads(DISPLAY_CACHE.get(('distribution_source',stamp),read))


def _distribution_bytes(job_id,query):
    def number(key,default):
        items=query.get(key,[str(default)])
        if len(items)!=1:raise ValueError('评价参数不能重复')
        return float(items[0])
    allowed={'radius_um','center_x_um','center_y_um','coverage_distance_um','resolution','download'}
    if set(query)-allowed:raise ValueError('未知电极分布评价参数')
    resolution=number('resolution',384)
    target=validate_target(number('radius_um',3000),
                           [number('center_x_um',0),number('center_y_um',0)],
                           number('coverage_distance_um',100),resolution)
    stamp,source=_distribution_source(job_id)
    key=('distribution',DISTRIBUTION_REVISION,stamp,target['radius_um'],
         *target['center_um'],target['coverage_distance_um'],target['resolution'])
    def calculate():
        report=evaluate_distribution(source['points_um'],radius_um=target['radius_um'],
                                     center_um=target['center_um'],
                                     coverage_distance_um=target['coverage_distance_um'],
                                     resolution=target['resolution'])
        return json.dumps({**source,**report,'points_um':None},ensure_ascii=False,
                          separators=(',',':'),allow_nan=False).encode('utf-8')
    return DISPLAY_CACHE.disk_bytes(RUNS_DIR,key,calculate)


def _vector_display(kind,identity,query,scene=False):
    """Cache unchanged exact vector views, independently of routing state."""
    if kind=='input':
        record=file_by_id(identity);path=Path(record['absolute_path'])
        details=_cached_gds_details(path)
        spec=query.get('spec',[details['suggested_support_layer']])[0]
        if spec not in details['layers']:raise ValueError('Selected layer does not exist in this GDS')
        stamp=(file_stamp(path),spec)
        def build():return _gds_svg(path,tuple(map(int,spec.split('/')))),None
    else:
        job=_display_job(identity)
        required='regions.json.gz' if kind=='feasible' else 'routing.gds'
        if required not in job.get('artifacts',[]):raise FileNotFoundError('This job has no requested vector view')
        path=RUNS_DIR/identity/required
        stamp=(file_stamp(RUNS_DIR/identity/'status.json'),file_stamp(path))
        if kind=='feasible':
            record=file_by_id(job['file_id']);stamp=(*stamp,file_stamp(record['absolute_path']))
        def build():
            focus=job.get('result',{}).get('routing',{}).get('preview_metrics',{}).get('first_route_bounds_um')
            if kind=='curve':
                full=_job_public(identity)
                return _curve_svg(full,focus_first_route=query.get('focus',[''])[0]=='first_route'),focus
            if kind=='feasible':return _feasible_svg(job),None
            if kind!='job':raise FileNotFoundError('Unknown vector view')
            layers=job['result']['routing']['gds_roundtrip_audit']['output_layers']
            highlight=query.get('highlight',[''])[0]=='metal'
            shown_focus=tuple(focus) if highlight and query.get('focus',[''])[0]=='first_route' and focus else None
            return _gds_svg(path,tuple(map(int,job.get('support_layer','10/0').split('/'))),
                layers.get('metal_layer'),layers.get('electrode_marker_layer'),
                highlight_metal=highlight,focus_bounds=shown_focus),focus
    normalized=tuple(sorted((key,tuple(values)) for key,values in query.items() if key in ('spec','highlight','focus')))
    key=('scene' if scene else 'vector',kind,stamp,normalized)
    def create():
        raw,focus=build()
        if scene:return json.dumps(svg_scene(raw,focus),ensure_ascii=False,separators=(',',':')).encode('utf-8')
        return raw
    return DISPLAY_CACHE.disk_bytes(RUNS_DIR,key,create)


def _compact_job(job):
    fields=('id','file_id','input_name','input_sha256','status','stage','progress',
            'message','created_at','started_at','finished_at','support_layer',
            'outlet_mode','method','artifacts','rules','pad_settings','batch_id')
    result=job.get('result') or {}
    return {**{key:job.get(key) for key in fields},
            'lower_bound':(result.get('capacity_interval') or {}).get('lower_bound',
                                                                      job.get('lower_bound'))}


def _queue_snapshot():
    with LOCK:
        jobs=list(JOBS.values())+list(RECENT_JOBS.values())
        latest={}
        for job in jobs:
            current=latest.get(job['file_id'])
            if current is None or job['created_at']>=current['created_at']:
                latest[job['file_id']]=job
        running=sum(job['status']=='running' for job in jobs)
        queued=sum(job['status']=='queued' for job in jobs)
        return {'max_workers':MAX_WORKERS,'max_pending_jobs':MAX_PENDING_JOBS,
                'running':running,'queued':queued,
                'available_slots':max(0,MAX_WORKERS-running),
                'jobs':{file_id:_compact_job(job) for file_id,job in latest.items()},
                'active_jobs':[_compact_job(job) for job in jobs
                               if job['status'] in ('queued','running')]}


def _validated_settings(data):
    outlet_mode=data.get('outlet_mode','auto_geometry')
    if outlet_mode not in ('collector','open_tips','auto_geometry','external_boundary'):
        raise ValueError('Unsupported outlet mode')
    method=data.get('method','attached_islands')
    if method not in ('attached_islands','four_side_pads','legacy_fixed_support'):
        raise ValueError('Unsupported routing method')
    if method=='four_side_pads':
        outlet_mode='external_boundary'
    elif outlet_mode=='external_boundary':
        raise ValueError('外边界出口仅适用于完整 Pad 布线模式')
    default=json.loads((FRONTEND_DIR/'process_rules.json').read_text(encoding='utf-8'))
    supplied=data.get('rules',{})
    if not isinstance(supplied,dict) or set(supplied)-set(default):
        raise ValueError('Unknown process rule')
    default.update(supplied)
    rules=ProcessRules(**default);rules.validate()
    if method=='legacy_fixed_support' and rules.minimum_center_spacing_um>0:
        raise ValueError('最小电极中心距仅支持圆岛联合布线模式')
    if method!='four_side_pads' and data.get('pad_settings') is not None:
        raise ValueError('Pad 参数仅适用于连接四边 Pad 模式')
    pad_settings=(_pad_settings_from_request(data.get('pad_settings'),rules)
                  if method=='four_side_pads' else PadSettings())
    return default,rules,pad_settings,method,outlet_mode


def _validated_file(file_id,support_layer=None):
    record=file_by_id(file_id)
    if not record['routable_input']:
        raise ValueError('请选择原始结构 GDS；已布线演示文件不会作为输入重新布线。')
    details=_read_gds_details(Path(record['absolute_path']))
    spec=details['suggested_support_layer'] if support_layer in (None,'auto') else support_layer
    if not isinstance(spec,str) or spec not in details['layers']:
        raise ValueError(f"{record['name']}: 请选择该 GDS 中存在的支撑层")
    return record,details,spec,tuple(int(x) for x in spec.split('/'))


def _submit_jobs(prepared,settings,batch_id=None):
    default,rules,pad_settings,method,outlet_mode=settings
    with LOCK:
        active=[job for job in JOBS.values() if job['status'] in ('queued','running')]
        active_files={job['file_id'] for job in active}
        conflicts=[record['name'] for record,_,_,_ in prepared if record['id'] in active_files]
        if conflicts:
            raise RuntimeError('这些 GDS 已在运行或排队：'+', '.join(conflicts))
        if len(active)+len(prepared)>MAX_PENDING_JOBS:
            raise OverflowError(f'队列最多容纳 {MAX_PENDING_JOBS} 个未完成任务')
        entries=[]
        for record,details,spec,support_spec in prepared:
            job_id=uuid.uuid4().hex[:16]
            output=RUNS_DIR/job_id;output.mkdir()
            JOBS[job_id]={'id':job_id,'file_id':record['id'],
                          'input_name':record['relative_path'],'input_path':record['absolute_path'],
                          'input_sha256':details['sha256'],
                          'solver_revision':PAD_LAYOUT_REVISION if method=='four_side_pads' else None,
                          'status':'queued','stage':'排队中','progress':0.0,
                          'message':'等待工作线程','created_at':utc_now(),
                          'rules':default,'pad_settings':asdict(pad_settings) if method=='four_side_pads' else None,
                          'outlet_mode':outlet_mode,'method':method,'artifacts':[],
                          'support_layer':spec,'batch_id':batch_id}
            entries.append((job_id,record,support_spec))
    for job_id,_,_ in entries:
        _save_job(job_id)
    for job_id,record,support_spec in entries:
        POOL.submit(_run_job,job_id,Path(record['absolute_path']),rules,outlet_mode,
                    support_spec,method,pad_settings)
    return [{'job_id':job_id,'file_id':record['id'],'status':'queued'}
            for job_id,record,_ in entries]


def _run_job(job_id,path:Path,rules:ProcessRules,outlet_mode:str,support_spec:tuple[int,int],method:str,
             pad_settings:PadSettings):
    output=RUNS_DIR/job_id
    try:
        if method in ('attached_islands','four_side_pads'):
            return _run_island_job(job_id,path,rules,outlet_mode,support_spec,
                                   pad_mode=method=='four_side_pads',pad_settings=pad_settings)
        _update_job(job_id,status='running',progress=.05,stage='读取 GDS 支撑层',message='正在读取并分析原始 GDS',
                    started_at=utc_now())
        result=extract(path,layer=support_spec[0],datatype=support_spec[1],
                       collector_clearance_um=rules.collector_clearance_um,
                       first_pitch=rules.first_pitch_um,minimum_pitch=rules.minimum_pitch_um)
        _update_job(job_id,progress=.30,stage='矢量可放区、走廊与节点窗口')
        process,regions=characterize(result,rules)
        with gzip.open(output/'regions.json.gz','wt',encoding='utf-8') as file:
            json.dump(regions,file,ensure_ascii=False,separators=(',',':'))
        save_graph(output/'graph.json.gz',result)
        with PLOT_LOCK:
            preview(output/'graph.png',result)
        _update_job(job_id,progress=.48,stage='选择电极、路径和候选出口')
        routing=route(result,electrode_diameter_um=rules.electrode_diameter_um,
                      electrode_region_radius_um=rules.electrode_region_radius_um,
                      wire_width_um=rules.wire_width_um,spacing_um=rules.spacing_um,
                      margin_um=rules.margin_um,outlet_mode=outlet_mode)
        if routing['status']=='checked_geometry_lower_bound':
            _update_job(job_id,progress=.82,stage='导出金属 GDS 并回读检查')
            output_layers=write_demo_gds(path,output/'routing.gds',routing)
            routing['gds_roundtrip_audit']=audit_demo_gds(output/'routing.gds',result,routing,output_layers)
            with PLOT_LOCK:
                routing_preview(output/'routing.png',result,routing)
            routing['gds_output']='routing.gds'
        else:
            _update_job(job_id,progress=.82,stage='未构造出合法金属路径')
        lower=(routing.get('retained_routes')
               if routing['status']=='checked_geometry_lower_bound' else None)
        summary={'input_path':str(path),'selected_support_layer':list(support_spec),
                 'geometry':result.summary,'process_geometry':process,
                 'routing':{k:v for k,v in routing.items() if k!='_chosen'},
                 'capacity_interval':{'lower_bound':lower,'upper_bound':None,
                                      'upper_bound_status':'连续几何容量上界未证明',
                                      'outlet_mode':outlet_mode,
                                      'scope':'有限候选的一轨几何下界；外环/开放端点不是已定义焊盘'}}
        (output/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,
                                                     indent=2,default=_plain),encoding='utf-8')
        files=['summary.json','graph.json.gz','regions.json.gz','graph.png']
        if (output/'routing.gds').exists():
            files+=['routing.gds','routing.png']
        incomplete_messages={
            'no_geometry_derived_outlet':'已完成几何分析；所选类型没有几何候选出口',
            'no_sampled_legal_electrode_candidate':'已完成几何分析；采样位置无满足电极直径及边缘余量的中心',
            'no_sampled_route_to_selected_outlet':'已完成几何分析；采样单轨图中无可达出口路线'
        }
        _update_job(job_id,status='complete',progress=1.0,stage='完成',
                    message=('提图与 GDS 试布线完成' if lower is not None
                             else incomplete_messages.get(routing['status'],'已完成几何分析；没有验证通过的布线')),
                    finished_at=utc_now(),result=summary,artifacts=files)
    except Exception as exc:
        (output/'error.log').write_text(traceback.format_exc(),encoding='utf-8')
        _update_job(job_id,status='error',stage='运行失败',finished_at=utc_now(),
                    message=f'{type(exc).__name__}: {exc}')


def _run_island_job(job_id,path,rules,outlet_mode,support_spec,*,pad_mode=False,
                    pad_settings=PadSettings()):
    output=RUNS_DIR/job_id
    def save_navigation_diagnostics(certificate):
        (output/'navigation_diagnostics.json').write_text(
            json.dumps(certificate,ensure_ascii=False,indent=2),encoding='utf-8')
        _update_job(job_id,artifacts=['navigation_diagnostics.json'])
    try:
        _update_job(job_id,status='running',progress=.05,stage='读取矢量支撑与构造拓扑图',
                    message='正在构造并认证完整网络；当前阶段见进度说明',started_at=utc_now())
        def progress(stage,value):
            _update_job(job_id,stage=stage,progress=value)
        # Use the task's requested gap for both metal and substrate islands.
        # Otherwise lowering the visible gap would retain a hidden 4 um rule.
        settings=IslandSettings(pad_gap_um=rules.spacing_um)
        exit_policy=None
        if pad_mode:
            source_support,_,source_meta=read_support(path,*support_spec)
            exit_policy=make_outer_exit_policy(source_support,
                wire_width_um=rules.wire_width_um,margin_um=rules.margin_um,
                numeric_guard_um=settings.numeric_guard_um,
                bridge_width_um=pad_settings.bridge_width_um,
                grid_um=source_meta['gds_native_precision_m']*1e6)
        def primary_progress(stage,value):
            progress(stage,value)
        nav=(load_navigation_snapshot(RUNS_DIR,path,rules,support_spec,settings,exit_policy,progress=primary_progress)
             if pad_mode else None)
        if nav is None:
            nav=build_navigation(path,rules,layer=support_spec[0],datatype=support_spec[1],
                             outlet_mode=outlet_mode,settings=settings,progress=primary_progress,
                             diagnostics_callback=save_navigation_diagnostics,
                             outer_exit_policy=exit_policy)
        else:
            save_navigation_diagnostics(nav.summary['navigation_geometry_certificate'])
        additional_navs=()
        auxiliary_skip_reason=None
        if pad_mode:
            auxiliary_skip_reason='Pad routes use one full-source navigation with mandatory outer-extremity ports'
        region={'rules':nav.summary,'navigation_settings':asdict(settings),
                'source_domain':'original vector support',
                'region_kind':('full_source_anchor_domain' if pad_mode else
                               'attached_island_anchor_domain'),
                'attachment_anchor_region_wkb_hex':nav.navigation_domain.wkb_hex,
                'support_wkb_hex':nav.support.wkb_hex,
                'auxiliary_geometry_summary':additional_navs[0].summary if additional_navs else None,
                'auxiliary_skip_reason':auxiliary_skip_reason}
        with gzip.open(output/'regions.json.gz','wt',encoding='utf-8') as file:
            json.dump(region,file,ensure_ascii=False,separators=(',',':'))
        save_navigation_graph(output/'graph.json.gz',nav)
        with PLOT_LOCK:
            preview_navigation(output/'graph.png',nav)
        routing=(solve_with_pads(nav,rules,settings=settings,pad_settings=pad_settings,
                                 progress=progress,additional_navigations=additional_navs,
                                 incumbent_records=discover_incumbents(RUNS_DIR,nav.summary['sha256'],support_spec,rules))
                 if pad_mode else solve_islands(nav,rules,settings=settings,progress=progress))
        files=['summary.json','navigation_diagnostics.json','graph.json.gz','regions.json.gz','graph.png']
        if routing['retained_routes']:
            if pad_mode:
                # Export checks use the actual uniform Pad geometry, rather
                # than the preferred sizes before automatic compaction.
                pad_settings=PadSettings(**routing['pad_settings'])
            progress('导出完整电极—Pad 网络并回读' if pad_mode else '导出衬底圆岛、电极及导线并回读',.85)
            layers=(write_pad_gds(path,output/'routing.gds',routing,support_spec) if pad_mode else
                    write_island_gds(path,output/'routing.gds',routing,support_spec))
            if pad_mode:
                files.append(Path(layers['wire_width_witness_path']).name)
            routing['gds_roundtrip_audit']=(audit_pad_gds(output/'routing.gds',nav,routing,layers,progress=progress) if pad_mode else
                                            audit_island_gds(output/'routing.gds',nav,routing,layers))
            if pad_mode:
                try:
                    integer_audit=audit_exported_gds(
                        path,output/'routing.gds',layers,
                        wire_spacing_um=rules.spacing_um,
                        metal_support_margin_um=rules.margin_um,
                        expected_nets=routing['retained_routes'],
                        minimum_electrode_diameter_um=rules.electrode_diameter_um,
                        maximum_electrode_center_radius_um=rules.electrode_region_radius_um,
                        minimum_electrode_center_distance_um=max(
                            rules.minimum_center_spacing_um,
                            rules.electrode_diameter_um+rules.spacing_um,
                            rules.electrode_diameter_um+2*rules.margin_um+
                            2*settings.numeric_guard_um+settings.pad_gap_um),
                        minimum_substrate_disk_radius_um=(
                            rules.electrode_diameter_um/2+rules.margin_um+
                            settings.numeric_guard_um),
                        maximum_substrate_disk_radius_um=(
                            rules.electrode_diameter_um/2+rules.margin_um+
                            settings.numeric_guard_um+.01),
                        minimum_island_spacing_um=settings.pad_gap_um,
                        minimum_pad_short_side_um=pad_settings.pad_width_um,
                        minimum_pad_long_side_um=pad_settings.pad_length_um,
                        minimum_wire_width_um=rules.wire_width_um,progress=progress)
                    integer_audit['status']='passed'
                except (ValueError,RuntimeError) as exc:
                    integer_audit={'status':'failed',
                                   'reason':f'{type(exc).__name__}: {exc}',
                                   'scope':'independent integer-grid polygon audit; geometric roundtrip result remains separate'}
                routing['integer_polygon_audit']=integer_audit
                (output/'integer_polygon_audit.json').write_text(
                    json.dumps(integer_audit,ensure_ascii=False,indent=2),
                    encoding='utf-8')
                files.append('integer_polygon_audit.json')
            with PLOT_LOCK:
                preview_navigation(output/'routing.png',nav,routing)
            routing['gds_output']='routing.gds'
            files.extend(['routing.gds','routing.png'])
        routing_public={k:v for k,v in routing.items() if not k.startswith('_')}
        region_summary={'rules':routing_public['rules'],
                        'electrode_region':routing_public['electrode_region'],
                        'attachment_anchor_region_area_um2':region_for_navigation(nav,rules).clip(nav.navigation_domain).area,
                        'wire_center_navigation_area_um2':nav.navigation_domain.area,
                        'electrode_center_region_area_um2':None,
                        'vector_legal_single_track_corridors':None,
                        'ordered_junction_windows':sum(len(d['ports_ccw'])>=3 for _,d in nav.graph.nodes(data=True)),
                        'region_meaning':('electrode anchors in the original-support navigation domain intersected with the declared disk; wires retain the full domain'
                                          if pad_mode else
                                          'wire-center anchor domain in original support, excluding geometry-derived collector; circles are added locally')}
        integer_audit=routing.get('integer_polygon_audit',{})
        strict_pad_lower=bool(
            integer_audit.get('passed') is True and
            integer_audit.get('electrode_disks_and_center_spacing_verified') is True and
            integer_audit.get('electrode_centers_in_original_support_verified') is True and
            integer_audit.get('electrode_region_verified') is True and
            integer_audit.get('maximum_electrode_center_radius_um') == rules.electrode_region_radius_um and
            integer_audit.get('substrate_disks_verified') is True and
            integer_audit.get('island_envelopes_verified') is True and
            integer_audit.get('support_provenance_verified') is True and
            integer_audit.get('source_island_topology_verified') is True and
            integer_audit.get('bridge_attachment_and_hole_exclusion_verified') is True and
            integer_audit.get('outer_exit_policy_verified') is True and
            integer_audit.get('island_island_spacing_verified') is True and
            integer_audit.get('functional_wire_width_verified') is True and
            integer_audit.get('pad_dimensions_verified') is True and
            integer_audit.get('minimum_electrode_diameter_um') == rules.electrode_diameter_um and
            integer_audit.get('minimum_electrode_center_distance_um',0) >= max(
                rules.minimum_center_spacing_um,
                rules.electrode_diameter_um+rules.spacing_um,
                rules.electrode_diameter_um+2*rules.margin_um+
                2*settings.numeric_guard_um+settings.pad_gap_um))
        lower=(routing['retained_routes'] if routing.get('pad_connection_verified') and
               strict_pad_lower else None) if pad_mode else (
            routing['retained_routes'] if routing['retained_routes']>0 else None)
        model_proven=bool(
            pad_mode and strict_pad_lower and lower is not None and
            routing['continuous_upper_bound'] is not None and
            lower==routing['continuous_upper_bound'] and
            routing.get('gds_roundtrip_audit',{}).get('passed') is True)
        summary={'input_path':str(path),'selected_support_layer':list(support_spec),
                 'method':'four_side_pads' if pad_mode else 'attached_islands',
                 'geometry':nav.summary,'process_geometry':region_summary,
                 'routing':routing_public,
                 'capacity_interval':{'lower_bound':lower,'upper_bound':routing['continuous_upper_bound'],
                                      'integer_polygon_lower_verified':strict_pad_lower,
                                      'finite_library_upper_bound':routing['optimization']['finite_upper_bound'],
                                      'initial_finite_library_upper_bound':routing.get('initial_finite_library_optimization',routing['optimization']).get('finite_upper_bound'),
                                      'initial_finite_library_bound_scope':'only the initial pruned library before dynamic route generation; not an upper bound on the final layout or continuous capacity',
                                      'bound_gap_closed':bool(lower is not None and
                                                               routing['continuous_upper_bound'] is not None and
                                                               lower==routing['continuous_upper_bound']),
                                      'declared_geometric_model_optimality_proven':model_proven,
                                      'global_optimality_proven':False,
                                      'global_optimality_note':('当前 GDS 与规则的连续容量上界和导出下界相等，已证明声明的几何模型内电极数量最多；完整制造 DRC 与模型外的桥、外框限制仍需单独认证' if model_proven else
                                                                '整数 GDS 审计核对原结构锚点、电极尺寸与中心距、衬底圆盘及来源、岛间距、Pad 矩形、功能性全宽走线、网络连通和间距；桥及外框的额外形状规则、未用金属毛刺和完整制造 DRC 尚未独立精确认证，且上下界未闭合'),
                                      'upper_bound_status':(('原始 GDS 的连续几何上界已与审计构造下界闭合，证明当前声明模型内最多' if model_proven else
                                                              '原始 GDS 的整数网格装填界、精确小网格冲突界与外侧 Pad 割界，独立于候选库；上下界尚未相等，不能证明最多') if pad_mode else
                                                            '原始 GDS 支撑多边形的整数网格覆盖与圆岛面积上界；缺少紧布线割上界，尚不能证明最多'),
                                      'electrode_domain':'original GDS support intersected with the declared electrode-center disk',
                                      'electrode_region':routing_public['electrode_region'],
                                      'minimum_center_spacing_um':rules.minimum_center_spacing_um,
                                      'outlet_mode':outlet_mode,
                                       'scope':('导出 GDS 回读验证的电极—专属实体 Pad 完整网络数；候选选择及连续几何增广的构造下界，不代表最优数量'
                                                if pad_mode else
                                                '有限圆岛、轨道和路由列到几何候选终端的构造性几何下界；不代表到实体 Pad 的完整连接')}}
        (output/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2,default=_plain),encoding='utf-8')
        _update_job(job_id,status='complete',progress=1.0,stage='完成',finished_at=utc_now(),
                     message=(('已回读验证 %d 个电极分别连到实体 Pad' % lower) if pad_mode and lower else
                              '本次有限候选未构造出电极到 Pad 的完整网络' if pad_mode else
                              '电极到几何候选终端的金属 GDS 回读通过；尚未接到实体 Pad' if lower is not None else
                              '已完成矢量拓扑与候选搜索；未构造出完整合法网络'),
                    result=summary,artifacts=files)
    except Exception as exc:
        (output/'error.log').write_text(traceback.format_exc(),encoding='utf-8')
        diagnostics=getattr(exc,'diagnostics',None)
        if diagnostics is not None:
            save_navigation_diagnostics(diagnostics)
        _update_job(job_id,status='error',stage='运行失败',finished_at=utc_now(),
                    message=f'{type(exc).__name__}: {exc}')


class Handler(BaseHTTPRequestHandler):
    server_version='GDSWorkbench/0.1'

    def log_message(self,format,*args):
        sys.stdout.write('%s - %s\n'%(self.address_string(),format%args))

    def _send_json(self,data,status=200):
        raw=json.dumps(data,ensure_ascii=False,default=_plain).encode('utf-8')
        compressed='gzip' in self.headers.get('Accept-Encoding','') and len(raw)>2048
        if compressed:raw=gzip.compress(raw,compresslevel=3,mtime=0)
        self.send_response(status)
        self.send_header('Content-Type','application/json; charset=utf-8')
        self.send_header('Content-Length',str(len(raw)))
        self.send_header('Cache-Control','no-store')
        self.send_header('Vary','Accept-Encoding')
        if compressed:self.send_header('Content-Encoding','gzip')
        self.end_headers();self.wfile.write(raw)

    def _send_file(self,path:Path,mime=None,download=False):
        raw=path.read_bytes()
        self.send_response(200)
        self.send_header('Content-Type',mime or ARTIFACT_TYPES.get(path.suffix,'application/octet-stream'))
        self.send_header('Content-Length',str(len(raw)))
        self.send_header('Cache-Control','no-store')
        self.send_header('X-Content-Type-Options','nosniff')
        if download:
            self.send_header('Content-Disposition',f'attachment; filename="{path.name}"')
        self.end_headers();self.wfile.write(raw)

    def _send_svg(self,raw):
        self.send_response(200)
        self.send_header('Content-Type','image/svg+xml; charset=utf-8')
        self.send_header('Content-Length',str(len(raw)))
        self.send_header('Cache-Control','no-store')
        self.send_header('X-Content-Type-Options','nosniff')
        self.send_header('Content-Security-Policy',"default-src 'none'; style-src 'unsafe-inline'; sandbox")
        self.end_headers();self.wfile.write(raw)

    def _send_cached(self,packed,mime,svg=False):
        etag='W/"'+sha256(packed).hexdigest()+'"'
        matches=[value.strip().removeprefix('W/') for value in self.headers.get('If-None-Match','').split(',')]
        if '*' in matches or etag.removeprefix('W/') in matches:
            self.send_response(304);self.send_header('ETag',etag)
            self.send_header('Cache-Control','no-cache');self.send_header('Vary','Accept-Encoding')
            self.end_headers();return
        compressed='gzip' in self.headers.get('Accept-Encoding','')
        raw=packed if compressed else gzip.decompress(packed)
        self.send_response(200)
        self.send_header('Content-Type',mime)
        self.send_header('Content-Length',str(len(raw)))
        self.send_header('Cache-Control','no-cache')
        self.send_header('ETag',etag)
        self.send_header('Vary','Accept-Encoding')
        self.send_header('X-Content-Type-Options','nosniff')
        if compressed:self.send_header('Content-Encoding','gzip')
        if svg:self.send_header('Content-Security-Policy',"default-src 'none'; style-src 'unsafe-inline'; sandbox")
        self.end_headers();self.wfile.write(raw)

    def _error(self,status,message):
        self._send_json({'error':message},status)

    def do_GET(self):
        parts=[unquote(p) for p in urlsplit(self.path).path.split('/') if p]
        query=parse_qs(urlsplit(self.path).query)
        try:
            if len(parts)==3 and parts[:2]==['api','distribution']:
                packed=_distribution_bytes(parts[2],query)
                return self._send_cached(packed,'application/json; charset=utf-8')
            if parts in (['distribution.js'],['distribution.css']):
                return self._send_file(HERE/parts[0],
                    'text/javascript; charset=utf-8' if parts[0].endswith('.js') else 'text/css; charset=utf-8')
            if len(parts)>=3 and parts[:2]==['vendor','katex']:
                root=(HERE/'vendor'/'katex').resolve();path=root.joinpath(*parts[2:]).resolve()
                if not path.is_relative_to(root) or not path.is_file() or path.suffix not in ('.js','.css','.woff2','.woff','.ttf'):
                    raise FileNotFoundError('Unknown formula resource')
                mime={'.js':'text/javascript; charset=utf-8','.css':'text/css; charset=utf-8',
                      '.woff2':'font/woff2','.woff':'font/woff','.ttf':'font/ttf'}[path.suffix]
                return self._send_file(path,mime)
            if parts in (['vector-viewer.js'],['vector-worker.js']):
                return self._send_file(HERE/parts[0],'text/javascript; charset=utf-8')
            if len(parts)==4 and parts[:2] in (['api','scene'],['api','vector']):
                if parts[2] not in ('input','job','curve','feasible'):raise FileNotFoundError('Unknown vector view')
                scene=parts[1]=='scene'
                return self._send_cached(_vector_display(parts[2],parts[3],query,scene),
                    'application/json; charset=utf-8' if scene else 'image/svg+xml; charset=utf-8',svg=not scene)
            if not parts:
                return self._send_file(HERE/'index.html','text/html; charset=utf-8')
            if parts==['admin']:
                return self._send_file(HERE/'admin.html','text/html; charset=utf-8')
            if parts==['admin.css']:
                return self._send_file(HERE/'admin.css','text/css; charset=utf-8')
            if parts==['admin.js']:
                return self._send_file(HERE/'admin.js','text/javascript; charset=utf-8')
            if parts==['app.css']:
                return self._send_file(HERE/'app.css','text/css; charset=utf-8')
            if parts==['app.js']:
                return self._send_file(HERE/'app.js','text/javascript; charset=utf-8')
            if parts==['api','health']:
                return self._send_json({'ok':True,'input_root':str(INPUT_DIR),
                                        'time':utc_now()})
            if parts==['api','queue']:
                return self._send_json(_queue_snapshot())
            if parts==['api','files']:
                return self._send_json({'root':str(INPUT_DIR),'files':discover_files(),
                                        'tasks':latest_tasks()})
            if len(parts)==3 and parts[:2]==['api','file']:
                record=file_by_id(parts[2])
                details=_cached_gds_details(Path(record['absolute_path']))
                spec=parse_qs(urlsplit(self.path).query).get('spec',[details['suggested_support_layer']])[0]
                if spec not in details['layers']:
                    return self._error(422,'Selected layer does not exist in this GDS')
                if spec!=details['suggested_support_layer']:
                    layer,datatype=(int(x) for x in spec.split('/'))
                    selected_support,_,details['geometry']=read_support(Path(record['absolute_path']),layer,datatype)
                    details['geometry']['electrode_region_reference_um']=region_for(selected_support,{},
                        details['geometry']['gds_native_precision_m']*1e6).record()['reference_um']
                return self._send_json({**record,**details,
                                        'selected_support_layer':spec,
                                        'reference':None})
            if len(parts)==3 and parts[:2]==['api','preview']:
                record=file_by_id(parts[2]);path=Path(record['absolute_path'])
                digest=sha256(path.read_bytes()).hexdigest()
                query=parse_qs(urlsplit(self.path).query)
                spec=query.get('spec',[None])[0]
                if spec is None:
                    spec=_read_gds_details(path)['suggested_support_layer']
                details=_read_gds_details(path)
                if spec not in details['layers']:
                    return self._error(422,'Selected layer does not exist in this GDS')
                layer,datatype=(int(x) for x in spec.split('/'))
                image=RUNS_DIR/'previews'/f'{digest}_{layer}_{datatype}.png'
                if not image.exists():
                    image.parent.mkdir(exist_ok=True)
                    with PLOT_LOCK:
                        if not image.exists():
                            import matplotlib
                            matplotlib.use('Agg')
                            import matplotlib.pyplot as plt
                            support,_,_=read_support(path,layer,datatype)
                            fig,ax=plt.subplots(figsize=(8,8))
                            try:
                                fig.patch.set_facecolor('#f7f8f5');ax.set_facecolor('#f7f8f5')
                                for polygon in _polygonal(support):
                                    ax.fill(*polygon.exterior.xy,color='#2d6f71',lw=0)
                                    for hole in polygon.interiors:
                                        ax.fill(*hole.xy,color='#f7f8f5',lw=0)
                                ax.set_aspect('equal');ax.axis('off')
                                fig.tight_layout(pad=.2);fig.savefig(image,dpi=160,facecolor='#f7f8f5')
                            finally:
                                plt.close(fig)
                return self._send_file(image,'image/png')
            if len(parts)==3 and parts[:2]==['api','job']:
                if query.get('detail',[''])[0]=='preview':
                    return self._send_cached(_preview_job_bytes(parts[2]),'application/json; charset=utf-8')
                return self._send_json(_job_public(parts[2]))
            if len(parts)==4 and parts[:2]==['api','artifact']:
                job=_display_job(parts[2])
                name=parts[3]
                if name not in job.get('artifacts',[]):
                    raise FileNotFoundError('Artifact not listed for this job')
                path=RUNS_DIR/parts[2]/name
                return self._send_file(path,download=path.suffix in ('.gds','.gz','.json'))
            self._error(404,'Unknown endpoint')
        except (BrokenPipeError,ConnectionResetError,ConnectionAbortedError):
            # A view/case switch deliberately cancels obsolete preview reads.
            return
        except ValueError as exc:
            self._error(422,str(exc))
        except FileNotFoundError as exc:
            self._error(404,str(exc))
        except Exception as exc:
            self._error(500,f'{type(exc).__name__}: {exc}')

    def do_POST(self):
        parts=[unquote(p) for p in urlsplit(self.path).path.split('/') if p]
        if parts not in (['api','jobs'],['api','jobs','batch']):
            return self._error(404,'Unknown endpoint')
        try:
            length=int(self.headers.get('Content-Length','0'))
            if length<=0 or length>INPUT_LIMIT:
                return self._error(413,'Invalid request size')
            data=json.loads(self.rfile.read(length))
            if not isinstance(data,dict):
                raise ValueError('Request body must be a JSON object')
            settings=_validated_settings(data)
            if parts==['api','jobs','batch']:
                file_ids=data.get('file_ids')
                if (not isinstance(file_ids,list) or not file_ids or
                        len(file_ids)>MAX_PENDING_JOBS or
                        any(not isinstance(item,str) for item in file_ids) or
                        len(set(file_ids))!=len(file_ids)):
                    raise ValueError('file_ids 必须是非空、无重复的 GDS ID 列表')
                layers=data.get('support_layers',{})
                if not isinstance(layers,dict) or set(layers)-set(file_ids):
                    raise ValueError('support_layers 只能包含本批次文件的图层选择')
                prepared=[_validated_file(file_id,layers.get(file_id,data.get('support_layer')))
                          for file_id in file_ids]
                batch_id=uuid.uuid4().hex[:16]
                jobs=_submit_jobs(prepared,settings,batch_id)
                return self._send_json({'batch_id':batch_id,'jobs':jobs,
                                        'accepted':len(jobs)},202)
            prepared=[_validated_file(str(data.get('file_id','')),data.get('support_layer'))]
            job=_submit_jobs(prepared,settings)[0]
            self._send_json(job,202)
        except RuntimeError as exc:
            self._error(409,str(exc))
        except OverflowError as exc:
            self._error(429,str(exc))
        except FileNotFoundError as exc:
            self._error(404,str(exc))
        except (ValueError,TypeError,KeyError) as exc:
            self._error(422,str(exc))
        except Exception as exc:
            self._error(500,f'{type(exc).__name__}: {exc}')


def main():
    host=os.environ.get('GDS_WORKBENCH_HOST','127.0.0.1').strip() or '127.0.0.1'
    port=int(os.environ.get('GDS_WORKBENCH_PORT','8769'))
    http=ThreadingHTTPServer((host,port),Handler)
    print(f'GDS Workbench listening on {host}:{port}',flush=True)
    print(f'Input root: {INPUT_DIR}',flush=True)
    print(f'Output root: {RUNS_DIR}',flush=True)
    print(f'Parallel workers: {MAX_WORKERS}',flush=True)
    try:
        http.serve_forever(poll_interval=.25)
    except KeyboardInterrupt:
        print('Stopping HTTP server; waiting for submitted jobs to finish.',flush=True)
    finally:
        http.server_close()
        POOL.shutdown(wait=True)


if __name__=='__main__':
    main()
