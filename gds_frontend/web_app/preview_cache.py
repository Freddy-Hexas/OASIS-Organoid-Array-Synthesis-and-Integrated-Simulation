"""Bounded, single-flight caches for presentation only; no solver state."""
from collections import OrderedDict
from hashlib import sha256
from pathlib import Path
from threading import Event, RLock
import gzip
import json
import os
import uuid
import xml.etree.ElementTree as ET
import re

REVISION='exact-vector-display-v3-electrode-region'


def file_stamp(path):
    path=Path(path);stat=path.stat()
    return str(path.resolve()),stat.st_size,stat.st_mtime_ns,stat.st_ctime_ns


class PresentationCache:
    def __init__(self,max_bytes=96*1024*1024,max_entries=128,max_disk_bytes=512*1024*1024):
        self.max_bytes=max_bytes;self.max_entries=max_entries
        self.max_disk_bytes=max_disk_bytes
        self.values=OrderedDict();self.pending={};self.bytes=0;self.lock=RLock()

    def get(self,key,factory,cost=None):
        while True:
            with self.lock:
                if key in self.values:
                    value,size=self.values.pop(key);self.values[key]=(value,size)
                    return value
                event=self.pending.get(key)
                if event is None:
                    event=Event();self.pending[key]=event;break
            event.wait()
        try:
            value=factory();size=cost(value) if cost else (len(value) if isinstance(value,bytes) else 2048)
            with self.lock:
                if size<=self.max_bytes:
                    self.values[key]=(value,size);self.bytes+=size
                    while self.bytes>self.max_bytes or len(self.values)>self.max_entries:
                        _,(_,old_size)=self.values.popitem(last=False);self.bytes-=old_size
            return value
        finally:
            with self.lock:
                self.pending.pop(key,None);event.set()

    def disk_bytes(self,root,key,factory):
        digest=sha256(repr((REVISION,key)).encode('utf-8')).hexdigest()
        path=Path(root)/'display_cache'/f'{digest}.gz'
        def create():
            if path.exists():
                try:
                    raw=path.read_bytes();os.utime(path,None);return raw
                except OSError:pass
            raw=gzip.compress(factory(),compresslevel=3,mtime=0)
            path.parent.mkdir(exist_ok=True)
            temp=path.with_name(f'{digest}.{uuid.uuid4().hex}.tmp')
            try:
                temp.write_bytes(raw);os.replace(temp,path)
                self._trim_disk(path.parent)
            finally:temp.unlink(missing_ok=True)
            return raw
        return self.get(('gzip',str(root),digest),create)

    def _trim_disk(self,directory):
        """Evict generated previews only; never touch GDS or job artifacts."""
        directory=directory.resolve()
        with self.lock:
            entries=[]
            for path in directory.glob('*.gz'):
                if not re.fullmatch(r'[0-9a-f]{64}\.gz',path.name):continue
                try:
                    if path.resolve().parent!=directory:continue
                    stat=path.stat();entries.append((stat.st_mtime_ns,stat.st_size,path))
                except OSError:continue
            total=sum(size for _,size,_ in entries)
            for _,size,path in sorted(entries):
                if total<=self.max_disk_bytes:break
                try:path.unlink();total-=size
                except OSError:pass


DISPLAY_CACHE=PresentationCache()


def svg_scene(raw,focus_bounds=None):
    """Reuse every exact SVG command and style; bound curves by control hulls."""
    root=ET.fromstring(raw);ns={'s':'http://www.w3.org/2000/svg'}
    group=root.find('s:g',ns)
    flip=float(re.search(r'translate\(0 ([^\)]+)\)',group.get('transform'))[1])
    number=re.compile(r'[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?')
    paths=[]
    for path in group.findall('s:path',ns):
        d=path.get('d');coords=list(map(float,number.findall(d)))
        xs=coords[::2];ys=coords[1::2];stroke=float(path.get('stroke-width',0))
        pad=stroke/2
        paths.append({'d':d,'bounds':[min(xs)-pad,min(ys)-pad,max(xs)+pad,max(ys)+pad],
            'fill':path.get('fill','none'),'opacity':float(path.get('fill-opacity',1)),
            'stroke':path.get('stroke'),'width':stroke})
    focus=None
    if focus_bounds:
        x1,y1,x2,y2=focus_bounds;pad=max(x2-x1,y2-y1)*.14+35
        focus=[x1-pad,flip-y2-pad,x2-x1+2*pad,y2-y1+2*pad]
    return {'viewBox':list(map(float,root.get('viewBox').split())),
            'flipY':flip,'paths':paths,'focusBox':focus,
            'pathCount':len(paths),'geometry':'original SVG commands; no vertex simplification'}


def preview_job(job):
    """Keep dashboard facts, remove geometry witnesses served separately."""
    result=job.get('result')
    if not result:return dict(job)
    routing=result.get('routing',{})
    routes=routing.get('routes',[])
    curves=[route.get('curve',{}) for route in routes if route.get('curve')]
    radii=[curve['minimum_bend_radius_um'] for curve in curves if curve.get('minimum_bend_radius_um') is not None]
    compact_routes=[{key:route[key] for key in ('terminal_kind','outlet_node') if key in route}
                    for route in routes]
    compact={key:value for key,value in routing.items() if key not in ('routes',)}
    # These detailed search trajectories contain millions of coordinates.
    # The UI uses their aggregate measures, never their individual witnesses.
    for key in ('center_compaction','initial_center_compaction','port_augmentation'):
        if isinstance(compact.get(key),dict):
            compact[key]={name:value for name,value in compact[key].items()
                          if not isinstance(value,(list,tuple))}
            if 'accepted_updates' in routing[key]:
                compact[key]['accepted_update_count']=len(routing[key]['accepted_updates'])
    compact['routes']=compact_routes
    first_bounds=None
    if routes and routes[0].get('points_um'):
        points=routes[0]['points_um']
        if routes[0].get('pad_target_um'):points=[*points,routes[0]['pad_target_um']]
        first_bounds=[min(p[0] for p in points),min(p[1] for p in points),
                      max(p[0] for p in points),max(p[1] for p in points)]
    compact['preview_metrics']={'curve_count':len(curves),
        'curved_corner_count':sum(curve.get('curved_corner_count',0) for curve in curves),
        'minimum_bend_radius_um':min(radii,default=None),
        'has_curves':any(curve.get('curve_segments') for curve in curves),
        'first_route_bounds_um':first_bounds}
    return {**job,'result':{**result,'routing':compact},
            'detail_level':'preview; full geometry and certificates in summary.json'}


def compact_status(path):
    """Read the small top-level header and last capacity facts of old jobs.

    The standard status writer puts result last. For other layouts, use the
    full JSON instead of guessing. This avoids parsing a 100 MB history file
    just to list its title, state and electrode count.
    """
    path=Path(path)
    with path.open(encoding='utf-8') as stream:
        head=stream.read(32768)
    match=re.search(r'\n  "result":',head)
    if match:
        try:
            job=json.loads(head[:match.start()].rstrip().rstrip(',')+'\n}')
            with path.open('rb') as stream:
                stream.seek(max(0,path.stat().st_size-131072));tail=stream.read().decode('utf-8',errors='replace')
            bound=re.search(r'\n    "capacity_interval":\s*',tail)
            capacity=json.JSONDecoder().raw_decode(tail[bound.end():])[0] if bound else None
            if capacity is not None and all(key in job for key in ('file_id','id','created_at')):
                job['lower_bound']=capacity.get('lower_bound')
                return job
        except (ValueError,TypeError):pass
    job=json.loads(path.read_text(encoding='utf-8'))
    job['lower_bound']=(job.pop('result',{}) or {}).get('capacity_interval',{}).get('lower_bound')
    return job
