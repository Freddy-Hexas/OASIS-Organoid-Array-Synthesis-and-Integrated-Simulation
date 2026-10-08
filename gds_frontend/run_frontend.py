"""Run the same GDS-only geometry frontend on every GDS in a directory."""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path
from time import perf_counter

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np

from frontend import extract, _polygonal
from demo_router import route, write_demo_gds, audit_demo_gds
from process_geometry import ProcessRules, characterize
from workspace_paths import DATA_DIR, GEOMETRY_DIR


def _plain(value):
    if isinstance(value,np.ndarray):
        return value.tolist()
    if isinstance(value,(np.integer,np.floating)):
        return value.item()
    raise TypeError(type(value).__name__)


def save_graph(path,result):
    g=result.graph
    record={'input_sha256':result.summary['sha256'],
            'pitch_um':result.pitch,
            'nodes':[{'id':n,**d} for n,d in g.nodes(data=True)],
            'corridors':[{'start':u,'end':v,'key':k,**d}
                         for u,v,k,d in g.edges(keys=True,data=True)]}
    with gzip.open(path,'wt',encoding='utf-8') as f:
        json.dump(record,f,ensure_ascii=False,default=_plain,separators=(',',':'))


def preview(path,result):
    g=result.graph
    segments=[d['points_um'] for _,_,d in g.edges(data=True)]
    fig,axes=plt.subplots(1,2,figsize=(14,7))
    junction=np.asarray([d['xy_um'] for _,d in g.nodes(data=True)
                         if d['kind']=='junction'])
    outlets=np.asarray([d['xy_um'] for _,d in g.nodes(data=True)
                        if d['kind']=='collector_interface'])
    for ax in axes:
        left,top=result.origin
        ax.imshow(result.mask,cmap='Greys',vmin=0,vmax=4,origin='upper',
                  extent=[left-result.pitch/2,left+(result.mask.shape[1]-.5)*result.pitch,
                          top-(result.mask.shape[0]-.5)*result.pitch,top+result.pitch/2],
                  interpolation='nearest',alpha=.45)
        ax.add_collection(LineCollection(segments,linewidths=.7,colors='#50646e'))
        if len(junction):
            ax.scatter(junction[:,0],junction[:,1],s=9,color='#2b7499',zorder=3)
        if len(outlets):
            ax.scatter(outlets[:,0],outlets[:,1],s=20,color='#d04d2b',zorder=4)
        ax.set_aspect('equal')
        ax.set_xlabel('x (um)');ax.set_ylabel('y (um)')
    left,bottom,right,top=result.support.bounds
    axes[0].set_xlim(left-100,right+100);axes[0].set_ylim(bottom-100,top+100)
    axes[0].set_title('Recovered corridors and collector interfaces')
    if len(junction):
        lo=junction.min(axis=0);hi=junction.max(axis=0)
        center=(lo+hi)/2;span=max((hi-lo).max()*1.15,100)
        axes[1].set_xlim(center[0]-span/2,center[0]+span/2)
        axes[1].set_ylim(center[1]-span/2,center[1]+span/2)
    else:
        axes[1].set_xlim(left-100,right+100);axes[1].set_ylim(bottom-100,top+100)
    axes[1].set_title('Junction and loop detail')
    fig.tight_layout();fig.savefig(path,dpi=160);plt.close(fig)


def routing_preview(path,result,routing):
    fig,axes=plt.subplots(1,2,figsize=(14,7))
    for ax in axes:
        for polygon in _polygonal(result.support):
            ax.fill(*polygon.exterior.xy,color='#d9e2e6',lw=0)
            for hole in polygon.interiors:
                ax.fill(*hole.xy,color='white',lw=0)
        for i,r in enumerate(routing['_chosen']):
            pts=r['metal']
            color=plt.cm.turbo((i+.5)/max(1,len(routing['_chosen'])))
            for polygon in _polygonal(pts):
                ax.fill(*polygon.exterior.xy,color=color,lw=0)
            ax.scatter(*r['source_um'],s=14,color=color,edgecolor='black',linewidth=.25,zorder=5)
        ax.set_aspect('equal')
        ax.set_xlabel('x (um)');ax.set_ylabel('y (um)')
    left,bottom,right,top=result.support.bounds
    axes[0].set_xlim(left-100,right+100);axes[0].set_ylim(bottom-100,top+100)
    axes[0].set_title(f'{len(routing["_chosen"])} GDS-only metal paths: full device')
    junction=np.asarray([d['xy_um'] for _,d in result.graph.nodes(data=True)
                         if d['kind']=='junction'])
    if len(junction):
        lo=junction.min(axis=0);hi=junction.max(axis=0)
        center=(lo+hi)/2;span=max((hi-lo).max()*1.15,100)
        axes[1].set_xlim(center[0]-span/2,center[0]+span/2)
        axes[1].set_ylim(center[1]-span/2,center[1]+span/2)
    axes[1].set_title('Core detail (exact drawn metal width)')
    fig.tight_layout();fig.savefig(path,dpi=160);plt.close(fig)


def main():
    here=Path(__file__).resolve().parent
    parser=argparse.ArgumentParser()
    parser.add_argument('--input-dir',type=Path,default=DATA_DIR)
    parser.add_argument('--output-dir',type=Path,default=GEOMETRY_DIR)
    parser.add_argument('--rules',type=Path,default=here/'process_rules.json')
    parser.add_argument('--no-route',action='store_true')
    args=parser.parse_args()
    rules=ProcessRules(**json.loads(args.rules.read_text(encoding='utf-8')))
    rules.validate()
    files=sorted(args.input_dir.glob('*.gds'))
    if not files:
        parser.error('input directory contains no GDS files')
    args.output_dir.mkdir(parents=True,exist_ok=True)
    summary_path=args.output_dir/'summary.json'
    summary=json.loads(summary_path.read_text(encoding='utf-8')) if args.no_route and summary_path.exists() else {}
    for path in files:
        began=perf_counter()
        print('GDS',path.name,flush=True)
        result=extract(path,collector_clearance_um=rules.collector_clearance_um,
                       first_pitch=rules.first_pitch_um,minimum_pitch=rules.minimum_pitch_um)
        process,regions=characterize(result,rules)
        with gzip.open(args.output_dir/(path.stem+'_regions.json.gz'),'wt',encoding='utf-8') as f:
            json.dump(regions,f,ensure_ascii=False,separators=(',',':'))
        save_graph(args.output_dir/(path.stem+'_graph.json.gz'),result)
        preview(args.output_dir/(path.stem+'_graph.png'),result)
        item={'geometry':result.summary,'process_geometry':process}
        if args.no_route and path.name in summary:
            for key in ('routing','capacity_interval'):
                if key in summary[path.name]:
                    item[key]=summary[path.name][key]
        if not args.no_route:
            routing=route(result,electrode_diameter_um=rules.electrode_diameter_um,
                          wire_width_um=rules.wire_width_um,spacing_um=rules.spacing_um,
                          margin_um=rules.margin_um)
            if routing['status']=='checked_geometry_lower_bound':
                out=args.output_dir/(path.stem+'_gds_only_routing.gds')
                output_layers=write_demo_gds(path,out,routing)
                routing['gds_roundtrip_audit']=audit_demo_gds(out,result,routing,output_layers)
                routing['gds_output']=out.name
                routing_preview(args.output_dir/(path.stem+'_gds_only_routing.png'),result,routing)
            item['routing']={k:v for k,v in routing.items() if k!='_chosen'}
            item['capacity_interval']={
                'problem':'geometry-only electrodes to detected collector interfaces; no connector/pad assignment',
                'lower_bound':routing.get('retained_routes') if routing['status']=='checked_geometry_lower_bound' else None,
                'upper_bound':None,
                'upper_bound_status':'not established; sampled widths and finite graph optima are not continuous upper bounds',
                'external_pad_capacity':None,
                'source_domain':'finite junction/corridor candidates in automatically detected cycle-bearing narrow core',
                'routing_is_multi_track':False}
        item['elapsed_seconds']=round(perf_counter()-began,2)
        summary[path.name]=item
        summary_path.write_text(json.dumps(summary,ensure_ascii=False,
                                           indent=2,default=_plain),encoding='utf-8')
        print(json.dumps({'file':path.name,'pitch_um':result.pitch,
                          'nodes':result.graph.number_of_nodes(),
                          'edges':result.graph.number_of_edges(),
                          'full_cycles':result.summary['full_graph_cycle_rank'],
                          'collector_interfaces':result.summary['collector_interface_windows'],
                          'routed':item.get('routing',{}).get('retained_routes'),
                          'seconds':item['elapsed_seconds']},ensure_ascii=False),flush=True)
    print('SUMMARY',args.output_dir/'summary.json',flush=True)


if __name__=='__main__':
    main()
