"""Vector feasible regions, sampled corridor widths, and ordered node ports.

Sampled width budgets are local diagnostics, never continuous capacity bounds.
Node windows expose cyclic order for the next multi-track solver; no legal
multi-net turns are inferred merely from a window's existence.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import math

import numpy as np
import shapely
from shapely.geometry import LineString, Point
from shapely.strtree import STRtree

from frontend import FrontendResult, _polygonal
from electrode_region import region_for


@dataclass(frozen=True)
class ProcessRules:
    electrode_diameter_um: float = 30.0
    minimum_center_spacing_um: float = 70.0
    electrode_region_radius_um: float = 3000.0
    wire_width_um: float = 5.0
    spacing_um: float = 4.0
    margin_um: float = 4.0
    collector_clearance_um: float = 25.0
    first_pitch_um: float = 2.0
    minimum_pitch_um: float = 0.5
    width_sample_step_um: float = 20.0
    node_window_radius_um: float = 12.0

    def validate(self):
        values=asdict(self)
        for name,value in values.items():
            if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value):
                raise ValueError(f'Rule {name} must be finite')
            if value<0 or (value==0 and name not in ('spacing_um','margin_um','minimum_center_spacing_um')):
                raise ValueError(f'Rule {name} must be positive (spacing/margin/center spacing may be zero)')
        if self.minimum_pitch_um>self.first_pitch_um:
            raise ValueError('minimum_pitch_um exceeds first_pitch_um')


class BoundaryIndex:
    """Indexed vector segments, avoiding distance to the whole giant boundary."""
    def __init__(self,support):
        records=[]
        for polygon in _polygonal(support):
            for ring in (polygon.exterior,*polygon.interiors):
                p=np.asarray(ring.coords)
                valid=np.linalg.norm(np.diff(p,axis=0),axis=1)>1e-10
                records.append(np.stack([p[:-1][valid],p[1:][valid]],axis=1))
        self.segments=np.concatenate(records,axis=0)
        self.tree=STRtree(shapely.linestrings(self.segments))

    def widths(self,point,tangent,half_length):
        normal=np.asarray([-tangent[1],tangent[0]])
        query=LineString([point-half_length*normal,point+half_length*normal])
        ids=self.tree.query(query)
        segments=self.segments[ids]
        aa=segments[:,0]-point
        dv=segments[:,1]-segments[:,0]
        denom=dv@tangent
        use=np.abs(denom)>1e-10
        aa=aa[use];dv=dv[use];denom=denom[use]
        frac=-(aa@tangent)/denom
        positions=(aa+frac[:,None]*dv)@normal
        hits=positions[(frac>=-1e-9)&(frac<=1+1e-9)&(np.abs(positions)<=half_length+1e-6)]
        neg=hits[hits<=-1e-9];pos=hits[hits>=1e-9]
        if not len(neg) or not len(pos):
            return None
        low=float(neg.max());high=float(pos.min())
        return low,high,normal


def characterize(result:FrontendResult,rules:ProcessRules):
    rules.validate()
    support=result.support
    shapely.prepare(support)
    electrode_region=support.buffer(-(rules.electrode_diameter_um/2+rules.margin_um),quad_segs=32)
    placement_region=region_for(support,rules,result.summary.get('gds_native_precision_m',1e-9)*1e6)
    electrode_region=placement_region.clip(electrode_region)
    wire_region=support.buffer(-(rules.wire_width_um/2+rules.margin_um),quad_segs=32)
    # The regions are approximate polygonal erosions. Actual metal containment
    # and actual boundary distance remain the acceptance tests.
    index=BoundaryIndex(support)
    legal_edges=0
    nominal_lanes=[]
    ordered_windows=0
    for eid,(u,v,k,data) in enumerate(result.graph.edges(keys=True,data=True)):
        data['corridor_id']=eid
        pts=data['points_um']
        cumulative=np.r_[0,np.cumsum(np.linalg.norm(np.diff(pts,axis=0),axis=1))]
        samples=np.arange(0,cumulative[-1]+1e-7,rules.width_sample_step_um)
        samples=np.r_[samples,cumulative[-1]]
        profiles=[]
        for length in np.unique(samples):
            i=int(np.clip(np.searchsorted(cumulative,length),1,len(pts)-1))
            denom=cumulative[i]-cumulative[i-1]
            fraction=(length-cumulative[i-1])/denom if denom>0 else 0
            point=(1-fraction)*pts[i-1]+fraction*pts[i]
            a=max(0,i-4);b=min(len(pts)-1,i+4)
            tangent=pts[b]-pts[a]
            norm=np.linalg.norm(tangent)
            if norm<1e-10 or not support.covers(Point(point)):
                continue
            tangent=tangent/norm
            cut=index.widths(point,tangent,max(30,4*data['median_raster_clearance_um']))
            if cut is None:
                continue
            lo,hi,normal=cut
            width=hi-lo
            capacity=max(0,math.floor((width-2*rules.margin_um+rules.spacing_um+1e-8)/
                                      (rules.wire_width_um+rules.spacing_um)))
            profiles.append({'arc_length_um':float(length),'point_um':point.tolist(),
                             'normal':normal.tolist(),'left_extent_um':lo,'right_extent_um':hi,
                             'support_width_um':width,'nominal_tracks':capacity})
        line=LineString(pts).simplify(result.pitch*.25,preserve_topology=False)
        valid=bool(support.covers(line.buffer(rules.wire_width_um/2+rules.margin_um,quad_segs=8)))
        legal_edges+=int(valid)
        data['vector_single_track_legal']=valid
        data['sampled_width_profile']=profiles
        data['minimum_sampled_width_um']=min((p['support_width_um'] for p in profiles),default=None)
        data['sampled_track_budget']=min((p['nominal_tracks'] for p in profiles),default=None)
        data['track_budget_is_continuous_upper_bound']=False
        if data['sampled_track_budget'] is not None:
            nominal_lanes.append(data['sampled_track_budget'])
    for node,data in result.graph.nodes(data=True):
        center=np.asarray(data['xy_um'])
        window=support.intersection(Point(center).buffer(rules.node_window_radius_um,quad_segs=32))
        ports=[]
        for u,v,k,edge in result.graph.edges(node,keys=True,data=True):
            pts=edge['points_um']
            if np.linalg.norm(pts[-1]-center)<np.linalg.norm(pts[0]-center):
                pts=pts[::-1]
            radii=np.linalg.norm(pts-center,axis=1)
            crossing=np.flatnonzero(radii>=rules.node_window_radius_um)
            point=pts[crossing[0]] if len(crossing) else pts[-1]
            direction=point-center
            ports.append({'corridor_id':edge['corridor_id'],'other_node':v if u==node else u,
                          'point_um':point.tolist(),
                          'angle_rad':float(math.atan2(direction[1],direction[0])),
                          'reaches_window_boundary':bool(len(crossing)),
                          'nominal_tracks':edge['sampled_track_budget']})
        ports.sort(key=lambda p:p['angle_rad'])
        data['window_radius_um']=rules.node_window_radius_um
        data['window_area_um2']=float(window.area)
        data['window_wkb_hex']=window.wkb_hex
        data['ports_ccw']=ports
        data['turn_states_solved']=False
        ordered_windows+=int(len(ports)>=3)
    counts={str(k):nominal_lanes.count(k) for k in sorted(set(nominal_lanes))}
    summary={'rules':asdict(rules),
             'electrode_region':placement_region.record(),
             'electrode_center_region_area_um2':float(electrode_region.area),
             'wire_center_region_area_um2':float(wire_region.area),
             'vector_legal_single_track_corridors':legal_edges,
             'vector_rejected_single_track_corridors':result.graph.number_of_edges()-legal_edges,
             'ordered_junction_windows':ordered_windows,
             'sampled_nominal_track_budget_histogram':counts,
             'regions_use_full_support_no_structure_specific_roi':True,
             'sampled_widths_are_capacity_bounds':False,
             'node_turn_combinations_solved':False}
    region_record={'rules':asdict(rules),
                   'support_wkb_hex':support.wkb_hex,
                   'electrode_center_region_wkb_hex':electrode_region.wkb_hex,
                   'wire_center_region_wkb_hex':wire_region.wkb_hex,
                   'polygonal_erosion_quad_segs':32,
                   'source_domain':'full support; routing candidates use the cycle-bearing narrow core'}
    return summary,region_record
