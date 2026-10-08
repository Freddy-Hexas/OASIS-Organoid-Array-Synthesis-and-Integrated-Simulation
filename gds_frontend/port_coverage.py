"""Geometric outer-window inventory; coverage never means capacity."""
from collections import Counter
import math

import numpy as np
from shapely.geometry import LineString,Point
from shapely.strtree import STRtree

from frontend import _polygonal
from outer_exit_policy import outside_circle_edge_intervals


def outer_port_groups(nav,policy):
    """Group portals by connected exterior arcs of the guarded source domain.

    No angular bins, arm-count preset, labels, or generator connections are
    used. Arcs are real exterior portions outside the policy search circle.
    """
    arcs=[]
    for part in _polygonal(nav.navigation_domain):
        vertices=np.asarray(part.exterior.coords)
        groups=[]
        for a,b in zip(vertices[:-1],vertices[1:]):
            for aa,bb in policy.edge_intervals(a,b,search=True):
                if np.linalg.norm(bb-aa)<1e-9:continue
                if groups and np.linalg.norm(groups[-1][-1]-aa)<1e-8:
                    groups[-1].append(bb)
                else:groups.append([aa,bb])
        if len(groups)>1 and np.linalg.norm(groups[-1][-1]-groups[0][0])<1e-8:
            groups[0]=groups[-1]+groups[0][1:];groups.pop()
        arcs.extend(LineString(group) for group in groups)
    arcs.sort(key=lambda line:(math.atan2(line.centroid.y-policy.origin[1],
                                          line.centroid.x-policy.origin[0]),*line.bounds))
    tree=STRtree(arcs)
    mapping={};records=[]
    for outlet in nav.outlets:
        xy=Point(nav.graph.nodes[outlet]['xy_um'])
        group=int(tree.nearest(xy)) if arcs else None
        if group is not None and arcs[group].distance(xy)<=max(8*policy.grid_um,1e-7):
            mapping[outlet]=group
    for i,line in enumerate(arcs):
        middle=line.interpolate(.5,normalized=True)
        records.append({'group':i,'arc_length_um':line.length,'representative_um':list(middle.coords[0]),
                        'angle_deg':math.degrees(math.atan2(middle.y-policy.origin[1],middle.x-policy.origin[0])),
                        'outlets':[node for node in nav.outlets if mapping.get(node)==i]})
    return mapping,records


def coverage_record(chosen,mapping,groups):
    counts=Counter(mapping[c['outlet_node']] for c in chosen if c.get('outlet_node') in mapping)
    return {'method':'connected_outer_arcs_of_guarded_original_support',
            'geometric_window_count':len(groups),'used_window_count':len(counts),
            'unused_window_count':len(groups)-len(counts),
            'groups':[{**record,'selected_networks':counts.get(record['group'],0)} for record in groups],
            'scope':'exit diversity diagnostic; one window may carry multiple tracks; unused window is not an infeasibility certificate and window count is not a capacity upper bound'}
