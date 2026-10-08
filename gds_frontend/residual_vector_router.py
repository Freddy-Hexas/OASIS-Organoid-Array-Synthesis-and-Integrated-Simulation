"""Find additional full-width paths in the residual continuous wire domain.

The original navigation graph is intentionally not reused: occupied metal
changes the available topology, and a second track may lie on either side of
the old centerline. The constrained triangle dual only proposes a path.
Every returned polyline and its buffered metal are checked against the
original support and the occupied-metal spacing.
"""
from __future__ import annotations

from collections import defaultdict
from heapq import heappop, heappush
import math

import numpy as np
import shapely
from shapely.geometry import LineString, Point, Polygon
from shapely.ops import unary_union
from shapely.strtree import STRtree
from navigation_simplification import simplify_inside,check_deadline


class ResidualVectorRouter:
    def __init__(self,support,occupied_metals,*,wire_width_um=5.0,
                 spacing_um=4.0,margin_um=4.0,numeric_guard_um=.05,
                 target_points=None,deadline=None):
        if (wire_width_um<=0 or spacing_um<0 or margin_um<0 or
                numeric_guard_um<0):
            raise ValueError('Invalid residual routing clearances')
        self.support=support
        self.occupied=unary_union(list(occupied_metals)) if occupied_metals else Polygon()
        self.wire_width_um=float(wire_width_um)
        self.spacing_um=float(spacing_um)
        self.margin_um=float(margin_um)
        self.numeric_guard_um=float(numeric_guard_um)
        self._target_cache={}
        self.deadline=deadline
        check_deadline(deadline,'residual domain construction')
        self.wire_radius=self.wire_width_um/2
        center=support.buffer(-(self.wire_radius+self.margin_um+
                                self.numeric_guard_um),quad_segs=32)
        if not self.occupied.is_empty:
            forbidden=self.occupied.buffer(
                self.wire_radius+self.spacing_um+self.numeric_guard_um,
                quad_segs=32)
            center=center.difference(forbidden)
        # A query with fixed exits only needs their connected components.
        # This is an exact restriction of this residual query, not pruning
        # of the original support or a global infeasibility certificate.
        self.component_filter_certificate=None
        if target_points is not None:
            targets=[Point(point) for point in target_points]
            pieces=(list(center.geoms) if center.geom_type=='MultiPolygon'
                    else [center] if center.geom_type=='Polygon' else [])
            keep=[any(piece.covers(point) for point in targets) for piece in pieces]
            self.component_filter_certificate={
                'method':'fixed_target_connected_component_restriction',
                'target_points_um':[list(point.coords[0]) for point in targets],
                'residual_component_count':len(pieces),'retained_component_count':sum(keep),
                'excluded_components':[{'area_um2':piece.area,'bounds_um':list(piece.bounds),
                                        'reason':'component_contains_no_fixed_target'}
                                       for piece,allowed in zip(pieces,keep) if not allowed],
                'scope':'conditional reachability in guarded residual polygon geometry with all other nets fixed; not physical impossibility in the original or joint routing problem'}
            center=unary_union([piece for piece,allowed in zip(pieces,keep) if allowed])
        self.center_domain=center
        self.physical_center_domain=center
        center,self.simplification_certificate=simplify_inside(
            center,self.wire_width_um/4,protected_points=target_points or (),deadline=deadline)
        self.center_domain=center
        check_deadline(deadline,'residual triangulation')
        self.triangles=(list(shapely.constrained_delaunay_triangles(center).geoms)
                        if not center.is_empty else [])
        check_deadline(deadline,'residual triangulation completed')
        if not self.triangles:
            self.tree=None
            self.centers=np.empty((0,2))
            self.adjacency=[]
            return
        union=shapely.coverage_union_all(self.triangles)
        if union.symmetric_difference(center).area>1e-8*max(1,center.area):
            raise RuntimeError('Residual triangulation does not cover center domain')
        self.tree=STRtree(self.triangles)
        self.centers=np.asarray([
            np.asarray(triangle.exterior.coords)[:3].mean(axis=0)
            for triangle in self.triangles])
        ownership={}
        self.adjacency=[[] for _ in self.triangles]
        for i,triangle in enumerate(self.triangles):
            if i%1024==0:check_deadline(deadline,'residual mesh adjacency')
            vertices=np.asarray(triangle.exterior.coords)[:3]
            for a,b in zip(vertices,np.roll(vertices,-1,axis=0)):
                key=tuple(sorted((tuple(a),tuple(b))))
                if key not in ownership:
                    ownership[key]=i
                    continue
                other=ownership.pop(key)
                portal=(a+b)/2
                weight=float(np.linalg.norm(self.centers[i]-portal)+
                             np.linalg.norm(self.centers[other]-portal))
                self.adjacency[i].append((other,weight,portal))
                self.adjacency[other].append((i,weight,portal))

    def _triangle_ids(self,point):
        if self.tree is None:
            return []
        query=Point(point)
        return [int(i) for i in self.tree.query(query,predicate='intersects')
                if self.triangles[int(i)].covers(query)]

    def _target_tree(self,targets):
        """One reverse Dijkstra tree shared by every source using these Pads."""
        key=tuple(tuple(map(float,target)) for target in targets)
        cached=self._target_cache.get(key)
        if cached is not None:
            return cached
        queue=[]
        distance={}
        successor={}
        target_for={}
        for target_id,point in enumerate(targets):
            for triangle in self._triangle_ids(point):
                length=float(np.linalg.norm(self.centers[triangle]-point))
                if length<distance.get(triangle,math.inf):
                    distance[triangle]=length
                    target_for[triangle]=target_id
                    successor.pop(triangle,None)
                    heappush(queue,(length,triangle))
        while queue:
            if len(distance)%1024==0:check_deadline(self.deadline,'residual shortest paths')
            length,triangle=heappop(queue)
            if length>distance[triangle]+1e-9:
                continue
            for neighbor,weight,portal in self.adjacency[triangle]:
                candidate=length+weight
                if candidate<distance.get(neighbor,math.inf)-1e-9:
                    distance[neighbor]=candidate
                    target_for[neighbor]=target_for[triangle]
                    successor[neighbor]=(triangle,portal)
                    heappush(queue,(candidate,neighbor))
        result=(distance,successor,target_for)
        self._target_cache[key]=result
        return result

    def route(self,source,targets):
        """Return (polyline, target_index, audit) or None.

        Target order has no semantics. A route is accepted only after its
        entire metal-width buffer passes continuous support and spacing
        checks; the caller must still certify any electrode, bridge and Pad.
        """
        check_deadline(self.deadline,'residual route proposal')
        source=np.asarray(source,dtype=float)
        targets=[np.asarray(target,dtype=float) for target in targets]
        origins=self._triangle_ids(source)
        if not origins or not targets:
            return None
        distance,successor,target_for=self._target_tree(targets)
        if not distance:
            return None
        available=[(distance[triangle]+float(np.linalg.norm(
                    source-self.centers[triangle])),triangle)
                   for triangle in origins if triangle in distance]
        if not available:
            return None
        _,first=min(available)
        target_id=target_for[first]
        target=targets[target_id]
        chain=[first]
        while chain[-1] in successor:
            chain.append(successor[chain[-1]][0])
        # Consecutive portal midpoints lie in the same convex triangle.  The
        # midpoint chain is therefore inside the triangle corridor without
        # visiting triangle centroids.  Alternating centroid/portal points
        # created enormous zigzags in long thin triangles (a 400 um detour
        # could become several millimeters of metal).
        points=[source]
        for triangle in chain[:-1]:
            portal=successor[triangle][1]
            if np.linalg.norm(points[-1]-portal)>1e-8:
                points.append(portal)
        if np.linalg.norm(points[-1]-target)>1e-8:
            points.append(target)
        line=LineString(points)
        if not self.center_domain.buffer(1e-6).covers(line):
            raise RuntimeError('Triangle-dual residual path left the center domain')
        shortcut=line.simplify(self.wire_radius,preserve_topology=False)
        if (shortcut.geom_type=='LineString' and
                self.center_domain.buffer(1e-6).covers(shortcut)):
            line=shortcut
        metal=line.buffer(self.wire_radius,quad_segs=32,
                          cap_style='round',join_style='round')
        if not self.support.covers(metal.buffer(
                self.margin_um,quad_segs=32)):
            return None
        if (not self.occupied.is_empty and
                metal.distance(self.occupied)<
                self.spacing_um+self.numeric_guard_um-1e-6):
            return None
        return line,target_id,{
            'triangles':len(self.triangles),
            'crossed_triangles':len(chain),
            'length_um':line.length,
            'minimum_existing_metal_gap_um':(
                None if self.occupied.is_empty else
                float(metal.distance(self.occupied))),
            'minimum_support_margin_um':
                float(metal.distance(self.support.boundary)),
            'support_containment_verified':True,
            'metal_spacing_verified':True}
