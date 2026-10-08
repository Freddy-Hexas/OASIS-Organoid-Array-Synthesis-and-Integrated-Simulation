"""Propose curved source-to-shell bridges in the original exterior face.

The bridge centerline travels outside a bridge-radius expansion of the
original support.  Existing metal is an obstacle to the wire centerline.
This is a constructive search only: every returned bridge and complete metal
path still needs exported-GDS integer readback before it counts as a lower
bound.  Search failure never implies an upper bound.
"""
from __future__ import annotations

import math

import numpy as np
from shapely.geometry import LineString, Point, Polygon, box
from shapely.ops import unary_union

from residual_vector_router import ResidualVectorRouter
from outer_exit_policy import make_outer_exit_policy


def _prefix_to_radius(line,origin,radius):
    """First outward crossing of a circle, retaining exact polyline order."""
    points=[np.asarray(point,dtype=float) for point in line.coords]
    origin=np.asarray(origin,dtype=float)
    if np.linalg.norm(points[0]-origin)>=radius:
        return None
    for index,(a,b) in enumerate(zip(points,points[1:])):
        vector=b-a
        relative=a-origin
        qa=float(vector@vector)
        if qa<1e-16:
            continue
        qb=2*float(relative@vector)
        qc=float(relative@relative)-radius*radius
        discriminant=qb*qb-4*qa*qc
        if discriminant<0:
            continue
        root=(-qb+math.sqrt(discriminant))/(2*qa)
        if -1e-10<=root<=1+1e-10:
            crossing=a+min(1,max(0,root))*vector
            prefix=[*points[:index+1],crossing]
            if len(prefix)>=2:
                return LineString(prefix)
    return None


class ExteriorBridgeRouter:
    def __init__(self,source,frame,occupied_metals,*,bridge_width_um=40.0,
                 wire_width_um=5.0,spacing_um=4.0,numeric_guard_um=.05,
                 bridge_shell_overlap_um=100.0,angular_samples=144,
                 exit_policy=None):
        if (bridge_width_um<=0 or wire_width_um<=0 or spacing_um<0 or
                numeric_guard_um<0 or bridge_shell_overlap_um<=0 or
                angular_samples<4):
            raise ValueError('Invalid exterior bridge routing dimensions')
        self.source=source
        self.exit_policy=exit_policy or make_outer_exit_policy(source,
            wire_width_um=wire_width_um,numeric_guard_um=numeric_guard_um,
            bridge_width_um=bridge_width_um)
        self.frame=frame
        self.bridge_width_um=float(bridge_width_um)
        self.bridge_shell_overlap_um=float(bridge_shell_overlap_um)
        self.origin=np.asarray(frame['origin_um'],dtype=float)
        self.inner_radius=float(frame['inner_radius_um'])
        overlap=float(bridge_shell_overlap_um)
        self.target_radius=self.inner_radius+overlap
        outer=Point(self.origin).buffer(self.inner_radius+2*overlap,
                                         quad_segs=64)
        forbidden_source=source.buffer(bridge_width_um/2+numeric_guard_um,
                                       quad_segs=8)
        free=outer.difference(forbidden_source)
        parts=list(free.geoms) if hasattr(free,'geoms') else [free]
        exterior=[part for part in parts if isinstance(part,Polygon) and
                  part.intersects(outer.boundary)]
        if not exterior:
            raise ValueError('No original-support exterior face reaches the outer domain')
        metals=[part.intersection(outer) for part in occupied_metals
                if part.intersects(outer)]
        self.router=ResidualVectorRouter(
            unary_union(exterior),metals,wire_width_um=wire_width_um,
            spacing_um=spacing_um,margin_um=0,
            numeric_guard_um=numeric_guard_um)
        self.targets=[self.origin+self.target_radius*np.asarray((
            math.cos(2*math.pi*i/angular_samples),
            math.sin(2*math.pi*i/angular_samples)))
            for i in range(angular_samples)]
        source_parts=(list(source.geoms) if hasattr(source,'geoms') else
                      [source])
        self.protected_holes=[Polygon(ring) for part in source_parts
                              if isinstance(part,Polygon)
                              for ring in part.interiors]
    def propose(self,exit_point,interior_point,inside_line):
        """Yield validated geometry proposals from one source exterior point."""
        exit_point=np.asarray(exit_point,dtype=float)
        if not self.exit_policy.allows_exit(exit_point):
            return
        interior_point=np.asarray(interior_point,dtype=float)
        outward=exit_point-interior_point
        length=float(np.linalg.norm(outward))
        if length<1e-8:
            return
        outward/=length
        seen=set()
        for distance in (self.bridge_width_um*.6,
                         self.bridge_width_um,
                         2*self.bridge_width_um):
            start=exit_point+distance*outward
            if not self.router.center_domain.covers(Point(start)):
                continue
            result=self.router.route(start,self.targets)
            if result is None:
                continue
            routed,_,path_audit=result
            points=[exit_point,start,*list(routed.coords)[1:]]
            bridge_line=LineString(points)
            shell_point=np.asarray(bridge_line.coords[-1])
            key=(round(shell_point[0],3),round(shell_point[1],3))
            if key in seen:
                continue
            seen.add(key)
            if bridge_line.intersection(self.source).length>.002:
                continue
            bridge=bridge_line.buffer(self.bridge_width_um/2,
                                      quad_segs=32,cap_style='round',
                                      join_style='round')
            if (bridge.intersection(self.source).area<=1 or
                    bridge.intersection(self.frame['shell']).area<=1):
                continue
            added=bridge.difference(self.source)
            if any(added.intersects(hole) for hole in self.protected_holes):
                continue
            allowed,_=self.exit_policy.check_bridge(self.source,bridge,exit_point=exit_point)
            if not allowed:continue
            yield {'inside_line':inside_line,
                   'exit_point':exit_point,
                   'shell_point':shell_point,
                   'bridge_line':bridge_line,
                   'bridge':bridge,
                   'external_path_audit':path_audit,
                   'external_strategy':'continuous_exterior_residual_triangles'}


class DirectExteriorPadRouter:
    """Join an exposed source boundary directly to an available Pad slot.

    The complete path is proposed in one continuous residual domain, so a
    source-to-shell split cannot discard paths whose first shell arrival is
    in a different metal component from the assigned Pad.
    """
    def __init__(self,source,frame,occupied_metals,*,bridge_width_um=40.0,
                 wire_width_um=5.0,spacing_um=4.0,margin_um=4.0,
                 numeric_guard_um=.05,bridge_shell_overlap_um=100.0,
                 exit_policy=None):
        self.source=source
        self.exit_policy=exit_policy or make_outer_exit_policy(source,
            wire_width_um=wire_width_um,margin_um=margin_um,
            numeric_guard_um=numeric_guard_um,bridge_width_um=bridge_width_um)
        self.frame=frame
        self.bridge_width_um=float(bridge_width_um)
        self.bridge_shell_overlap_um=float(bridge_shell_overlap_um)
        self.wire_width_um=float(wire_width_um)
        self.spacing_um=float(spacing_um)
        self.margin_um=float(margin_um)
        self.guard=float(numeric_guard_um)
        ox,oy=frame['origin_um']
        self.origin=np.asarray((ox,oy),dtype=float)
        half=frame['square_side_um']/2
        self.square=box(ox-half,oy-half,ox+half,oy+half)
        # A full 32 mm square, re-triangulated for each residual search, is
        # needlessly large.  Each Pad search only needs a corridor from the
        # complete source envelope to that Pad.  This is a lower-bound search:
        # restricting the proposal domain cannot strengthen a capacity proof.
        self.occupied=(unary_union(list(occupied_metals)) if occupied_metals
                       else Polygon())
        corners=[(x,y) for x in (source.bounds[0],source.bounds[2])
                       for y in (source.bounds[1],source.bounds[3])]
        self.corridor_radius=max(math.dist((ox,oy),point) for point in corners)+(
            2*bridge_width_um+wire_width_um+spacing_um+margin_um)
        self._source_exclusion=source.buffer(
            bridge_width_um/2+numeric_guard_um,quad_segs=8)
        self._router_cache={}
        source_parts=(list(source.geoms) if hasattr(source,'geoms') else
                      [source])
        self.protected_holes=[Polygon(ring) for part in source_parts
                              if isinstance(part,Polygon)
                              for ring in part.interiors]
        self.diagnostics={'starts_outside_domain':0,
                          'route_not_found':0,
                          'route_reenters_source':0,
                          'bridge_invalid':0,
                          'bridge_enters_hole':0,
                          'pad_marker_existing_gap':0,
                          'wire_existing_gap':0,
                          'wire_pad_contact_failed':0,
                          'wire_or_metal_gap':0,
                          'support_margin':0,
                          'local_triangulations':0,
                          'largest_local_triangle_count':0,
                          'largest_local_corridor_area_um2':0.0,
                          'outer_exit_policy_rejected':0,
                          'accepted':0}

    def _router_for_slot(self,slot):
        key=(slot['pad_id'],tuple(map(float,slot['target_um'])))
        cached=self._router_cache.get(key)
        if cached is not None:
            return cached
        target=tuple(map(float,slot['target_um']))
        corridor=LineString([self.origin,target]).buffer(
            self.corridor_radius,quad_segs=8).intersection(self.square)
        free=corridor.difference(self._source_exclusion)
        target_point=Point(target)
        parts=list(free.geoms) if hasattr(free,'geoms') else [free]
        exterior=[part for part in parts if isinstance(part,Polygon) and
                  part.covers(target_point)]
        if not exterior:
            return None
        domain=unary_union(exterior)
        obstacle_window=corridor.buffer(
            self.wire_width_um+self.spacing_um+self.guard,quad_segs=4)
        local_occupied=(self.occupied.intersection(obstacle_window)
                        if not self.occupied.is_empty else Polygon())
        router=ResidualVectorRouter(
            domain,[local_occupied] if not local_occupied.is_empty else [],
            wire_width_um=self.wire_width_um,spacing_um=self.spacing_um,
            margin_um=0,numeric_guard_um=self.guard)
        self.diagnostics['local_triangulations']+=1
        self.diagnostics['largest_local_triangle_count']=max(
            self.diagnostics['largest_local_triangle_count'],
            len(router.triangles))
        self.diagnostics['largest_local_corridor_area_um2']=max(
            self.diagnostics['largest_local_corridor_area_um2'],
            float(corridor.area))
        # The usual residual call has only a handful of adjacent Pad slots.
        # A bounded cache also works when a larger bank is passed in.
        if len(self._router_cache)>=8:
            self._router_cache.pop(next(iter(self._router_cache)))
        self._router_cache[key]=router
        return router

    def propose(self,exit_point,interior_point,inside_line,slots,
                *,pad_target_trials=3):
        if not slots:
            return
        exit_point=np.asarray(exit_point,dtype=float)
        if not self.exit_policy.allows_exit(exit_point):
            self.diagnostics['outer_exit_policy_rejected']+=1
            return
        interior_point=np.asarray(interior_point,dtype=float)
        outward=exit_point-interior_point
        length=float(np.linalg.norm(outward))
        if length<1e-8:
            return
        outward/=length
        if pad_target_trials<1:
            raise ValueError('pad_target_trials must be positive')
        diagnostics=self.diagnostics
        ordered=sorted(slots,key=lambda slot:math.dist(
            exit_point,slot['target_um']))[:pad_target_trials]
        for distance in (self.bridge_width_um*.6,
                         self.bridge_width_um,
                         2*self.bridge_width_um):
            for side_shift in (0,-self.bridge_width_um,
                               self.bridge_width_um):
                tangent=np.asarray((-outward[1],outward[0]))
                start=exit_point+distance*outward+side_shift*tangent
                for slot in ordered:
                    router=self._router_for_slot(slot)
                    if router is None or not router.center_domain.covers(Point(start)):
                        diagnostics['starts_outside_domain']+=1
                        continue
                    result=router.route(start,[slot['target_um']])
                    if result is None:
                        diagnostics['route_not_found']+=1
                        continue
                    routed,_,path_audit=result
                    exterior_line=LineString([exit_point,start,*list(routed.coords)[1:]])
                    if exterior_line.intersection(self.source).length>.002:
                        diagnostics['route_reenters_source']+=1
                        continue
                    bridge_line=_prefix_to_radius(
                        exterior_line,self.frame['origin_um'],
                        self.frame['inner_radius_um']+
                        self.bridge_shell_overlap_um)
                    if bridge_line is None:
                        diagnostics['bridge_invalid']+=1
                        continue
                    bridge=bridge_line.buffer(self.bridge_width_um/2,
                                              quad_segs=32,cap_style='round',
                                              join_style='round')
                    if (not self.square.covers(bridge) or
                            bridge.intersection(self.source).area<=1 or
                            bridge.intersection(self.frame['shell']).area<=1):
                        diagnostics['bridge_invalid']+=1
                        continue
                    added=bridge.difference(self.source)
                    if any(added.intersects(hole) for hole in self.protected_holes):
                        diagnostics['bridge_enters_hole']+=1
                        continue
                    allowed,_=self.exit_policy.check_bridge(self.source,bridge,
                                                            exit_point=exit_point)
                    if not allowed:
                        diagnostics['outer_exit_policy_rejected']+=1;continue
                    full_line=LineString([*inside_line.coords,
                                          *list(exterior_line.coords)[1:]])
                    allowed,_=self.exit_policy.check_exterior_line(self.source,full_line)
                    if not allowed:
                        diagnostics['outer_exit_policy_rejected']+=1;continue
                    wire=full_line.buffer(self.wire_width_um/2+.01,
                                          quad_segs=32,cap_style='round',
                                          join_style='round')
                    if not wire.intersects(slot['polygon']):
                        diagnostics['wire_pad_contact_failed']+=1
                        diagnostics['wire_or_metal_gap']+=1
                        continue
                    if wire.distance(self.occupied)<self.spacing_um+self.guard:
                        diagnostics['wire_existing_gap']+=1
                        diagnostics['wire_or_metal_gap']+=1
                        continue
                    if slot['polygon'].distance(self.occupied)<self.spacing_um+self.guard:
                        diagnostics['pad_marker_existing_gap']+=1
                        diagnostics['wire_or_metal_gap']+=1
                        continue
                    outer_metal=unary_union([wire,slot['polygon']])
                    support=unary_union([self.source,bridge,self.frame['shell']])
                    if not support.covers(outer_metal.buffer(
                            self.margin_um+.005,quad_segs=32)):
                        diagnostics['support_margin']+=1
                        continue
                    diagnostics['accepted']+=1
                    yield {'pad_id':slot['pad_id'],'pad_side':slot['side'],
                           'pad_index':slot['index'],'pad':slot['polygon'],
                           'pad_target_um':slot['target_um'],
                           'bridge':bridge,'outer_wire':wire,
                           'outer_line':full_line,
                           'outer_curve':{'method':'joint_exterior_to_pad_residual_triangles',
                                          'path_audit':path_audit},
                           'outer_length_um':full_line.length,
                           'portal_escape_um':exit_point.tolist(),
                           'outer_route_strategy':'joint_exterior_to_pad_residual_triangles'}
