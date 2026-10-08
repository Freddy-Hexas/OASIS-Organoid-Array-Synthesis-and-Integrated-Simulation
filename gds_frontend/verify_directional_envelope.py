"""Unknown rectangular lattice: side exits, joint routes and independent GDS DRC."""
from pathlib import Path
from tempfile import TemporaryDirectory
import json
from copy import deepcopy
from dataclasses import replace

import gdstk
import numpy as np
from shapely.geometry import LineString
from shapely.ops import unary_union
from shapely import affinity

from exact_gds_audit import audit_exported_gds
from island_router import IslandSettings,build_navigation
from outer_exit_policy import make_outer_exit_policy,OuterExitPolicy
from pad_router import PadSettings,solve_with_pads,write_pad_gds,audit_pad_gds
from process_geometry import ProcessRules
from route_incumbent import recertify_incumbent


def main():
    rules=ProcessRules();settings=IslandSettings();pads=PadSettings()
    xs=np.arange(-520,521,130);ys=np.arange(-390,391,130)
    support=unary_union([LineString([(x,-480),(x,480)]).buffer(7.5,cap_style='flat') for x in xs]+
                        [LineString([(-610,y),(610,y)]).buffer(7.5,cap_style='flat') for y in ys])
    policy=make_outer_exit_policy(support)
    legacy=OuterExitPolicy(policy.grid_um,policy.center_ticks,policy.radius_ticks,
                           policy.maximum_vertex_radius_squared_ticks,5,4,.05,40)
    assert policy.allows_exit((610,0)) and not legacy.allows_exit((610,0))
    assert not policy.allows_exit((130,7.5))
    transformed=make_outer_exit_policy(affinity.rotate(affinity.translate(support,133,-77),37,origin=(133,-77)))
    point=affinity.rotate(affinity.translate(LineString([(610,0),(611,0)]),133,-77),37,origin=(133,-77))
    assert transformed.allows_exit(point.coords[0])
    with TemporaryDirectory(prefix='unknown_directional_lattice_') as tmp:
        source=Path(tmp)/'unclassified.gds';output=Path(tmp)/'routed.gds'
        lib=gdstk.Library(unit=1e-6,precision=1e-9);cell=lib.new_cell('UNKNOWN')
        for p in gdstk.boolean([gdstk.Polygon(support.exterior.coords)],
                 [gdstk.Polygon(r.coords) for r in support.interiors],'not',precision=.001,layer=42):cell.add(p)
        lib.write_gds(str(source))
        nav=build_navigation(source,rules,layer=42,outlet_mode='external_boundary',settings=settings,outer_exit_policy=policy)
        assert len(nav.outlets)>=20
        result=solve_with_pads(nav,rules,settings=settings,pad_settings=pads)
        assert result['retained_routes']>=8,result['joint_path_flow']
        layers=write_pad_gds(source,output,result,(42,0))
        audit_pad_gds(output,nav,result,layers)
        audit=audit_exported_gds(source,output,layers,wire_spacing_um=4,metal_support_margin_um=4,
            expected_nets=result['retained_routes'],minimum_electrode_diameter_um=30,
            minimum_electrode_center_distance_um=70,minimum_wire_width_um=5,
            minimum_pad_short_side_um=500,minimum_pad_long_side_um=3000)
        assert audit['outer_exit_policy_verified']
        assert result['optimization']['finite_upper_bound'] is None
        report={'geometry':nav.summary,'selected_support_layer':[42,0],'routing':result}
        from pad_router import make_pad_frame
        frame=make_pad_frame(nav.support,pads,spacing_um=4,margin_um=4,exit_policy=policy)
        warm,_,warm_diagnostic=recertify_incumbent([{'job_id':'unknown_previous','report':report}],
            nav,rules,settings,pads,frame,policy)
        assert len(warm)==result['retained_routes'],warm_diagnostic
        sparse_pads=replace(pads,pad_pitch_um=4000)
        sparse_frame=make_pad_frame(nav.support,sparse_pads,spacing_um=4,margin_um=4,exit_policy=policy)
        grown,grown_frame,grown_diagnostic=recertify_incumbent([{'job_id':'unknown_previous','report':report}],
            nav,rules,settings,sparse_pads,sparse_frame,policy)
        assert len(grown)==len(warm) and grown_frame['square_side_um']>sparse_frame['square_side_um'],grown_diagnostic
        compact_pads=replace(sparse_pads,sizing_mode='compact_then_expand',
                             minimum_pad_width_um=300,minimum_pad_length_um=1800)
        compact=solve_with_pads(nav,rules,settings=settings,pad_settings=compact_pads)
        assert compact['retained_routes']>=len(warm)
        actual=compact['pad_settings']
        assert 300<=actual['pad_width_um']<=500 and 1800<=actual['pad_length_um']<=3000,actual
        assert compact['pad_sizing']['pad_dimensions_reduced'],compact['pad_sizing']
        compact_output=Path(tmp)/'compact.gds'
        compact_layers=write_pad_gds(source,compact_output,compact,(42,0))
        audit_pad_gds(compact_output,nav,compact,compact_layers)
        compact_audit=audit_exported_gds(source,compact_output,compact_layers,wire_spacing_um=4,
            metal_support_margin_um=4,expected_nets=compact['retained_routes'],
            minimum_electrode_diameter_um=30,minimum_electrode_center_distance_um=70,
            minimum_wire_width_um=5,minimum_pad_short_side_um=actual['pad_width_um'],
            minimum_pad_long_side_um=actual['pad_length_um'])
        assert compact_audit['passed']
        tampered=deepcopy(report)
        tampered['routing']['routes'][0]['points_um'][0][0]+=10000
        invalid,_,invalid_diagnostic=recertify_incumbent([{'job_id':'tampered','report':tampered}],
            nav,rules,settings,pads,frame,policy)
        assert not invalid,invalid_diagnostic
        print(json.dumps({'passed':True,'portals':len(nav.outlets),'complete_nets':result['retained_routes'],
                          'joint':result['joint_path_flow'],'integer_audit':True,
                          'rechecked_warm_start_count':len(warm),'warm_start_pad_frame_expands':True,
                          'adaptive_pad_dimensions_um':[actual['pad_width_um'],actual['pad_length_um']],
                          'adaptive_pad_integer_audit':compact_audit['passed'],
                          'tampered_warm_start_rejected':True},ensure_ascii=False,default=str))


if __name__=='__main__':main()
