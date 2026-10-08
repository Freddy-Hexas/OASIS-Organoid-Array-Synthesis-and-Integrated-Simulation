"""End-to-end physical Pad witness on an unfamiliar support-only GDS."""
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import patch
import json

import gdstk
from shapely.geometry import Point, box

from island_router import IslandSettings,build_navigation,_select
from pad_router import (PadSettings,_pad_occupancy,make_pad_frame,solve_with_pads,
                        write_pad_gds,audit_pad_gds,_diverse_pad_columns)
from process_geometry import ProcessRules
import pad_router
import island_router
from exact_gds_audit import audit_exported_gds,audit_support_polygon_provenance
from ordered_pad_fanout import ordered_pad_assignment


def main():
    ranked=[{'source_node':i,'complete_length_um':i+1.0,'depth_um':i+1.0}
            for i in range(20)]
    diversified=_diverse_pad_columns(ranked,12,6)
    assert len(diversified)==18
    assert [item['source_node'] for item in diversified[:12]]==list(range(12))
    assert {item['source_node'] for item in diversified[12:]}==set(range(14,20))
    support=Point(0,0).buffer(3500)
    reference=make_pad_frame(support,PadSettings())
    expanded=make_pad_frame(support,PadSettings(square_side_um=40000))
    assert reference['pad_slots_per_side']==26
    assert expanded['pad_slots_per_side']==34 and len(expanded['slots'])==136
    assert expanded['corner_keepout_um']>=3104
    assert all(expanded['shell'].buffer(-4).covers(slot['polygon'])
               for slot in expanded['slots'])
    tight_frame=make_pad_frame(support,PadSettings(pad_pitch_um=504))
    tight_top=[slot for slot in tight_frame['slots'] if slot['side']=='top']
    assert tight_top[0]['polygon'].distance(tight_top[1]['polygon'])>=4.05
    odd_frame=make_pad_frame(support,PadSettings(pad_width_um=500.001,
        pad_length_um=3000.001,pad_pitch_um=1000.001))
    assert all(abs(odd_frame['effective_pad_dimensions'][key]-value)<1e-8 for key,value in
        {'pad_width_um':500.002,'pad_length_um':3000.002,'pad_pitch_um':1000.002}.items())
    sparse_settings=PadSettings(pad_pitch_um=12000)
    sparse_frame=make_pad_frame(support,sparse_settings)
    exits=[(x,1000) for x in (-150,-50,50,150)]+[(1000,y) for y in (-200,-100,0,100,200)]
    assignment,growth=ordered_pad_assignment(exits,sparse_frame)
    assert assignment is None and growth['required']==5,growth
    assert growth['required_slots_by_side']['top']==4,growth
    grow_side=2*sparse_frame['corner_keepout_um']+sparse_settings.pad_width_um+4*sparse_settings.pad_pitch_um
    grown=make_pad_frame(support,sparse_settings,minimum_side_um=grow_side)
    assignment,growth=ordered_pad_assignment(exits,grown)
    assert assignment is not None and len(assignment)==9,growth
    from dataclasses import replace
    compact_settings=replace(PadSettings(),sizing_mode='compact_then_expand',
                             minimum_pad_width_um=300,minimum_pad_length_um=1800)
    compact_frame=make_pad_frame(support,compact_settings,required_slots_per_side=40)
    assert compact_frame['pad_slots_per_side']>=40
    assert compact_frame['square_side_um']==32000
    actual=compact_frame['effective_pad_dimensions']
    assert 300<=actual['pad_width_um']<500 and 1800<=actual['pad_length_um']<3000
    huge_frame=make_pad_frame(support,compact_settings,required_slots_per_side=100)
    assert huge_frame['pad_slots_per_side']>=100 and huge_frame['square_side_um']>32000
    assert huge_frame['effective_pad_dimensions']['pad_width_um']==300
    assert huge_frame['effective_pad_dimensions']['pad_length_um']==1800
    # A wider process wire must also raise the Pad's contact-size floor in
    # direct Python calls, rather than only being checked by the HTTP parser.
    too_small=replace(compact_settings,minimum_pad_width_um=20,minimum_pad_length_um=120)
    try:make_pad_frame(support,too_small,wire_width_um=20,required_slots_per_side=40)
    except ValueError as exc:assert 'lower dimensions' in str(exc),exc
    else:raise AssertionError('Pad floor smaller than the actual wire contact was accepted')
    try:solve_with_pads(SimpleNamespace(),ProcessRules(wire_width_um=20),
                        settings=IslandSettings(),pad_settings=too_small)
    except ValueError as exc:assert 'lower dimensions' in str(exc),exc
    else:raise AssertionError('Direct solver did not validate the actual wire contact')
    try:make_pad_frame(support,PadSettings(),wire_width_um=38)
    except ValueError as exc:assert 'Bridge is narrower' in str(exc),exc
    else:raise AssertionError('Bridge smaller than the actual wire plus margins was accepted')
    # Exercise automatic growth after a bank fills.  The dummy route selector
    # isolates frame policy from the expensive physical routing fixture below.
    fake_columns=[{'source_node':i} for i in range(30)]
    def fake_frame_selection(nav,rules,settings,pad_settings,frame,central,info,bound,notify):
        count=min(frame['pad_slots_per_side'],sum(map(len,central)))
        banks={side:{'connected_pad_count':count if side=='top' else 0}
               for side in ('top','right','bottom','left')}
        return {'retained_routes':count,'pad_bank_occupancy':banks,'_chosen':[]}
    with (patch.object(pad_router,'_route_columns',return_value=(fake_columns,{})),
          patch.object(pad_router,'continuous_center_upper_bound',return_value={'value':100}),
          patch.object(pad_router,'_solve_on_pad_frame',side_effect=fake_frame_selection)):
        adaptive=solve_with_pads(SimpleNamespace(support=support,navigation_domain=support,
                                                  summary={'outlet_mode':'external_boundary'}),
                                 ProcessRules(),settings=IslandSettings())
    assert [trial['available_slots_per_side'] for trial in adaptive['pad_frame_trials']]==[26,30]
    assert adaptive['retained_routes']==30
    # The independent-column construction must apply the same compact policy
    # as the joint-flow construction, before deciding to increase the frame.
    with (patch.object(pad_router,'_route_columns',return_value=(fake_columns,{})),
          patch.object(pad_router,'continuous_center_upper_bound',return_value={'value':100}),
          patch.object(pad_router,'_solve_on_pad_frame',side_effect=fake_frame_selection)):
        compact_adaptive=solve_with_pads(SimpleNamespace(support=support,navigation_domain=support,
                                                          summary={'outlet_mode':'external_boundary'}),
                                         ProcessRules(),settings=IslandSettings(),pad_settings=compact_settings)
    assert compact_adaptive['retained_routes']==30
    assert compact_adaptive['pad_frame_trials'][-1]['square_side_um']==32000
    # A finite Pad library with an unavailable middle slot used to leave a
    # visible hole.  Contiguity is part of the MILP, before GDS export.
    columns=[]
    for index in (1,3,4,5):
        x=index*1000.0
        columns.append({'source_node':index,'source_um':(x,0),
                        'metal':box(x-10,-10,x+10,10),
                        'island':box(x-20,-20,x+20,20),'depth_um':float(index),
                        'pad_id':f'top-{index:02d}','pad_side':'top','pad_index':index})
    order={'top':[f'top-{index:02d}' for index in range(1,6)]}
    selected,info=_select(columns,ProcessRules(),IslandSettings(milp_time_limit_s=5),
                          lambda *args:None,pad_order=order)
    assert len(selected)==3 and [c['pad_index'] for c in selected]==[3,4,5]
    assert info['pad_contiguous_per_side']
    assert _pad_occupancy(selected,PadSettings())['top']['internal_empty_slots']==0
    four=[]
    for index,side in enumerate(('top','right','bottom','left')):
        x=index*6000.0
        four.append({'source_node':index,'source_um':(x,0),
                     'metal':box(x-10,-10,x+10,10),
                     'island':box(x-20,-20,x+20,20),
                     'depth_um':float(index+1),
                     'pad_id':f'{side}-01','pad_side':side,'pad_index':1})
    four_order={side:[f'{side}-01'] for side in ('top','right','bottom','left')}
    four_selected,four_info=_select(four,ProcessRules(),IslandSettings(),
                                    lambda *args:None,center_upper=4,
                                    pad_order=four_order)
    assert len(four_selected)==4 and four_info['finite_model_proven_optimal']
    with patch.object(island_router,'milp',return_value=SimpleNamespace(
            x=None,status=1,message='time limit')):
        fallback,fallback_info=_select(four[:1],ProcessRules(),IslandSettings(),
                                       lambda *args:None,center_upper=2)
    assert len(fallback)==1 and fallback_info['finite_lower_bound']==1
    many=[]
    for index in range(1,35):
        x=index*1000.0
        many.append({'source_node':index,'source_um':(x,0),
                     'metal':box(x-10,-10,x+10,10),
                     'island':box(x-20,-20,x+20,20),'depth_um':float(index),
                     'pad_id':f'top-{index:02d}','pad_side':'top','pad_index':index})
    many_order={'top':[f'top-{index:02d}' for index in range(1,35)]}
    many_selected,many_info=_select(many,ProcessRules(),IslandSettings(milp_time_limit_s=5),
                                    lambda *args:None,pad_order=many_order,
                                    baseline_candidate_ids=list(range(20)))
    assert len(many_selected)==34
    assert many_info['baseline_finite_lower_bound']==20
    assert _pad_occupancy(many_selected,PadSettings())['top']['internal_empty_slots']==0
    with TemporaryDirectory(prefix='verify_pad_') as tmp:
        source=Path(tmp)/'unknown_source.gds'
        output=Path(tmp)/'electrode_to_pad.gds'
        lib=gdstk.Library(unit=1e-6,precision=1e-9)
        cell=lib.new_cell('UNKNOWN_STRUCTURE')
        cell.add(gdstk.rectangle((-700,-20),(700,20),layer=42),
                 gdstk.rectangle((600,-150),(900,150),layer=42))
        lib.write_gds(str(source))
        rules=ProcessRules(minimum_center_spacing_um=5000)
        settings=IslandSettings(candidate_step_um=140,lane_half_steps=0,milp_time_limit_s=10)
        nav=build_navigation(source,rules,layer=42,outlet_mode='external_boundary',settings=settings)
        high_threshold=build_navigation(source,ProcessRules(
            minimum_center_spacing_um=5000,collector_clearance_um=80),
            layer=42,outlet_mode='external_boundary',settings=settings)
        assert nav.collector.is_empty and nav.summary['external_boundary_terminals']>0
        assert nav.summary['collector_interface_windows']==0
        assert nav.summary['open_tip_terminals']==0
        assert nav.navigation_domain.equals(high_threshold.navigation_domain)
        assert nav.summary['external_boundary_terminals']==high_threshold.summary['external_boundary_terminals']
        assert nav.summary['collector_clearance_threshold_um'] is None
        assert all(nav.graph.nodes[n]['terminal_kind']=='external_boundary' for n in nav.outlets)
        routing=solve_with_pads(nav,rules,settings=settings)
        assert routing['retained_routes']>=1,routing['candidate_library']
        layers=write_pad_gds(source,output,routing,(42,0))
        audit=audit_pad_gds(output,nav,routing,layers)
        assert audit['passed'] and audit['connected_pad_count']==routing['retained_routes']
        integer_audit=audit_exported_gds(
            source,output,layers,wire_spacing_um=rules.spacing_um,
            metal_support_margin_um=rules.margin_um,
            expected_nets=routing['retained_routes'],
            minimum_electrode_diameter_um=rules.electrode_diameter_um,
            minimum_electrode_center_distance_um=rules.minimum_center_spacing_um,
            minimum_substrate_disk_radius_um=(rules.electrode_diameter_um/2+
                                              rules.margin_um+settings.numeric_guard_um),
            maximum_substrate_disk_radius_um=(rules.electrode_diameter_um/2+
                                              rules.margin_um+settings.numeric_guard_um+.01),
            minimum_island_spacing_um=settings.pad_gap_um,
            minimum_pad_short_side_um=PadSettings().pad_width_um,
            minimum_pad_long_side_um=PadSettings().pad_length_um,
            minimum_wire_width_um=rules.wire_width_um)
        assert integer_audit['passed'] and integer_audit['original_support_exactly_preserved']
        assert integer_audit['support_provenance_verified']
        assert integer_audit['island_island_spacing_verified']
        assert integer_audit['island_envelopes_verified']
        assert integer_audit['functional_wire_width_verified']
        assert integer_audit['outer_exit_policy_verified']
        witness_file=Path(layers['wire_width_witness_path'])
        intact_witness=witness_file.read_bytes()
        altered=json.loads(intact_witness)
        altered['networks'][0]['points_grid_ticks'][0]=[40000000,40000000]
        witness_file.write_text(json.dumps(altered),encoding='utf-8')
        try:
            audit_exported_gds(source,output,layers,wire_spacing_um=4,
                               metal_support_margin_um=4,
                               expected_nets=routing['retained_routes'],
                               minimum_wire_width_um=5)
        except ValueError as exc:
            assert 'witness does not start' in str(exc),exc
        else:
            raise AssertionError('Invalid width witness passed')
        finally:
            witness_file.write_bytes(intact_witness)
        tampered=Path(tmp)/'unmarked_support.gds'
        copied=gdstk.read_gds(str(output),unit=1e-6)
        copied.top_level()[0].add(gdstk.rectangle((40000,40000),
                                                   (40010,40010),layer=42))
        copied.write_gds(str(tampered),max_points=4000)
        assert audit_support_polygon_provenance(source,output,layers)['passed']
        try:audit_support_polygon_provenance(source,tampered,layers)
        except ValueError as exc:assert 'unmarked or altered geometry' in str(exc),exc
        else:raise AssertionError('Exact polygon provenance accepted unmarked support')
        try:
            audit_exported_gds(source,tampered,layers,wire_spacing_um=4,
                               metal_support_margin_um=4,
                               expected_nets=routing['retained_routes'])
        except ValueError as exc:
            assert 'unmarked or altered geometry' in str(exc),exc
        else:
            raise AssertionError('Unmarked added support passed provenance audit')
        assert all(side['internal_empty_slots']==0 for side in audit['pad_bank_occupancy'].values())
        print(json.dumps({'routes':routing['retained_routes'],'pads':audit['pad_ids'],
                          'minimum_inter_net_gap_um':audit['minimum_inter_net_gap_um'],
                          'support_topology':audit['output_support_topology']},indent=2))


if __name__=='__main__':main()
