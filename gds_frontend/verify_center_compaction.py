"""Residual center placement: two tracks, blocked passage, rigid transforms."""
from types import SimpleNamespace
from unittest.mock import patch
import math

import numpy as np
from shapely.affinity import rotate, translate
from shapely.geometry import GeometryCollection, LineString, Point, box

from center_compaction import CompactionSettings, compact_centers, _try_replacement
from island_router import IslandSettings, _select, disk
from process_geometry import ProcessRules
from navigation_simplification import SearchDeadlineExceeded
import island_router
import center_compaction


def simple_net(x,y,index):
    source=(x,y);target=(485,y);terminal=(490,y);pad_target=(650,y)
    central=LineString([source,target]).buffer(2.52,quad_segs=16)
    join=LineString([target,terminal]).buffer(2.51,quad_segs=16)
    outer=LineString([terminal,pad_target]);outer_wire=outer.buffer(2.51,quad_segs=32)
    electrode=disk(source,15);island=disk(source,19.05)
    pad=box(630,y-15,670,y+15)
    return {'source_node':index,'source_um':source,'points_um':[source,target],
        'outlet_um':target,'declared_terminal_um':terminal,
        'pad_id':f'right-{index:02d}','pad_side':'right','pad_index':index,
        'pad_target_um':pad_target,'portal_escape_um':(500,y),
        'bridge':outer.buffer(20),'pad':pad,'island':island,'electrode':electrode,
        'central_wire':central,'wire':central,'join_wire':join,
        'outer_wire':outer_wire,'outer_line':outer,'outer_length_um':outer.length,
        'metal':GeometryCollection([central.union(electrode),join,outer_wire,pad])}


def transformed(net,angle,offset):
    radians=math.radians(angle)
    rotation=np.asarray([[math.cos(radians),-math.sin(radians)],
                         [math.sin(radians),math.cos(radians)]])
    result=dict(net)
    for name in ['source_um','outlet_um','declared_terminal_um','pad_target_um','portal_escape_um']:
        result[name]=(rotation@np.asarray(net[name])+offset).tolist()
    result['points_um']=(np.asarray(net['points_um'])@rotation.T+offset).tolist()
    for name in ['bridge','pad','island','electrode','central_wire','wire','join_wire','outer_wire','outer_line','metal']:
        result[name]=translate(rotate(net[name],angle,origin=(0,0)),*offset)
    return result


def main():
    rules=ProcessRules();assert rules.minimum_center_spacing_um==70
    settings=IslandSettings();search=CompactionSettings(max_sweeps=2,time_limit_s=30)
    base=box(-500,-60,500,60)
    fixture=[simple_net(300,35,1),simple_net(300,-35,2)]
    for angle,offset in [(0,(0,0)),(37,(172,-291))]:
        support=translate(rotate(base,angle,origin=(0,0)),*offset)
        chosen=[transformed(c,angle,offset) for c in fixture]
        result,certificate=compact_centers(support,chosen,rules,settings,offset,search=search)
        assert len(result)==2
        assert certificate['after']['maximum_radius_um']<80,certificate
        assert certificate['after']['sum_squared_radius_um2']<certificate['before']['sum_squared_radius_um2']
        assert math.dist(result[0]['source_um'],result[1]['source_um'])>=70
        for old,new in zip(chosen,result):
            for key in ['bridge','pad','outer_line','outer_wire','join_wire']:
                assert new[key].equals_exact(old[key],0)
            assert new['declared_terminal_um']==old['declared_terminal_um']
            assert support.covers(new['central_wire'].buffer(4,quad_segs=32))
        assert result[0]['metal'].distance(result[1]['metal'])>=4
        assert not certificate['global_center_optimality_proven']

    # A full transverse metal obstacle separates the original structure.
    # The source must remain in the exit's component despite a closer region.
    old=simple_net(300,0,1)
    wall={'source_um':(0,0),'island':disk((0,0),19.05),
          'metal':LineString([(0,-60),(0,60)]).buffer(2.5)}
    moved,detail=_try_replacement(base,old,[wall],rules,settings,np.asarray((0,0)),.001,32)
    assert moved is not None,detail
    assert moved['source_um'][0]>0 and 70<=math.dist(moved['source_um'],(0,0))<100
    assert min(x for x,y in moved['points_um'])>0
    # Disconnected source islands cannot be bridged by compaction.
    disconnected=box(-500,-60,-100,60).union(box(100,-60,500,60))
    moved,detail=_try_replacement(disconnected,old,[],rules,settings,np.asarray((0,0)),.001,32)
    assert moved is not None and moved['source_um'][0]>=106
    assert disconnected.covers(moved['central_wire'].buffer(4,quad_segs=32))

    # Exhausting the search after a fully checked move must retain that move.
    # This models a computational deadline, not a failed geometric condition.
    real_next=center_compaction.AttachmentAwarePlacementSearch.next_proposal
    calls=0
    def timeout_after_one(self):
        nonlocal calls
        calls+=1
        if calls>1:raise SearchDeadlineExceeded('fixture after certified move')
        return real_next(self)
    with patch.object(center_compaction.AttachmentAwarePlacementSearch,'next_proposal',timeout_after_one):
        moved,detail=_try_replacement(base,old,[],rules,settings,np.asarray((0,0)),.001,32)
    assert moved is not None and math.dist(moved['source_um'],(0,0))<10
    assert detail['search_deadline_hit'] and detail['placement_search']['deadline_hit']

    # A time-limited secondary MILP incumbent must not worsen centrality.
    columns=[]
    for i,x in enumerate((100,1000,200,2000)):
        columns.append({'source_node':i,'source_um':(x,0),'depth_um':1,
                        'metal':box(x-2,-2,x+2,2),'island':box(x-5,-5,x+5,5),
                        'pad_id':f'top-{i//2+1:02d}','pad_side':'top','pad_index':i//2+1})
    optimum=SimpleNamespace(x=np.asarray([1,0,1,0,1,0]),status=0,message='fixture',mip_dual_bound=-2)
    worse=SimpleNamespace(x=np.asarray([0,1,0,1,1,0]),status=1,message='time limit')
    with patch.object(island_router,'milp',side_effect=[optimum,worse]):
        selected,report=_select(columns,rules,settings,lambda *args:None,
            pad_order={'top':['top-01','top-02']},center_origin=(0,0))
    assert [c['source_node'] for c in selected]==[0,2]
    assert not report['secondary_incumbent_accepted']
    print('PASS: two residual tracks, rigid transforms, 70 um spacing, blocked components, fixed exits, monotone timeout handling')


if __name__=='__main__':main()
