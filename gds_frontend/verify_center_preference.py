"""Regression checks for center-aware Pad candidates and lexicographic choice."""
from random import Random

import numpy as np
from shapely.geometry import box

from island_router import IslandSettings, _physical_conflicts, _select
from pad_router import _PadColumnReservoir
from process_geometry import ProcessRules


def column(index, x, pad_id, length, depth):
    return {'source_node':index, 'source_um':(float(x),0.0),
            'complete_length_um':float(length), 'depth_um':float(depth),
            'metal':box(x-5,-5,x+5,5), 'island':box(x-10,-10,x+10,10),
            'pad_id':pad_id, 'pad_side':'top',
            'pad_index':int(pad_id.split('-')[1])}


def main():
    random=Random(20261003)
    shapes=[]
    for _ in range(100):
        x=random.randrange(40);y=random.randrange(40)
        shapes.append(box(x,y,x+random.randrange(1,6),y+random.randrange(1,6)))
    source_codes=np.asarray([i//3 for i in range(100)],dtype=np.int32)
    pad_codes=np.asarray([i//7 for i in range(100)],dtype=np.int32)
    expected={(j,i) for i in range(100) for j in range(i)
              if source_codes[i]!=source_codes[j] and pad_codes[i]!=pad_codes[j]
              and shapes[i].distance(shapes[j])<=4-1e-6}
    assert _physical_conflicts(shapes,4,source_codes,pad_codes,16)==expected

    for size in (0,1,50,1000):
        reservoir=_PadColumnReservoir(4,0,4,(0,0))
        choices=[]
        for index in range(size):
            source=random.randrange(max(1,size//3))
            candidate={'id':index,'source_node':source,
                       'source_um':(float(source*17),float(source*3)),
                       'complete_length_um':float(random.randrange(1,1000)),
                       'depth_um':float(random.randrange(1,1000))}
            choices.append(candidate)
            reservoir.add(candidate)
        expected=[]
        seen=set()
        for candidate in sorted(choices,key=lambda c:(c['complete_length_um'],
                                                     -c['depth_um'],c['id'])):
            if len(expected)==4:break
            if candidate['source_node'] not in seen:
                expected.append(candidate);seen.add(candidate['source_node'])
        for candidate in sorted(choices,key=lambda c:(c['source_um'][0]**2+
                                                     c['source_um'][1]**2,
                                                     c['complete_length_um'],c['id'])):
            if len(expected)==8:break
            if candidate['source_node'] not in seen:
                expected.append(candidate);seen.add(candidate['source_node'])
        assert [c['id'] for c in reservoir.selected()]==[c['id'] for c in expected]

    reservoir=_PadColumnReservoir(4,0,4,(0,0))
    for index in range(12):
        reservoir.add(column(index,index*100,f'top-{index+1:02d}',
                             100-index,100-index))
    retained=reservoir.selected()
    assert {c['source_node'] for c in retained[:4]}=={8,9,10,11}
    assert {c['source_node'] for c in retained[4:]}=={0,1,2,3}

    routes=[column(1,100,'top-01',200,1),
            column(2,1000,'top-01',10,100),
            column(3,200,'top-02',200,1),
            column(4,2000,'top-02',10,100)]
    order={'top':['top-01','top-02']}
    chosen,report=_select(routes,ProcessRules(),
                          IslandSettings(milp_time_limit_s=10),
                          lambda *args:None,pad_order=order,
                          center_origin=(0,0))
    assert report['finite_model_proven_optimal']
    assert report['secondary_proven_optimal']
    assert {c['source_node'] for c in chosen}=={1,3}

    routes.append(column(5,3000,'top-03',200,1))
    order={'top':['top-01','top-02','top-03']}
    chosen,report=_select(routes,ProcessRules(),
                          IslandSettings(milp_time_limit_s=10),
                          lambda *args:None,pad_order=order,
                          center_origin=(0,0))
    assert report['finite_lower_bound']==3
    assert {c['source_node'] for c in chosen}=={1,3,5}
    print('center candidates retained; count remains first; center cost breaks ties')


if __name__=='__main__':
    main()
