"""Geometric regression: a second track and a node-like detour."""
from shapely.geometry import LineString, box

from residual_vector_router import ResidualVectorRouter


def main():
    support=box(0,0,300,100)
    existing=LineString([(10,40),(290,40)]).buffer(2.5,cap_style='round')
    residual=ResidualVectorRouter(support,[existing])
    route=residual.route((15,65),[(285,65)])
    assert route is not None
    line,target,audit=route
    assert target==0 and audit['support_containment_verified']
    assert audit['minimum_existing_metal_gap_um']>=4
    assert line.length>=270
    second=residual.route((15,85),[(285,65)])
    assert second is not None and second[1]==0
    assert len(residual._target_cache)==1
    assert residual.route((15,40),[(285,40)]) is None

    obstacle=LineString([(150,10),(150,80)]).buffer(2.5,cap_style='round')
    detour=ResidualVectorRouter(support,[obstacle])
    turned=detour.route((40,50),[(260,50)])
    assert turned is not None
    bend=turned[0]
    assert 220<bend.length<350 and max(y for x,y in bend.coords)>85
    assert turned[2]['minimum_existing_metal_gap_um']>=4

    narrow=box(0,0,300,20)
    blocked=LineString([(0,10),(300,10)]).buffer(2.5)
    assert ResidualVectorRouter(narrow,[blocked]).route(
        (10,10),[(290,10)]) is None
    print('residual vector routing: parallel track, certified bend, blocked strip passed')


if __name__=='__main__':
    main()
