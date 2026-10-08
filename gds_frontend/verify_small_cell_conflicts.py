"""Counterexamples and independent geometry checks for small-cell bounds."""
from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory
from random import Random

import gdstk
from shapely.geometry import MultiPoint, Polygon, box

from small_cell_conflicts import (_clip,_simple_integer_polygon,
                                  gds_small_cell_conflict_upper_bound,
                                  verify_small_cell_conflict_certificate)


def main():
    assert _simple_integer_polygon([(0,0),(10,0),(10,10),(0,10)])
    assert not _simple_integer_polygon([(0,0),(10,10),(0,10),(10,0)])
    assert not _simple_integer_polygon([(0,0),(10,0),(5,0),(5,10),(0,10)])
    rng=Random(2031)
    for _ in range(100):
        shape=MultiPoint([(rng.randrange(-20,21),rng.randrange(-20,21))
                          for _ in range(9)]).convex_hull
        if not isinstance(shape,Polygon):continue
        vertices=[tuple(map(int,p)) for p in shape.exterior.coords[:-1]]
        square=(-7,-5,8,9)
        clipped=_clip(vertices,square)
        actual=shape.intersection(box(*square))
        assert bool(clipped)==(not actual.is_empty)
        if not actual.is_empty:
            xs=[p[0] for p in clipped];ys=[p[1] for p in clipped]
            bounds=(float(min(xs)),float(min(ys)),
                    float(max(xs)),float(max(ys)))
            assert all(abs(a-b)<1e-9 for a,b in zip(bounds,actual.bounds))
    with TemporaryDirectory(prefix='small_cells_') as folder:
        source=Path(folder)/'unknown.gds'
        lib=gdstk.Library(unit=1e-6,precision=1e-9)
        lib.new_cell('UNKNOWN').add(gdstk.rectangle((-700,-20),(700,20),
                                                    layer=42))
        lib.write_gds(str(source))
        cert=gds_small_cell_conflict_upper_bound(source,42,0,5000)
        assert cert['value']==1
        assert verify_small_cell_conflict_certificate(source,cert)
        modified={**cert,'value':0}
        assert not verify_small_cell_conflict_certificate(source,modified)
        modified={**cert,'conflict_edges':[[0,0]]}
        assert not verify_small_cell_conflict_certificate(source,modified)
    print('exact rational cell conflict tests passed')


if __name__=='__main__':main()
