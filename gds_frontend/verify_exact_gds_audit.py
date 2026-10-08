"""Integer-grid counterexamples for the independent export checker."""
from pathlib import Path
from tempfile import TemporaryDirectory

import gdstk

from exact_gds_audit import (_center_in_original_paths, _certified_marker_center,
                             _certified_pad_bbox, _read_layer_paths, _union,
                             _HullGate, _audit_bridge_outer_window,
                             audit_exported_gds)


LAYERS={'support_layer':[10,0],'metal_layer':20,
        'electrode_marker_layer':30,'pad_contact_layer':31,
        'bridge_marker_layer':32}


def _fixture(folder,second_y_um=19.0,first_y_um=10.0):
    source=folder/'source.gds'
    output=folder/'output.gds'
    lib=gdstk.Library(unit=1e-6,precision=1e-9)
    cell=lib.new_cell('UNFAMILIAR')
    cell.add(gdstk.rectangle((0,0),(100,50),layer=10))
    lib.write_gds(str(source),max_points=4000)
    for network,low in ((1,first_y_um),(2,second_y_um)):
        cell.add(gdstk.rectangle((10,low),(90,low+5),layer=20,datatype=network),
                 gdstk.rectangle((10,low),(20,low+5),layer=30,datatype=network),
                 gdstk.rectangle((80,low),(90,low+5),layer=31,datatype=network),
                 gdstk.rectangle((70,low),(80,low+5),layer=32,datatype=network))
    lib.write_gds(str(output),max_points=4000)
    return source,output


def main():
    # An acute hull's eroded core need not contain the bounding-circle center.
    # A bridge enclosing the whole core must still fail the independent audit.
    gate=_HullGate(((0,0),(1000,0),(200,100)),20)
    assert not gate.segment_enters((500,0),(500,0))
    assert gate.core_reference is not None and gate.segment_enters(gate.core_reference,gate.core_reference)
    enclosing=_union([[(-100,-100),(1100,-100),(1100,200),(-100,200)]])
    try:
        _audit_bridge_outer_window(enclosing,enclosing,(500,0),gate)
    except ValueError as exc:
        assert 'forbidden interior core' in str(exc),exc
    else:
        raise AssertionError('A bridge containing an off-center hull core passed')
    original=[[(0,0),(100000,0),(100000,100000),(0,100000)]]
    assert _center_in_original_paths((100000,100000),original)
    assert not _center_in_original_paths((230000,100000),original)
    rectangle=_union([[(0,0),(500000,0),(500000,3000000),(0,3000000)]])
    diagonal=_union([[(0,0),(500000,0),(500000,3000000),(400000,3000000)]])
    assert _certified_pad_bbox(rectangle,500000,3000000)
    assert not _certified_pad_bbox(diagonal,500000,3000000)
    assert _certified_marker_center(_union([[(-15000,-15000),(15000,-15000),
                                             (15000,15000),(-15000,15000)]]),
                                    30000)==(0,0)
    try:
        _certified_marker_center(_union([[(-14999,-15000),(14999,-15000),
                                         (14999,15000),(-14999,15000)]]),30000)
    except ValueError as exc:
        assert 'inscribed diameter' in str(exc)
    else:
        raise AssertionError('29.998 um inscribed diameter passed a 30 um rule')
    from island_router import disk
    with TemporaryDirectory(prefix='verify_disk_quantization_') as temporary:
        output=Path(temporary)/'disk.gds'
        lib=gdstk.Library(unit=1e-6,precision=1e-9)
        cell=lib.new_cell('DISK')
        for index,center in enumerate(((0,0),(0.137,0.281),(37.333,-19.779)),1):
            shape=disk(center,15.0)
            cell.add(gdstk.Polygon(shape.exterior.coords[:-1],layer=30,
                                   datatype=index))
        lib.write_gds(str(output),max_points=4000)
        paths,_,tick=_read_layer_paths(output,.001)
        for index in range(1,4):
            assert _certified_marker_center(_union(paths[(30,index)]),30000)
    with TemporaryDirectory(prefix='verify_exact_gds_') as temporary:
        folder=Path(temporary)
        source,output=_fixture(folder)
        result=audit_exported_gds(source,output,LAYERS,
                                  wire_spacing_um=4,metal_support_margin_um=4,
                                  expected_nets=2)
        assert result['passed'] and result['original_support_exactly_preserved']
        assert result['network_pairs_checked']==1
        source,output=_fixture(folder,second_y_um=5.001,first_y_um=0)
        result=audit_exported_gds(source,output,LAYERS,
                                  wire_spacing_um=0,metal_support_margin_um=0,
                                  expected_nets=2)
        assert result['passed'] and result['minimum_inter_net_spacing_um']==0
        assert result['minimum_metal_support_margin_um']==0
        # Zero means no positive distance rule, never permission for a short.
        for second_y in (5.0,4.999):
            source,output=_fixture(folder,second_y_um=second_y,first_y_um=0)
            try:
                audit_exported_gds(source,output,LAYERS,
                                   wire_spacing_um=0,metal_support_margin_um=0,
                                   expected_nets=2)
            except ValueError as exc:
                assert 'violate spacing' in str(exc),exc
            else:
                raise AssertionError('Touching or overlapping metal passed zero-gap audit')
        source,output=_fixture(folder,second_y_um=18.999)
        try:
            audit_exported_gds(source,output,LAYERS,
                               wire_spacing_um=4,metal_support_margin_um=4,
                               expected_nets=2)
        except ValueError as exc:
            assert 'violate spacing' in str(exc),exc
        else:
            raise AssertionError('A 3.999 um inter-net gap passed')
        source,output=_fixture(folder,first_y_um=3.999)
        try:
            audit_exported_gds(source,output,LAYERS,
                               wire_spacing_um=4,metal_support_margin_um=4,
                               expected_nets=2)
        except ValueError as exc:
            assert 'support margin' in str(exc),exc
        else:
            raise AssertionError('A 3.999 um metal-support margin passed')
    print('integer-grid audit: default and zero-gap rules pass; undersized gaps, margins and shorts are rejected')


if __name__=='__main__':main()
