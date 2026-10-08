"""Meaningful geometry checks for source independence and topology stability."""
from pathlib import Path
from tempfile import TemporaryDirectory
import json
import sys
import argparse

import gdstk

from frontend import extract
from workspace_paths import DATA_DIR, PACKAGE_DIR

HERE=Path(__file__).resolve().parent
INPUT=DATA_DIR


def topology(result):
    s=result.summary
    return {'support_components':s['support_components'],
            'support_holes':s['support_holes'],
            'full_graph_cycle_rank':s['full_graph_cycle_rank'],
            'collector_interfaces':s['collector_interface_windows'],
            'graph_nodes':s['graph_nodes'],'graph_edges':s['graph_edges'],
            'selected_pitch_um':s['selected_pitch_um']}


def main():
    global INPUT
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir',type=Path,default=DATA_DIR)
    parser.add_argument('--output',type=Path,default=PACKAGE_DIR/'outputs'/'verification'/'frontend.json')
    args=parser.parse_args()
    INPUT=args.input_dir.resolve()
    files=sorted(INPUT.glob('*.gds'))
    if not files:
        parser.error('input directory contains no GDS files')
    results={}
    denied_source_reads=[]
    def forbid_input_source(event,args):
        if event!='open' or not isinstance(args[0],(str,bytes)):
            return
        path=str(args[0]).lower()
        if str(INPUT).lower() in path and path.endswith(('.py','.pyc')):
            denied_source_reads.append(path)
            raise RuntimeError('Input generator source access is forbidden in this verification')
    sys.addaudithook(forbid_input_source)
    with TemporaryDirectory(prefix='blind_gds_verification_') as temp:
        temp=Path(temp)
        for i,path in enumerate(files):
            print('VERIFY',path.name,flush=True)
            baseline=extract(path)
            source=temp/'read.gds';source.write_bytes(path.read_bytes())
            lib=gdstk.read_gds(str(source))
            for j,cell in enumerate(lib.cells):
                cell.name=f'unknown_cell_{j}'
            anonymous=temp/f'unknown_input_{i}.gds'
            lib.write_gds(str(anonymous))
            blind=extract(anonymous)
            a=topology(baseline);b=topology(blind)
            if a!=b:
                raise RuntimeError(f'Renaming changed extraction: {a} vs {b}')
            # Refine the accepted grid once; this checks cycles and interfaces,
            # not equality of every raster-generated branch/spur.
            refined=extract(anonymous,first_pitch=baseline.pitch/2,
                            minimum_pitch=baseline.pitch/2)
            c=topology(refined)
            stable_keys=['support_components','support_holes','full_graph_cycle_rank','collector_interfaces']
            if any(a[k]!=c[k] for k in stable_keys):
                raise RuntimeError(f'Refinement changed structural invariants: {a} vs {c}')
            # Flatten and translate/rotate while using a fresh unrelated cell
            # name. No type/shape family is passed to the extractor.
            lib=gdstk.read_gds(str(source),unit=1e-6)
            transformed_lib=gdstk.Library(unit=1e-6,precision=1e-9)
            cell=transformed_lib.new_cell('arbitrary_rotated_geometry')
            for top in lib.top_level():
                for polygon in top.get_polygons():
                    polygon.rotate(1.5707963267948966)
                    polygon.translate(100,200)
                    cell.add(polygon)
            transformed=temp/f'transformed_{i}.gds'
            transformed_lib.write_gds(str(transformed))
            rotated=extract(transformed)
            d=topology(rotated)
            if any(a[k]!=d[k] for k in stable_keys):
                raise RuntimeError(f'Rigid transform changed invariants: {a} vs {d}')
            results[path.name]={'rename_and_source_denial_passed':True,
                                'baseline':a,'refined':c,'rotated_translated':d,
                                'refinement_invariants_passed':True,
                                'rigid_transform_invariants_passed':True,
                                'raw_node_counts_equal_on_refinement':a['graph_nodes']==c['graph_nodes']}
        # A loop without a branch must survive compression as a self-loop.
        lib=gdstk.Library(unit=1e-6,precision=1e-9)
        cell=lib.new_cell('isolated_loop')
        outer=gdstk.rectangle((0,0),(100,100))
        inner=gdstk.rectangle((10,10),(90,90))
        for polygon in gdstk.boolean([outer],[inner],'not',layer=10):
            cell.add(polygon)
        path=temp/'loop.gds';lib.write_gds(str(path))
        loop=extract(path,collector_clearance_um=1000)
        if loop.summary['full_graph_cycle_rank']!=1:
            raise RuntimeError('A branch-free loop was lost')
        # GDS native units are arbitrary; input in mm should give the same um geometry.
        mm=gdstk.Library(unit=1e-3,precision=1e-9)
        mmcell=mm.new_cell('millimetre_units')
        for p in cell.get_polygons():
            p.scale(.001);mmcell.add(p)
        mmfile=temp/'mm.gds';mm.write_gds(str(mmfile))
        mmresult=extract(mmfile,collector_clearance_um=1000)
        if topology(mmresult)!=topology(loop):
            raise RuntimeError('GDS unit conversion changed extraction')
        results['synthetic_closed_loop']={'passed':True,'topology':topology(loop),
                                          'native_mm_unit_conversion_passed':True}
        # Two almost-touching components merge on all allowed rasters. The
        # frontend must refuse the graph rather than manufacture a connection.
        lib=gdstk.Library(unit=1e-6,precision=1e-9)
        cell=lib.new_cell('unresolved_separation')
        cell.add(gdstk.rectangle((0,0),(20,20),layer=10))
        cell.add(gdstk.rectangle((20.1,0),(40.1,20),layer=10))
        path=temp/'too_close.gds';lib.write_gds(str(path))
        try:
            extract(path)
        except RuntimeError as exc:
            if 'Raster topology did not match' not in str(exc):
                raise
            results['unresolved_0_1_um_gap']={'refused_as_expected':True,'reason':str(exc)}
        else:
            raise RuntimeError('Unresolved subpixel separation was silently accepted')
    if denied_source_reads:
        raise RuntimeError('The algorithm attempted to read an input source file')
    payload={'all_required_checks_passed':True,'input_source_reads':0,'checks':results,
             'scope':'Invariant checks do not certify every raster branch, planar turn, or arbitrary GDS.'}
    out=args.output
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(payload,ensure_ascii=False,indent=2),encoding='utf-8')
    print('VERIFIED',out,flush=True)


if __name__=='__main__':
    main()
