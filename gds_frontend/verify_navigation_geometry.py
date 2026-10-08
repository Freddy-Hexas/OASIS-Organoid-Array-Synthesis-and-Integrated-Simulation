"""Counterexamples for component certificates and topology-aware navigation."""
from __future__ import annotations

from shapely.geometry import box
from shapely.ops import unary_union

from island_router import (IslandSettings,NavigationGeometryError,
                           _navigation_domain,betti)
from process_geometry import ProcessRules


def main():
    rules=ProcessRules()
    settings=IslandSettings()
    # A 0.04 um neck is useful: its large connected component admits a
    # distant terminal, so the default 0.05 um target must shrink locally.
    neck=unary_union([box(-100,-20,0,20),box(0,-.02,100,.02),
                      box(100,-20,200,20)])
    # The 0.3 um island is genuinely too small for the model's required
    # anchor-to-terminal separation, independent of graph sampling.
    tiny=box(300,0,300.3,.3)
    center=unary_union([neck,tiny])
    nav,certificate=_navigation_domain(center,rules,settings,.001)
    assert certificate['excluded_component_count']==1
    excluded=certificate['excluded_components'][0]
    assert excluded['bbox_diameter_upper_um']+excluded['precision_allowance_um'] < excluded['required_center_to_terminal_distance_um']
    assert certificate['retained_center_topology']==certificate['navigation_topology']=={'components':1,'holes':0}
    assert betti(nav)=={'components':1,'holes':0}
    assert 0 < certificate['component_insets'][0]['selected_inset_um'] < settings.numeric_guard_um

    # Small area alone is not a certificate: a very long, narrow component
    # can have widely separated anchors and terminals and must be retained.
    thin=box(0,0,100,.01)
    _,certificate=_navigation_domain(thin,rules,settings,.001)
    assert certificate['excluded_component_count']==0
    assert certificate['retained_center_topology']=={'components':1,'holes':0}

    try:
        _navigation_domain(tiny,rules,settings,.001)
    except NavigationGeometryError as exc:
        assert exc.diagnostics['excluded_component_count']==1
        assert exc.diagnostics['excluded_components'][0]['reason']=='no_anchor_terminal_pair_can_meet_clearance'
    else:
        raise AssertionError('An impossible component was accepted')
    # The same tiny component must not be called physically unusable when
    # short island-to-exterior connections are allowed by the search model.
    unguarded=IslandSettings(navigation_component_exclusion_distance_um=0.0,
                             navigation_prune_leaf_length_um=0.0)
    nav,certificate=_navigation_domain(tiny,rules,unguarded,.001)
    assert certificate['excluded_component_count']==0
    assert betti(nav)=={'components':1,'holes':0}
    print('navigation geometry certificates passed',flush=True)


if __name__=='__main__':main()
