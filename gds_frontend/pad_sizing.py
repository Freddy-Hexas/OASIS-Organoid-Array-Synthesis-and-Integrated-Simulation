"""Uniform Pad sizing coupled to the most populated continuous side bank."""
from dataclasses import replace
import math


def validate_pad_sizing(settings,*,wire_width_um=5,margin_um=4):
    if settings.sizing_mode not in ('expand_frame','compact_then_expand'):
        raise ValueError('Unknown Pad sizing mode')
    minimum_width=settings.minimum_pad_width_um if settings.minimum_pad_width_um is not None else settings.pad_width_um
    minimum_length=settings.minimum_pad_length_um if settings.minimum_pad_length_um is not None else settings.pad_length_um
    contact=wire_width_um+2*margin_um
    if (not all(math.isfinite(value) for value in (minimum_width,minimum_length)) or
            not contact<=minimum_width<=settings.pad_width_um or
            not contact<=minimum_length<=settings.pad_length_um):
        raise ValueError('Pad lower dimensions must fit the wire and margins and not exceed preferred dimensions')
    return minimum_width,minimum_length


def size_pad_bank(settings,opening_um,required_slots,current_side_um,*,spacing_um,margin_um,grid_um,
                  wire_width_um=5.0,clearance_guard_um=0.0):
    """Maximize a common dimension factor in the current footprint, then grow.

    This optimizes a declared uniform-sizing policy, not arbitrary rectangle
    arrangements. Bank occupancy and full metal geometry are still audited.
    """
    minimum_width,minimum_length=validate_pad_sizing(settings,wire_width_um=wire_width_um,margin_um=margin_um)
    if isinstance(required_slots,bool) or not isinstance(required_slots,int) or required_slots<1:
        raise ValueError('Requested Pad slots must be a positive integer')
    floor=max(minimum_width/settings.pad_width_um,minimum_length/settings.pad_length_um)
    def geometry(scale):
        # Even grid ticks keep both rectangle edges on the native grid when
        # the square/bank is centered; odd half-grid rounding can shrink a Pad.
        def up(value):return math.ceil(value/(2*grid_um)-1e-10)*(2*grid_um)
        width=up(max(minimum_width,settings.pad_width_um*scale))
        length=up(max(minimum_length,settings.pad_length_um*scale))
        pitch=up(max(settings.pad_pitch_um*scale,width+spacing_um+clearance_guard_um))
        corner=settings.pad_outer_setback_um+length+spacing_um
        side=max(settings.square_side_um,2*(opening_um+length+max(1000.,margin_um)),
                 2*corner+width+(required_slots-1)*pitch)
        return width,length,pitch,up(side)
    scale=1.
    if settings.sizing_mode=='compact_then_expand' and geometry(1.)[3]>current_side_um+1e-8:
        if geometry(floor)[3]>current_side_um+1e-8:scale=floor
        else:
            lo,hi=floor,1.
            for _ in range(52):
                mid=(lo+hi)/2
                if geometry(mid)[3]<=current_side_um+1e-8:lo=mid
                else:hi=mid
            scale=lo
    width,length,pitch,side=geometry(scale)
    actual=replace(settings,pad_width_um=width,pad_length_um=length,pad_pitch_um=pitch)
    return actual,max(current_side_um,side),{
        'policy':settings.sizing_mode,'required_slots_per_side':required_slots,
        'preferred_dimensions_um':[settings.pad_width_um,settings.pad_length_um],
        'minimum_dimensions_um':[minimum_width,minimum_length],
        'actual_dimensions_um':[width,length],'actual_pitch_um':pitch,
        'common_scale_factor':scale,'pad_dimensions_reduced':width<settings.pad_width_um or length<settings.pad_length_um,
        'previous_side_um':current_side_um,'actual_side_um':max(current_side_um,side),
        'frame_expanded':side>current_side_um+1e-8,
        'scope':'uniform dimension factor; prefer current footprint then grow square; full fanout certification remains mandatory'}


def pad_settings_for_frame(settings,frame):
    actual=frame.get('effective_pad_dimensions')
    return replace(settings,**actual) if actual else settings
