"""No default wall-clock limit for layout construction or optimization.

None means no deadline and no solver time_limit option. Explicit finite
limits remain compatible with offline experiments and old test fixtures;
the web workbench creates its settings with no limit.
"""
from __future__ import annotations

import math
import time


def validate_time_limit(seconds):
    if seconds is not None and (isinstance(seconds, bool) or
            not isinstance(seconds, (int, float)) or
            not math.isfinite(seconds) or seconds <= 0):
        raise ValueError('Explicit time limit must be finite and positive, or None for no limit')
    return seconds


def deadline_after(seconds, *, started=None):
    validate_time_limit(seconds)
    if seconds is None:
        return None
    return (time.perf_counter() if started is None else started) + seconds


def deadline_expired(deadline):
    return deadline is not None and time.perf_counter() >= deadline


def solver_options(seconds=None, **options):
    validate_time_limit(seconds)
    # Omitting the option invokes the solver's unlimited wall-clock policy.
    # Passing None, infinity or a large placeholder number would not be the
    # same API contract across LP and MILP solvers.
    if seconds is not None:
        options['time_limit'] = float(seconds)
    return options
