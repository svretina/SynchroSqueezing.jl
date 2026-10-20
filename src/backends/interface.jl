"""
    plan_cwt(n; kwargs...)

Create a reusable, backend-neutral CWT plan. Backends implement this function
for their own configuration types.
"""
function plan_cwt end

"""
    cwt(plan, signal)

Compute a CWT and its physical-time derivative. Results always use
scale/frequency × time orientation, independent of backend storage.
"""
function cwt end

"""Return a full-signal inverse CWT using the result's recorded metadata."""
function inverse_cwt end

"""Return a stable phase-based instantaneous-frequency estimate in Hz."""
function instantaneous_frequency end
