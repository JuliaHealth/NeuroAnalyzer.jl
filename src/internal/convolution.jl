"""
Return the group delay (in samples) of a symmetric FIR kernel.
For an odd-length kernel of length K, the group delay is exactly `(K-1)/2`.
For an even-length kernel, the value is rounded to the nearest integer.
"""
_group_delay(kernel::AbstractVector)::Int64 = round(Int64, (length(kernel) - 1) / 2)

"""
Trim edge artifacts from a signal after kernel-based convolution.

Removes `group_delay` samples from each end of `s`, where `group_delay = round((length(kernel) - 1) / 2)`. This recovers the central, artifact-free portion of the convolved signal.

# Arguments

- `s::AbstractVector`: convolved signal
- `kernel::AbstractVector`: the kernel that was applied

# Returns

- `AbstractVector`: trimmed signal with both edges removed
"""
function _remove_kernel(s::AbstractVector, kernel::AbstractVector)::AbstractVector
    gd = _group_delay(kernel)
    # nothing to trim for a length-1 kernel (gd = 0)
    gd == 0 && return s
    # remove exactly gd samples from each end
    return s[(gd + 1):(end - gd)]
end
