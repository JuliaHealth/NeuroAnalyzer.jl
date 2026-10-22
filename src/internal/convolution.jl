"""
Return the group delay (in samples) of a symmetric FIR kernel.
For an odd-length kernel of length K, the group delay is exactly `(K-1)/2`.
For an even-length kernel, the value is rounded to the nearest integer.
"""
_group_delay(kernel::AbstractVector)::Int64 = round(Int64, (length(kernel) - 1) / 2)

"""
Trim the full convolution of a signal with a kernel to the length of the signal ("same" convolution).

The full convolution has `length(signal) + length(kernel) - 1` samples. The first `(length(kernel) - 1) ÷ 2` samples (the group delay of a symmetric kernel) are removed and the next `length(signal)` samples are kept, so that the output is aligned with the input (no delay for a symmetric kernel). The edges of the output are still affected by the kernel overlapping the ends of the signal (edge effects).

# Arguments

- `s::AbstractVector`: full convolution of a signal with `kernel`
- `kernel::AbstractVector`: the kernel that was applied

# Returns

- `AbstractVector`: trimmed signal, `length(s) - length(kernel) + 1` samples
"""
function _remove_kernel(s::AbstractVector, kernel::AbstractVector)::AbstractVector
    n = length(s) - length(kernel) + 1
    n >= 1 || throw(ArgumentError("s must be the full convolution with kernel."))
    gd = (length(kernel) - 1) ÷ 2
    return s[(gd + 1):(gd + n)]
end
