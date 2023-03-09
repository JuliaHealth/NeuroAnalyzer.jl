export erp_peaks

"""
    erp_peaks(obj)

Detect a pair of positive and negative peaks of ERP.

# Arguments

- `obj::NeuroAnalyzer.NEURO`:

# Returns
 
- `p::Array{Int64, 2}`: peaks: channels × positive peak position, negative peak position
"""
function erp_peaks(obj::NeuroAnalyzer.NEURO)

    channels = signal_channels(obj)
    erp = erp(obj).signals[channels, :]

    ch_n = size(erp, 1)
    p = zeros(Int64, ch_n, 2)
    @inbounds @simd for ch_idx in 1:ch_n
        pp_pos = @views maximum(erp[ch_idx, :])
        pp_neg = @views minimum(erp[ch_idx, :])
        p[ch_idx, :] = @views [vsearch(pp_pos, erp[ch_idx, :]), vsearch(pp_neg, erp[ch_idx, :])]
    end

    return p
end
