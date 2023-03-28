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

    ch = signal_channels(obj)
    erp_obj = erp(obj).data[ch, :]

    ch_n = size(erp_obj, 1)
    p = zeros(Int64, ch_n, 2)
    @inbounds @simd for ch_idx in 1:ch_n
        pp_pos = @views maximum(erp_obj[ch_idx, :])
        pp_neg = @views minimum(erp_obj[ch_idx, :])
        p[ch_idx, :] = @views [vsearch(pp_pos, erp_obj[ch_idx, :]), vsearch(pp_neg, erp_obj[ch_idx, :])]
    end

    return p
    
end
