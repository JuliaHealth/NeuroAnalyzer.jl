export mep_peaks

"""
    mep_peaks(obj)

Detect a pair of positive and negative peaks of MEP.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns
 
- `p::Array{Int64, 2}`: peaks: channels × positive peak position, negative peak position
"""
function mep_peaks(obj::NeuroAnalyzer.NEURO)

    _check_datatype(obj, "mep")

    ch_n = size(obj)[1]
    p = zeros(Int64, ch_n, 2)
    @inbounds @simd for ch_idx in 1:ch_n
        pp_pos = @views maximum(obj.data[ch_idx, obj.header.recording[:stimulation_sample][ch_idx] + 20:end, 1])
        pp_neg = @views minimum(obj.data[ch_idx, obj.header.recording[:stimulation_sample][ch_idx] + 20:end, 1])
        p[ch_idx, :] = @views [vsearch(pp_pos, obj.data[ch_idx, obj.header.recording[:stimulation_sample][ch_idx] + 20:end, 1]), vsearch(pp_neg, obj.data[ch_idx, obj.header.recording[:stimulation_sample][ch_idx] + 20:end, 1])]
        p[ch_idx, :] .+= (obj.header.recording[:stimulation_sample][ch_idx] + 20)
    end

    return p
    
end
