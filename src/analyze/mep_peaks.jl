export mep_peaks

"""
    mep_peaks(obj)

Detect a pair of positive and negative peaks of MEP.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `p::Matrix{Int64}`: peaks, shape (channels × positive peak position, negative peak position)
"""
function mep_peaks(obj::NeuroAnalyzer.NEURO)::Matrix{Int64}

    _check_datatype(obj, "mep")

    # number of channels
    ch_n = nchannels(obj)

    # pre-allocate output
    p = zeros(Int64, ch_n, 2)

    # calculate over channels
    ss = obj.header.recording[:stimulation_sample] .+ 20
    @inbounds Threads.@threads :dynamic for ch_idx in 1:ch_n

        p[ch_idx, :] = [
            argmax(@view(obj.data[ch_idx, ss[idx]:end, 1])),
            argmin(@view(obj.data[ch_idx, ss[idx]:end, 1]))
        ]

        p[ch_idx, :] .+= ss[idx]

    end

    return p

end
