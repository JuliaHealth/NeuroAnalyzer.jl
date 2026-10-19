export mep_peaks

"""
    mep_peaks(obj)

Detect a pair of positive and negative peaks of MEP.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object

# Returns

- `Matrix{Int64}`: peak positions (samples), shape (channels, 2): positive peak, negative peak (searched from 20 samples after the stimulation sample)
"""
function mep_peaks(obj::NeuroAnalyzer.NEURO)::Matrix{Int64}
    # validate
    _check_datatype(obj, "mep")

    # number of channels
    ch_n = nchannels(obj)

    # pre-allocate output
    p = zeros(Int64, ch_n, 2)

    # calculate over channels
    ss = obj.header.recording[:stimulation_sample] .+ 20
    @inbounds Threads.@threads :static for ch_idx = 1:ch_n
        p[ch_idx, :] = [
            argmax(@view(obj.data[ch_idx, ss[ch_idx]:end, 1])),
            argmin(@view(obj.data[ch_idx, ss[ch_idx]:end, 1])),
        ]

        # slice positions → sample positions
        p[ch_idx, :] .+= ss[ch_idx] - 1
    end

    return p
end
