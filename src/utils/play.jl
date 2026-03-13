export play

"""
    play(obj; <keyword arguments>)

Play a channel signal segment as audio.

The signal is normalised to `[−1, 1]` and scaled to `[−1000, +1000]` before playback. For best results the channel should contain an audio-range signal (e.g. speech, auditory ERP); arbitrary EEG/MEG data will be audible but may not be meaningful.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: channel name; must resolve to exactly one channel
- `seg::Tuple{Real, Real}`: time segment to play in seconds, as `(start, stop)`; both values must lie within the epoch time range
- `ep::Int64`: epoch index

# Returns

- `Nothing`

# Throws
- `ArgumentError`: if `ep` is out of range, `ch` does not resolve to exactly one channel, or `seg` boundaries fall outside the epoch time axis
"""
function play(
    obj::NeuroAnalyzer.NEURO;
    ch::String,
    seg::Tuple{Real, Real},
    ep::Int64
)::Nothing

    _check_epochs(obj, ep)
    ch = get_channel(obj, ch=ch)
    @assert length(ch) == 1 "ch must resolve to exactly one channel."
    ch = ch[1]
    _check_tuple(seg, (obj.epoch_time[1], obj.epoch_time[end]), "seg")

    # map time values to nearest sample indices on the epoch time axis
    s1 = vsearch(seg[1], obj.epoch_time)
    s2 = vsearch(seg[2], obj.epoch_time)

    # extract the signal slice and normalise to [−1000, +1000] for audio output
    s = normalize_minmax(@view obj.data[ch, s1:s2, ep]) .* 1000
    wavplay(s, sr(obj))

    return nothing

end
