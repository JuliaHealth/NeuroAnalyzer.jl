export play

"""
    play(obj; <keyword arguments>)

Interactive play channel signal as audio

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::String`: channel name
- `seg::Tuple{Real, Real}`: time segment to play
- `ep::Int64`: epoch number
"""
function play(obj::NeuroAnalyzer.NEURO; ch::String, seg::Tuple{Real, Real}, ep::Int64)

    _check_epochs(obj, ep)

    ch = get_channel(obj, ch=ch)[1]
    seg = (vsearch(seg[1], obj.epoch_time), vsearch(seg[2], obj.epoch_time))
    s = @views obj.data[ch, seg[1]:seg[2], ep]
    s = normalize_minmax(s) .* 1000
    fs = sr(obj)

    wavplay(s, fs)

    return nothing

end