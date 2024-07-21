export play

"""
    play(obj; <keyword arguments>)

Interactive play channel signal as audio

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::String`: channel name
- `ep::Int64`: epoch number
- `mono::Bool=true`: play mono or stereo
"""
function play(obj::NeuroAnalyzer.NEURO; ch::String, ep::Int64, mono::Bool=true)

    ch = _ch_idx(obj, ch)
    _check_epochs(obj, ep)

    s = @views obj.data[ch, :, ep]
    fs = sr(obj)

    !mono && (s = [s s])

    wavplay(s, fs)

    return nothing

end