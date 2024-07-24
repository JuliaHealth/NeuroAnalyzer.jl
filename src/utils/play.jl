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

    ch = get_channel(obj, ch=ch)[1]
    _check_epochs(obj, ep)

    s = @views obj.data[ch, :, ep]
    fs = sr(obj)

    !mono && (s = [s s])

    wavplay(s, fs)

    return nothing

end