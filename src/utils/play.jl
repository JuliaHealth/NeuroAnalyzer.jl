export play

"""
    play(obj; <keyword arguments>)

Interactive play channel signal as audio

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Int64`: channel number
- `ep::Int64`: epoch number
- `mono::Bool=true`: play mono or stereo
"""
function play(obj::NeuroAnalyzer.NEURO; ch::Int64, ep::Int64, mono::Bool=true, maxvol::Bool=false)

    _check_channels(obj, ch)
    _check_epochs(obj, ep)

    s = @views obj.data[ch, :, ep]
    fs = sr(obj)

    mono == false && (s = [s s])
    
    wavplay(s, fs)

    return nothing

end