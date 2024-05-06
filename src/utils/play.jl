export play

"""
    play(obj; <keyword arguments>)

Interactive play channel signal as audio

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Int64`: channel number
- `ep::Int64`: epoch number
- `mono::Bool=true`: play mono or stereo
- `maxvol::Bool=false`: play at maximum volume (scaled to unit amplitude)
"""
function play(obj::NeuroAnalyzer.NEURO; ch::Int64, ep::Int64, mono::Bool=true, maxvol::Bool=false)

    _check_channels(obj, ch)
    _check_epochs(obj, ep)

    s = @views obj.data[ch, :, ep]
    fs = sr(obj)

    if maxvol == true
        s = normalize_minmax(s)
    end

    if mono == false
        s = [s s]
    end
    
    wavplay(s, fs)

    return nothing

end