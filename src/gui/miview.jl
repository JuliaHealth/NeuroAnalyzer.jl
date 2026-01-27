export miview

"""
    miview(obj; <keyword arguments>)

Interactive view of continuous signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `mch::Bool=true`: draw multichannel signal (up to 20 channels in one plot)
- `zoom::Real=10`: how many seconds are displayed in one segment
- `bad::Bool=true`: list of bad channels; if not false - plot bad channels using this list
- `snap::Bool=true`: snap region markers to grid at 0.0, 0.25, 0.5 and 0.75 time points

# Returns

- `seg::Union{Nothing, Tuple{Float64, Float64}}`
"""
function miview(obj::NeuroAnalyzer.NEURO; mch::Bool=true, zoom::Real=10, bad::Bool=true, snap::Bool=true)::Union{Nothing, Tuple{Float64, Float64}}

    @assert nepochs(obj) == 1 "For epoched object iview_ep() must be used."

    obj.time_pts[end] < zoom && (zoom = round(obj.time_pts[end]) / 2)

    @assert zoom > 0 "zoom must be > 0."
    @assert zoom <= signal_len(obj) / sr(obj) "zoom must be ≤ $(signal_len(obj) / sr(obj))."

    ch_order = obj.header.recording[:channel_order]
    cl = labels(obj)[ch_order]
    ch = 1:nchannels(obj)

    mono = false
    scale = true
    ts1 = 0
    ts2 = 0
    k = nothing

    plot_type = 0

    if mch
        if length(ch) > 15
            ch_first = 1
            ch_last = ch_first + 14
        else
            ch_first = 1
            ch_last = ch[end]
        end
    else
        ch_first = 1
        ch_last = ch[end]
    end



end
