export plot_hb

"""
    plot_hb(obj; <keyword arguments>)

Plot the oxygenated (HbO, red) and deoxygenated (HbR, blue) hemoglobin concentration changes of NIRS source-detector pairs: one panel per pair, both chromophores in the same panel (optionally also HbT, green). For a continuous object the segment `seg` is shown; for an epoched object the epoch `ep` (for an object created by `average_epochs()`, epoch 1 is the averaged response); markers of the stimuli are drawn as shaded blocks.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object (with concentration channels, see `od2conc()`)
- `pair::Union{Nothing, String, Vector{String}}=nothing`: source-detector pair(s), e.g. `"S1_D1"`; `nothing`: all pairs with concentration channels
- `seg::Tuple{Real, Real}=(0, 60)`: time segment in seconds (continuous objects)
- `ep::Int64=1`: epoch number (epoched objects)
- `hbt::Bool=false`: if `true`, also plot HbT
- `title::String=""`: plot title

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_hb(
    obj::NeuroAnalyzer.NEURO;
    pair::Union{Nothing, String, Vector{String}} = nothing,
    seg::Tuple{Real, Real} = (0, 60),
    ep::Int64 = 1,
    hbt::Bool = false,
    title::String = "",
)::GLMakie.Figure
    # validate
    _check_chtype(obj, "nirs")
    rec = obj.header.recording
    all_pairs = [
        p for p in _nirs_pairs(obj) if
        !isempty(_nirs_pair_ch(obj, p, "nirs_hbo")) || !isempty(_nirs_pair_ch(obj, p, "nirs_hbr"))
    ]
    isempty(all_pairs) &&
        throw(ArgumentError("OBJ does not contain concentration channels, use od2conc() first."))
    sel = if isnothing(pair)
        all_pairs
    else
        pl = pair isa String ? [pair] : pair
        for p in pl
            p in [q.label for q in all_pairs] ||
                throw(ArgumentError("Pair $p (with concentration channels) not found."))
        end
        [q for q in all_pairs if q.label in pl]
    end
    _check_epochs(obj, ep)

    # time segment
    if nepochs(obj) == 1
        _check_tuple(seg, (obj.time_pts[1], obj.time_pts[end]), "seg")
        i1, i2 = vsearch(seg[1], obj.time_pts), vsearch(seg[2], obj.time_pts)
        t = obj.time_pts[i1:i2]
        ep = 1
    else
        i1, i2 = 1, epoch_len(obj)
        t = obj.epoch_time
    end

    # prepare plot
    GLMakie.activate!(; title = "plot_hb()")
    p_n = length(sel)
    fig = GLMakie.Figure(; size = (900, max(250, 50 + 200 * p_n)))
    axs = GLMakie.Axis[]
    for (idx, p) in enumerate(sel)
        ax = GLMakie.Axis(
            fig[idx, 1];
            title              = p.label,
            xlabel             = idx == p_n ? "Time [s]" : "",
            ylabel             = "Δ [" * rec[:unit][first(_nirs_pair_ch(obj, p, "nirs_hbo"))] * "]",
            xminorticksvisible = true,
            xminorticks        = IntervalsBetween(10),
            xautolimitmargin   = (0, 0),
            yautolimitmargin   = (0.1, 0.1),
            _AXIS_LOCK_KWARGS...,
        )
        _style_axis!(ax)
        push!(axs, ax)

        # stimulus blocks (markers) within the shown time range
        if nepochs(obj) == 1 && _has_markers(obj)
            for r in eachrow(obj.markers)
                m1 = r[:start]
                m2 = r[:start] + max(r[:length], 0)
                (m2 < t[1] || m1 > t[end]) && continue
                GLMakie.vspan!(ax, max(m1, t[1]), min(max(m2, m1 + 1 / sr(obj)), t[end]);
                               color = (:gray, 0.2))
            end
        end
        nepochs(obj) > 1 && GLMakie.vlines!(ax, 0; color = :gray, linestyle = :dash)
        GLMakie.hlines!(ax, 0; color = :black, linewidth = 0.5)

        for (type, col, lbl) in (("nirs_hbo", :red, "HbO"), ("nirs_hbr", :blue, "HbR"),
                                 ("nirs_hbt", :green, "HbT"))
            (type == "nirs_hbt" && !hbt) && continue
            c = _nirs_pair_ch(obj, p, type)
            isempty(c) && continue
            GLMakie.lines!(ax, t, obj.data[c[1], i1:i2, ep]; color = col, linewidth = 1.5,
                           label = lbl)
        end
        idx == 1 && GLMakie.Legend(fig[1, 2], ax; framevisible = false)
    end
    p_n > 1 && GLMakie.linkxaxes!(axs...)
    title != "" &&
        GLMakie.Label(fig[0, 1], title; fontsize = 18, font = :bold, tellwidth = false)

    return fig
end
