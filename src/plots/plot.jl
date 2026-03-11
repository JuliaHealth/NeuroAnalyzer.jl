export plot

"""
    plot(obj; <keyword arguments>)

Plot signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{String, Vector{String}, Regex}="all"`: channel name or list of channel names
- `ep::Int64=1`: first epoch to display
- `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
- `tm::Union{Nothing, Int64, Vector{Int64}}=nothing`: time markers (in milliseconds) to plot as vertical lines, useful for adding topoplots at these time points (for ERP/ERF)
- `rt::Union{Nothing, Real, AbstractVector}=nothing`: response time for each epoch; if provided, the response time line will be plotted over the `:stack` plot (for ERP/ERF)
- `xlabel::String="default"`: x-axis label
- `ylabel::String="default"`: y-axis label
- `title::String="default"`: plot title
- `markers::Bool=true`: draw markers if available
- `scale::Bool=true`: draw scales
- `group_ch::Bool=true`: group channels by type
- `type::Symbol=:normal`: plot type:
    - `:normal`
    - `:butterfly`
    - `:stack`: for ERP/ERF/MEP
    - `:topo`: for ERP/ERF
- `avg::Bool=false`: plot averaged channel in butterfly plot
- `ci95::Bool=false`: plot averaged channels and 95% CI in butterfly plot
- `n_channels::Int64=20`: number of visible channels
- `n_epochs::Int64=5`: number of visible epochs
- `cb::Bool=true`: plot color bar (for ERP/ERF/MEP)
- `cb_title::String="default"`: color bar title (for ERP/ERF/MEP)
- `peaks::Bool=true`: draw peaks (for ERP/ERF/MEP)
- `leg::Bool=true`: if true, add legend with channel labels (for ERP/ERF/MEP)
- `yrev::Bool=false`: reverse y-axis (for ERP/ERF/MEP)
- `smooth::Bool=false`: smooth the image using Gaussian blur (for ERP/ERF/MEP)
- `ks::Int64=3`: kernel size of the Gaussian blur (larger kernel means more smoothing) (for ERP/ERF/MEP)
- `zl::Bool=true`: draw line at t = 0 (for ERP/ERF/MEP)
- `mono::Bool=false`: use color or gray palette
- `res::Int64=1`: resampling factor (draw every res-nth sample)
- `snap::Bool=true`: snap to grid when placing markers
- `gui::Bool=true`: if true, keep window open and use it interactively

# Returns

- `p::GLMakie.Figure`
"""
function plot(
        obj::NeuroAnalyzer.NEURO;
        ch::Union{String, Vector{String}, Regex} = "all",
        ep::Int64 = 1,
        seg::Tuple{Real, Real} = (0, 10),
        tm::Union{Nothing, Int64, Vector{Int64}} = nothing,
        rt::Union{Nothing, Real, AbstractVector} = nothing,
        xlabel::String = "default",
        ylabel::String = "default",
        title::String = "default",
        markers::Bool = true,
        scale::Bool = true,
        group_ch::Bool = true,
        type::Symbol = :normal,
        avg::Bool = true,
        ci95::Bool = false,
        n_channels::Int64 = 20,
        n_epochs::Int64 = 5,
        cb::Bool = true,
        cb_title::String = "default",
        peaks::Bool = true,
        leg::Bool = true,
        yrev::Bool = false,
        smooth::Bool = false,
        ks::Int64 = 3,
        zl::Bool = true,
        mono::Bool = false,
        res::Int64 = 1,
        snap::Bool = true,
        gui::Bool = true,
    )::GLMakie.Figure

    n_channels > nchannels(obj) && (n_channels = nchannels(obj))
    n_epochs > nepochs(obj) && (n_epochs = nepochs(obj))

    if datatype(obj) in ["erp", "erf"]
        p = plot_erp(
            obj;
            ch = ch,
            tm = tm,
            rt = rt,
            xlabel = xlabel,
            ylabel = ylabel,
            title = title,
            type = type,
            avg = avg,
            ci95 = ci95,
            cb = cb,
            cb_title = cb_title,
            peaks = peaks,
            leg = leg,
            yrev = yrev,
            smooth = smooth,
            ks = ks,
            zl = zl,
            mono = mono,
            gui = gui,
        )
    elseif datatype(obj) == "mep"
        p = plot_mep(
            obj;
            ch = ch,
            xlabel = xlabel,
            ylabel = ylabel,
            title = title,
            type = type,
            avg = avg,
            ci95 = ci95,
            cb = cb,
            cb_title = cb_title,
            peaks = peaks,
            leg = leg,
            yrev = yrev,
            smooth = smooth,
            ks = ks,
            zl = zl,
            mono = mono,
            gui = gui,
        )
    else
        if nepochs(obj) == 1
            p = plot_cont(
                obj;
                ch = ch,
                seg = seg,
                xlabel = xlabel,
                ylabel = ylabel,
                title = title,
                markers = markers,
                scale = scale,
                group_ch = group_ch,
                type = type,
                avg = avg,
                ci95 = ci95,
                n_channels = n_channels,
                mono = mono,
                res = res,
                snap = snap,
                gui = gui,
            )
        else
            p = plot_ep(
                obj;
                ch = ch,
                ep = ep,
                xlabel = xlabel,
                ylabel = ylabel,
                title = title,
                markers = markers,
                scale = scale,
                group_ch = group_ch,
                type = type,
                avg = avg,
                ci95 = ci95,
                n_channels = n_channels,
                n_epochs = n_epochs,
                mono = mono,
                res = res,
                gui = gui,
            )
        end
    end

    return p

end


"""
    plot(obj1, obj2; <keyword arguments>)

Plot signal.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `obj2::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{String, Vector{String}, Regex}="all"`: channel name or list of channel names
- `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
- `xlabel::String="default"`: x-axis label
- `ylabel::String="default"`: y-axis label
- `title::String="default"`: plot title
- `scale::Bool=true`: draw scales
- `group_ch::Bool=true`: group channels by type
- `n_channels::Int64=20`: number of visible channels
- `n_epochs::Int64=5`: number of visible epochs
- `res::Int64=1`: resampling factor (draw every res-nth sample)
- `gui::Bool=true`: if true, keep window open and use it interactively

# Returns

- `p::GLMakie.Figure`
"""
function plot(
        obj1::NeuroAnalyzer.NEURO,
        obj2::NeuroAnalyzer.NEURO;
        ch::Union{String, Vector{String}, Regex} = "all",
        seg::Tuple{Real, Real} = (0, 10),
        xlabel::String = "default",
        ylabel::String = "default",
        title::String = "default",
        scale::Bool = true,
        group_ch::Bool = true,
        n_channels::Int64 = 20,
        res::Int64 = 1,
        gui::Bool = true,
    )::GLMakie.Figure

    @assert datatype(obj1) in ["eeg", "meg"] "This function works for continuous EEG and MEG objects."
    @assert datatype(obj2) in ["eeg", "meg"] "This function works for continuous EEG and MEG objects."
    @assert nepochs(obj1) == 1 "This function works for continuous EEG and MEG objects."

    p = plot_cont(
        obj1,
        obj2;
        ch = ch,
        seg = seg,
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        scale = scale,
        group_ch = group_ch,
        n_channels = n_channels,
        res = res,
        gui = gui,
    )

    return p

end

export plot

"""
    plot(t, s; <keyword arguments>)

Plot continuous signal.

# Arguments

- `t::AbstractVector`: time points
- `s::AbstractVector`: signal vector: signal data
- `xlabel::String="Time [s]"`: x-axis label
- `ylabel::String="Amplitude"`: y-axis label
- `title::String=""`: plot title

# Returns

- `p::GLMakie.Figure`
"""
function plot(
        t::AbstractVector, s::AbstractVector; xlabel::String = "Time [s]", ylabel::String = "Amplitude", title::String = ""
    )::GLMakie.Figure

    @assert length(t) == length(s) "Length of s must equal length of t."

    # prepare plot
    GLMakie.activate!(title = "plot()")
    plot_size = (900, 450)
    p = GLMakie.Figure(size = plot_size)
    ax = GLMakie.Axis(
        p[1, 1];
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        xticks = LinearTicks(10),
        xminorticksvisible = true,
        xminorticks = IntervalsBetween(10),
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
        xzoomlock = true,
        yzoomlock = true,
        xpanlock = true,
        ypanlock = true,
        xrectzoom = false,
        yrectzoom = false,
    )
    if minimum(s) == 0
        GLMakie.ylims!(ax, 0, extrema(s)[2] * 1.1)
    else
        GLMakie.ylims!(ax, extrema(s) .* 1.1)
    end
    ax.titlesize = 18
    ax.xlabelsize = 12
    ax.ylabelsize = 12
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    GLMakie.lines!(ax, t, s, color = :black)

    return p

end
