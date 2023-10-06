export iplot_icatopo

"""
    iplot_icatopo(obj; <keyword arguments>)

Interactive topographical plot of embedded ("ic" and "ic_mw") ICA components.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0`: component(s) to plot, default is all components
- `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
- `cb::Bool=false`: plot color bar
- `cb_label::String="[A.U.]"`: color bar label
- `amethod::Symbol=:mean`: averaging method:
    - `:mean`
    - `:median`
- `imethod::Symbol=:sh`: interpolation method:
    - `:sh`: Shepard
    - `:mq`: Multiquadratic
    - `:imq`: InverseMultiquadratic
    - `:tp`: ThinPlate
    - `:nn`: NearestNeighbour
    - `:ga`: Gaussian
- `nmethod::Symbol=:minmax`: method for normalization, see `normalize()`
- `plot_contours::Bools=true`: plot contours over topo plot
- `plot_electrodes::Bools=true`: plot electrodes over topo plot
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function iplot_icatopo(obj::NeuroAnalyzer.NEURO; ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0, seg::Tuple{Real, Real}=(0, 10), cb::Bool=false, cb_label::String="default", amethod::Symbol=:mean, imethod::Symbol=:sh, nmethod::Symbol=:minmax, plot_contours::Bool=true, plot_electrodes::Bool=true, kwargs...)

    _wip()

    @assert :ic in keys(obj.components) "OBJ does not contain :ic component. Perform ica_decompose() first."
    @assert :ic_mw in keys(obj.components) "OBJ does not contain :ic_mw component. Perform ica_decompose() first."
    
    ic = obj.components[:ic]
    ic_mw = obj.components[:ic_mw]

    # select component channels, default is all channels
    ic_idx == 0 && (ic_idx = _select_cidx(ic, ic_idx))
    _check_cidx(ic, ic_idx)

    p_topo = Vector{Plots.Plot{Plots.GRBackend}}()

    for idx in eachindex(ic_idx)
        obj_tmp = ica_reconstruct(obj, ic, ic_mw, ch=signal_channels(obj), ic_idx=ic_idx[idx], keep=true)
        p_tmp = plot_topo(obj_tmp, title="IC $(ic_idx[idx])", cb=cb, cb_label=cb_label, amethod=amethod, imethod=imethod, nmethod=nmethod, plot_contours=plot_contours, plot_electrodes=plot_electrodes, seg=seg, kwargs...)
        push!(p_topo, p_tmp)
    end
    
    if length(ic_idx) <= 4
        p = plot_compose(p_topo, layout=(1, 4))
    elseif length(ic_idx) <= 8
        p = plot_compose(p_topo, layout=(2, ceil(Int64, length(ic_idx) / 2)))
    elseif length(ic_idx) <= 12
        p = plot_compose(p_topo, layout=(3, ceil(Int64, length(ic_idx) / 3)))
    else
        p = plot_compose(p_topo, layout=(4, ceil(Int64, length(ic_idx) / 4)))
    end

    return p

end

"""
    iplot_icatopo(obj, ic, ic_mw; <keyword arguments>)

Interactive topographical plot of external ICA components.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ic::Matrix{Float64}`: components IC(1)..IC(n)
- `ic_mw::Matrix{Float64}`: weighting matrix IC(1)..IC(n)
- `ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0`: component(s) to plot, default is all components
- `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
- `cb::Bool=false`: plot color bar
- `cb_label::String="[A.U.]"`: color bar label
- `amethod::Symbol=:mean`: averaging method:
    - `:mean`
    - `:median`
- `imethod::Symbol=:sh`: interpolation method:
    - `:sh`: Shepard
    - `:mq`: Multiquadratic
    - `:imq`: InverseMultiquadratic
    - `:tp`: ThinPlate
    - `:nn`: NearestNeighbour
    - `:ga`: Gaussian
- `nmethod::Symbol=:minmax`: method for normalization, see `normalize()`
- `plot_contours::Bools=true`: plot contours over topo plot
- `plot_electrodes::Bools=true`: plot electrodes over topo plot
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function iplot_icatopo(obj::NeuroAnalyzer.NEURO, ic::Matrix{Float64}, ic_mw::Matrix{Float64}; ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0, seg::Tuple{Real, Real}=(0, 10), amethod::Symbol=:mean, imethod::Symbol=:sh, nmethod::Symbol=:minmax, plot_contours::Bool=true, plot_electrodes::Bool=true, kwargs...)

    _wip()

    # select component channels, default is all channels
    ic_idx == 0 && (ic_idx = _select_cidx(ic, ic_idx))
    _check_cidx(ic, ic_idx)

    cx_set = Vector{Cairo.CairoSurfaceBase{UInt32}}()
    @inbounds @simd for idx in ic_idx
        obj_tmp = ica_reconstruct(obj, ic, ic_mw, ch=signal_channels(obj), ic_idx=idx, keep=true)
        p_tmp = plot_topo(obj_tmp, cb=false, large=false)
        push!(cx_set, NeuroAnalyzer._p2c(p_tmp))
    end

    nc = 4
    nr = mod(length(ic_idx), 4) == 0 ? div(length(ic_idx), 4) : div((length(ic_idx) + (4 - mod(length(ic_idx), 4))), 4)

    w = Int64(cx_set[1].width) + 10
    h = Int64(cx_set[1].height) + 10

    win = GtkWindow("NeuroAnalyzer: iplot_icatopo()", w * 4 + 10, h * 3 - 10)
    set_gtk_property!(win, :border_width, 20)
    set_gtk_property!(win, :resizable, true)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    win_view = GtkScrolledWindow()
    set_gtk_property!(win_view, :min_content_width, w * 4 + 10)
    set_gtk_property!(win_view, :min_content_height, h * 3 - 10)

    g_opts = GtkGrid()
    set_gtk_property!(g_opts, :column_homogeneous, false)
    set_gtk_property!(g_opts, :row_spacing, 10)
    set_gtk_property!(g_opts, :column_spacing, 10)
    entry_ic = GtkSpinButton(1, length(ic_idx), 1)
    set_gtk_property!(entry_ic, :tooltip_text, "ICA component")
    
    g_opts[2, 1] = entry_ic
    g_opts[1, 1] = GtkLabel("#IC:")
    vbox = GtkBox(:v)
    push!(vbox, g_opts)

    can_set = Vector{Gtk.GtkCanvas}()
    for idx in 1:(nc * nr)
        push!(can_set, GtkCanvas(w, h))
    end

    cans_g = GtkGrid()
    idx = 1
    for idx1 in 1:nr, idx2 in 1:nc
        cans_g[idx2, idx1] = can_set[idx]
        idx += 1
    end
    push!(win_view, cans_g)
    
    g = GtkGrid()
    set_gtk_property!(g, :column_homogeneous, false)
    set_gtk_property!(g, :row_spacing, 10)
    set_gtk_property!(g, :column_spacing, 10)
    g[1, 1] = vbox
    g[2, 1] = win_view
    push!(win, g)
    showall(win)

    for idx in ic_idx
        can = can_set[idx]
        @guarded draw(can) do widget
            ctx = getgc(can)
            Cairo.set_source_surface(ctx, cx_set[idx], 0, 0)
            Cairo.paint(ctx)
        end
    end

end