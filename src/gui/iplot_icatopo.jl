export iplot_icatopo

"""
    iplot_icatopo(obj; <keyword arguments>)

Interactive topographical plot of embedded ("ic" and "ic_mw") ICA components.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0`: component(s) to plot, default is all components
- `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
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


# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function iplot_icatopo(obj::NeuroAnalyzer.NEURO; ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0, seg::Tuple{Real, Real}=(0, 10), amethod::Symbol=:mean, imethod::Symbol=:sh, nmethod::Symbol=:minmax)

    @assert :ic in keys(obj.components) "OBJ does not contain :ic component. Perform ica_decompose() first."
    @assert :ic_mw in keys(obj.components) "OBJ does not contain :ic_mw component. Perform ica_decompose() first."
    
    ic = obj.components[:ic]
    ic_mw = obj.components[:ic_mw]

    iplot_icatopo(obj, ic, ic_mw, ic_idx=ic_idx, seg=seg, amethod=amethod, imethod=imethod, nmethod=nmethod)

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

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function iplot_icatopo(obj::NeuroAnalyzer.NEURO, ic::Matrix{Float64}, ic_mw::Matrix{Float64}; ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0, seg::Tuple{Real, Real}=(0, 10), amethod::Symbol=:mean, imethod::Symbol=:sh, nmethod::Symbol=:minmax)

    obj_new = deepcopy(obj)
    obj_edited = false

    # select component channels, default is all channels
    ic_idx == 0 && (ic_idx = _select_cidx(ic, ic_idx))
    _check_cidx(ic, ic_idx)

    ic_remove_idx = zeros(Bool, length(ic_idx))
    ic_available_for_removal_idx = ones(Bool, length(ic_idx))

    _info("Preparing signal reconstructed from ICA components")
    cx_set = Vector{Cairo.CairoSurfaceBase{UInt32}}()
    @inbounds for idx in ic_idx
        obj_tmp = ica_reconstruct(obj, ic, ic_mw, ch=signal_channels(obj), ic_idx=idx, keep=true)
        p_tmp = plot_topo(obj_tmp, seg=seg, amethod=amethod, imethod=imethod, nmethod=nmethod, cb=true, large=false)
        cx_tmp = _p2c(p_tmp)
        push!(cx_set, cx_tmp)
    end

    nc = 4
    nr = mod(length(ic_idx), 4) == 0 ? div(length(ic_idx), 4) : div((length(ic_idx) + (4 - mod(length(ic_idx), 4))), 4)

    w = Int64(cx_set[1].width)
    h = Int64(cx_set[1].height)

    win = GtkWindow("NeuroAnalyzer: iplot_icatopo()", w * 4 + 40, h * 3 + 20)
    set_gtk_property!(win, :border_width, 20)
    set_gtk_property!(win, :resizable, true)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    set_gtk_property!(win, :startup_id, "org.neuroanalyzer")
    win_view = GtkScrolledWindow()
    set_gtk_property!(win_view, :min_content_width, w * 4 + 40)
    set_gtk_property!(win_view, :min_content_height, h * 3 + 20)

    g_opts = GtkGrid()
    set_gtk_property!(g_opts, :column_homogeneous, false)
    set_gtk_property!(g_opts, :row_spacing, 10)
    set_gtk_property!(g_opts, :column_spacing, 10)

    entry_ic = GtkSpinButton(1, length(ic_idx), 1)
    set_gtk_property!(entry_ic, :tooltip_text, "ICA component")
    current_ic = ic_idx[1]
    
    bt_details_amp = GtkButton("Show amplitude")
    set_gtk_property!(bt_details_amp, :tooltip_text, "Show component amplitude")
    bt_details_psd = GtkButton("Show PSD")
    set_gtk_property!(bt_details_psd, :tooltip_text, "Show component PSD")
    bt_details_spec = GtkButton("Show spectrogram")
    set_gtk_property!(bt_details_spec, :tooltip_text, "Show component spectrogram")
    bt_details_topo = GtkButton("Show topomap")
    set_gtk_property!(bt_details_topo, :tooltip_text, "Show component topomap")

    bt_signal_remove = GtkButton("Show signal\n  (remove)")
    set_gtk_property!(bt_signal_remove, :tooltip_text, "Show signal with the current component removed")
    bt_signal_reconstruct = GtkButton(" Show signal \n(reconstruct)")
    set_gtk_property!(bt_signal_reconstruct, :tooltip_text, "Show signal reconstructed without the current component")

    cb_mark = GtkCheckButton("")
    set_gtk_property!(cb_mark, :active, false)
    set_gtk_property!(cb_mark, :tooltip_text, "Mark this ICA component for removal")

    cb_preview = GtkCheckButton("")
    set_gtk_property!(cb_preview, :active, true)
    set_gtk_property!(cb_preview, :tooltip_text, "Preview signal after ICA component(s) removal/reconstruction")

    bt_remove = GtkButton("Remove")
    set_gtk_property!(bt_remove, :tooltip_text, "Remove marked ICA components")
    bt_reconstruct = GtkButton("Reconstruct")
    set_gtk_property!(bt_reconstruct, :tooltip_text, "Reconstruct the signal without marked ICA components")

    bt_apply = GtkButton("Apply")
    set_gtk_property!(bt_apply, :tooltip_text, "Apply changes and close this window")
    bt_cancel = GtkButton("Cancel")
    set_gtk_property!(bt_cancel, :tooltip_text, "Close this window and abandon changes")

    g_opts[1, 1] = GtkLabel("#IC:")
    g_opts[2, 1] = entry_ic
    g_opts[1, 2] = GtkLabel("Mark:")
    g_opts[2, 2] = cb_mark
    g_opts[1:2, 3] = GtkLabel("")
    g_opts[1:2, 4] = bt_details_amp
    g_opts[1:2, 5] = bt_details_psd
    g_opts[1:2, 6] = bt_details_spec
    g_opts[1:2, 7] = bt_details_topo
    g_opts[1:2, 8] = GtkLabel("")
    g_opts[1:2, 9] = bt_signal_reconstruct
    g_opts[1:2, 10] = bt_signal_remove
    g_opts[1:2, 11] = GtkLabel("")
    g_opts[1:2, 12] = bt_reconstruct
    g_opts[1:2, 13] = bt_remove
    g_opts[1, 14] = GtkLabel("Preview:")
    g_opts[2, 14] = cb_preview
    g_opts[1:2, 15] = GtkLabel("")
    g_opts[1:2, 16] = bt_apply
    g_opts[1:2, 17] = bt_cancel
    vbox = GtkBox(:v)
    push!(vbox, g_opts)

    can_set = Vector{Gtk.GtkCanvas}()
    for idx in 1:(nc * nr)
        push!(can_set, GtkCanvas(w, h))
    end

    g_cans = GtkGrid()
    set_gtk_property!(g_cans, :column_homogeneous, false)
    set_gtk_property!(g_cans, :row_spacing, 10)
    set_gtk_property!(g_cans, :column_spacing, 10)
    idx = 1
    for idx1 in 1:nr, idx2 in 1:nc
        g_cans[idx2, idx1] = can_set[idx]
        idx += 1
    end
    push!(win_view, g_cans)
    
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
            Cairo.move_to(ctx, 10.0, 12.0)
            Cairo.set_source_rgb(ctx, 0, 0, 0)
            Cairo.show_text(ctx, "#IC: $idx")
        end
    end

    signal_connect(cb_mark, "clicked") do widget
        current_ic = get_gtk_property(entry_ic, :value, Int64)
        ic_remove_idx[current_ic] = get_gtk_property(cb_mark, :active, Bool)
    end

    signal_connect(entry_ic, "value-changed") do widget
        current_ic = get_gtk_property(entry_ic, :value, Int64)
        Gtk.@sigatom begin
            set_gtk_property!(cb_mark, :sensitive, ic_available_for_removal_idx[current_ic])
            set_gtk_property!(cb_mark, :active, ic_remove_idx[current_ic])
            if ic_available_for_removal_idx[current_ic] == false
                set_gtk_property!(cb_mark, :opacity, 0.2)
            else
                set_gtk_property!(cb_mark, :opacity, 1)
            end
        end
    end

    signal_connect(bt_reconstruct, "clicked") do widget
        if length(ic_idx[ic_remove_idx]) > 0
            if ask_dialog("Reconstruct the signal without the ICA component$(_pl(ic_idx[ic_remove_idx])) $(_v2s(ic_idx[ic_remove_idx])) ?", "No", "Yes")
                current_ic = get_gtk_property(entry_ic, :value, Int64)
                _info("Reconstructing the signal without the ICA component$(_pl(ic_idx[ic_remove_idx])): $(_v2s(ic_idx[ic_remove_idx]))")
                ica_reconstruct!(obj_new, ic, ic_mw, ch=signal_channels(obj), ic_idx=ic_idx[ic_remove_idx])
                ic_available_for_removal_idx[ic_idx[ic_remove_idx]] .= false
                ic_remove_idx[ic_idx[ic_remove_idx]] .= false
                Gtk.@sigatom begin
                    set_gtk_property!(cb_mark, :sensitive, ic_available_for_removal_idx[current_ic])
                    if ic_available_for_removal_idx[current_ic] == false
                        set_gtk_property!(cb_mark, :opacity, 0.2)
                    else
                        set_gtk_property!(cb_mark, :opacity, 1)
                    end
                    set_gtk_property!(cb_mark, :active, ic_remove_idx[current_ic])
                end
                if get_gtk_property(cb_preview, :active, Bool)
                    iview(obj, obj_new)
                end
            end
        else
            warn_dialog("No ICA components marked for removal!")
        end
    end

    signal_connect(bt_remove, "clicked") do widget
        if length(ic_idx[ic_remove_idx]) > 0
            if ask_dialog("Remove ICA component$(_pl(ic_idx[ic_remove_idx])) $(_v2s(ic_idx[ic_remove_idx])) ?", "No", "Yes")
                current_ic = get_gtk_property(entry_ic, :value, Int64)
                _info("Removing ICA component$(_pl(ic_idx[ic_remove_idx])): $(_v2s(ic_idx[ic_remove_idx]))")
                ica_remove!(obj_new, ic, ic_mw, ch=signal_channels(obj), ic_idx=ic_idx[ic_remove_idx])
                ic_available_for_removal_idx[ic_idx[ic_remove_idx]] .= false
                ic_remove_idx[ic_idx[ic_remove_idx]] .= false
                Gtk.@sigatom begin
                    set_gtk_property!(cb_mark, :sensitive, ic_available_for_removal_idx[current_ic])
                    if ic_available_for_removal_idx[current_ic] == false
                        set_gtk_property!(cb_mark, :opacity, 0.2)
                    else
                        set_gtk_property!(cb_mark, :opacity, 1)
                    end
                    set_gtk_property!(cb_mark, :active, ic_remove_idx[current_ic])
                end
                if get_gtk_property(cb_preview, :active, Bool)
                    iview(obj, obj_new)
                end
            end
        else
            warn_dialog("No ICA components marked for removal!")
        end
    end

    signal_connect(bt_details_amp, "clicked") do widget
        current_ic = get_gtk_property(entry_ic, :value, Int64)
        iview(obj, ic, c_idx=current_ic)
    end

    signal_connect(bt_details_psd, "clicked") do widget
        current_ic = get_gtk_property(entry_ic, :value, Int64)
        iview_plot(plot_psd(obj, ic, c_idx=current_ic, seg=seg))
    end

    signal_connect(bt_details_spec, "clicked") do widget
        current_ic = get_gtk_property(entry_ic, :value, Int64)
        iview_plot(plot_spectrogram(obj, ic, c_idx=current_ic, seg=seg))
    end

    signal_connect(bt_details_topo, "clicked") do widget
        current_ic = get_gtk_property(entry_ic, :value, Int64)
        obj_tmp = ica_reconstruct(obj, ic, ic_mw, ch=signal_channels(obj), ic_idx=current_ic, keep=true)
        iview_plot(plot_topo(obj_tmp, cb=true, seg=seg, amethod=amethod, imethod=imethod, nmethod=nmethod, large=true))
    end

    signal_connect(bt_signal_reconstruct, "clicked") do widget
        current_ic = get_gtk_property(entry_ic, :value, Int64)
        ica_reconstruct!(obj_new, ic, ic_mw, ch=signal_channels(obj), ic_idx=current_ic, keep=true)
        iview(obj, obj_new)
        obj_edited = true
    end
    
    signal_connect(bt_signal_remove, "clicked") do widget
        current_ic = get_gtk_property(entry_ic, :value, Int64)
        ica_remove!(obj_new, ic, ic_mw, ch=signal_channels(obj), ic_idx=current_ic)
        iview(obj, obj_new)
        obj_edited = true
    end
    
    signal_connect(bt_apply, "clicked") do widget
        if obj_edited
            if ask_dialog("This operation will apply all changes.\nPlease confirm.", "No", "Yes")
                obj.header = obj_new.header
                obj.data = obj_new.data
                obj.locs = obj_new.locs
                obj.components = obj_new.components
                obj.history = obj_new.history
                obj.markers = obj_new.markers
                Gtk.destroy(win)
                return nothing
            end
        else
            warn_dialog("There are no changes to apply!")
        end
    end

    signal_connect(bt_cancel, "clicked") do widget
        Gtk.destroy(win)
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        s = event.state
        if s == 20
            if k == 113 # q
                Gtk.destroy(win)
            end
        end
    end

    c = Condition()
    signal_connect(win, :destroy) do widget
        notify(c)
    end
    @async Gtk.gtk_main()
    wait(c)

    return nothing

end
