export iview_ica

"""
    iview_ica(obj; <keyword arguments>)

Interactive view of embedded ("ic" and "ic_mw") ICA components.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{String, Vector{String}}`: channel name or list of channel names

# Returns

Nothing
"""
function iview_ica(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}})::Nothing

    @assert :ic in keys(obj.components) "OBJ does not contain :ic component. Perform ica_decompose() first."
    @assert :ic_mw in keys(obj.components) "OBJ does not contain :ic_mw component. Perform ica_decompose() first."

    iview_ica(obj, obj.components[:ic], obj.components[:ic_mw], ch=ch)

end

"""
    iview_ica(obj, ic, ic_mw; <keyword arguments>)

Interactive view of external ICA components.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ic::Matrix{Float64}`: components IC(1)..IC(n)
- `ic_mw::Matrix{Float64}`: weighting matrix IC(1)..IC(n)
- `ch::Union{String, Vector{String}}`: channel name or list of channel names

# Returns

Nothing
"""
function iview_ica(obj::NeuroAnalyzer.NEURO, ic::Matrix{Float64}, ic_mw::Matrix{Float64}; ch::Union{String, Vector{String}})::Nothing

    @assert size(ic_mw, 1) == nchannels(obj) "ICA weighting matrix size does not match number of OBJ channels."
    @assert size(ic_mw, 2) == size(ic, 1) "ICA weighting matrix size does not match number of ICA components."
    @assert size(ic, 2) == signal_len(obj) "ICA components length does not match OBJ signal length."
    @assert size(ic_mw, 1) >= length(ch) "ICA weighting matrix size does not match number of selected channels."

    plot_sig_type = 0
    plot_psd_type = 0

    # select all ICA components
    ic_idx = _select_cidx(ic, 1:size(ic, 1))

    zoom = 10
    obj.time_pts[end] < zoom && (zoom = obj.time_pts[end])
    seg = (obj.time_pts[1], obj.time_pts[1] + zoom)
    chn = get_channel(obj, ch=ch)
    cl = labels(obj)
    ch_idx = 1

    obj_new = deepcopy(obj)
    obj_edited = false

    ic_remove_idx = zeros(Bool, length(ic_idx))
    ic_available_for_removal_idx = ones(Bool, length(ic_idx))

    _info("Preparing reconstructed objects")
    obj_reconstructed = Vector{NeuroAnalyzer.NEURO}()
    obj_removed = Vector{NeuroAnalyzer.NEURO}()
    @inbounds for idx in ic_idx
        push!(obj_reconstructed, ica_reconstruct(obj, ic, ic_mw, ch=ch, ic_idx=idx, keep=true))
        push!(obj_removed, ica_reconstruct(obj, ic, ic_mw, ch=ch, ic_idx=idx))
    end

    # ICA topos
    ica_set = Vector{Cairo.CairoSurfaceBase{UInt32}}()
    for idx in ic_idx
        p_tmp = plot_topo(obj_reconstructed[idx], ch=datatype(obj_reconstructed[1]), seg=seg, amethod=:mean, imethod=:sh, nmethod=:minmax, cb=false, large=false)
        cx_tmp = plot2canvas(p_tmp)
        push!(ica_set, cx_tmp)
    end
    ica_can_set = Vector{Gtk.GtkCanvas}()
    for idx in 1:length(ic_idx)
        push!(ica_can_set, GtkCanvas(ica_set[1].width, ica_set[1].height))
    end
    for idx in ic_idx
        ica_can = ica_can_set[idx]
        @guarded draw(ica_can) do widget
            ctx_ica = getgc(ica_can)
            Cairo.set_source_surface(ctx_ica, ica_set[idx], 0, 0)
            Cairo.paint(ctx_ica)
            Cairo.move_to(ctx_ica, 10.0, 12.0)
            Cairo.set_source_rgb(ctx_ica, 0, 0, 0)
            Cairo.show_text(ctx_ica, "IC: $idx")
        end
    end
    ica_view = GtkScrolledWindow()
    g_cans = GtkGrid()
    set_gtk_property!(g_cans, :column_homogeneous, false)
    set_gtk_property!(g_cans, :row_spacing, 5)
    set_gtk_property!(g_cans, :column_spacing, 5)
    for idx in eachindex(ica_can_set)
        g_cans[1, idx] = ica_can_set[idx]
    end
    push!(ica_view, g_cans)

    p_sig = NeuroAnalyzer.plot(obj, ch=cl[chn[ch_idx]], title="Channel: $(cl[chn[ch_idx]]) (original)")
    p_psd = NeuroAnalyzer.plot_psd(obj, ch=cl[chn[ch_idx]], title="Channel: $(cl[chn[ch_idx]]) (original)")
    win = GtkWindow("NeuroAnalyzer: iview_ica()", ica_set[1].width + Int32(p_sig.attr[:size][1]) + 20, Int32(p_sig.attr[:size][2]) + Int32(p_psd.attr[:size][2]) + 20)
    set_gtk_property!(win, :border_width, 5)
    set_gtk_property!(win, :resizable, false)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    set_gtk_property!(win, :startup_id, "org.neuroanalyzer")
    signal_view = GtkCanvas(Int32(p_sig.attr[:size][1]), Int32(p_sig.attr[:size][2]))
    psd_view = GtkCanvas(Int32(p_psd.attr[:size][1]), Int32(p_psd.attr[:size][2]))
    set_gtk_property!(ica_view, :min_content_width, ica_set[1].width + 10)
    set_gtk_property!(ica_view, :max_content_height, Int32(p_sig.attr[:size][2]) + Int32(p_psd.attr[:size][2]) + 20)

    g_opts = GtkGrid()
    set_gtk_property!(g_opts, :column_homogeneous, false)
    set_gtk_property!(g_opts, :row_spacing, 5)
    set_gtk_property!(g_opts, :column_spacing, 5)

    entry_ic = GtkSpinButton(1, length(ic_idx), 1)
    set_gtk_property!(entry_ic, :tooltip_text, "ICA component")
    current_ic = ic_idx[1]

    cb_mark = GtkCheckButton("Mark component")
    set_gtk_property!(cb_mark, :active, false)
    set_gtk_property!(cb_mark, :tooltip_text, "Mark this ICA component for removal")

    bt_remove = GtkButton("Remove")
    set_gtk_property!(bt_remove, :tooltip_text, "Remove marked ICA components from the signal")
    bt_reconstruct = GtkButton("Reconstruct")
    set_gtk_property!(bt_reconstruct, :tooltip_text, "Reconstruct the signal from marked ICA components")

    bt_apply = GtkButton("Apply")
    set_gtk_property!(bt_apply, :tooltip_text, "Close this window and apply changes")
    bt_close = GtkButton("Close")
    set_gtk_property!(bt_close, :tooltip_text, "Close this window and abandon changes")
    bt_help = GtkButton("Help")
    set_gtk_property!(bt_help, :tooltip_text, "Show help")

    entry_time = GtkSpinButton(obj.time_pts[1], obj.time_pts[end] - zoom, 1)
    set_gtk_property!(entry_time, :digits, 2)
    set_gtk_property!(entry_time, :value, obj.time_pts[1])
    set_gtk_property!(entry_time, :tooltip_text, "Time position [s]")
    bt_start = GtkButton("⇤")
    set_gtk_property!(bt_start, :tooltip_text, "Go to the signal start")
    bt_prev = GtkButton("←")
    set_gtk_property!(bt_prev, :tooltip_text, "Go back by $(round(zoom)) seconds")
    bt_next = GtkButton("→")
    set_gtk_property!(bt_next, :tooltip_text, "Go forward by $(round(zoom)) seconds")
    bt_end = GtkButton("⇥")
    set_gtk_property!(bt_end, :tooltip_text, "Go to the signal end")

    ch_slider = GtkScale(false, 1:length(chn))
    set_gtk_property!(ch_slider, :draw_value, false)
    set_gtk_property!(ch_slider, :tooltip_text, "Scroll channels")
    set_gtk_property!(ch_slider, :vexpand, true)
    oc = GtkOrientable(ch_slider)
    set_gtk_property!(oc, :orientation, 1)

    signal_slider = GtkScale(false, obj.time_pts[1]:obj.time_pts[end] - zoom)
    set_gtk_property!(signal_slider, :draw_value, false)
    set_gtk_property!(signal_slider, :tooltip_text, "Time position")

    combo_sig = GtkComboBoxText()
    for idx in ["signal (original)", "signal (reconstructed from IC)", "signal (IC removed)", "IC"]
        push!(combo_sig, idx)
    end
    set_gtk_property!(combo_sig, :active, 0)
    set_gtk_property!(combo_sig, :tooltip_text, "Viewed signal")

    combo_psd = GtkComboBoxText()
    for idx in ["signal (original)", "signal (reconstructed from IC)", "signal (IC removed)", "IC"]
        push!(combo_psd, idx)
    end
    set_gtk_property!(combo_psd, :active, 0)
    set_gtk_property!(combo_psd, :tooltip_text, "Viewed PSD")

    lab_amp = GtkLabel("Amplitude:")
    set_gtk_property!(lab_amp, :halign, 2)
    lab_psd = GtkLabel("PSD:")
    set_gtk_property!(lab_psd, :halign, 2)

    g = GtkGrid()
    set_gtk_property!(g, :column_homogeneous, false)
    set_gtk_property!(g, :row_spacing, 5)
    set_gtk_property!(g, :column_spacing, 5)
    g[1, 1:2] = ica_view
    g[2:8, 1] = psd_view
    g[2:8, 2] = signal_view
    g[9, 2] = ch_slider
    g[1, 3] = entry_ic
    g[1, 4] = cb_mark
    g[1, 5] = bt_help
    g[2:8, 3] = signal_slider
    g[2, 4] = lab_psd
    g[3, 4] = combo_psd
    g[2, 5] = lab_amp
    g[3, 5] = combo_sig
    g[4, 4] = bt_prev
    g[4, 5] = bt_start
    g[5, 4] = entry_time
    g[6, 4] = bt_next
    g[6, 5] = bt_end
    g[7, 4] = bt_reconstruct
    g[7, 5] = bt_remove
    g[8, 4] = bt_apply
    g[8, 5] = bt_close

    push!(win, g)
    showall(win)

    @guarded draw(signal_view) do widget
        plot_sig_type = get_gtk_property(combo_sig, :active, Int64)
        time1 = get_gtk_property(entry_time, :value, Float64)
        time2 = time1 + zoom
        time2 > obj.time_pts[end] && (time2 = obj.time_pts[end])
        current_ic = get_gtk_property(entry_ic, :value, Int64)
        if plot_sig_type == 0
            p_sig = NeuroAnalyzer.plot(obj_new,
                                       ch=cl[chn[ch_idx]],
                                       seg=(time1, time2),
                                       mono=true,
                                       title="Channel: $(cl[chn[ch_idx]]) (original)")
        elseif plot_sig_type == 1
            p_sig = NeuroAnalyzer.plot(obj_new,
                                       obj_reconstructed[current_ic],
                                       ch=cl[chn[ch_idx]],
                                       seg=(time1, time2),
                                       title="Channel: $(cl[chn[ch_idx]]) (reconstructed from IC: $current_ic)")
        elseif plot_sig_type == 2
            p_sig = NeuroAnalyzer.plot(obj_new,
                                       obj_removed[current_ic],
                                       ch=cl[chn[ch_idx]],
                                       seg=(time1, time2),
                                       title="Channel: $(cl[chn[ch_idx]]) (removed IC: $current_ic)")
        elseif plot_sig_type == 3
            p_sig = NeuroAnalyzer.plot(obj_new,
                                       ic,
                                       ch=ch,
                                       c_idx=current_ic,
                                       seg=(time1, time2),
                                       mono=true,
                                       title="IC: $(current_ic)")
        end
        show(io, MIME("image/png"), p_sig)
        ctx_sig = getgc(signal_view)
        img_sig = read_from_png(io)
        set_source_surface(ctx_sig, img_sig, 0, 0)
        paint(ctx_sig)
    end

    @guarded draw(psd_view) do widget
        plot_psd_type = get_gtk_property(combo_psd, :active, Int64)
        time1 = get_gtk_property(entry_time, :value, Float64)
        time2 = time1 + zoom
        time2 > obj.time_pts[end] && (time2 = obj.time_pts[end])
        current_ic = get_gtk_property(entry_ic, :value, Int64)
        if plot_psd_type == 0
            p_psd = NeuroAnalyzer.plot_psd(obj_new,
                                           ch=cl[chn[ch_idx]],
                                           seg=(time1, time2),
                                           mono=true,
                                           title="Channel: $(cl[chn[ch_idx]]) (original)")
        elseif plot_psd_type == 1
            p_psd = NeuroAnalyzer.plot_psd(obj_reconstructed[current_ic],
                                           ch=cl[chn[ch_idx]],
                                           seg=(time1, time2),
                                           title="Channel: $(cl[chn[ch_idx]]) (reconstructed from IC: $current_ic)")
        elseif plot_psd_type == 2
            p_psd = NeuroAnalyzer.plot_psd(obj_removed[current_ic],
                                           ch=cl[chn[ch_idx]],
                                           seg=(time1, time2),
                                           title="Channel: $(cl[chn[ch_idx]]) (removed IC: $current_ic)")
        elseif plot_psd_type == 3
            p_psd = NeuroAnalyzer.plot_psd(obj_new,
                                           ic,
                                           ch=ch,
                                           c_idx=current_ic,
                                           seg=(time1, time2),
                                           mono=true,
                                           title="IC: $(current_ic)")
        end
        show(io, MIME("image/png"), p_psd)
        ctx_psd = getgc(psd_view)
        img_psd = read_from_png(io)
        set_source_surface(ctx_psd, img_psd, 0, 0)
        paint(ctx_psd)
    end

    signal_connect(signal_slider, "value-changed") do widget
        Gtk.@sigatom begin
            set_gtk_property!(entry_time, :value, round(GAccessor.value(signal_slider)))
        end
        draw(signal_view)
        draw(psd_view)
        time1 = get_gtk_property(entry_time, :value, Float64)
        time2 = time1 + zoom
        time2 > obj.time_pts[end] && (time2 = obj.time_pts[end])
        _refresh_ica_can_set(obj_reconstructed, ica_can_set, ic_idx, time1, time2)
    end

    signal_connect(ch_slider, "value-changed") do widget, others...
        ch_idx = round(Int64, GAccessor.value(ch_slider))
        draw(signal_view)
        draw(psd_view)
    end

    signal_view.mouse.scroll = @guarded (widget, event) -> begin
        plot_sig_type = get_gtk_property(combo_sig, :active, Int64)
        s = event.state
        if event.direction == 1 # down
            if s == 0x00000001
                time_current = get_gtk_property(entry_time, :value, Float64)
                if time_current < obj.time_pts[end] - zoom
                    time_current += 1
                else
                    time_current = obj.time_pts[end] - zoom
                end
                Gtk.@sigatom begin
                    set_gtk_property!(entry_time, :value, time_current)
                end
            elseif s == 0x00000004
                time_current = get_gtk_property(entry_time, :value, Float64)
                if time_current < obj.time_pts[end] - zoom
                    time_current += zoom
                else
                    time_current = obj.time_pts[end] - zoom
                end
                Gtk.@sigatom begin
                    set_gtk_property!(entry_time, :value, time_current)
                end
            else
                if plot_sig_type != 3
                    if ch_idx < length(chn)
                        ch_idx += 1
                        Gtk.@sigatom begin
                            GAccessor.value(ch_slider, chn[ch_idx])
                        end
                        draw(signal_view)
                        draw(psd_view)
                    end
                end
            end
        elseif event.direction == 0 # up
            if s == 0x00000001
                time_current = get_gtk_property(entry_time, :value, Float64)
                if time_current >= obj.time_pts[1] + 1
                    time_current -= 1
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_time, :value, time_current)
                    end
                end
            elseif s == 0x00000004
                time_current = get_gtk_property(entry_time, :value, Float64)
                if time_current >= obj.time_pts[1] + zoom
                    time_current = time_current - zoom
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_time, :value, time_current)
                    end
                end
            else
                if plot_sig_type != 3
                    if ch_idx > 1
                        ch_idx -= 1
                        Gtk.@sigatom begin
                            GAccessor.value(ch_slider, chn[ch_idx])
                        end
                        draw(signal_view)
                        draw(psd_view)
                    end
                end
            end
        end
    end

    signal_connect(entry_time, "value-changed") do widget
        Gtk.@sigatom begin
            GAccessor.value(signal_slider, get_gtk_property(entry_time, :value, Float64))
        end
        draw(signal_view)
        draw(psd_view)
    end

    signal_connect(combo_sig, "changed") do widget
        draw(signal_view)
    end

    signal_connect(combo_psd, "changed") do widget
        draw(psd_view)
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
            if !ic_available_for_removal_idx[current_ic]
                set_gtk_property!(cb_mark, :opacity, 0.2)
            else
                set_gtk_property!(cb_mark, :opacity, 1)
            end
        end
        draw(signal_view)
        draw(psd_view)
    end

    signal_connect(bt_reconstruct, "clicked") do widget
        if length(ic_idx[ic_remove_idx]) > 0
            if ask_dialog("Reconstruct the signal from marked ICA component$(_pl(ic_idx[ic_remove_idx])): $(_v2s(ic_idx[ic_remove_idx])) ?", "No", "Yes")
                current_ic = get_gtk_property(entry_ic, :value, Int64)
                _info("Reconstructing the signal using the IC$(_pl(ic_idx[ic_remove_idx])): $(_v2s(ic_idx[ic_remove_idx]))")
                ica_reconstruct!(obj_new, ic, ic_mw, ch=ch, ic_idx=ic_idx[ic_remove_idx], keep=true)
                ic_available_for_removal_idx[ic_idx[ic_remove_idx]] .= false
                ic_remove_idx[ic_idx[ic_remove_idx]] .= false
                Gtk.@sigatom begin
                    set_gtk_property!(cb_mark, :sensitive, ic_available_for_removal_idx[current_ic])
                    if !ic_available_for_removal_idx[current_ic]
                        set_gtk_property!(cb_mark, :opacity, 0.2)
                    else
                        set_gtk_property!(cb_mark, :opacity, 1)
                    end
                    set_gtk_property!(cb_mark, :active, ic_remove_idx[current_ic])
                    draw(signal_view)
                    draw(psd_view)
                    obj_edited = true
                end
            end
        else
            warn_dialog("No ICA components marked for reconstruction!")
        end
    end

    signal_connect(bt_remove, "clicked") do widget
        if length(ic_idx[ic_remove_idx]) > 0
            if ask_dialog("Remove marked ICA component$(_pl(ic_idx[ic_remove_idx])) from the signal: $(_v2s(ic_idx[ic_remove_idx])) ?", "No", "Yes")
                current_ic = get_gtk_property(entry_ic, :value, Int64)
                _info("Removing IC$(_pl(ic_idx[ic_remove_idx])): $(_v2s(ic_idx[ic_remove_idx]))")
                ica_reconstruct!(obj_new, ic, ic_mw, ch=ch, ic_idx=ic_idx[ic_remove_idx])
                ic_available_for_removal_idx[ic_idx[ic_remove_idx]] .= false
                ic_remove_idx[ic_idx[ic_remove_idx]] .= false
                Gtk.@sigatom begin
                    set_gtk_property!(cb_mark, :sensitive, ic_available_for_removal_idx[current_ic])
                    if !ic_available_for_removal_idx[current_ic]
                        set_gtk_property!(cb_mark, :opacity, 0.2)
                    else
                        set_gtk_property!(cb_mark, :opacity, 1)
                    end
                    set_gtk_property!(cb_mark, :active, ic_remove_idx[current_ic])
                end
                draw(signal_view)
                draw(psd_view)
                obj_edited = true
            end
        else
            warn_dialog("No ICA components marked for removal!")
        end
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

    signal_connect(bt_prev, "clicked") do widget
        time_current = get_gtk_property(entry_time, :value, Float64)
        if time_current >= obj.time_pts[1] + zoom
            time_current = time_current - zoom
            Gtk.@sigatom begin
                set_gtk_property!(entry_time, :value, time_current)
            end
        end
    end

    signal_connect(bt_next, "clicked") do widget
        time_current = get_gtk_property(entry_time, :value, Float64)
        if time_current < obj.time_pts[end] - zoom
            time_current += zoom
        else
            time_current = obj.time_pts[end] - zoom
        end
        Gtk.@sigatom begin
            set_gtk_property!(entry_time, :value, time_current)
        end
    end

    signal_connect(bt_start, "clicked") do widget
        Gtk.@sigatom begin
            set_gtk_property!(entry_time, :value, obj.time_pts[1])
        end
    end

    signal_connect(bt_end, "clicked") do widget
        time_current = obj.time_pts[end] - zoom
        Gtk.@sigatom begin
            set_gtk_property!(entry_time, :value, time_current)
        end
    end

    signal_connect(bt_close, "clicked") do widget
        Gtk.destroy(win)
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        s = event.state
        if s == 0x00000004 || s == 0x00000014 # ctrl
            if k == 113 # q
                Gtk.destroy(win)
            end
        end

    end

    cnd = Condition()
    signal_connect(win, :destroy) do widget
        notify(cnd)
    end
    @async Gtk.gtk_main()
    wait(cnd)

    return nothing

end
