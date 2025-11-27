export iview_ica

"""
    iview_ica(obj; <keyword arguments>)

Interactive view of embedded ("ic" and "ic_mw") ICA components.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names

# Returns

Nothing
"""
function iview_ica(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::Nothing

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
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names

# Returns

Nothing
"""
function iview_ica(obj::NeuroAnalyzer.NEURO, ic::Matrix{Float64}, ic_mw::Matrix{Float64}; ch::Union{String, Vector{String}, Regex})::Nothing

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

    @assert size(ic_mw, 1) == length(chn) "ICA weighting matrix size does not match number of OBJ channels."
    @assert size(ic_mw, 2) == size(ic, 1) "ICA weighting matrix size does not match number of ICA components."
    @assert size(ic, 2) == signal_len(obj) "ICA components length does not match OBJ signal length."
    @assert size(ic_mw, 1) >= length(ch) "ICA weighting matrix size does not match number of selected channels."

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
    ica_can_set = Vector{Gtk4.GtkCanvas}()
    for idx in 1:length(ic_idx)
        push!(ica_can_set, GtkCanvas(Int64(ica_set[1].width), Int64(ica_set[1].height)))
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

    p_sig = NeuroAnalyzer.plot(obj, ch=cl[chn[ch_idx]], title="Channel: $(cl[chn[ch_idx]]) (original)")
    p_psd = NeuroAnalyzer.plot_psd(obj, ch=cl[chn[ch_idx]], title="Channel: $(cl[chn[ch_idx]]) (original)")

    k = nothing
    scaled_sig = false
    scaled_psd = false

    function _activate(app)

        win = GtkApplicationWindow(app, "NeuroAnalyzer: iview_ica()")
        win.width_request = Int64(ica_set[1].width) + round(Int64, p_sig.attr[:size][1] * 0.75) + 20
        win.height_request = round(Int64, p_sig.attr[:size][2] * 0.75) + round(Int64, p_psd.attr[:size][2] * 0.75) + 20

        ica_view = GtkScrolledWindow()
        g_cans = GtkGrid()
        g_cans.column_homogeneous = false
        g_cans.row_spacing = 5
        g_cans.column_spacing = 5
        g_cans.margin_start = 5
        g_cans.margin_end = 5
        g_cans.margin_top = 5
        g_cans.margin_bottom = 5

        for idx in eachindex(ica_can_set)
            g_cans[1, idx] = ica_can_set[idx]
        end
        ica_view.child = g_cans

        signal_view = GtkCanvas(round(Int64, p_sig.attr[:size][1] * 0.75), round(Int64, p_sig.attr[:size][2] * 0.75))
        psd_view = GtkCanvas(round(Int64, p_psd.attr[:size][1] * 0.75), round(Int64, p_psd.attr[:size][2] * 0.75))
        ica_view.min_content_width = Int64(ica_set[1].width) + 10
        ica_view.max_content_height = round(Int64, p_sig.attr[:size][2] * 0.75) + round(Int64, p_psd.attr[:size][2] * 0.75) + 20

        g_opts = GtkGrid()
        g_opts.column_homogeneous = false
        g_opts.row_spacing = 5
        g_opts.column_spacing = 5
        g_opts.margin_start = 5
        g_opts.margin_end = 5
        g_opts.margin_top = 5
        g_opts.margin_bottom = 5

        entry_ic = GtkSpinButton(1, length(ic_idx), 1)
        entry_ic.tooltip_text = "ICA component"
        current_ic = ic_idx[1]

        cb_mark = GtkCheckButton("Mark component")
        cb_mark.active = false
        cb_mark.tooltip_text = "Mark this ICA component for removal"

        bt_remove = GtkButton("Remove")
        bt_remove.tooltip_text = "Remove marked ICA components from the signal"
        bt_reconstruct = GtkButton("Reconstruct")
        bt_reconstruct.tooltip_text = "Reconstruct the signal from marked ICA components"

        bt_apply = GtkButton("Apply")
        bt_apply.tooltip_text = "Close this window and apply changes"
        bt_close = GtkButton("Close")
        bt_close.tooltip_text = "Close this window and abandon changes"
        bt_help = GtkButton("Help")
        bt_help.tooltip_text = "Show help"

        entry_time = GtkSpinButton(obj.time_pts[1], obj.time_pts[end] - zoom, 1)
        entry_time.digits = 2
        entry_time.value = obj.time_pts[1]
        entry_time.tooltip_text = "Time position [s]"
        bt_start = GtkButton("⇤")
        bt_start.tooltip_text = "Go to the signal start"
        bt_prev = GtkButton("←")
        bt_prev.tooltip_text = "Go back by $(round(zoom)) seconds"
        bt_next = GtkButton("→")
        bt_next.tooltip_text = "Go forward by $(round(zoom)) seconds"
        bt_end = GtkButton("⇥")
        bt_end.tooltip_text = "Go to the signal end"

        ch_slider = GtkScale(:h, 1:length(chn))
        ch_slider.draw_value = false
        ch_slider.tooltip_text = "Scroll channels"
        ch_slider.vexpand = true
        oc = GtkOrientable(ch_slider)
        oc.orientation = 1

        signal_slider = GtkScale(:h, obj.time_pts[1]:obj.time_pts[end] - zoom)
        signal_slider.draw_value = false
        signal_slider.tooltip_text = "Time position"

        combo_sig = GtkComboBoxText()
        for idx in ["signal (original)", "signal (reconstructed from IC)", "signal (IC removed)", "IC"]
            push!(combo_sig, idx)
        end
        combo_sig.active = 0
        combo_sig.tooltip_text = "Viewed signal"

        combo_psd = GtkComboBoxText()
        for idx in ["signal (original)", "signal (reconstructed from IC)", "signal (IC removed)", "IC"]
            push!(combo_psd, idx)
        end
        combo_psd.active = 0
        combo_psd.tooltip_text = "Viewed PSD"

        lab_amp = GtkLabel("Amplitude:")
        lab_amp.halign = 2
        lab_psd = GtkLabel("PSD:")
        lab_psd.halign = 2

        g = GtkGrid()
        g.column_homogeneous = false
        g.row_spacing = 5
        g.column_spacing = 5
        g.margin_start = 5
        g.margin_end = 5
        g.margin_top = 5
        g.margin_bottom = 5

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

        Gtk4.show(win)

        @guarded draw(signal_view) do widget
            plot_sig_type = combo_sig.active
            time1 = entry_time.value
            time2 = time1 + zoom
            time2 > obj.time_pts[end] && (time2 = obj.time_pts[end])
            current_ic = Int64(entry_ic.value)
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
            if !scaled_sig
                Cairo.scale(ctx_sig, 0.75, 0.75)
                scaled_sig = true
            end
            img_sig = read_from_png(io)
            set_source_surface(ctx_sig, img_sig, 0, 0)
            paint(ctx_sig)
        end

        @guarded draw(psd_view) do widget
            plot_psd_type = combo_psd.active
            time1 = get_gtk_property(entry_time, :value, Float64)
            time2 = time1 + zoom
            time2 > obj.time_pts[end] && (time2 = obj.time_pts[end])
            current_ic = Int64(entry_ic.value)
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
            if !scaled_psd
                Cairo.scale(ctx_psd, 0.75, 0.75)
                scaled_psd = true
            end
            img_psd = read_from_png(io)
            set_source_surface(ctx_psd, img_psd, 0, 0)
            paint(ctx_psd)
        end

        function _mwheel_scroll(_, dx, dy)
            if dy > 0.5
                if k == 0x0000ffe1 # shift
                    time_current = entry_time.value
                    if time_current < obj.time_pts[end] - zoom
                        time_current += 1
                    else
                        time_current = obj.time_pts[end] - zoom
                    end
                    @idle_add entry_time.value = time_current
                elseif k == 0x0000ffe9 # alt
                    time_current = entry_time.value
                    if time_current < obj.time_pts[end] - zoom
                        time_current += zoom
                    else
                        time_current = obj.time_pts[end] - zoom
                    end
                    @idle_add entry_time.value = time_current
                elseif k == 0x0000ffe3 # ctrl
                    if plot_sig_type != 3
                        if ch_idx < length(chn)
                            ch_idx += 1
                            @idle_add Gtk4.value(ch_slider, chn[ch_idx])
                            draw(signal_view)
                            draw(psd_view)
                        end
                    end
                end
            elseif dy < -0.5
                if k == 0x0000ffe1 # shift
                    time_current = entry_time.value
                    if time_current >= obj.time_pts[1] + 1
                        time_current -= 1
                        @idle_add entry_time.value = time_current
                    end
                elseif k == 0x0000ffe9 # alt
                    time_current = entry_time.value
                    if time_current >= obj.time_pts[1] + zoom
                        time_current = time_current - zoom
                        @idle_add entry_time.value = time_current
                    end
                elseif k == 0x0000ffe3 # ctrl
                    if plot_sig_type != 3
                        if ch_idx > 1
                            ch_idx -= 1
                            @idle_add Gtk4.value(ch_slider, chn[ch_idx])
                            draw(signal_view)
                            draw(psd_view)
                        end
                    end
                end
            end
        end
        ecsf = Gtk4.GtkEventControllerScroll(Gtk4.EventControllerScrollFlags_VERTICAL)
        push!(signal_view, ecsf)
        push!(psd_view, ecsf)
        signal_connect(_mwheel_scroll, ecsf, "scroll")

        signal_connect(signal_slider, "value-changed") do widget
            @idle_add entry_time.value = round(Gtk4.value(signal_slider))
        end

        signal_connect(entry_time, "value-changed") do widget
            @idle_add Gtk4.value(signal_slider, entry_time.value)
            draw(signal_view)
            draw(psd_view)
            time1 = entry_time.value
            time2 = time1 + zoom
            time2 > obj.time_pts[end] && (time2 = obj.time_pts[end])
            _refresh_ica_can_set(obj_reconstructed, ica_can_set, ic_idx, time1, time2)
        end

        signal_connect(ch_slider, "value-changed") do widget, others...
            ch_idx = round(Int64, Gtk4.value(ch_slider))
            draw(signal_view)
            draw(psd_view)
        end

        signal_connect(combo_sig, "changed") do widget
            draw(signal_view)
        end

        signal_connect(combo_psd, "changed") do widget
            draw(psd_view)
        end

        signal_connect(cb_mark, "toggled") do widget
            current_ic = Int64(entry_ic.value)
            ic_remove_idx[current_ic] = cb_mark.active
        end

        signal_connect(entry_ic, "value-changed") do widget
            current_ic = Int64(entry_ic.value)
            cb_mark.sensitive = ic_available_for_removal_idx[current_ic]
            cb_mark.active = ic_remove_idx[current_ic]
            if !ic_available_for_removal_idx[current_ic]
                cb_mark.opacity = 0.2
            else
                cb_mark.opacity = 1
            end
            draw(signal_view)
            draw(psd_view)
        end

        signal_connect(bt_reconstruct, "clicked") do widget
            if length(ic_idx[ic_remove_idx]) > 0
                ask_dialog("Reconstruct the signal from marked ICA component$(_pl(ic_idx[ic_remove_idx])): $(_v2s(ic_idx[ic_remove_idx])) ?", win) do ans
                    if ans
                        current_ic = Int64(entry_ic.value)
                        _info("Reconstructing the signal using the IC$(_pl(ic_idx[ic_remove_idx])): $(_v2s(ic_idx[ic_remove_idx]))")
                        ica_reconstruct!(obj_new, ic, ic_mw, ch=ch, ic_idx=ic_idx[ic_remove_idx], keep=true)
                        ic_available_for_removal_idx[ic_idx[ic_remove_idx]] .= false
                        ic_remove_idx[ic_idx[ic_remove_idx]] .= false
                        cb_mark.sensitive = ic_available_for_removal_idx[current_ic]
                        if !ic_available_for_removal_idx[current_ic]
                            cb_mark.opacity = 0.2
                        else
                            cb_mark.opacity = 1
                        end
                        cb_mark.active = ic_remove_idx[current_ic]
                        draw(signal_view)
                        draw(psd_view)
                        obj_edited = true
                    end
                end
            else
                warn_dialog(_nill, "No ICA components marked for reconstruction!", win)
            end
        end

        signal_connect(bt_remove, "clicked") do widget
            if length(ic_idx[ic_remove_idx]) > 0
                ask_dialog("Remove marked ICA component$(_pl(ic_idx[ic_remove_idx])) from the signal: $(_v2s(ic_idx[ic_remove_idx])) ?", win) do ans
                    if ans
                        current_ic = Int64(entry_ic.value)
                        _info("Removing IC$(_pl(ic_idx[ic_remove_idx])): $(_v2s(ic_idx[ic_remove_idx]))")
                        ica_reconstruct!(obj_new, ic, ic_mw, ch=ch, ic_idx=ic_idx[ic_remove_idx])
                        ic_available_for_removal_idx[ic_idx[ic_remove_idx]] .= false
                        ic_remove_idx[ic_idx[ic_remove_idx]] .= false
                        cb_mark.sensitive = ic_available_for_removal_idx[current_ic]
                        if !ic_available_for_removal_idx[current_ic]
                            cb_mark.opacity = 0.2
                        else
                            cb_mark.opacity = 1
                        end
                        cb_mark.active = ic_remove_idx[current_ic]
                        draw(signal_view)
                        draw(psd_view)
                        obj_edited = true
                    end
                end
            else
                warn_dialog(_nill, "No ICA components marked for removal!", win)
            end
        end

        signal_connect(bt_apply, "clicked") do widget
            if obj_edited
                ask_dialog("This operation will apply all changes.\nPlease confirm.", win) do ans
                    if ans
                        obj.header = obj_new.header
                        obj.data = obj_new.data
                        obj.locs = obj_new.locs
                        obj.components = obj_new.components
                        obj.history = obj_new.history
                        obj.markers = obj_new.markers
                        close(win)
                        return nothing
                    end
                end
            else
                warn_dialog(_nill, "There are no changes to apply!", win)
            end
        end

        signal_connect(bt_prev, "clicked") do widget
            time_current = entry_time.value
            if time_current >= obj.time_pts[1] + zoom
                time_current = time_current - zoom
                @idle_add Gtk4.value(entry_time, time_current)
            end
        end

        signal_connect(bt_next, "clicked") do widget
            time_current = entry_time.value
            if time_current < obj.time_pts[end] - zoom
                time_current += zoom
            else
                time_current = obj.time_pts[end] - zoom
            end
            @idle_add Gtk4.value(entry_time, time_current)
        end

        signal_connect(bt_start, "clicked") do widget
            @idle_add Gtk4.value(entry_time, obj.time_pts[1])
        end

        signal_connect(bt_end, "clicked") do widget
            time_current = obj.time_pts[end] - zoom
            @idle_add Gtk4.value(entry_time, time_current)
        end

        signal_connect(bt_close, "clicked") do widget
            close(win)
        end

        win_key = Gtk4.GtkEventControllerKey(win)

        signal_connect(win_key, "key-pressed") do widget, keyval, keycode, state
            k = keyval
            # CONTROL
            if ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('q'))
                close(win)
            end
        end
    end

    app = GtkApplication("org.neuroanalyzer.iview")
    Gtk4.signal_connect(_activate, app, :activate)
    Gtk4.GLib.stop_main_loop()
    Gtk4.run(app)

    return nothing

end
