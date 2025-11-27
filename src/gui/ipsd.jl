export ipsd
export ipsd_ep

"""
    ipsd(obj; <keyword arguments>)

Interactive PSD of continuous signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::String`: channel name
- `zoom::Real=10`: how many seconds are displayed in one segment

# Returns

Nothing
"""
function ipsd(obj::NeuroAnalyzer.NEURO; ch::String, zoom::Real=10)::Nothing

    @assert nepochs(obj) == 1 "For epoched object ipsd_ep() must be used."

    obj.time_pts[end] < zoom && (zoom = round(obj.time_pts[end]) / 2)

    @assert zoom > 0 "zoom must be > 0."
    @assert zoom <= signal_len(obj) / sr(obj) "zoom must be ≤ $(signal_len(obj) / sr(obj))."

    ch_init = ch
    ch = get_channel(obj, ch=ch)
    clabels = labels(obj)

    k = nothing
    mono = false

    p = NeuroAnalyzer.plot_psd(obj, ch=clabels[ch])

    function _activate(app)

        win = GtkApplicationWindow(app, "NeuroAnalyzer: ipsd()")
        win.width_request = 1200
        win.height_request = 850

        can = GtkCanvas(p.attr[:size][1], p.attr[:size][2])
        win_view = GtkScrolledWindow()
        win_view.min_content_width = 1200
        win_view.min_content_height = 800
        win_view.child = can

        g = GtkGrid()
        g.column_homogeneous = false
        g.column_spacing = 5
        g.row_spacing = 5
        g.margin_start = 5
        g.margin_end = 5
        g.margin_top = 5
        g.margin_bottom = 5

        g_opts = GtkGrid()
        g_opts.column_homogeneous = false
        g_opts.column_spacing = 5
        g_opts.row_spacing = 5
        g_opts.margin_start = 5
        g_opts.margin_end = 5
        g_opts.margin_top = 5
        g_opts.margin_bottom = 5

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
        bt_help = GtkButton("Help")
        bt_help.tooltip_text = "Show help"
        bt_close = GtkButton("Close")
        bt_close.tooltip_text = "Close this window"

        entry_title = GtkEntry()
        entry_title.text = "default"
        entry_title.tooltip_text = "Plot title"
        entry_xlab = GtkEntry()
        entry_xlab.text = "default"
        entry_xlab.tooltip_text = "X axis label"
        entry_ylab = GtkEntry()
        entry_ylab.text = "default"
        entry_ylab.tooltip_text = "Y axis label"

        combo_ch = GtkComboBoxText()
        ctypes = uppercase.(unique(obj.header.recording[:channel_type]))
        ch_types = [ctypes; clabels]
        for idx in ch_types
            push!(combo_ch, idx)
        end
        if ch_init in lowercase.(ctypes)
            n = findfirst(isequal(ch_init), lowercase.(ctypes))
            combo_ch.active = n - 1
        else
            combo_ch.active = (length(ctypes) + ch[1] - 1)
        end
        combo_ch.tooltip_text = "Channels"

        cb_mono = GtkCheckButton()
        cb_mono.tooltip_text = "Use color or gray palette"

        cb_db = GtkCheckButton()
        cb_db.tooltip_text = "Normalize powers to dB"
        cb_db.active = true

        cb_hw = GtkCheckButton()
        cb_hw.tooltip_text = "Apply Hanning window"
        cb_hw.active = true

        combo_method = GtkComboBoxText()
        psd_methods = ["Welch's periodogram", "fast Fourier transform", "short-time Fourier transform", "multi-taper", "Morlet wavelet", "Gaussian-Hilbert transform"]
        for idx in psd_methods
            push!(combo_method, idx)
        end
        combo_method.active = 0
        combo_method.tooltip_text = "PSD method"

        combo_type = GtkComboBoxText()
        psd_types = ["normal", "butterfly", "mean", "w3d", "s3d", "topo"]
        for idx in psd_types
            push!(combo_type, idx)
        end
        combo_type.active = 0
        combo_type.tooltip_text = "Plot type"
        if length(ch) > 1
            combo_type.sensitive = true
        else
            combo_type.sensitive = false
        end

        combo_ref = GtkComboBoxText()
        ref_types = ["absolute", "total power", "delta", "theta", "alpha", "alpha lower", "alpha higher", "beta", "beta lower", "beta higher", "gamma", "gamma 1", "gamma 2", "gamma lower", "gamma higher"]
        for idx in ref_types
            push!(combo_ref, idx)
        end
        combo_ref.active = 0
        combo_ref.tooltip_text = "Powers referenced to"

        cb_frq = GtkCheckButton()
        cb_frq.tooltip_text = "Linear or logarithmic frequencies"
        cb_frq.active = true

        entry_nt = GtkSpinButton(1, 128, 1)
        entry_nt.value = 8
        entry_nt.tooltip_text = "Number of Slepian tapers (for MT)"

        entry_ncyc = GtkSpinButton(1, 256, 1)
        entry_ncyc.value = 32
        entry_ncyc.tooltip_text = "Number of Morlet wavelet cycles"

        entry_frq1 = GtkSpinButton(0.0, (sr(obj) / 2) - 0.5, 0.5)
        entry_frq1.value = 0.0
        entry_frq1.tooltip_text = "Start frequency"

        entry_frq2 = GtkSpinButton(0.0, sr(obj) / 2, 0.5)
        entry_frq2.value = sr(obj) / 2
        entry_frq2.tooltip_text = "End frequency"

        entry_wlen = GtkSpinButton(2, zoom * sr(obj) + 1, 1)
        entry_wlen.value = sr(obj)
        entry_wlen.tooltip_text = "Window length (samples)"

        entry_woverlap = GtkSpinButton(0, zoom * sr(obj), 1)
        entry_woverlap.value = round(Int64, sr(obj) * 0.97)
        entry_woverlap.tooltip_text = "Window overlap (samples)"

        entry_gw = GtkSpinButton(1, 128, 1)
        entry_gw.value = 6
        entry_gw.tooltip_text = "Gaussian width in Hz"

        bt_refresh = GtkButton("Refresh")
        bt_refresh.tooltip_text = "Refresh the plot"

        lab_type = GtkLabel("Plot type:")
        lab_type.halign = 2
        lab_method = GtkLabel("PSD method:")
        lab_method.halign = 2
        lab_ch = GtkLabel("Channels:")
        lab_ch.halign = 2
        lab_t = GtkLabel("Title:")
        lab_t.halign = 2
        lab_x = GtkLabel("X lab:")
        lab_x.halign = 2
        lab_y = GtkLabel("Y lab:")
        lab_y.halign = 2
        lab_mono = GtkLabel("Grayscale:")
        lab_mono.halign = 2
        lab_norm = GtkLabel("Normalize:")
        lab_norm.halign = 2
        lab_nt = GtkLabel("Slepians:")
        lab_nt.halign = 2
        lab_wlen = GtkLabel("Window length:")
        lab_wlen.halign = 2
        lab_woverlap = GtkLabel("Window overlap:")
        lab_woverlap.halign = 2
        lab_frq1 = GtkLabel("Start frequency:")
        lab_frq1.halign = 2
        lab_frq2 = GtkLabel("End frequency:")
        lab_frq2.halign = 2
        lab_nc = GtkLabel("Cycles:")
        lab_nc.halign = 2
        lab_ref = GtkLabel("PSD reference:")
        lab_ref.halign = 2
        lab_frq = GtkLabel("Linear frequencies:")
        lab_frq.halign = 2
        lab_hw = GtkLabel("Hanning:")
        lab_hw.halign = 2
        lab_gw = GtkLabel("Gaussian width:")
        lab_gw.halign = 2

        signal_slider = GtkScale(:h, obj.time_pts[1]:obj.time_pts[end] - zoom)
        signal_slider.draw_value = false
        signal_slider.tooltip_text = "Time position"

        g_opts[1, 1] = lab_ch
        g_opts[1, 2] = lab_type
        g_opts[1, 3] = lab_method
        g_opts[1, 4] = lab_ref
        g_opts[1, 5] = lab_t
        g_opts[1, 6] = lab_x
        g_opts[1, 7] = lab_y
        g_opts[1, 8] = lab_frq1
        g_opts[1, 9] = lab_frq2
        g_opts[1, 10] = lab_nc
        g_opts[1, 11] = lab_nt
        g_opts[1, 12] = lab_wlen
        g_opts[1, 13] = lab_woverlap
        g_opts[1, 14] = lab_gw
        g_opts[1, 15] = lab_frq
        g_opts[1, 16] = lab_norm
        g_opts[1, 17] = lab_mono
        g_opts[1, 18] = lab_hw
        g_opts[2, 1] = combo_ch
        g_opts[2, 2] = combo_type
        g_opts[2, 3] = combo_method
        g_opts[2, 4] = combo_ref
        g_opts[2, 5] = entry_title
        g_opts[2, 6] = entry_xlab
        g_opts[2, 7] = entry_ylab
        g_opts[2, 8] = entry_frq1
        g_opts[2, 9] = entry_frq2
        g_opts[2, 10] = entry_ncyc
        g_opts[2, 11] = entry_nt
        g_opts[2, 12] = entry_wlen
        g_opts[2, 13] = entry_woverlap
        g_opts[2, 14] = entry_gw
        g_opts[2, 15] = cb_frq
        g_opts[2, 16] = cb_db
        g_opts[2, 17] = cb_mono
        g_opts[2, 18] = cb_hw
        g_opts[1:2, 19] = bt_refresh
        vbox = GtkBox(:v)
        push!(vbox, g_opts)

        g[1, 1] = vbox
        g[2:9, 1] = win_view
        g[2:9, 2] = signal_slider
        g[2, 3] = bt_start
        g[3, 3] = bt_prev
        g[4, 3] = entry_time
        g[5, 3] = bt_next
        g[6, 3] = bt_end
        g[7, 3] = GtkLabel("")
        g[8, 3] = bt_help
        g[9, 3] = bt_close
        push!(win, g)

        Gtk4.show(win)

        @guarded draw(can) do widget
            ch = ch_types[Int64(combo_ch.active) + 1]
            ch in ctypes && (ch = lowercase(ch))
            title = entry_title.text
            xlab = entry_xlab.text
            ylab = entry_ylab.text
            mono = cb_mono.active
            db = cb_db.active
            hw = cb_hw.active
            method = combo_method.active
            method == 0 && (method = :welch)
            method == 1 && (method = :fft)
            method == 2 && (method = :stft)
            method == 3 && (method = :mt)
            method == 4 && (method = :mw)
            method == 5 && (method = :gh)
            type = combo_type.active
            type == 0 && (type = :normal)
            type == 1 && (type = :butterfly)
            type == 2 && (type = :mean)
            type == 3 && (type = :w3d)
            type == 4 && (type = :s3d)
            type == 5 && (type = :topo)
            ref = combo_ref.active
            ref == 0 && (ref = :abs)
            ref == 1 && (ref = :total)
            ref == 2 && (ref = :delta)
            ref == 3 && (ref = :theta)
            ref == 4 && (ref = :alpha)
            ref == 5 && (ref = :alpha_lower)
            ref == 6 && (ref = :alpha_higher)
            ref == 7 && (ref = :beta)
            ref == 8 && (ref = :beta_lower)
            ref == 9 && (ref = :beta_higher)
            ref == 10 && (ref = :gamma)
            ref == 11 && (ref = :gamma_1)
            ref == 12 && (ref = :gamma_2)
            ref == 13 && (ref = :gamma_lower)
            ref == 14 && (ref = :gamma_higher)
            frq = cb_frq.active ? :lin : :log
            gw = entry_gw.value
            frq1 = entry_frq1.value
            frq2 = entry_frq2.value
            nt = Int64(entry_nt.value)
            ncyc = Int64(entry_ncyc.value)
            wlen = Int64(entry_wlen.value)
            woverlap = Int64(entry_woverlap.value)

            no_error = true
            if frq1 == frq2
                warn_dialog(_nill, "Start and end frequencies must be different.", win)
                no_error = false
            elseif frq1 > frq2
                warn_dialog(_nill, "Start frequency must be < end frequency.", win)
                no_error = false
            elseif woverlap >= wlen
                warn_dialog(_nill, "Window overlap must be < window length.", win)
                no_error = false
            elseif length(get_channel(obj, ch=ch)) < 2 && type === :butterfly
                warn_dialog(_nill, "For butterfly plot, the signal must contain ≥ 2 channels.", win)
                no_error = false
            elseif length(get_channel(obj, ch=ch)) < 2 && type === :mean
                warn_dialog(_nill, "For mean plot, the signal must contain ≥ 2 channels.", win)
                no_error = false
            elseif length(get_channel(obj, ch=ch)) < 2 && type === :w3d
                warn_dialog(_nill, "For w3d plot, the signal must contain ≥ 2 channels.", win)
                no_error = false
            elseif length(get_channel(obj, ch=ch)) < 2 && type === :s3d
                warn_dialog(_nill, "For s3d plot, the signal must contain ≥ 2 channels.", win)
                no_error = false
            elseif DataFrames.nrow(obj.locs) == 0 && type === :topo
                warn_dialog(_nill, "Electrode locations not available.", win)
                no_error = false
            elseif length(unique(obj.header.recording[:channel_type][get_channel(obj, ch=ch)])) > 1 && (type in [:butterfly, :mean, :w3d, :s3d, :topo] || ch == "all")
                warn_dialog(_nill, "For multi-channel $(string(type)) plot all channels must be of the same type.", win)
                no_error = false
            end

            if no_error
                time1 = entry_time.value
                time2 = time1 + zoom
                time2 > obj.time_pts[end] && (time2 = obj.time_pts[end])
                p = NeuroAnalyzer.plot_psd(obj,
                                           ch=ch,
                                           seg=(time1, time2),
                                           mono=mono,
                                           title=title,
                                           xlabel=xlab,
                                           ylabel=ylab,
                                           db=db,
                                           method=method,
                                           type=type,
                                           frq=frq,
                                           ref=ref,
                                           frq_lim=(frq1, frq2),
                                           ncyc=ncyc,
                                           nt=nt,
                                           wlen=wlen,
                                           woverlap=woverlap,
                                           w=hw,
                                           gw=gw)
                img = read_from_png(io)
                can.width_request = p.attr[:size][1]
                can.height_request = p.attr[:size][2]
                ctx = getgc(can)
                show(io, MIME("image/png"), p)
                img = read_from_png(io)
                set_source_surface(ctx, img, 0, 0)
                paint(ctx)
            end
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
                end
            end
        end
        ecsf = Gtk4.GtkEventControllerScroll(Gtk4.EventControllerScrollFlags_VERTICAL)
        push!(can, ecsf)
        signal_connect(_mwheel_scroll, ecsf, "scroll")

        signal_connect(combo_ch, "changed") do widget
            ch = Int64(combo_ch.active) + 1
            if ch in 1:length(ctypes)
                ch = lowercase(ctypes[ch])
                if length(get_channel(obj, type=ch)) > 1
                    combo_type.sensitive = true
                else
                    combo_type.active = 0
                    combo_type.sensitive = false
                end
            else
                combo_type.active = 0
                combo_type.sensitive = false
            end
            draw(can)
        end

        signal_connect(bt_refresh, "clicked") do widget
            draw(can)
        end
        signal_connect(signal_slider, "value-changed") do widget
            @idle_add entry_time.value = round(Gtk4.value(signal_slider))
        end
        signal_connect(entry_time, "value-changed") do widget
            @idle_add Gtk4.value(signal_slider, entry_time.value)
            draw(can)
        end
        signal_connect(combo_type, "changed") do widget
            draw(can)
        end
        signal_connect(combo_method, "changed") do widget
            draw(can)
        end
        signal_connect(combo_ref, "changed") do widget
            draw(can)
        end
        signal_connect(cb_frq, "toggled") do widget
            draw(can)
        end
        signal_connect(entry_frq1, "value-changed") do widget
            draw(can)
        end
        signal_connect(entry_frq2, "value-changed") do widget
            draw(can)
        end
        signal_connect(entry_ncyc, "value-changed") do widget
            draw(can)
        end
        signal_connect(entry_nt, "value-changed") do widget
            draw(can)
        end
        signal_connect(entry_wlen, "value-changed") do widget
            draw(can)
        end
        signal_connect(entry_woverlap, "value-changed") do widget
            draw(can)
        end
        signal_connect(cb_db, "toggled") do widget
            draw(can)
        end
        signal_connect(cb_mono, "toggled") do widget
            draw(can)
        end
        signal_connect(cb_hw, "toggled") do widget
            draw(can)
        end
        signal_connect(entry_gw, "value-changed") do widget
            draw(can)
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

        help = "Keyboard shortcuts:\n\nCtrl + ,\t\t\tGo back by 1 second\nCtrl + .\t\t\tGo forward by 1 second\nAlt + ,\t\t\tGo back by $(round(zoom)) seconds\nAlt + .\t\t\tGo forward by $(round(zoom)) seconds\n\n[\t\t\t\tZoom in\n]\t\t\t\tZoom out\n\nCtrl + s\t\t\tSave as PNG\nAlt + m\t\t\tToggle monochromatic mode\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tClose\n"

        signal_connect(bt_help, "clicked") do widget
            info_dialog(_nill, help, win)
        end

        win_key = Gtk4.GtkEventControllerKey(win)

        signal_connect(win_key, "key-released") do widget, keyval, keycode, state
            k = nothing
        end

        signal_connect(win_key, "key-pressed") do widget, keyval, keycode, state
            k = keyval
            if keyval == UInt('[')
                if zoom > 1
                    zoom -= 1
                    bt_next.tooltip_text = "Go forward by $(round(zoom)) seconds"
                    bt_prev.tooltip_text = "Go back by $(round(zoom)) seconds"

                    signal_adj = GtkAdjustment(signal_slider)
                    signal_adj.lower = obj.time_pts[1]
                    signal_adj.upper = obj.time_pts[end] - zoom
                    @idle_add Gtk4.adjustment(signal_slider, signal_adj)

                    time_adj = GtkAdjustment(entry_time)
                    time_adj.lower = obj.time_pts[1]
                    time_adj.upper = obj.time_pts[end] - zoom
                    @idle_add Gtk4.adjustment(entry_time, time_adj)
                    draw(can)
                end
                help = "Keyboard shortcuts:\n\nCtrl + ,\t\t\tGo back by 1 second\nCtrl + .\t\t\tGo forward by 1 second\nAlt + ,\t\t\tGo back by $(round(zoom)) seconds\nAlt + .\t\t\tGo forward by $(round(zoom)) seconds\n\n[\t\t\t\tZoom in\n]\t\t\t\tZoom out\n\nCtrl + s\t\t\tSave as PNG\nAlt + m\t\t\tToggle monochromatic mode\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tClose\n"
            elseif keyval == UInt(']')
                if zoom < 30 && zoom < obj.time_pts[end] - 1
                    zoom += 1
                    bt_next.tooltip_text = "Go forward by $(round(zoom)) seconds"
                    bt_prev.tooltip_text = "Go back by $(round(zoom)) seconds"

                    signal_adj = GtkAdjustment(signal_slider)
                    signal_adj.lower = obj.time_pts[1]
                    signal_adj.upper = obj.time_pts[end] - zoom
                    @idle_add Gtk4.adjustment(signal_slider, signal_adj)

                    time_adj = GtkAdjustment(entry_time)
                    time_adj.lower = obj.time_pts[1]
                    time_adj.upper = obj.time_pts[end] - zoom
                    @idle_add Gtk4.adjustment(entry_time, time_adj)
                    draw(can)
                else
                    zoom = obj.time_pts[end]
                    bt_next.tooltip_text = "Go forward by $(round(zoom)) seconds"
                    bt_prev.tooltip_text = "Go back by $(round(zoom)) seconds"

                    signal_adj = GtkAdjustment(signal_slider)
                    signal_adj.lower = obj.time_pts[1]
                    signal_adj.upper = obj.time_pts[end] - zoom
                    @idle_add Gtk4.adjustment(signal_slider, signal_adj)

                    time_adj = GtkAdjustment(entry_time)
                    time_adj.lower = obj.time_pts[1]
                    time_adj.upper = obj.time_pts[end] - zoom
                    @idle_add Gtk4.adjustment(entry_time, time_adj)
                    draw(can)
                end
                help = "Keyboard shortcuts:\n\nCtrl + ,\t\t\tGo back by 1 second\nCtrl + .\t\t\tGo forward by 1 second\nAlt + ,\t\t\tGo back by $(round(zoom)) seconds\nAlt + .\t\t\tGo forward by $(round(zoom)) seconds\n\n[\t\t\t\tZoom in\n]\t\t\t\tZoom out\n\nCtrl + s\t\t\tSave as PNG\nAlt + m\t\t\tToggle monochromatic mode\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tClose\n"
            end

            # ALT
            if ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_alt == mask_alt) && keyval == UInt(','))
                time_current = entry_time.value
                if time_current >= obj.time_pts[1] + zoom
                    time_current = time_current - zoom
                    @idle_add entry_time.value = time_current
                end
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_alt == mask_alt) && keyval == UInt('.'))
                time_current = entry_time.value
                if time_current < obj.time_pts[end] - zoom
                    time_current += zoom
                    @idle_add entry_time.value = time_current
                else
                    time_current = obj.time_pts[end] - zoom
                    @idle_add entry_time.value = time_current
                end
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_alt == mask_alt) && keyval == UInt('m'))
                mono = !mono
                cb_mono.active = mono
            end

            # CONTROL
            if ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('q'))
                close(win)
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('h'))
                info_dialog(_nill, help, win)
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('s'))
                save_dialog("Pick an image file", win, ["*.png"]) do file_name
                    if file_name != ""
                        surface_buf = Gtk4.cairo_surface(can)
                        if Cairo.write_to_png(surface_buf, file_name) == Cairo.STATUS_SUCCESS
                            _info("Plot saved as: $file_name")
                        else
                            warn_dialog(_nill, "File cannot be saved!", win)
                        end
                    end
                end
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt(','))
                time_current = entry_time.value
                if time_current >= obj.time_pts[1] + 1
                    time_current = time_current - 1
                    @idle_add entry_time.value = time_current
                end
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('.'))
                time_current = entry_time.value
                if time_current < obj.time_pts[end] - 1
                    time_current += 1
                    @idle_add entry_time.value = time_current
                else
                    time_current = obj.time_pts[end] - 1
                    @idle_add entry_time.value = time_current
                end
            end
        end
    end

    app = GtkApplication("org.neuroanalyzer.ipsd")
    Gtk4.signal_connect(_activate, app, :activate)
    Gtk4.GLib.stop_main_loop()
    Gtk4.run(app)

    return nothing

end

"""
    ipsd_ep(obj, ch)

Interactive PSD of epoched signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::String`: channel name

# Returns

Nothing
"""
function ipsd_ep(obj::NeuroAnalyzer.NEURO; ch::String)::Nothing

    @assert nepochs(obj) > 1 "For continuous object ipsd() must be used."

    ch_init = ch
    ch = get_channel(obj, ch=ch)
    clabels = labels(obj)

    p = NeuroAnalyzer.plot_psd(obj, ch=clabels[ch], ep=1)

    function _activate(app)

        win = GtkApplicationWindow(app, "NeuroAnalyzer: ipsd_ep()")
        win.width_request = 1200
        win.height_request = 850

        can = GtkCanvas(p.attr[:size][1], p.attr[:size][2])
        win_view = GtkScrolledWindow()
        win_view.min_content_width = 1200
        win_view.min_content_height = 800
        win_view.child = can

        g = GtkGrid()
        g.column_homogeneous = false
        g.column_spacing = 5
        g.row_spacing = 5
        g.margin_start = 5
        g.margin_end = 5
        g.margin_top = 5
        g.margin_bottom = 5

        g_opts = GtkGrid()
        g_opts.column_homogeneous = false
        g_opts.column_spacing = 5
        g_opts.row_spacing = 5
        g_opts.margin_start = 5
        g_opts.margin_end = 5
        g_opts.margin_top = 5
        g_opts.margin_bottom = 5

        entry_epoch = GtkSpinButton(1, nepochs(obj), 1)
        entry_epoch.tooltip_text = "Epoch"
        bt_start = GtkButton("⇤")
        bt_start.tooltip_text = "Go to the first epoch"
        bt_end = GtkButton("⇥")
        bt_end.tooltip_text = "Go to the last epoch"
        bt_help = GtkButton("Help")
        bt_help.tooltip_text = "Show help"
        bt_close = GtkButton("Close")
        bt_close.tooltip_text = "Close this window"

        entry_title = GtkEntry()
        entry_title.text = "default"
        entry_title.tooltip_text = "Plot title"
        entry_xlab = GtkEntry()
        entry_xlab.text = "default"
        entry_xlab.tooltip_text = "X axis label"
        entry_ylab = GtkEntry()
        entry_ylab.text = "default"
        entry_ylab.tooltip_text = "Y axis label"

        combo_ch = GtkComboBoxText()
        ctypes = uppercase.(unique(obj.header.recording[:channel_type]))
        ch_types = [ctypes; clabels]
        for idx in ch_types
            push!(combo_ch, idx)
        end
        if ch_init in lowercase.(ctypes)
            n = findfirst(isequal(ch_init), lowercase.(ctypes))
            combo_ch.active = n - 1
        else
            combo_ch.active = (length(ctypes) + ch[1] - 1)
        end
        combo_ch.tooltip_text = "Channels"

        cb_mono = GtkCheckButton()
        cb_mono.tooltip_text = "Use color or gray palette"

        cb_db = GtkCheckButton()
        cb_db.tooltip_text = "Normalize powers to dB"
        cb_db.active = true

        cb_hw = GtkCheckButton()
        cb_hw.tooltip_text = "Apply Hanning window"
        cb_hw.active = true

        combo_method = GtkComboBoxText()
        psd_methods = ["Welch's periodogram", "fast Fourier transform", "short-time Fourier transform", "multi-taper", "Morlet wavelet", "Gaussian-Hilbert transform"]
        for idx in psd_methods
            push!(combo_method, idx)
        end
        combo_method.active = 0
        combo_method.tooltip_text = "PSD method"

        combo_type = GtkComboBoxText()
        psd_types = ["normal", "butterfly", "mean", "w3d", "s3d", "topo"]
        for idx in psd_types
            push!(combo_type, idx)
        end
        combo_type.active = 0
        combo_type.tooltip_text = "Plot type"
        if length(ch) > 1
            combo_type.sensitive = true
        else
            combo_type.sensitive = false
        end

        combo_ref = GtkComboBoxText()
        ref_types = ["absolute", "total power", "delta", "theta", "alpha", "alpha lower", "alpha higher", "beta", "beta lower", "beta higher", "gamma", "gamma 1", "gamma 2", "gamma lower", "gamma higher"]
        for idx in ref_types
            push!(combo_ref, idx)
        end
        combo_ref.active = 0
        combo_ref.tooltip_text = "Powers referenced to"

        cb_frq = GtkCheckButton()
        cb_frq.tooltip_text = "Linear or logarithmic frequencies"
        cb_frq.active = true

        entry_nt = GtkSpinButton(1, 128, 1)
        entry_nt.value = 8
        entry_nt.tooltip_text = "Number of Slepian tapers (for MT)"

        entry_ncyc = GtkSpinButton(1, 256, 1)
        entry_ncyc.value = 32
        entry_ncyc.tooltip_text = "Number of Morlet wavelet cycles"

        entry_frq1 = GtkSpinButton(0.0, (sr(obj) / 2) - 0.5, 0.5)
        entry_frq1.value = 0.0
        entry_frq1.tooltip_text = "Start frequency"

        entry_frq2 = GtkSpinButton(0.0, sr(obj) / 2, 0.5)
        entry_frq2.value = sr(obj) / 2
        entry_frq2.tooltip_text = "End frequency"

        entry_wlen = GtkSpinButton(2, epoch_len(obj) * sr(obj) + 1, 1)
        entry_wlen.value = sr(obj)
        entry_wlen.tooltip_text = "Window length (samples)"

        entry_woverlap = GtkSpinButton(0, epoch_len(obj) * sr(obj), 1)
        entry_woverlap.value = round(Int64, sr(obj) * 0.97)
        entry_woverlap.tooltip_text = "Window overlap (samples)"

        entry_gw = GtkSpinButton(1, 128, 1)
        entry_gw.value = 6
        entry_gw.tooltip_text = "Gaussian width in Hz"

        bt_refresh = GtkButton("Refresh")
        bt_refresh.tooltip_text = "Refresh the plot"

        lab_type = GtkLabel("Plot type:")
        lab_type.halign = 2
        lab_method = GtkLabel("PSD method:")
        lab_method.halign = 2
        lab_ch = GtkLabel("Channels:")
        lab_ch.halign = 2
        lab_t = GtkLabel("Title:")
        lab_t.halign = 2
        lab_x = GtkLabel("X lab:")
        lab_x.halign = 2
        lab_y = GtkLabel("Y lab:")
        lab_y.halign = 2
        lab_mono = GtkLabel("Grayscale:")
        lab_mono.halign = 2
        lab_norm = GtkLabel("Normalize:")
        lab_norm.halign = 2
        lab_nt = GtkLabel("Slepians:")
        lab_nt.halign = 2
        lab_wlen = GtkLabel("Window length:")
        lab_wlen.halign = 2
        lab_woverlap = GtkLabel("Window overlap:")
        lab_woverlap.halign = 2
        lab_frq1 = GtkLabel("Start frequency:")
        lab_frq1.halign = 2
        lab_frq2 = GtkLabel("End frequency:")
        lab_frq2.halign = 2
        lab_nc = GtkLabel("Cycles:")
        lab_nc.halign = 2
        lab_ref = GtkLabel("PSD reference:")
        lab_ref.halign = 2
        lab_frq = GtkLabel("Linear frequencies:")
        lab_frq.halign = 2
        lab_hw = GtkLabel("Hanning:")
        lab_hw.halign = 2
        lab_gw = GtkLabel("Gaussian width:")
        lab_gw.halign = 2

        signal_slider = GtkScale(:h, 1:nepochs(obj))
        signal_slider.draw_value = false
        signal_slider.tooltip_text = "Current epoch"

        g_opts[1, 1] = lab_ch
        g_opts[1, 2] = lab_type
        g_opts[1, 3] = lab_method
        g_opts[1, 4] = lab_ref
        g_opts[1, 5] = lab_t
        g_opts[1, 6] = lab_x
        g_opts[1, 7] = lab_y
        g_opts[1, 8] = lab_frq1
        g_opts[1, 9] = lab_frq2
        g_opts[1, 10] = lab_nc
        g_opts[1, 11] = lab_nt
        g_opts[1, 12] = lab_wlen
        g_opts[1, 13] = lab_woverlap
        g_opts[1, 14] = lab_gw
        g_opts[1, 15] = lab_frq
        g_opts[1, 16] = lab_norm
        g_opts[1, 17] = lab_mono
        g_opts[1, 18] = lab_hw
        g_opts[2, 1] = combo_ch
        g_opts[2, 2] = combo_type
        g_opts[2, 3] = combo_method
        g_opts[2, 4] = combo_ref
        g_opts[2, 5] = entry_title
        g_opts[2, 6] = entry_xlab
        g_opts[2, 7] = entry_ylab
        g_opts[2, 8] = entry_frq1
        g_opts[2, 9] = entry_frq2
        g_opts[2, 10] = entry_ncyc
        g_opts[2, 11] = entry_nt
        g_opts[2, 12] = entry_wlen
        g_opts[2, 13] = entry_woverlap
        g_opts[2, 14] = entry_gw
        g_opts[2, 15] = cb_frq
        g_opts[2, 16] = cb_db
        g_opts[2, 17] = cb_mono
        g_opts[2, 18] = cb_hw
        g_opts[1:2, 19] = bt_refresh
        vbox = GtkBox(:v)
        push!(vbox, g_opts)

        g[1, 1] = vbox
        g[2:7, 1] = win_view
        g[2:7, 2] = signal_slider
        g[2, 3] = bt_start
        g[3, 3] = entry_epoch
        g[4, 3] = bt_end
        g[5, 3] = GtkLabel("")
        g[6, 3] = bt_help
        g[7, 3] = bt_close
        push!(win, g)

        Gtk4.show(win)

        @guarded draw(can) do widget
            ch = ch_types[Int64(combo_ch.active) + 1]
            ch in ctypes && (ch = lowercase(ch))
            title = entry_title.text
            xlab = entry_xlab.text
            ylab = entry_ylab.text
            mono = cb_mono.active
            db = cb_db.active
            hw = cb_hw.active
            method = combo_method.active
            method == 0 && (method = :welch)
            method == 1 && (method = :fft)
            method == 2 && (method = :stft)
            method == 3 && (method = :mt)
            method == 4 && (method = :mw)
            method == 5 && (method = :gh)
            type = combo_type.active
            type == 0 && (type = :normal)
            type == 1 && (type = :butterfly)
            type == 2 && (type = :mean)
            type == 3 && (type = :w3d)
            type == 4 && (type = :s3d)
            type == 5 && (type = :topo)
            ref = combo_ref.active
            ref == 0 && (ref = :abs)
            ref == 1 && (ref = :total)
            ref == 2 && (ref = :delta)
            ref == 3 && (ref = :theta)
            ref == 4 && (ref = :alpha)
            ref == 5 && (ref = :alpha_lower)
            ref == 6 && (ref = :alpha_higher)
            ref == 7 && (ref = :beta)
            ref == 8 && (ref = :beta_lower)
            ref == 9 && (ref = :beta_higher)
            ref == 10 && (ref = :gamma)
            ref == 11 && (ref = :gamma_1)
            ref == 12 && (ref = :gamma_2)
            ref == 13 && (ref = :gamma_lower)
            ref == 14 && (ref = :gamma_higher)
            frq = cb_frq.active ? :lin : :log
            gw = entry_gw.value
            frq1 = entry_frq1.value
            frq2 = entry_frq2.value
            nt = Int64(entry_nt.value)
            ncyc = Int64(entry_ncyc.value)
            wlen = Int64(entry_wlen.value)
            woverlap = Int64(entry_woverlap.value)

            no_error = true
            if frq1 == frq2
                warn_dialog(_nill, "Start and end frequencies must be different.", win)
                no_error = false
            elseif frq1 > frq2
                warn_dialog(_nill, "Start frequency must be < end frequency.", win)
                no_error = false
            elseif woverlap >= wlen
                warn_dialog(_nill, "Window overlap must be < window length.", win)
                no_error = false
            elseif length(get_channel(obj, ch=ch)) < 2 && type === :butterfly
                warn_dialog(_nill, "For butterfly plot, the signal must contain ≥ 2 channels.", win)
                no_error = false
            elseif length(get_channel(obj, ch=ch)) < 2 && type === :mean
                warn_dialog(_nill, "For mean plot, the signal must contain ≥ 2 channels.", win)
                no_error = false
            elseif length(get_channel(obj, ch=ch)) < 2 && type === :w3d
                warn_dialog(_nill, "For w3d plot, the signal must contain ≥ 2 channels.", win)
                no_error = false
            elseif length(get_channel(obj, ch=ch)) < 2 && type === :s3d
                warn_dialog(_nill, "For s3d plot, the signal must contain ≥ 2 channels.", win)
                no_error = false
            elseif DataFrames.nrow(obj.locs) == 0 && type === :topo
                warn_dialog(_nill, "Electrode locations not available.", win)
                no_error = false
            elseif length(unique(obj.header.recording[:channel_type][get_channel(obj, ch=ch)])) > 1 && (type in [:butterfly, :mean, :w3d, :s3d, :topo] || ch == "all")
                warn_dialog(_nill, "For multi-channel $(string(type)) plot all channels must be of the same type.", win)
                no_error = false
            end

            if no_error
                ep = Int64(entry_epoch.value)
                p = NeuroAnalyzer.plot_psd(obj,
                                           ch=ch,
                                           ep=ep,
                                           mono=mono,
                                           title=title,
                                           xlabel=xlab,
                                           ylabel=ylab,
                                           db=db,
                                           method=method,
                                           type=type,
                                           frq=frq,
                                           ref=ref,
                                           frq_lim=(frq1, frq2),
                                           ncyc=ncyc,
                                           nt=nt,
                                           wlen=wlen,
                                           woverlap=woverlap,
                                           w=hw,
                                           gw=gw)
                img = read_from_png(io)
                set_gtk_property!(can, :width_request, p.attr[:size][1])
                set_gtk_property!(can, :height_request, p.attr[:size][2])
                ctx = getgc(can)
                show(io, MIME("image/png"), p)
                img = read_from_png(io)
                set_source_surface(ctx, img, 0, 0)
                paint(ctx)
            end
        end

        function _mwheel_scroll(_, dx, dy)
            if dy > 0.5
                if k == 0x0000ffe1 # shift
                    ep = Int64(entry_epoch.value)
                    if ep < nepochs(obj1)
                        ep += 1
                        @idle_add entry_epoch.value = ep
                    end
                end
            elseif dy < 0.5
                if k == 0x0000ffe1 # shift
                    ep = Int64(entry_epoch.value)
                    if ep > 1
                        ep -= 1
                        @idle_add entry_epoch.value = ep
                    end
                end
            end
        end
        ecsf = Gtk4.GtkEventControllerScroll(Gtk4.EventControllerScrollFlags_VERTICAL)
        push!(can, ecsf)
        signal_connect(_mwheel_scroll, ecsf, "scroll")

        signal_connect(combo_ch, "changed") do widget
            ch = Int64(combo_ch.active) + 1
            if ch in 1:length(ctypes)
                ch = lowercase(ctypes[ch])
                if length(get_channel(obj, type=ch)) > 1
                    combo_type.sensitive = true
                else
                    combo_type.active = 0
                    combo_type.sensitive = false
                end
            else
                combo_type.active = 0
                combo_type.sensitive = false
            end
            draw(can)
        end

        signal_connect(bt_refresh, "clicked") do widget
            draw(can)
        end
        signal_connect(combo_type, "changed") do widget
            draw(can)
        end
        signal_connect(combo_method, "changed") do widget
            draw(can)
        end
        signal_connect(combo_ref, "changed") do widget
            draw(can)
        end
        signal_connect(cb_frq, "toggled") do widget
            draw(can)
        end
        signal_connect(entry_frq1, "value-changed") do widget
            draw(can)
        end
        signal_connect(entry_frq2, "value-changed") do widget
            draw(can)
        end
        signal_connect(entry_ncyc, "value-changed") do widget
            draw(can)
        end
        signal_connect(entry_nt, "value-changed") do widget
            draw(can)
        end
        signal_connect(entry_wlen, "value-changed") do widget
            draw(can)
        end
        signal_connect(entry_woverlap, "value-changed") do widget
            draw(can)
        end
        signal_connect(cb_db, "toggled") do widget
            draw(can)
        end
        signal_connect(cb_mono, "toggled") do widget
            draw(can)
        end
        signal_connect(cb_hw, "toggled") do widget
            draw(can)
        end
        signal_connect(entry_gw, "value-changed") do widget
            draw(can)
        end

        signal_connect(signal_slider, "value-changed") do widget
            @idle_add entry_epoch.value = round(Int64, Gtk4.value(signal_slider))
        end
        signal_connect(entry_epoch, "value-changed") do widget
            @idle_add Gtk4.value(signal_slider, entry_epoch.value)
            draw(can)
        end

        signal_connect(bt_start, "clicked") do widget
            @idle_add Gtk4.value(entry_epoch, 1)
        end

        signal_connect(bt_end, "clicked") do widget
            @idle_add Gtk4.value(entry_epoch, nepochs(obj))
        end

        signal_connect(bt_close, "clicked") do widget
            close(win)
        end

        help = "Keyboard shortcuts:\n\nCtrl + ,\t\t\tPrevious epoch\nCtrl + .\t\t\tNext epoch\n\nCtrl + s\t\t\tSave as PNG\nAlt + m\t\t\tToggle monochromatic mode\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tClose\n"

        signal_connect(bt_help, "clicked") do widget
            info_dialog(_nill, help, win)
        end

        win_key = Gtk4.GtkEventControllerKey(win)

        signal_connect(win_key, "key-released") do widget, keyval, keycode, state
            k = nothing
        end

        signal_connect(win_key, "key-pressed") do widget, keyval, keycode, state
            # ALT
            if ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_alt == mask_alt) && keyval == UInt('m'))
                mono = !mono
                cb_mono.active = mono
            end
            # CONTROL
            if ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('q'))
                close(win)
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('h'))
                info_dialog(_nill, help, win)
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('s'))
                save_dialog("Pick an image file", win, ["*.png"]) do file_name
                    if file_name != ""
                        surface_buf = Gtk4.cairo_surface(can)
                        if Cairo.write_to_png(surface_buf, file_name) == Cairo.STATUS_SUCCESS
                            _info("Plot saved as: $file_name")
                        else
                            warn_dialog(_nill, "File cannot be saved!", win)
                        end
                    end
                end
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt(','))
                ep = Int64(entry_epoch.value)
                if ep > 1
                    ep -= 1
                    @idle_add entry_epoch.value = ep
                end
                draw(can)
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('.'))
                ep = Int64(entry_epoch.value)
                if ep > 1
                    ep -= 1
                    @idle_add entry_epoch.value = ep
                end
            end
        end
    end

    app = GtkApplication("org.neuroanalyzer.ipsd_ep")
    Gtk4.signal_connect(_activate, app, :activate)
    Gtk4.GLib.stop_main_loop()
    Gtk4.run(app)

    return nothing

end
