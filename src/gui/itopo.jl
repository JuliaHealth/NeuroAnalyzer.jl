export itopo
export itopo_ep

"""
    itopo(obj; <keyword arguments>)

Interactive topographical map of continuous signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channels to plot

# Returns

Nothing
"""
function itopo(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::Nothing

    @assert nepochs(obj) == 1 "For epoched object itopo_ep() must be used."

    _check_datatype(obj, ["eeg", "meg", "erp"])

    p = NeuroAnalyzer.plot_topo(obj, ch=ch)

    function _activate(app)

        if p.attr[:size][1] > 900
            win = GtkApplicationWindow(app, "NeuroAnalyzer: itopo()")
            win.content_width = p.attr[:size][1] + 100
            win.content_height = 1000
            can = GtkCanvas(1000, 1000)
        else
            win = GtkApplicationWindow(app, "NeuroAnalyzer: itopo()")
            win.content_width = p.attr[:size][1] + 100
            win.content_height = 800
            can = GtkCanvas(800, 800)
        end

        g_opts = GtkGrid()
        g_opts.column_homogeneous = true
        g_opts.row_spacing = 5
        g_opts.column_spacing = 5

        entry_ts1 = GtkSpinButton(obj.time_pts[1], obj.time_pts[end], 0.01)
        entry_ts1.value = obj.time_pts[1]
        entry_ts1.digits = 3
        entry_ts1.tooltip_text = "Segment start [s]"

        entry_ts2 = GtkSpinButton(obj.time_pts[1], obj.time_pts[end], 0.01)
        entry_ts2.value = obj.time_pts[1]
        entry_ts2.digits = 3
        entry_ts2.tooltip_text = "Segment end [s]"

        bt_help = GtkButton("Help")
        bt_help.tooltip_text = "Show help"
        bt_close = GtkButton("Close")
        bt_close.tooltip_text = "Close this window"

        entry_title = GtkEntry()
        entry_title.text = "default"
        entry_title.tooltip_text = "Plot title"

        entry_cblab = GtkEntry()
        entry_cblab.text = "default"
        entry_cblab.tooltip_text = "Color bar label"

        cb_cb = GtkCheckButton()
        cb_cb.tooltip_text = "Draw color bar"
        cb_cb.active = true

        cb_cart = GtkCheckButton()
        cb_cart.tooltip_text = "Use Cartesian coordinates of electrodes"
        cb_cart.active = false

        cb_large = GtkCheckButton()
        cb_large.tooltip_text = "Draw large plot"
        cb_large.active = true

        cb_elec = GtkCheckButton()
        cb_elec.tooltip_text = "Draw electrodes"
        cb_elec.active = true

        cb_contour = GtkCheckButton()
        cb_contour.tooltip_text = "Draw contour lines"
        cb_contour.active = true

        combo_imethod = GtkComboBoxText()
        imethod_types = ["shepard", "multiquadratic", "inv multiquadratic", "thin plate", "nearest neighbour", "gaussian"]
        for idx in imethod_types
            push!(combo_imethod, idx)
        end
        combo_imethod.active = 0
        combo_imethod.tooltip_text = "Interpolate type"

        combo_amethod = GtkComboBoxText()
        amethod_types = ["mean", "median"]
        for idx in amethod_types
            push!(combo_amethod, idx)
        end
        combo_amethod.active = 0
        combo_amethod.tooltip_text = "Averaging method"
        combo_amethod.sensitive = false

        combo_nmethod = GtkComboBoxText()
        nmethod_types = ["zscore", "gauss", "invroot", "log", "minmax", "neg", "neglog", "neglog10", "perc", "pos", "softmax", "none"]
        for idx in nmethod_types
            push!(combo_nmethod, idx)
        end
        combo_nmethod.active = 0
        combo_nmethod.tooltip_text = "Normalization method"

        bt_refresh = GtkButton("Refresh")
        bt_refresh.tooltip_text = "Refresh the plot"

        lab_type = GtkLabel("Interpolate:")
        lab_type.halign = 2
        lab_amethod = GtkLabel("Averaging:")
        lab_amethod.halign = 2
        lab_nmethod = GtkLabel("Normalization:")
        lab_nmethod.halign = 2
        lab_t = GtkLabel("Title:")
        lab_t.halign = 2
        lab_cb = GtkLabel("Color bar title:")
        lab_cb.halign = 2
        lab_cb_draw = GtkLabel("Draw color bar:")
        lab_cb_draw.halign = 2
        lab_cart = GtkLabel("Cartesian:")
        lab_cart.halign = 2
        lab_large = GtkLabel("Large plot:")
        lab_large.halign = 2
        lab_elec = GtkLabel("Draw electrodes:")
        lab_elec.halign = 2
        lab_contour = GtkLabel("Draw contours:")
        lab_contour.halign = 2
        lab_ts = GtkLabel("Time segment:")
        g_opts[1, 2] = entry_ts1
        g_opts[2, 2] = entry_ts2
        g_opts[1, 3] = lab_t
        g_opts[1, 4] = lab_type
        g_opts[1, 5] = lab_amethod
        g_opts[1, 6] = lab_nmethod
        g_opts[1, 7] = lab_cb
        g_opts[1, 8] = lab_cb_draw
        g_opts[1, 9] = lab_cart
        g_opts[1, 10] = lab_large
        g_opts[1, 11] = lab_elec
        g_opts[1, 12] = lab_contour
        g_opts[1:2, 1] = lab_ts
        g_opts[2, 3] = entry_title
        g_opts[2, 4] = combo_imethod
        g_opts[2, 5] = combo_amethod
        g_opts[2, 6] = combo_nmethod
        g_opts[2, 7] = entry_cblab
        g_opts[2, 8] = cb_cb
        g_opts[2, 9] = cb_cart
        g_opts[2, 10] = cb_large
        g_opts[2, 11] = cb_elec
        g_opts[2, 12] = cb_contour
        g_opts[1:2, 14] = bt_refresh
        g_opts[1, 15] = bt_help
        g_opts[2, 15] = bt_close

        g = GtkGrid()
        g.column_homogeneous = false
        g.row_spacing = 5
        g.column_spacing = 5
        g.margin_start = 5
        g.margin_end = 5
        g.margin_top = 5
        g.margin_bottom = 5

        vbox = GtkBox(:v)
        push!(vbox, g_opts)
        g[1, 1] = vbox
        g[2:11, 1] = can
        push!(win, g)

        Gtk4.show(win)

        @guarded draw(can) do widget
            seg = round.((entry_ts1.value, entry_ts2.value), digits=3)
            title = entry_title.text
            cblab = entry_cblab.text
            cb = cb_cb.active
            imethod = combo_imethod.active
            imethod == 0 && (imethod = :sh)
            imethod == 1 && (imethod = :mq)
            imethod == 2 && (imethod = :imq)
            imethod == 3 && (imethod = :tp)
            imethod == 4 && (imethod = :nn)
            imethod == 5 && (imethod = :ga)
            amethod = combo_amethod.active
            amethod == 0 && (amethod = :mean)
            amethod == 1 && (amethod = :median)
            nmethod = combo_nmethod.active
            nmethod == 0 && (nmethod = :zscore)
            nmethod == 1 && (nmethod = :gauss)
            nmethod == 2 && (nmethod = :invroot)
            nmethod == 3 && (nmethod = :log)
            nmethod == 4 && (nmethod = :minmax)
            nmethod == 5 && (nmethod = :neg)
            nmethod == 6 && (nmethod = :neglog)
            nmethod == 7 && (nmethod = :neglog10)
            nmethod == 8 && (nmethod = :perc)
            nmethod == 9 && (nmethod = :pos)
            nmethod == 10 && (nmethod = :softmax)
            nmethod == 11 && (nmethod = :none)
            plot_contours = cb_contour.active
            plot_electrodes = cb_elec.active
            cart = cb_cart.active
            large = cb_large.active
            p = NeuroAnalyzer.plot_topo(obj,
                                        ch=ch,
                                        seg=seg,
                                        title=title,
                                        cb=cb,
                                        cb_label=cblab,
                                        amethod=amethod,
                                        imethod=imethod,
                                        nmethod=nmethod,
                                        large=large,
                                        plot_contours=plot_contours,
                                        plot_electrodes=plot_electrodes,
                                        cart=cart)
            ctx = getgc(can)
            if p.attr[:size][1] > 900
                Gtk4.rectangle(ctx, 0, 0, 999, 999)
                Cairo.set_source_rgb(ctx, 255, 255, 255)
                Gtk4.fill(ctx)
                withenv("GKSwstype" => "100") do
                    png(p, io)
                end
                img = read_from_png(io)
                set_source_surface(ctx, img, 500 - (p.attr[:size][1] ÷ 2) - 1, 500 - (p.attr[:size][1] ÷ 2) - 1)
            else
                Gtk4.rectangle(ctx, 0, 0, 799, 799)
                Cairo.set_source_rgb(ctx, 255, 255, 255)
                Gtk4.fill(ctx)
                withenv("GKSwstype" => "100") do
                    png(p, io)
                end
                img = read_from_png(io)
                set_source_surface(ctx, img, 400 - (p.attr[:size][1] ÷ 2) - 1, 400 - (p.attr[:size][1] ÷ 2) - 1)
            end
            paint(ctx)
        end

        signal_connect(entry_ts1, "value-changed") do widget
            seg = round.((entry_ts1.value, entry_ts2.value), digits=3)
            if seg[1] > seg[2]
                warn_dialog(_nill, "Cannot plot!\nSegment start is larger than segment end.", win)
            else
                draw(can)
            end
            seg[1] == seg[2] && (combo_amethod.sensitive = false)
            seg[1] < seg[2] && (combo_amethod.sensitive = true)
        end

        signal_connect(entry_ts2, "value-changed") do widget
            seg = round.((entry_ts1.value, entry_ts2.value), digits=3)
            if seg[1] > seg[2]
                warn_dialog(_nill, "Cannot plot!\nSegment start is larger than segment end.", win)
            else
                draw(can)
            end
            seg[1] == seg[2] && (combo_amethod.sensitive = false)
            seg[1] < seg[2] && (combo_amethod.sensitive = true)
        end

        signal_connect(bt_refresh, "clicked") do widget
            draw(can)
        end
        signal_connect(combo_imethod, "changed") do widget
            draw(can)
        end
        signal_connect(combo_amethod, "changed") do widget
            draw(can)
        end
        signal_connect(combo_nmethod, "changed") do widget
            draw(can)
        end
        signal_connect(cb_cb, "toggled") do widget
            draw(can)
        end
        signal_connect(cb_cart, "toggled") do widget
            draw(can)
        end
        signal_connect(cb_large, "toggled") do widget
            draw(can)
        end
        signal_connect(cb_elec, "toggled") do widget
            draw(can)
        end
        signal_connect(cb_contour, "toggled") do widget
            draw(can)
        end
        signal_connect(bt_close, "clicked") do widget
            close(win)
        end

        help = "Keyboard shortcuts:\n\nCtrl + s\t\t\tSave as PNG\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tExit\n"

        signal_connect(bt_help, "clicked") do widget
            info_dialog(_nill, help, win)
        end

        win_key = Gtk4.GtkEventControllerKey(win)

        signal_connect(win_key, "key-pressed") do widget, keyval, keycode, state
            if ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('s'))
                save_dialog("Pick an image file", win, ["*.png"]) do file_name
                    if file_name != ""
                        surface_buf = Gtk4.cairo_surface(can)
                        if Cairo.write_to_png(surface_buf, file_name) == Cairo.STATUS_SUCCESS
                            _info("Plot saved as: $file_name")
                        else
                            warn_dialog(_nill, "File $file_name cannot be written!", win)
                        end
                    end
                end
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('h'))
                info_dialog(_nill, help, win)
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('q'))
                close(win)
            end
        end
    end

    app = GtkApplication("org.neuroanalyzer.itopo")
    Gtk4.signal_connect(_activate, app, :activate)
    Gtk4.GLib.stop_main_loop()
    Gtk4.run(app)

    return nothing

end

"""
    itopo_ep(obj; <keyword arguments>)

Interactive topographical map of epoched signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channels to plot

# Returns

Nothing
"""
function itopo_ep(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::Nothing

    @assert nepochs(obj) > 1 "For continuous object itopo() must be used."

    _check_datatype(obj, ["eeg", "meg", "erp"])

    p = NeuroAnalyzer.plot_topo(obj, ch=ch)

    function _activate(app)

        if p.attr[:size][1] > 900
            win = GtkApplicationWindow(app, "NeuroAnalyzer: itopo_ep()")
            Gtk4.default_size(win, p.attr[:size][1] + 100, 1000)
            can = GtkCanvas()
            can.content_width = 1000
            can.content_height = 1000
        else
            win = GtkApplicationWindow(app, "NeuroAnalyzer: itopo_ep()")
        Gtk4.default_size(win, p.attr[:size][1] + 100, 800)
            can = GtkCanvas()
            can.content_width = 800
            can.content_height = 800
        end

        g_opts = GtkGrid()
        g_opts.column_homogeneous = true
        g_opts.row_spacing = 5
        g_opts.column_spacing = 5

        entry_epoch = GtkSpinButton(1, nepochs(obj), 1)
        entry_epoch.tooltip_text = "Epoch"

        entry_ts1 = GtkSpinButton(obj.time_pts[1], obj.time_pts[end], 0.01)
        entry_ts1.value = obj.time_pts[1]
        entry_ts1.digits = 3
        entry_ts1.tooltip_text = "Segment start [s]"

        entry_ts2 = GtkSpinButton(obj.time_pts[1], obj.time_pts[end], 0.01)
        entry_ts2.value = obj.time_pts[1]
        entry_ts2.digits = 3
        entry_ts2.tooltip_text = "Segment end [s]"

        bt_help = GtkButton("Help")
        bt_help.tooltip_text = "Show help"
        bt_close = GtkButton("Close")
        bt_close.tooltip_text = "Close this window"

        entry_title = GtkEntry()
        entry_title.text = "default"
        entry_title.tooltip_text = "Plot title"

        entry_cblab = GtkEntry()
        entry_cblab.text = "default"
        entry_cblab.tooltip_text = "Color bar label"

        cb_cb = GtkCheckButton()
        cb_cb.tooltip_text = "Draw color bar"
        cb_cb.active = true

        cb_cart = GtkCheckButton()
        cb_cart.tooltip_text = "Use Cartesian coordinates of electrodes"
        cb_cart.active = false

        cb_large = GtkCheckButton()
        cb_large.tooltip_text = "Draw large plot"
        cb_large.active = true

        cb_elec = GtkCheckButton()
        cb_elec.tooltip_text = "Draw electrodes"
        cb_elec.active = true

        cb_contour = GtkCheckButton()
        cb_contour.tooltip_text = "Draw contour lines"
        cb_contour.active = true

        combo_imethod = GtkComboBoxText()
        imethod_types = ["shepard", "multiquadratic", "inv multiquadratic", "thin plate", "nearest neighbour", "gaussian"]
        for idx in imethod_types
            push!(combo_imethod, idx)
        end
        combo_imethod.active = 0
        combo_imethod.tooltip_text = "Interpolate type"

        combo_amethod = GtkComboBoxText()
        amethod_types = ["mean", "median"]
        for idx in amethod_types
            push!(combo_amethod, idx)
        end
        combo_amethod.active = 0
        combo_amethod.tooltip_text = "Averaging method"
        combo_amethod.sensitive = false

        combo_nmethod = GtkComboBoxText()
        nmethod_types = ["zscore", "gauss", "invroot", "log", "minmax", "neg", "neglog", "neglog10", "perc", "pos", "softmax", "none"]
        for idx in nmethod_types
            push!(combo_nmethod, idx)
        end
        combo_nmethod.active = 0
        combo_nmethod.tooltip_text = "Normalization method"

        bt_refresh = GtkButton("Refresh")
        bt_refresh.tooltip_text = "Refresh the plot"

        lab_type = GtkLabel("Interpolate:")
        lab_type.halign = 2
        lab_amethod = GtkLabel("Averaging:")
        lab_amethod.halign = 2
        lab_nmethod = GtkLabel("Normalization:")
        lab_nmethod.halign = 2
        lab_t = GtkLabel("Title:")
        lab_t.halign = 2
        lab_cb = GtkLabel("Color bar title:")
        lab_cb.halign = 2
        lab_cb_draw = GtkLabel("Draw color bar:")
        lab_cb_draw.halign = 2
        lab_cart = GtkLabel("Cartesian:")
        lab_cart.halign = 2
        lab_large = GtkLabel("Large plot:")
        lab_large.halign = 2
        lab_elec = GtkLabel("Draw electrodes:")
        lab_elec.halign = 2
        lab_contour = GtkLabel("Draw contours:")
        lab_contour.halign = 2
        lab_ts = GtkLabel("Time segment:")
        lab_epoch = GtkLabel("Epoch:")
        lab_epoch.halign = 2
        g_opts[1:2, 2] = lab_ts
        g_opts[1, 1] = lab_epoch
        g_opts[1, 4] = lab_t
        g_opts[1, 5] = lab_type
        g_opts[1, 6] = lab_amethod
        g_opts[1, 7] = lab_nmethod
        g_opts[1, 8] = lab_cb
        g_opts[1, 9] = lab_cb_draw
        g_opts[1, 10] = lab_cart
        g_opts[1, 11] = lab_large
        g_opts[1, 12] = lab_elec
        g_opts[1, 13] = lab_contour
        g_opts[1, 3] = entry_ts1
        g_opts[2, 3] = entry_ts2
        g_opts[2, 1] = entry_epoch
        g_opts[2, 4] = entry_title
        g_opts[2, 5] = combo_imethod
        g_opts[2, 6] = combo_amethod
        g_opts[2, 7] = combo_nmethod
        g_opts[2, 8] = entry_cblab
        g_opts[2, 9] = cb_cb
        g_opts[2, 10] = cb_cart
        g_opts[2, 11] = cb_large
        g_opts[2, 12] = cb_elec
        g_opts[2, 13] = cb_contour
        g_opts[1:2, 14] = bt_refresh
        g_opts[1, 15] = bt_help
        g_opts[2, 15] = bt_close

        g = GtkGrid()
        g.column_homogeneous = false
        g.row_spacing = 5
        g.column_spacing = 5
        g.margin_start = 5
        g.margin_end = 5
        g.margin_top = 5
        g.margin_bottom = 5

        vbox = GtkBox(:v)
        push!(vbox, g_opts)
        g[1, 1] = vbox
        g[2:11, 1] = can
        push!(win, g)

        Gtk4.show(win)

        @guarded draw(can) do widget
            ep = Int64(entry_epoch.value)
            ts1 = entry_ts1.value
            ts2 = entry_ts2.value
            ts1 = (epoch_len(obj) / sr(obj) * (ep - 1)) + obj.epoch_time[1] + ts1
            ts2 = (epoch_len(obj) / sr(obj) * (ep - 1)) + obj.epoch_time[1] + ts2
            seg = round.((ts1, ts2), digits=3)
            title = entry_title.text
            cblab = entry_cblab.text
            cb = cb_cb.active
            imethod = combo_imethod.active
            imethod == 0 && (imethod = :sh)
            imethod == 1 && (imethod = :mq)
            imethod == 2 && (imethod = :imq)
            imethod == 3 && (imethod = :tp)
            imethod == 4 && (imethod = :nn)
            imethod == 5 && (imethod = :ga)
            amethod = combo_amethod.active
            amethod == 0 && (amethod = :mean)
            amethod == 1 && (amethod = :median)
            nmethod = combo_nmethod.active
            nmethod == 0 && (nmethod = :zscore)
            nmethod == 1 && (nmethod = :gauss)
            nmethod == 2 && (nmethod = :invroot)
            nmethod == 3 && (nmethod = :log)
            nmethod == 4 && (nmethod = :minmax)
            nmethod == 5 && (nmethod = :neg)
            nmethod == 6 && (nmethod = :neglog)
            nmethod == 7 && (nmethod = :neglog10)
            nmethod == 8 && (nmethod = :perc)
            nmethod == 9 && (nmethod = :pos)
            nmethod == 10 && (nmethod = :softmax)
            nmethod == 11 && (nmethod = :none)
            plot_contours = cb_contour.active
            plot_electrodes = cb_elec.active
            cart = cb_cart.active
            large = cb_large.active
            p = NeuroAnalyzer.plot_topo(obj,
                                        ch=ch,
                                        seg=seg,
                                        title=title,
                                        cb=cb,
                                        cb_label=cblab,
                                        amethod=amethod,
                                        imethod=imethod,
                                        nmethod=nmethod,
                                        large=large,
                                        plot_contours=plot_contours,
                                        plot_electrodes=plot_electrodes,
                                        cart=cart)
            ctx = getgc(can)
            if p.attr[:size][1] > 900
                Gtk4.rectangle(ctx, 0, 0, 999, 999)
                Cairo.set_source_rgb(ctx, 255, 255, 255)
                Gtk4.fill(ctx)
                withenv("GKSwstype" => "100") do
                    png(p, io)
                end
                img = read_from_png(io)
                set_source_surface(ctx, img, 500 - (p.attr[:size][1] ÷ 2) - 1, 500 - (p.attr[:size][1] ÷ 2) - 1)
            else
                Gtk4.rectangle(ctx, 0, 0, 799, 799)
                Cairo.set_source_rgb(ctx, 255, 255, 255)
                Gtk4.fill(ctx)
                withenv("GKSwstype" => "100") do
                    png(p, io)
                end
                img = read_from_png(io)
                set_source_surface(ctx, img, 400 - (p.attr[:size][1] ÷ 2) - 1, 400 - (p.attr[:size][1] ÷ 2) - 1)
            end
            paint(ctx)
        end

        signal_connect(entry_ts1, "value-changed") do widget
            seg = round.((entry_ts1.value, entry_ts2.value), digits=3)
            if seg[1] > seg[2]
                warn_dialog(_nill, "Cannot plot!\nSegment start is larger than segment end.", win)
            else
                draw(can)
            end
            seg[1] == seg[2] && (combo_amethod.sensitive = false)
            seg[1] < seg[2] && (combo_amethod.sensitive = true)
        end

        signal_connect(entry_ts2, "value-changed") do widget
            seg = round.((entry_ts1.value, entry_ts2.value), digits=3)
            if seg[1] > seg[2]
                warn_dialog(_nill, "Cannot plot!\nSegment start is larger than segment end.", win)
            else
                draw(can)
            end
            seg[1] == seg[2] && (combo_amethod.sensitive = false)
            seg[1] < seg[2] && (combo_amethod.sensitive = true)
        end

        signal_connect(entry_epoch, "value-changed") do widget
            draw(can)
        end
        signal_connect(bt_refresh, "clicked") do widget
            draw(can)
        end
        signal_connect(combo_imethod, "changed") do widget
            draw(can)
        end
        signal_connect(combo_amethod, "changed") do widget
            draw(can)
        end
        signal_connect(combo_nmethod, "changed") do widget
            draw(can)
        end
        signal_connect(cb_cb, "toggled") do widget
            draw(can)
        end
        signal_connect(cb_cart, "toggled") do widget
            draw(can)
        end
        signal_connect(cb_large, "toggled") do widget
            draw(can)
        end
        signal_connect(cb_elec, "toggled") do widget
            draw(can)
        end
        signal_connect(cb_contour, "toggled") do widget
            draw(can)
        end
        signal_connect(bt_close, "clicked") do widget
            close(win)
        end

        help = "Keyboard shortcuts:\n\nCtrl + s\t\t\tSave as PNG\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tExit\n"

        signal_connect(bt_help, "clicked") do widget
            info_dialog(_nill, help, win)
        end

        win_key = Gtk4.GtkEventControllerKey(win)

        signal_connect(win_key, "key-pressed") do widget, keyval, keycode, state
            if ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('s'))
                save_dialog("Pick an image file", win, ["*.png"]) do file_name
                    if file_name != ""
                        surface_buf = Gtk4.cairo_surface(can)
                        if Cairo.write_to_png(surface_buf, file_name) == Cairo.STATUS_SUCCESS
                            _info("Plot saved as: $file_name")
                        else
                            warn_dialog(_nill, "File $file_name cannot be written!", win)
                        end
                    end
                end
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('h'))
                info_dialog(_nill, help, win)
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('q'))
                close(win)
            end
        end
    end

    app = GtkApplication("org.neuroanalyzer.itopo_ep")
    Gtk4.signal_connect(_activate, app, :activate)
    Gtk4.GLib.stop_main_loop()
    Gtk4.run(app)

    return nothing

end
