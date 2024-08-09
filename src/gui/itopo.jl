export itopo
export itopo_cont
export itopo_ep

"""
    itopo(obj; <keyword arguments>)

Interactive topographical map.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{String, Vector{String}}`: channels to plot
"""
function itopo(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}})

    if nepochs(obj) == 1
        itopo_cont(obj, ch=ch)
    else
        itopo_ep(obj, ch=ch)
    end

    return nothing

end

"""
    itopo_cont(obj; <keyword arguments>)

Interactive topographical map of continuous signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{String, Vector{String}}`: channels to plot
"""
function itopo_cont(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}})

    redraw = true

    _check_datatype(obj, ["eeg", "meg", "erp"])

    p = NeuroAnalyzer.plot_topo(obj, ch=ch)
    size = p.attr[:size]
    if size[1] > 900
        win = GtkWindow("NeuroAnalyzer: itopo_cont()", p.attr[:size][1] + 100, 1000)
        can = GtkCanvas(1000, 1000)
    else
        win = GtkWindow("NeuroAnalyzer: itopo_cont()", p.attr[:size][1] + 100, 800)
        can = GtkCanvas(800, 800)
    end
    set_gtk_property!(win, :border_width, 5)
    set_gtk_property!(win, :resizable, false)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    set_gtk_property!(win, :startup_id, "org.neuroanalyzer")

    g_opts = GtkGrid()
    set_gtk_property!(g_opts, :column_homogeneous, true)
    set_gtk_property!(g_opts, :row_spacing, 5)
    set_gtk_property!(g_opts, :column_spacing, 5)

    entry_ts1 = GtkSpinButton(obj.time_pts[1], obj.time_pts[end], 0.01)
    set_gtk_property!(entry_ts1, :value, obj.time_pts[1])
    set_gtk_property!(entry_ts1, :digits, 3)
    set_gtk_property!(entry_ts1, :tooltip_text, "Segment start [s]")

    entry_ts2 = GtkSpinButton(obj.time_pts[1], obj.time_pts[end], 0.01)
    set_gtk_property!(entry_ts2, :value, obj.time_pts[1])
    set_gtk_property!(entry_ts2, :digits, 3)
    set_gtk_property!(entry_ts2, :tooltip_text, "Segment end [s]")

    bt_help = GtkButton("Help")
    set_gtk_property!(bt_help, :tooltip_text, "Show help")
    bt_close = GtkButton("Close")
    set_gtk_property!(bt_close, :tooltip_text, "Close this window")

    entry_title = GtkEntry()
    set_gtk_property!(entry_title, :text, "default")
    set_gtk_property!(entry_title, :tooltip_text, "Plot title")

    entry_cblab = GtkEntry()
    set_gtk_property!(entry_cblab, :text, "default")
    set_gtk_property!(entry_cblab, :tooltip_text, "Color bar label")

    cb_cb = GtkCheckButton()
    set_gtk_property!(cb_cb, :tooltip_text, "Draw color bar")
    set_gtk_property!(cb_cb, :active, true)

    cb_cart = GtkCheckButton()
    set_gtk_property!(cb_cart, :tooltip_text, "Use Cartesian coordinates of electrodes")
    set_gtk_property!(cb_cart, :active, false)

    cb_large = GtkCheckButton()
    set_gtk_property!(cb_large, :tooltip_text, "Draw large plot")
    set_gtk_property!(cb_large, :active, true)

    cb_elec = GtkCheckButton()
    set_gtk_property!(cb_elec, :tooltip_text, "Draw electrodes")
    set_gtk_property!(cb_elec, :active, true)

    cb_contour = GtkCheckButton()
    set_gtk_property!(cb_contour, :tooltip_text, "Draw contour lines")
    set_gtk_property!(cb_contour, :active, true)

    combo_imethod = GtkComboBoxText()
    imethod_types = ["shepard", "multiquadratic", "inv multiquadratic", "thin plate", "nearest neighbour", "gaussian"]
    for idx in imethod_types
        push!(combo_imethod, idx)
    end
    set_gtk_property!(combo_imethod, :active, 0)
    set_gtk_property!(combo_imethod, :tooltip_text, "Interpolate type")

    combo_amethod = GtkComboBoxText()
    amethod_types = ["mean", "median"]
    for idx in amethod_types
        push!(combo_amethod, idx)
    end
    set_gtk_property!(combo_amethod, :active, 0)
    set_gtk_property!(combo_amethod, :tooltip_text, "Averaging method")
    set_gtk_property!(combo_amethod, :sensitive, false)
    combo_nmethod = GtkComboBoxText()
    nmethod_types = ["zscore", "gauss", "invroot", "log", "minmax", "neg", "neglog", "neglog10", "perc", "pos", "softmax", "none"]
    for idx in nmethod_types
        push!(combo_nmethod, idx)
    end
    set_gtk_property!(combo_nmethod, :active, 0)
    set_gtk_property!(combo_nmethod, :tooltip_text, "Normalization method")

    bt_refresh = GtkButton("Refresh")
    set_gtk_property!(bt_refresh, :tooltip_text, "Refresh the plot")

    lab_type = GtkLabel("Interpolate:")
    set_gtk_property!(lab_type, :halign, 2)
    lab_amethod = GtkLabel("Averaging:")
    set_gtk_property!(lab_amethod, :halign, 2)
    lab_nmethod = GtkLabel("Normalization:")
    set_gtk_property!(lab_nmethod, :halign, 2)
    lab_t = GtkLabel("Title:")
    set_gtk_property!(lab_t, :halign, 2)
    lab_cb = GtkLabel("Color bar title:")
    set_gtk_property!(lab_cb, :halign, 2)
    lab_cb_draw = GtkLabel("Draw color bar:")
    set_gtk_property!(lab_cb_draw, :halign, 2)
    lab_cart = GtkLabel("Cartesian:")
    set_gtk_property!(lab_cart, :halign, 2)
    lab_large = GtkLabel("Large plot:")
    set_gtk_property!(lab_large, :halign, 2)
    lab_elec = GtkLabel("Draw electrodes:")
    set_gtk_property!(lab_elec, :halign, 2)
    lab_contour = GtkLabel("Draw contours:")
    set_gtk_property!(lab_contour, :halign, 2)
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
    set_gtk_property!(g, :column_homogeneous, false)
    set_gtk_property!(g, :row_spacing, 5)
    set_gtk_property!(g, :column_spacing, 5)
    vbox = GtkBox(:v)
    push!(vbox, g_opts)
    g[1, 1] = vbox
    g[2:11, 1] = can
    push!(win, g)

    showall(win)

    @guarded draw(can) do widget
        if redraw
            seg = round.((get_gtk_property(entry_ts1, :value, Float64), get_gtk_property(entry_ts2, :value, Float64)), digits=3)
            title = get_gtk_property(entry_title, :text, String)
            cblab = get_gtk_property(entry_cblab, :text, String)
            cb = get_gtk_property(cb_cb, :active, Bool)
            imethod = get_gtk_property(combo_imethod, :active, String)
            imethod == "0" && (imethod = :sh)
            imethod == "1" && (imethod = :mq)
            imethod == "2" && (imethod = :imq)
            imethod == "3" && (imethod = :tp)
            imethod == "4" && (imethod = :nn)
            imethod == "5" && (imethod = :ga)
            amethod = get_gtk_property(combo_amethod, :active, String)
            amethod == "0" && (amethod = :mean)
            amethod == "1" && (amethod = :median)
            nmethod = get_gtk_property(combo_nmethod, :active, String)
            nmethod == "0" && (nmethod = :zscore)
            nmethod == "1" && (nmethod = :gauss)
            nmethod == "2" && (nmethod = :invroot)
            nmethod == "3" && (nmethod = :log)
            nmethod == "4" && (nmethod = :minmax)
            nmethod == "5" && (nmethod = :neg)
            nmethod == "6" && (nmethod = :neglog)
            nmethod == "7" && (nmethod = :neglog10)
            nmethod == "8" && (nmethod = :perc)
            nmethod == "9" && (nmethod = :pos)
            nmethod == "10" && (nmethod = :softmax)
            nmethod == "11" && (nmethod = :none)
            plot_contours = get_gtk_property(cb_contour, :active, Bool)
            plot_electrodes = get_gtk_property(cb_elec, :active, Bool)
            cart = get_gtk_property(cb_cart, :active, Bool)
            large = get_gtk_property(cb_large, :active, Bool)
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
            if size[1] > 900
                Gtk.rectangle(ctx, 0, 0, 999, 999)
                Cairo.set_source_rgb(ctx, 255, 255, 255)
                Gtk.fill(ctx)
                show(io, MIME("image/png"), p)
                img = read_from_png(io)
                set_source_surface(ctx, img, 500 - (p.attr[:size][1] ÷ 2) - 1, 500 - (p.attr[:size][1] ÷ 2) - 1)
            else
                Gtk.rectangle(ctx, 0, 0, 799, 799)
                Cairo.set_source_rgb(ctx, 255, 255, 255)
                Gtk.fill(ctx)
                show(io, MIME("image/png"), p)
                img = read_from_png(io)
                set_source_surface(ctx, img, 400 - (p.attr[:size][1] ÷ 2) - 1, 400 - (p.attr[:size][1] ÷ 2) - 1)
            end
            paint(ctx)
        end
    end

    signal_connect(entry_ts1, "value-changed") do widget
        seg = round.((get_gtk_property(entry_ts1, :value, Float64), get_gtk_property(entry_ts2, :value, Float64)), digits=3)
        if seg[1] > seg[2]
            warn_dialog("Cannot plot!\nSegment start is larger than segment end.")
        else
            draw(can)
        end
        Gtk.@sigatom begin
            seg[1] == seg[2] && set_gtk_property!(combo_amethod, :sensitive, false)
            seg[1] != seg[2] && set_gtk_property!(combo_amethod, :sensitive, true)
        end
    end
    signal_connect(entry_ts2, "value-changed") do widget
        seg = round.((get_gtk_property(entry_ts1, :value, Float64), get_gtk_property(entry_ts2, :value, Float64)), digits=3)
        if seg[1] > seg[2]
            warn_dialog("Cannot plot!\nSegment start is larger than segment end.")
        else
            draw(can)
        end
        Gtk.@sigatom begin
            seg[1] == seg[2] && set_gtk_property!(combo_amethod, :sensitive, false)
            seg[1] != seg[2] && set_gtk_property!(combo_amethod, :sensitive, true)
        end
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
    signal_connect(cb_cb, "clicked") do widget
        draw(can)
    end
    signal_connect(cb_cart, "clicked") do widget
        draw(can)
    end
    signal_connect(cb_large, "clicked") do widget
        draw(can)
    end
    signal_connect(cb_elec, "clicked") do widget
        draw(can)
    end
    signal_connect(cb_contour, "clicked") do widget
        draw(can)
    end
    signal_connect(bt_close, "clicked") do widget
        Gtk.destroy(win)
    end

    help = "Keyboard shortcuts:\n\nCtrl + s\t\t\tSave as PNG\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tExit\n"

    signal_connect(bt_help, "clicked") do widgete
        info_dialog(help)
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        s = event.state
        if s == 0x00000004 # ctrl
            if k == 0x00000071 # q
                Gtk.destroy(win)
            elseif k == 0x00000068 # h
                info_dialog(help)
            elseif k == 0x00000073 # s
                file_name = save_dialog_native("Save as PNG", GtkNullContainer(), (GtkFileFilter("*.png", name="All supported formats"), "*.png"));
                if file_name != ""
                    splitext(file_name)[2] == "" && (file_name *= ".png")
                    if splitext(file_name)[2] == ".png"
                        if plot_save(p, file_name=file_name) == -1
                            warn_dialog("File $file_name cannot be written!")
                        else
                            _info("Plot saved as: $file_name")
                        end
                    else
                        warn_dialog("Incorrect filename!")
                    end
                end
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

"""
    itopo_ep(obj; <keyword arguments>)

Interactive topographical map of epoched signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{String, Vector{String}}`: channels to plot
"""
function itopo_ep(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}})

    @assert nepochs(obj) > 1 "itopo_cont() must be used for continuous object."

    _check_datatype(obj, ["eeg", "meg", "erp"])

    p = NeuroAnalyzer.plot_topo(obj, ch=ch)
    size = p.attr[:size]
    if size[1] > 900
        win = GtkWindow("NeuroAnalyzer: itopo_ep()", p.attr[:size][1] + 100, 1000)
        can = GtkCanvas(1000, 1000)
    else
        win = GtkWindow("NeuroAnalyzer: itopo_ep()", p.attr[:size][1] + 100, 800)
        can = GtkCanvas(800, 800)
    end
    set_gtk_property!(win, :border_width, 5)
    set_gtk_property!(win, :resizable, false)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    set_gtk_property!(win, :startup_id, "org.neuroanalyzer")

    g_opts = GtkGrid()
    set_gtk_property!(g_opts, :column_homogeneous, true)
    set_gtk_property!(g_opts, :row_spacing, 5)
    set_gtk_property!(g_opts, :column_spacing, 5)

    entry_epoch = GtkSpinButton(1, nepochs(obj), 1)
    set_gtk_property!(entry_epoch, :tooltip_text, "Epoch")

    entry_ts1 = GtkSpinButton(obj.epoch_time[1], obj.epoch_time[end], 0.01)
    set_gtk_property!(entry_ts1, :value, obj.epoch_time[1])
    set_gtk_property!(entry_ts1, :digits, 3)
    set_gtk_property!(entry_ts1, :tooltip_text, "Segment start [s]")

    entry_ts2 = GtkSpinButton(obj.epoch_time[1], obj.epoch_time[end], 0.01)
    set_gtk_property!(entry_ts2, :value, obj.epoch_time[end])
    set_gtk_property!(entry_ts2, :digits, 3)
    set_gtk_property!(entry_ts2, :tooltip_text, "Segment end [s]")

    bt_help = GtkButton("Help")
    set_gtk_property!(bt_help, :tooltip_text, "Show help")
    bt_close = GtkButton("Close")
    set_gtk_property!(bt_close, :tooltip_text, "Close this window")

    entry_title = GtkEntry()
    set_gtk_property!(entry_title, :text, "default")
    set_gtk_property!(entry_title, :tooltip_text, "Plot title")

    entry_cblab = GtkEntry()
    set_gtk_property!(entry_cblab, :text, "default")
    set_gtk_property!(entry_cblab, :tooltip_text, "Color bar label")

    cb_cb = GtkCheckButton()
    set_gtk_property!(cb_cb, :tooltip_text, "Draw color bar")
    set_gtk_property!(cb_cb, :active, true)

    cb_cart = GtkCheckButton()
    set_gtk_property!(cb_cart, :tooltip_text, "Use Cartesian coordinates of electrodes")
    set_gtk_property!(cb_cart, :active, false)

    cb_large = GtkCheckButton()
    set_gtk_property!(cb_large, :tooltip_text, "Draw large plot")
    set_gtk_property!(cb_large, :active, true)

    cb_elec = GtkCheckButton()
    set_gtk_property!(cb_elec, :tooltip_text, "Draw electrodes")
    set_gtk_property!(cb_elec, :active, true)

    cb_contour = GtkCheckButton()
    set_gtk_property!(cb_contour, :tooltip_text, "Draw contour lines")
    set_gtk_property!(cb_contour, :active, true)

    combo_imethod = GtkComboBoxText()
    imethod_types = ["shepard", "multiquadratic", "inv multiquadratic", "thin plate", "nearest neighbour", "gaussian"]
    for idx in imethod_types
        push!(combo_imethod, idx)
    end
    set_gtk_property!(combo_imethod, :active, 0)
    set_gtk_property!(combo_imethod, :tooltip_text, "Interpolate type")

    combo_amethod = GtkComboBoxText()
    amethod_types = ["mean", "median"]
    for idx in amethod_types
        push!(combo_amethod, idx)
    end
    set_gtk_property!(combo_amethod, :active, 0)
    set_gtk_property!(combo_amethod, :tooltip_text, "Averaging method")

    combo_nmethod = GtkComboBoxText()
    nmethod_types = ["zscore", "gauss", "invroot", "log", "minmax", "neg", "neglog", "neglog10", "perc", "pos", "softmax", "none"]
    for idx in nmethod_types
        push!(combo_nmethod, idx)
    end
    set_gtk_property!(combo_nmethod, :active, 0)
    set_gtk_property!(combo_nmethod, :tooltip_text, "Normalization method")

    bt_refresh = GtkButton("Refresh")
    set_gtk_property!(bt_refresh, :tooltip_text, "Refresh the plot")

    lab_type = GtkLabel("Interpolate:")
    set_gtk_property!(lab_type, :halign, 2)
    lab_amethod = GtkLabel("Averaging:")
    set_gtk_property!(lab_amethod, :halign, 2)
    lab_nmethod = GtkLabel("Normalization:")
    set_gtk_property!(lab_nmethod, :halign, 2)
    lab_t = GtkLabel("Title:")
    set_gtk_property!(lab_t, :halign, 2)
    lab_cb = GtkLabel("Color bar title:")
    set_gtk_property!(lab_cb, :halign, 2)
    lab_cb_draw = GtkLabel("Draw color bar:")
    set_gtk_property!(lab_cb_draw, :halign, 2)
    lab_cart = GtkLabel("Cartesian:")
    set_gtk_property!(lab_cart, :halign, 2)
    lab_large = GtkLabel("Large plot:")
    set_gtk_property!(lab_large, :halign, 2)
    lab_elec = GtkLabel("Draw electrodes:")
    set_gtk_property!(lab_elec, :halign, 2)
    lab_contour = GtkLabel("Draw contours:")
    set_gtk_property!(lab_contour, :halign, 2)
    lab_epoch = GtkLabel("Epoch:")
    set_gtk_property!(lab_epoch, :halign, 2)
    lab_ts = GtkLabel("Time segment:")
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
    set_gtk_property!(g, :column_homogeneous, false)
    set_gtk_property!(g, :row_spacing, 5)
    set_gtk_property!(g, :column_spacing, 5)
    vbox = GtkBox(:v)
    push!(vbox, g_opts)
    g[1, 1] = vbox
    g[2:11, 1] = can
    push!(win, g)

    showall(win)

    @guarded draw(can) do widget
        ep = get_gtk_property(entry_epoch, :value, Int64)
        ts1 = get_gtk_property(entry_ts1, :value, Float64)
        ts2 = get_gtk_property(entry_ts2, :value, Float64)
        ts1 = (epoch_len(obj) / sr(obj) * (ep - 1)) + obj.epoch_time[1] + ts1
        ts2 = (epoch_len(obj) / sr(obj) * (ep - 1)) + obj.epoch_time[1] + ts2
        title = get_gtk_property(entry_title, :text, String)
        cblab = get_gtk_property(entry_cblab, :text, String)
        cb = get_gtk_property(cb_cb, :active, Bool)
        imethod = get_gtk_property(combo_imethod, :active, String)
        imethod == "0" && (imethod = :sh)
        imethod == "1" && (imethod = :mq)
        imethod == "2" && (imethod = :imq)
        imethod == "3" && (imethod = :tp)
        imethod == "4" && (imethod = :nn)
        imethod == "5" && (imethod = :ga)
        amethod = get_gtk_property(combo_amethod, :active, String)
        amethod == "0" && (amethod = :mean)
        amethod == "1" && (amethod = :median)
        nmethod = get_gtk_property(combo_nmethod, :active, String)
        nmethod == "0" && (nmethod = :zscore)
        nmethod == "1" && (nmethod = :gauss)
        nmethod == "2" && (nmethod = :invroot)
        nmethod == "3" && (nmethod = :log)
        nmethod == "4" && (nmethod = :minmax)
        nmethod == "5" && (nmethod = :neg)
        nmethod == "6" && (nmethod = :neglog)
        nmethod == "7" && (nmethod = :neglog10)
        nmethod == "8" && (nmethod = :perc)
        nmethod == "9" && (nmethod = :pos)
        nmethod == "10" && (nmethod = :softmax)
        nmethod == "11" && (nmethod = :none)
        plot_contours = get_gtk_property(cb_contour, :active, Bool)
        plot_electrodes = get_gtk_property(cb_elec, :active, Bool)
        cart = get_gtk_property(cb_cart, :active, Bool)
        large = get_gtk_property(cb_large, :active, Bool)
        p = NeuroAnalyzer.plot_topo(obj,
                                    ch=ch,
                                    seg=(ts1, ts2),
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
        if size[1] > 900
            Gtk.rectangle(ctx, 0, 0, 999, 999)
            Cairo.set_source_rgb(ctx, 255, 255, 255)
            Gtk.fill(ctx)
            show(io, MIME("image/png"), p)
            img = read_from_png(io)
            set_source_surface(ctx, img, 500 - (p.attr[:size][1] ÷ 2) - 1, 500 - (p.attr[:size][1] ÷ 2) - 1)
        else
            Gtk.rectangle(ctx, 0, 0, 799, 799)
            Cairo.set_source_rgb(ctx, 255, 255, 255)
            Gtk.fill(ctx)
            show(io, MIME("image/png"), p)
            img = read_from_png(io)
            set_source_surface(ctx, img, 400 - (p.attr[:size][1] ÷ 2) - 1, 400 - (p.attr[:size][1] ÷ 2) - 1)
        end
        paint(ctx)
    end

    signal_connect(entry_ts1, "value-changed") do widget
        draw(can)
    end
    signal_connect(entry_ts2, "value-changed") do widget
        draw(can)
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
    signal_connect(cb_cb, "clicked") do widget
        draw(can)
    end
    signal_connect(cb_cart, "clicked") do widget
        draw(can)
    end
    signal_connect(cb_large, "clicked") do widget
        draw(can)
    end
    signal_connect(cb_elec, "clicked") do widget
        draw(can)
    end
    signal_connect(cb_contour, "clicked") do widget
        draw(can)
    end
    signal_connect(bt_close, "clicked") do widget
        Gtk.destroy(win)
    end

    help = "Keyboard shortcuts:\n\nCtrl + s\t\t\tSave as PNG\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tExit\n"

    signal_connect(bt_help, "clicked") do widgete
        info_dialog(help)
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        s = event.state
        if s == 0x00000004 # ctrl
            if k == 0x00000071 # q
                Gtk.destroy(win)
            elseif k == 0x00000068 # h
                info_dialog(help)
            elseif k == 0x00000073 # s
                file_name = save_dialog_native("Save as PNG", GtkNullContainer(), (GtkFileFilter("*.png", name="All supported formats"), "*.png"));
                if file_name != ""
                    splitext(file_name)[2] == "" && (file_name *= ".png")
                    if splitext(file_name)[2] == ".png"
                        if plot_save(p, file_name=file_name) == -1
                            warn_dialog("File $file_name cannot be written!")
                        else
                            _info("Plot saved as: $file_name")
                        end
                    else
                        warn_dialog("Incorrect filename!")
                    end
                end
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