export itopo

"""
    itopo(obj; <keyword arguments>)

Interactive topographical map.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{String, Vector{String}}`: channels to plot
- `seg::Tuple{Real, Real}=(0, 0)`: segment (from, to) in seconds to display
"""
function itopo(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}}, seg::Tuple{Real, Real}=(0, 0))

    redraw = true

    _check_datatype(obj, ["eeg", "meg", "erp"])

    ch = get_channel(obj, ch=ch)

    p = NeuroAnalyzer.plot_topo(obj, ch=ch)

    win = GtkWindow("NeuroAnalyzer: itopo()", p.attr[:size][1] + 100, p.attr[:size][2] + 40)
    can = GtkCanvas(p.attr[:size][1], p.attr[:size][2])
    set_gtk_property!(win, :border_width, 20)
    set_gtk_property!(win, :resizable, true)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    set_gtk_property!(win, :startup_id, "org.neuroanalyzer")

    g_opts = GtkGrid()
    set_gtk_property!(g_opts, :column_homogeneous, false)
    set_gtk_property!(g_opts, :row_spacing, 10)
    set_gtk_property!(g_opts, :column_spacing, 10)

    entry_ts1 = GtkSpinButton(obj.time_pts[1], obj.time_pts[end], 0.01)
    set_gtk_property!(entry_ts1, :value, seg[1])
    set_gtk_property!(entry_ts1, :digits, 3)
    set_gtk_property!(entry_ts1, :tooltip_text, "Segment start [s]")

    entry_ts2 = GtkSpinButton(obj.time_pts[1], obj.time_pts[end], 0.01)
    set_gtk_property!(entry_ts2, :value, seg[2])
    set_gtk_property!(entry_ts2, :digits, 3)
    set_gtk_property!(entry_ts2, :tooltip_text, "Segment end [s]")

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
    seg[1] == seg[2] && set_gtk_property!(combo_amethod, :sensitive, false)
    combo_nmethod = GtkComboBoxText()
    nmethod_types = ["zscore", "gauss", "invroot", "log", "minmax", "neg", "neglog", "neglog10", "perc", "pos", "softmax", "none"]
    for idx in nmethod_types
        push!(combo_nmethod, idx)
    end
    set_gtk_property!(combo_nmethod, :active, 0)
    set_gtk_property!(combo_nmethod, :tooltip_text, "Normalization method")

    combo_save = GtkComboBoxText()
    file_types = ["PNG", "PDF"]
    for idx in file_types
        push!(combo_save, idx)
    end
    set_gtk_property!(combo_save, :active, 0)
    bt_save = GtkButton("Save as:")

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
    lab_ts1 = GtkLabel("Segment start")
    lab_ts2 = GtkLabel("Segment end")
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
    g_opts[1, 13] = bt_save
    g_opts[1, 1] = lab_ts1
    g_opts[2, 1] = lab_ts2
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
    g_opts[2, 13] = combo_save
    g_opts[1:2, 14] = bt_refresh
    g_opts[1:2, 15] = bt_close

    g = GtkGrid()
    set_gtk_property!(g, :column_homogeneous, false)
    set_gtk_property!(g, :row_spacing, 10)
    set_gtk_property!(g, :column_spacing, 10)
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
            Gtk.resize!(win, p.attr[:size][2] + 100, p.attr[:size][2] + 40)
            set_gtk_property!(can, :width_request, Int32(p.attr[:size][1]))
            set_gtk_property!(can, :height_request, Int32(p.attr[:size][2]))
            ctx = getgc(can)
            # Gtk.rectangle(ctx, 0, 0, 705, 705)
            # Cairo.set_source_rgb(ctx, 255, 255, 255)
            # Gtk.fill(ctx)
            show(io, MIME("image/png"), p)
            img = read_from_png(io)
            set_source_surface(ctx, img, 0, 0)
            paint(ctx)
        end
    end

    signal_connect(entry_ts1, "value-changed") do widget
        Gtk.@sigatom begin
            seg = round.((get_gtk_property(entry_ts1, :value, Float64), get_gtk_property(entry_ts2, :value, Float64)), digits=3)
        end
        if seg[1] > seg[2]
            warn_dialog("Cannot plot!\nSegment start is larger than segment end.")
        else
            draw(can)
        end
    end
    signal_connect(entry_ts2, "value-changed") do widget
        Gtk.@sigatom begin
            seg = round.((get_gtk_property(entry_ts1, :value, Float64), get_gtk_property(entry_ts2, :value, Float64)), digits=3)
        end
        if seg[1] > seg[2]
            warn_dialog("Cannot plot!\nSegment start is larger than segment end.")
        else
            draw(can)
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

    signal_connect(bt_save, "clicked") do widget
        format = get_gtk_property(combo_save, :active, String)
        if format == "0"
            file_name = save_dialog("Save as PNG", GtkNullContainer(), (GtkFileFilter("*.png", name="All supported formats"), "*.png"))
            if file_name != ""
                splitext(file_name)[2] == "" && (file_name *= ".png")
                if splitext(file_name)[2] == ".png"
                    plot_save(p, file_name=file_name)
                    _info("Plot saved as: $file_name")
                else
                    warn_dialog("Incorrect filename!")
                end
            end
        else
            file_name = save_dialog("Save as PDF", GtkNullContainer(), (GtkFileFilter("*.pdf", name="All supported formats"), "*.pdf"))
            if file_name != ""
                splitext(file_name)[2] == "" && (file_name *= ".pdf")
                if splitext(file_name)[2] == ".pdf"
                    plot_save(p, file_name=file_name)
                    _info("Plot saved as: $file_name")
                else
                    warn_dialog("Incorrect filename!")
                end
            end
        end
    end

    signal_connect(bt_close, "clicked") do widget
        Gtk.destroy(win)
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        s = event.state
        if s == 4
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
