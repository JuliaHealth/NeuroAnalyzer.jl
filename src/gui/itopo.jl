export itopo
export itopo_cont
export itopo_ep

"""
    itopo(obj, ch, zoom)

Interactive topoplot of continuous or epoched signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=get_channel_bytype(obj, type="eeg")`: channel(s) to plot, default is EEG channels
- `zoom::Real=5`: how many seconds are displayed in one segment
"""
function itopo(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=get_channel_bytype(obj, type="eeg"), zoom::Real=5)

    _wip()

    _check_datatype(obj, "eeg")

    if nepochs(obj) == 1
        itopo_cont(obj, ch=ch, zoom=zoom)
    else
        itopo_ep(obj, ch=ch)
    end

    return nothing

end

"""
    itopo_cont(obj, ch, zoom)

Interactive topoplot of continuous signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=get_channel_bytype(obj, type="eeg")`: channel(s) to plot, default is all channels
- `zoom::Real=5`: how many seconds are displayed in one segment
"""
function itopo_cont(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=get_channel_bytype(obj, type="eeg"), zoom::Real=5)

    _check_datatype(obj, "eeg")

    (signal_len(obj) / sr(obj)) < zoom && (zoom = signal_len(obj) / sr(obj))

    @assert zoom > 0 "zoom must be > 0."
    @assert zoom <= signal_len(obj) / sr(obj) "zoom must be ≤ $(signal_len(obj) / sr(obj))."
    @assert nepochs(obj) == 1 "itopo_ep() should be used for epoched object."
    _check_channels(obj, ch)
    ch_init = ch

    p = NeuroAnalyzer.plot_topo(obj, ch=ch)

    win = GtkWindow("NeuroAnalyzer: itopo_cont()", 700, 700)
    win_view = GtkScrolledWindow()
    Gtk.resize!(win, p.attr[:size][2] + 100, p.attr[:size][2] + 40)
    can = GtkCanvas(Int32(p.attr[:size][1]), Int32(p.attr[:size][2]))
    set_gtk_property!(win_view, :min_content_width, 700)
    set_gtk_property!(win_view, :min_content_height, 700)
    set_gtk_property!(win, :border_width, 20)
    set_gtk_property!(win, :resizable, true)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    set_gtk_property!(win, :startup_id, "org.neuroanalyzer")
    push!(win_view, can)
    g = GtkGrid()
    g_opts = GtkGrid()
    set_gtk_property!(g, :column_homogeneous, false)
    set_gtk_property!(g_opts, :column_homogeneous, false)
    set_gtk_property!(g, :column_spacing, 10)
    set_gtk_property!(g, :row_spacing, 10)
    set_gtk_property!(g_opts, :row_spacing, 10)
    set_gtk_property!(g_opts, :column_spacing, 10)
    entry_time = GtkSpinButton(obj.time_pts[1], obj.time_pts[end] - zoom, zoom)
    set_gtk_property!(entry_time, :digits, 2)
    set_gtk_property!(entry_time, :value, obj.time_pts[1])
    set_gtk_property!(entry_time, :tooltip_text, "Time position [s]")
    bt_start = GtkButton("⇤")
    set_gtk_property!(bt_start, :tooltip_text, "Go to the signal beginning")
    bt_prev5 = GtkButton("↞")
    set_gtk_property!(bt_prev5, :tooltip_text, "Go back by $zoom seconds")
    bt_prev = GtkButton("←")
    set_gtk_property!(bt_prev, :tooltip_text, "Go back by 1 second")
    bt_next = GtkButton("→")
    set_gtk_property!(bt_next, :tooltip_text, "Go forward by 1 second")
    bt_next5 = GtkButton("↠")
    set_gtk_property!(bt_next5, :tooltip_text, "Go forward by $zoom seconds")
    bt_end = GtkButton("⇥")
    set_gtk_property!(bt_end, :tooltip_text, "Go to the signal end")
    bt_help = GtkButton("🛈")
    set_gtk_property!(bt_help, :tooltip_text, "Show keyboard shortcuts")
    bt_close = GtkButton("✖")
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

    combo_imethod = GtkComboBoxText()
    imethod_types = ["shepard", "multiquadratic", "inv multiquadratic", "thin plate", "nearest neighbour", "gaussian"]
    for idx in imethod_types
        push!(combo_imethod, idx)
    end
    set_gtk_property!(combo_imethod, :active, 0)
    set_gtk_property!(combo_imethod, :tooltip_text, "Interpolate type")

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
    lab_t = GtkLabel("Title:")
    set_gtk_property!(lab_t, :halign, 2)
    lab_cb = GtkLabel("CB lab:")
    set_gtk_property!(lab_cb, :halign, 2)
    lab_cb = GtkLabel("Draw CB:")
    set_gtk_property!(lab_cb, :halign, 2)
    g_opts[1, 1] = lab_type
    g_opts[1, 2] = lab_t
    g_opts[1, 3] = lab_cb
    g_opts[1, 4] = lab_cb
    g_opts[1, 5] = bt_save
    g_opts[2, 1] = combo_imethod
    g_opts[2, 2] = entry_title
    g_opts[2, 3] = entry_cblab
    g_opts[2, 4] = cb_cb
    g_opts[2, 5] = combo_save
    g_opts[1:2, 6] = bt_refresh
    vbox = GtkBox(:v)
    push!(vbox, g_opts)

    g[1, 1] = vbox
    g[2:11, 1] = win_view
    g[1, 2] = GtkLabel("")
    g[2, 2] = bt_start
    g[3, 2] = bt_prev5
    g[4, 2] = bt_prev
    g[5, 2] = entry_time
    g[6, 2] = bt_next
    g[7, 2] = bt_next5
    g[8, 2] = bt_end
    g[9, 2] = GtkLabel("")
    g[10, 2] = bt_help
    g[11, 2] = bt_close
    push!(win, g)

    showall(win)

    @guarded draw(can) do widget
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
        amethod = :mean
        nmethod = :minmax
        plot_contours = true
        plot_electrodes = true
        cart = false
        time1 = get_gtk_property(entry_time, :value, Float64)
        time2 = time1 + zoom
        time2 > obj.time_pts[end] && (time2 = obj.time_pts[end])
        p = NeuroAnalyzer.plot_topo(obj,
                                    ch=ch,
                                    seg=(time1, time2),
                                    title=title,
                                    cb=cb,
                                    cb_label=cblab,
                                    amethod=amethod,
                                    imethod=imethod,
                                    nmethod=nmethod,
                                    plot_contours=plot_contours,
                                    plot_electrodes=plot_electrodes,
                                    cart=cart)
        Gtk.resize!(win, p.attr[:size][2] + 100, p.attr[:size][2] + 40)
        set_gtk_property!(win_view, :min_content_width, p.attr[:size][1])
        set_gtk_property!(win_view, :min_content_height, p.attr[:size][2])
        ctx = getgc(can)
        show(io, MIME("image/png"), p)
        img = read_from_png(io)
        set_source_surface(ctx, img, 0, 0)
        paint(ctx)
    end

    signal_connect(bt_refresh, "clicked") do widget
        draw(can)
    end
    signal_connect(entry_time, "value-changed") do widget
        draw(can)
    end
    signal_connect(combo_imethod, "changed") do widget
        draw(can)
    end
    signal_connect(cb_cb, "clicked") do widget
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
                    warn_dialog("Incorrect file name!")
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
                    warn_dialog("Incorrect file name!")
                end
            end
        end
    end

    signal_connect(bt_prev, "clicked") do widget
        time_current = get_gtk_property(entry_time, :value, Float64)
        if time_current >= obj.time_pts[1] + 1
            time_current -= 1
            Gtk.@sigatom begin
                set_gtk_property!(entry_time, :value, time_current)
            end
        end
    end

    signal_connect(bt_prev5, "clicked") do widget
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
            time_current += 1
            Gtk.@sigatom begin
                set_gtk_property!(entry_time, :value, time_current)
            end
        else
            time_current = obj.time_pts[end] - zoom
            Gtk.@sigatom begin
                set_gtk_property!(entry_time, :value, time_current)
            end
        end
    end

    signal_connect(bt_next5, "clicked") do widget
        time_current = get_gtk_property(entry_time, :value, Float64)
        if time_current < obj.time_pts[end] - zoom
            time_current += zoom
            Gtk.@sigatom begin
                set_gtk_property!(entry_time, :value, time_current)
            end
        else
            time_current = obj.time_pts[end] - zoom
            Gtk.@sigatom begin
                set_gtk_property!(entry_time, :value, time_current)
            end
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

    signal_connect(bt_help, "clicked") do widgete
        info_dialog("Keyboard shortcuts:\n\nctrl-a\tgo to the signal beginning\nctrl-s\tgo to the signal end\nctrl-z\tgo back by 1 second\nctrl-x\tgo forward by 1 second\nctrl-c\tgo back by $zoom seconds\nctrl-v\tgo forward by $zoom seconds\n\nctrl-h\tthis info\nctrl-q\texit\n")
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        s = event.state
        if s == 20
            if k == 113 # q
                Gtk.destroy(win)
            elseif k == 104 # h
                info_dialog("Keyboard shortcuts:\n\nctrl-a\tgo to the signal beginning\nctrl-s\tgo to the signal end\nctrl-z\tgo back by 1 second\nctrl-x\tgo forward by 1 second\nctrl-c\tgo back by $zoom seconds\nctrl-v\tgo forward by $zoom seconds\n\nctrl-h\tthis info\nctrl-q\texit\n")
            elseif k == 97 # a
                Gtk.@sigatom begin
                    set_gtk_property!(entry_time, :value, obj.time_pts[1])
                end
            elseif k == 115 # s
                time_current = obj.time_pts[end] - zoom
                Gtk.@sigatom begin
                    set_gtk_property!(entry_time, :value, time_current)
                end
            elseif k == 122 # z
                time_current = get_gtk_property(entry_time, :value, Float64)
                if time_current >= obj.time_pts[1] + 1
                    time_current -= 1
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_time, :value, time_current)
                    end
                end
            elseif k == 99 # c
                time_current = get_gtk_property(entry_time, :value, Float64)
                if time_current >= obj.time_pts[1] + zoom
                    time_current = time_current - zoom
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_time, :value, time_current)
                    end
                end
            elseif k == 120 # x
                time_current = get_gtk_property(entry_time, :value, Float64)
                if time_current < obj.time_pts[end] - zoom
                    time_current += 1
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_time, :value, time_current)
                    end
                else
                    time_current = obj.time_pts[end] - zoom
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_time, :value, time_current)
                    end
                end
            elseif k == 118 # v
                time_current = get_gtk_property(entry_time, :value, Float64)
                if time_current < obj.time_pts[end] - zoom
                    time_current += zoom
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_time, :value, time_current)
                    end
                else
                    time_current = obj.time_pts[end] - zoom
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_time, :value, time_current)
                    end
                end
            end
        end
    end

    return nothing

end

"""
    itopo_ep(obj, ch)

Interactive topoplot of epoched signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=get_channel_bytype(obj, type="eeg")`: channel(s) to plot, default is all channels
"""
function itopo_ep(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=get_channel_bytype(obj, type="eeg"))

    _check_datatype(obj, "eeg")

    @assert nepochs(obj) > 1 "itopo_cont() should be used for continuous object."
    _check_channels(obj, ch)

    p = NeuroAnalyzer.plot(obj, ch=ch, ep=1)

    win = GtkWindow("NeuroAnalyzer: itopo_ep()", 1200, 800)
    win_view = GtkScrolledWindow()
    set_gtk_property!(win_view, :min_content_width, 1200)
    set_gtk_property!(win_view, :min_content_height, 800)
    set_gtk_property!(win, :border_width, 20)
    set_gtk_property!(win, :resizable, true)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    set_gtk_property!(win, :startup_id, "org.neuroanalyzer")
    can = GtkCanvas(Int32(p.attr[:size][1]), Int32(p.attr[:size][2]))
    push!(win_view, can)
    g = GtkGrid()
    g_opts = GtkGrid()
    set_gtk_property!(g, :column_homogeneous, false)
    set_gtk_property!(g_opts, :column_homogeneous, false)
    set_gtk_property!(g, :column_spacing, 10)
    set_gtk_property!(g, :row_spacing, 10)
    set_gtk_property!(g_opts, :row_spacing, 10)
    set_gtk_property!(g_opts, :column_spacing, 10)
    entry_epoch = GtkSpinButton(1, nepochs(obj), 1)
    set_gtk_property!(entry_epoch, :tooltip_text, "Epoch")
    bt_start = GtkButton("⇤")
    set_gtk_property!(bt_start, :tooltip_text, "Go to the signal beginning")
    bt_prev = GtkButton("←")
    set_gtk_property!(bt_prev, :tooltip_text, "Go back by 1 epoch")
    bt_next = GtkButton("→")
    set_gtk_property!(bt_next, :tooltip_text, "Go forward by 1 epoch")
    bt_end = GtkButton("⇥")
    set_gtk_property!(bt_end, :tooltip_text, "Go to the signal end")
    bt_help = GtkButton("🛈")
    set_gtk_property!(bt_help, :tooltip_text, "Show keyboard shortcuts")
    bt_close = GtkButton("✖")
    set_gtk_property!(bt_close, :tooltip_text, "Close this window")

    entry_title = GtkEntry()
    set_gtk_property!(entry_title, :text, "default")
    set_gtk_property!(entry_title, :tooltip_text, "Plot title")
    entry_xlab = GtkEntry()
    set_gtk_property!(entry_xlab, :text, "default")
    set_gtk_property!(entry_xlab, :tooltip_text, "X axis label")
    entry_ylab = GtkEntry()
    set_gtk_property!(entry_ylab, :text, "default")
    set_gtk_property!(entry_ylab, :tooltip_text, "Y axis label")

    entry_ch = GtkEntry()
    ch = _i2s(ch)
    ch_init = ch
    set_gtk_property!(entry_ch, :text, string(ch))
    set_gtk_property!(entry_ch, :tooltip_text, "Channels")

    cb_mono = GtkCheckButton()
    set_gtk_property!(cb_mono, :tooltip_text, "Use color or gray palette")

    cb_cb = GtkCheckButton()
    set_gtk_property!(cb_cb, :tooltip_text, "Add color bar to the plot")
    set_gtk_property!(cb_cb, :active, true)

    cb_markers = GtkCheckButton()
    set_gtk_property!(cb_markers, :tooltip_text, "Draw event markers")
    set_gtk_property!(cb_markers, :active, true)

    cb_norm = GtkCheckButton()
    set_gtk_property!(cb_norm, :tooltip_text, "Normalize signal for butterfly and averaged plots")
    set_gtk_property!(cb_norm, :active, false)

    combo_imethod = GtkComboBoxText()
    imethod_types = ["normal", "mean", "butterfly"]
    for idx in imethod_types
        push!(combo_imethod, idx)
    end
    set_gtk_property!(combo_imethod, :active, 0)
    set_gtk_property!(combo_imethod, :tooltip_text, "Plot type")

    combo_save = GtkComboBoxText()
    file_types = ["PNG", "PDF"]
    for idx in file_types
        push!(combo_save, idx)
    end
    set_gtk_property!(combo_save, :active, 0)
    bt_save = GtkButton("Save as:")

    bt_refresh = GtkButton("Refresh")
    set_gtk_property!(bt_refresh, :tooltip_text, "Refresh the plot")

    lab_type = GtkLabel("Plot type:")
    set_gtk_property!(lab_type, :halign, 2)
    lab_ch = GtkLabel("Channels:")
    set_gtk_property!(lab_ch, :halign, 2)
    lab_t = GtkLabel("Title:")
    set_gtk_property!(lab_t, :halign, 2)
    lab_x = GtkLabel("X lab:")
    set_gtk_property!(lab_x, :halign, 2)
    lab_y = GtkLabel("Y lab:")
    set_gtk_property!(lab_y, :halign, 2)
    lab_cb = GtkLabel("Draw scale:")
    set_gtk_property!(lab_cb, :halign, 2)
    lab_mono = GtkLabel("Grayscale:")
    set_gtk_property!(lab_mono, :halign, 2)
    lab_cb = GtkLabel("Draw scale:")
    set_gtk_property!(lab_cb, :halign, 2)
    lab_norm = GtkLabel("Normalize:")
    set_gtk_property!(lab_norm, :halign, 2)
    lab_markers = GtkLabel("Draw markers:")
    set_gtk_property!(lab_markers, :halign, 2)
    g_opts[1, 1] = lab_type
    g_opts[1, 2] = lab_ch
    g_opts[1, 3] = lab_t
    g_opts[1, 4] = lab_x
    g_opts[1, 5] = lab_y
    g_opts[1, 6] = lab_mono
    g_opts[1, 7] = lab_cb
    g_opts[1, 8] = lab_norm
    g_opts[1, 9] = lab_markers
    g_opts[1, 10] = bt_save
    g_opts[2, 1] = combo_imethod
    g_opts[2, 2] = entry_ch
    g_opts[2, 3] = entry_title
    g_opts[2, 4] = entry_xlab
    g_opts[2, 5] = entry_ylab
    g_opts[2, 6] = cb_mono
    g_opts[2, 7] = cb_cb
    g_opts[2, 8] = cb_norm
    g_opts[2, 9] = cb_markers
    g_opts[2, 10] = combo_save
    g_opts[1:2, 11] = bt_refresh
    vbox = GtkBox(:v)
    push!(vbox, g_opts)

    g[1, 1] = vbox
    g[2:9, 1] = win_view
    g[1, 2] = GtkLabel("")
    g[2, 2] = bt_start
    g[3, 2] = bt_prev
    g[4, 2] = entry_epoch
    g[5, 2] = bt_next
    g[6, 2] = bt_end
    g[7, 2] = GtkLabel("")
    g[8, 2] = bt_help
    g[9, 2] = bt_close
    push!(win, g)

    showall(win)

    @guarded draw(can) do widget
        ch = get_gtk_property(entry_ch, :text, String)
        if length(ch) < 1
            warn_dialog("Incorrect list of channels.")
            ch = ch_init
            set_gtk_property!(entry_ch, :text, string(ch))
        end
        if (occursin(", ", ch) && _check_svec(ch)) || (occursin(":", ch) && _check_srange(ch)) || _check_sint(ch)
            ch = _s2i(ch)
            if ch isa Int64 && !in(ch, get_channel_bytype(obj))
                warn_dialog("Incorrect list of channels.")
                ch = ch_init
                set_gtk_property!(entry_ch, :text, string(ch))
            elseif ch isa Int64 == false && intersect(ch, get_channel_bytype(obj)) != ch
                warn_dialog("Incorrect list of channels.")
                ch = ch_init
                set_gtk_property!(entry_ch, :text, string(ch))
            end
            title = get_gtk_property(entry_title, :text, String)
            xlab = get_gtk_property(entry_xlab, :text, String)
            ylab = get_gtk_property(entry_ylab, :text, String)
            mono = get_gtk_property(cb_mono, :active, Bool)
            norm = get_gtk_property(cb_norm, :active, Bool)
            markers = get_gtk_property(cb_markers, :active, Bool)
            scale = get_gtk_property(cb_cb, :active, Bool)
            type = get_gtk_property(combo_imethod, :active, String)
            type == "0" && (type = :normal)
            type == "1" && (type = :mean)
            type == "2" && (type = :butterfly)
            if length(unique(obj.header.recording[:channel_type][ch])) > 1 && type === :butterfly
                warn_dialog("For butterfly plot, all channels should be of the same type.")
            elseif length(unique(obj.header.recording[:channel_type][ch])) > 1 && type === :mean
                warn_dialog("For mean plot all, channels should be of the same type.")
            elseif length(ch) < 2 && type === :butterfly
                warn_dialog("For butterfly plot, the signal must contain ≥ 2 channels.")
            elseif length(ch) < 2 && type === :mean
                warn_dialog("For mean plot, the signal must contain ≥ 2 channels.")
            else
                ep = get_gtk_property(entry_epoch, :value, Int64)
                p = NeuroAnalyzer.plot(obj,
                                       ep=ep,
                                       ch=ch,
                                       type=type,
                                       scale=scale,
                                       mono=mono,
                                       title=title,
                                       xlabel=xlab,
                                       ylabel=ylab,
                                       norm=norm,
                                       markers=markers)
                if p.attr[:size][2] + 40 < 800
                    Gtk.resize!(win, 1200, p.attr[:size][2] + 40)
                    set_gtk_property!(win_view, :min_content_height, p.attr[:size][2])
                else
                    Gtk.resize!(win, 1200, 800)
                    set_gtk_property!(win_view, :min_content_height, 800)
                end
                set_gtk_property!(can, :width_request, Int32(p.attr[:size][1]))
                set_gtk_property!(can, :height_request, Int32(p.attr[:size][2]))
                ctx = getgc(can)
                show(io, MIME("image/png"), p)
                img = read_from_png(io)
                set_source_surface(ctx, img, 0, 0)
                paint(ctx)
            end
        else
            warn_dialog("Incorrect channels!")
            set_gtk_property!(entry_ch, :text, string(ch_init))
        end
    end

    signal_connect(bt_refresh, "clicked") do widget
        draw(can)
    end
    signal_connect(entry_epoch, "value-changed") do widget
        draw(can)
    end
    signal_connect(combo_imethod, "changed") do widget
        draw(can)
    end
    signal_connect(cb_norm, "clicked") do widget
        draw(can)
    end
    signal_connect(cb_mono, "clicked") do widget
        draw(can)
    end
    signal_connect(cb_cb, "clicked") do widget
        draw(can)
    end
    signal_connect(cb_markers, "clicked") do widget
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
                    warn_dialog("Incorrect file name!")
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
                    warn_dialog("Incorrect file name!")
                end
            end
        end
    end

    signal_connect(bt_prev, "clicked") do widget
        ep = get_gtk_property(entry_epoch, :value, Int64)
        if ep >= 2
            ep -= 1
            Gtk.@sigatom begin
                set_gtk_property!(entry_epoch, :value, ep)
            end
        end
    end

    signal_connect(bt_next, "clicked") do widget
        ep = get_gtk_property(entry_epoch, :value, Int64)
        if ep < nepochs(obj)
            ep += 1
            Gtk.@sigatom begin
                set_gtk_property!(entry_epoch, :value, ep)
            end
        end
    end

    signal_connect(bt_start, "clicked") do widget
        Gtk.@sigatom begin
            set_gtk_property!(entry_epoch, :value, 1)
        end
    end

    signal_connect(bt_end, "clicked") do widget
        Gtk.@sigatom begin
            set_gtk_property!(entry_epoch, :value, nepochs(obj))
        end
    end

    signal_connect(bt_close, "clicked") do widget
        Gtk.destroy(win)
    end

    signal_connect(bt_help, "clicked") do widgete
        info_dialog("Keyboard shortcuts:\nctrl-a\tgo to first epoch\nctrl-s\tgo to last epoch\nctrl-z\tprevious epoch\nctrl-x\tnext epoch\n\nctrl-h\tthis info\nctrl-q\texit\n")
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        s = event.state
        if s == 20
            if k == 113 # q
                Gtk.destroy(win)
            elseif k == 104 # h
                info_dialog("Keyboard shortcuts:\nctrl-a\tgo to first epoch\nctrl-s\tgo to last epoch\nctrl-z\tprevious epoch\nctrl-x\tnext epoch\n\nctrl-h\tthis info\nctrl-q\texit\n")
            elseif k == 97 # a
                Gtk.@sigatom begin
                    set_gtk_property!(entry_epoch, :value, 1)
                end
            elseif k == 115 # a
                Gtk.@sigatom begin
                    set_gtk_property!(entry_epoch, :value, nepochs(obj))
                end
            elseif k == 122 # z
                ep = get_gtk_property(entry_epoch, :value, Int64)
                if ep >= 2
                    ep -= 1
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_epoch, :value, ep)
                    end
                end
            elseif k == 120 # x
                ep = get_gtk_property(entry_epoch, :value, Int64)
                if ep < nepochs(obj)
                    ep += 1
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_epoch, :value, ep)
                    end
                end
            end
        end
    end

    return nothing

end
