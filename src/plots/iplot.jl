export iplot
export iplot_cont
export iplot_ep

"""
    iplot(obj, ch, zoom)

Interactive plot of continuous or epoched signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: channel(s) to plot, default is all channels
- `zoom::Real=5`: how many seconds are displayed in one segment

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function iplot(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=NeuroAnalyzer._c(nchannels(obj)), zoom::Real=5)

    if nepochs(obj) == 1
        iplot_cont(obj, ch=ch, zoom=zoom)
    else
        iplot_ep(obj, ch=ch)
    end

end

"""
    iplot_cont(obj, ch, zoom)

Interactive plot of continuous signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: channel(s) to plot, default is all channels
- `zoom::Real=5`: how many seconds are displayed in one segment

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function iplot_cont(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=NeuroAnalyzer._c(nchannels(obj)), zoom::Real=5)

    @assert zoom > 0 "zoom must be > 0."
    @assert zoom <= signal_len(obj) / sr(obj) "zoom must be ≤ $(signal_len(obj) / sr(obj))."
    @assert nepochs(obj) == 1 "iplot_ep() should be used for epoched object."
    _check_channels(obj, ch)
    ch_init = ch

    p = NeuroAnalyzer.plot(obj, ch=ch)

    win = GtkWindow("NeuroAnalyzer: iplot_cont()", 1200, 800)
    win_view = GtkScrolledWindow()
    set_gtk_property!(win_view, :min_content_width, 1200)
    set_gtk_property!(win_view, :min_content_height, 800)
    set_gtk_property!(win, :border_width, 20)
    set_gtk_property!(win, :resizable, true)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
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

    cb_scale = GtkCheckButton()
    set_gtk_property!(cb_scale, :tooltip_text, "Add scale to the plot")
    set_gtk_property!(cb_scale, :active, true)

    cb_markers = GtkCheckButton()
    set_gtk_property!(cb_markers, :tooltip_text, "Draw event markers")
    set_gtk_property!(cb_markers, :active, true)

    cb_norm = GtkCheckButton()
    set_gtk_property!(cb_norm, :tooltip_text, "Normalize signal for butterfly and averaged plots")
    set_gtk_property!(cb_norm, :active, false)

    combo_type = GtkComboBoxText()
    plot_types = ["normal", "mean", "butterfly"]
    for idx in plot_types
        push!(combo_type, idx)
    end
    set_gtk_property!(combo_type, :active, 0)
    set_gtk_property!(combo_type, :tooltip_text, "Plot type")

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
    lab_scale = GtkLabel("Draw scale:")
    set_gtk_property!(lab_scale, :halign, 2)
    lab_mono = GtkLabel("Grayscale:")
    set_gtk_property!(lab_mono, :halign, 2)
    lab_scale = GtkLabel("Draw scale:")
    set_gtk_property!(lab_scale, :halign, 2)
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
    g_opts[1, 7] = lab_scale
    g_opts[1, 8] = lab_norm
    g_opts[1, 9] = lab_markers
    g_opts[1, 10] = bt_save
    g_opts[2, 1] = combo_type
    g_opts[2, 2] = entry_ch
    g_opts[2, 3] = entry_title
    g_opts[2, 4] = entry_xlab
    g_opts[2, 5] = entry_ylab
    g_opts[2, 6] = cb_mono
    g_opts[2, 7] = cb_scale
    g_opts[2, 8] = cb_norm
    g_opts[2, 9] = cb_markers
    g_opts[2, 10] = combo_save
    g_opts[1:2, 11] = bt_refresh
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
        ch = get_gtk_property(entry_ch, :text, String)
        if length(ch) < 1
            warn_dialog("Incorrect list of channels.")
            ch = ch_init
            set_gtk_property!(entry_ch, :text, string(ch))
        end
        if (occursin(", ", ch) && _check_svec(ch)) || (occursin(":", ch) && _check_srange(ch)) || _check_sint(ch)
            ch = _s2i(ch)
            if length(ch) == 1 && !in(ch, get_channel_bytype(obj))
                warn_dialog("Incorrect list of channels.")
                ch = ch_init
                set_gtk_property!(entry_ch, :text, string(ch))
            elseif length(ch) > 1 && intersect(ch, get_channel_bytype(obj)) != ch
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
            scale = get_gtk_property(cb_scale, :active, Bool)
            type = get_gtk_property(combo_type, :active, String)
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
                time1 = get_gtk_property(entry_time, :value, Float64)
                time2 = time1 + zoom
                time2 > obj.time_pts[end] && (time2 = obj.time_pts[end])
                p = NeuroAnalyzer.plot(obj,
                                       ch=ch,
                                       type=type,
                                       seg=(time1, time2),
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
    signal_connect(entry_time, "value-changed") do widget
        draw(can)
    end
    signal_connect(combo_type, "changed") do widget
        draw(can)
    end
    signal_connect(cb_norm, "clicked") do widget
        draw(can)
    end
    signal_connect(cb_mono, "clicked") do widget
        draw(can)
    end
    signal_connect(cb_scale, "clicked") do widget
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
        info_dialog("Keyboard shortcuts:\n\na\tgo to the signal beginning\ns\tgo to the signal end\nz\tgo back by 1 second\nx\tgo forward by 1 second\nc\tgo back by $zoom seconds\nv\tgo forward by $zoom seconds\n\nh\tthis info\nq\texit\n")
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        if k == 113 # q
            Gtk.destroy(win)
        elseif k == 104 # h
            info_dialog("Keyboard shortcuts:\n\na\tgo to the signal beginning\ns\tgo to the signal end\nz\tgo back by 1 second\nx\tgo forward by 1 second\nc\tgo back by $zoom seconds\nv\tgo forward by $zoom seconds\n\nh\tthis info\nq\texit\n")
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

    return nothing

end

"""
    iplot_ep(obj, ch)

Interactive plot of epoched signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: channel(s) to plot, default is all channels

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function iplot_ep(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=NeuroAnalyzer._c(nchannels(obj)))

    @assert nepochs(obj) > 1 "iplot_cont() should be used for continuous object."
    _check_channels(obj, ch)

    p = NeuroAnalyzer.plot(obj, ch=ch, ep=1)

    win = GtkWindow("NeuroAnalyzer: iplot_ep()", 1200, 800)
    win_view = GtkScrolledWindow()
    set_gtk_property!(win_view, :min_content_width, 1200)
    set_gtk_property!(win_view, :min_content_height, 800)
    set_gtk_property!(win, :border_width, 20)
    set_gtk_property!(win, :resizable, true)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
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

    cb_scale = GtkCheckButton()
    set_gtk_property!(cb_scale, :tooltip_text, "Add scale to the plot")
    set_gtk_property!(cb_scale, :active, true)

    cb_markers = GtkCheckButton()
    set_gtk_property!(cb_markers, :tooltip_text, "Draw event markers")
    set_gtk_property!(cb_markers, :active, true)

    cb_norm = GtkCheckButton()
    set_gtk_property!(cb_norm, :tooltip_text, "Normalize signal for butterfly and averaged plots")
    set_gtk_property!(cb_norm, :active, false)

    combo_type = GtkComboBoxText()
    plot_types = ["normal", "mean", "butterfly"]
    for idx in plot_types
        push!(combo_type, idx)
    end
    set_gtk_property!(combo_type, :active, 0)
    set_gtk_property!(combo_type, :tooltip_text, "Plot type")

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
    lab_scale = GtkLabel("Draw scale:")
    set_gtk_property!(lab_scale, :halign, 2)
    lab_mono = GtkLabel("Grayscale:")
    set_gtk_property!(lab_mono, :halign, 2)
    lab_scale = GtkLabel("Draw scale:")
    set_gtk_property!(lab_scale, :halign, 2)
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
    g_opts[1, 7] = lab_scale
    g_opts[1, 8] = lab_norm
    g_opts[1, 9] = lab_markers
    g_opts[1, 10] = bt_save
    g_opts[2, 1] = combo_type
    g_opts[2, 2] = entry_ch
    g_opts[2, 3] = entry_title
    g_opts[2, 4] = entry_xlab
    g_opts[2, 5] = entry_ylab
    g_opts[2, 6] = cb_mono
    g_opts[2, 7] = cb_scale
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
            if length(ch) == 1 && !in(ch, get_channel_bytype(obj))
                warn_dialog("Incorrect list of channels.")
                ch = ch_init
                set_gtk_property!(entry_ch, :text, string(ch))
            elseif length(ch) > 1 && intersect(ch, get_channel_bytype(obj)) != ch
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
            scale = get_gtk_property(cb_scale, :active, Bool)
            type = get_gtk_property(combo_type, :active, String)
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
    signal_connect(combo_type, "changed") do widget
        draw(can)
    end
    signal_connect(cb_norm, "clicked") do widget
        draw(can)
    end
    signal_connect(cb_mono, "clicked") do widget
        draw(can)
    end
    signal_connect(cb_scale, "clicked") do widget
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
        info_dialog("Keyboard shortcuts:\na\tgo to first epoch\ns\tgo to last epoch\nz\tprevious epoch\nx\tnext epoch\n\nh\tthis info\nq\texit\n")
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        if k == 113 # q
            Gtk.destroy(win)
        elseif k == 104 # h
            info_dialog("Keyboard shortcuts:\na\tgo to first epoch\ns\tgo to last epoch\nz\tprevious epoch\nx\tnext epoch\n\nh\tthis info\nq\texit\n")
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

    return nothing

end

"""
    iplot(obj1, obj2, ch, zoom)

Interactive plot of two continuous or epoched signals.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `obj2::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: channel(s) to plot, default is all channels
- `zoom::Real=5`: how many seconds are displayed in one segment

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function iplot(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=NeuroAnalyzer._c(nchannels(obj1)), zoom::Real=5)

    if nepochs(obj1) == 1 && nepochs(obj2) == 1
        iplot_cont(obj1, obj2, ch=ch, zoom=zoom)
    else
        iplot_ep(obj1, obj2, ch=ch)
    end

end

"""
    iplot_cont(obj1, obj2, ch, zoom)

Interactive plot of two continuous signals.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `obj2::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: channel(s) to plot, default is all channels
- `zoom::Real=5`: how many seconds are displayed in one segment

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function iplot_cont(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=NeuroAnalyzer._c(nchannels(obj1)), zoom::Real=5)

    @assert zoom > 0 "zoom must be > 0."
    @assert zoom <= signal_len(obj1) / sr(obj1) "zoom must be ≤ $(signal_len(obj1) / sr(obj1))."
    @assert zoom <= signal_len(obj2) / sr(obj2) "zoom must be ≤ $(signal_len(obj2) / sr(obj2))."
    @assert nepochs(obj1) == 1 "iplot_ep() should be used for epoched object."
    @assert nepochs(obj2) == 1 "iplot_ep() should be used for epoched object."
    _check_channels(obj1, ch)
    _check_channels(obj2, ch)
    ch_init = ch

    p = NeuroAnalyzer.plot(obj1, obj2, ch=ch)
    win = GtkWindow("NeuroAnalyzer: iplot_cont()", 1200, 800)
    win_view = GtkScrolledWindow()
    set_gtk_property!(win_view, :min_content_width, 1200)
    set_gtk_property!(win_view, :min_content_height, 800)
    set_gtk_property!(win, :border_width, 20)
    set_gtk_property!(win, :resizable, true)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
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
    entry_time = GtkSpinButton(obj1.time_pts[1], obj1.time_pts[end] - zoom, zoom)
    set_gtk_property!(entry_time, :digits, 2)
    set_gtk_property!(entry_time, :value, obj1.time_pts[1])
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

    cb_scale = GtkCheckButton()
    set_gtk_property!(cb_scale, :tooltip_text, "Add scale to the plot")
    set_gtk_property!(cb_scale, :active, true)

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
    lab_scale = GtkLabel("Draw scale:")
    set_gtk_property!(lab_scale, :halign, 2)
    lab_mono = GtkLabel("Grayscale:")
    set_gtk_property!(lab_mono, :halign, 2)
    lab_norm = GtkLabel("Normalize:")
    set_gtk_property!(lab_norm, :halign, 2)
    g_opts[1, 1] = lab_ch
    g_opts[1, 2] = lab_t
    g_opts[1, 3] = lab_x
    g_opts[1, 4] = lab_y
    g_opts[1, 5] = lab_scale
    g_opts[1, 6] = bt_save
    g_opts[2, 1] = entry_ch
    g_opts[2, 2] = entry_title
    g_opts[2, 3] = entry_xlab
    g_opts[2, 4] = entry_ylab
    g_opts[2, 5] = cb_scale
    g_opts[2, 6] = combo_save
    g_opts[1:2, 7] = bt_refresh
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
        ch = get_gtk_property(entry_ch, :text, String)
        if length(ch) < 1
            warn_dialog("Incorrect list of channels.")
            ch = ch_init
            set_gtk_property!(entry_ch, :text, string(ch))
        end
        if (occursin(", ", ch) && _check_svec(ch)) || (occursin(":", ch) && _check_srange(ch)) || _check_sint(ch)
            ch = _s2i(ch)
            if length(ch) == 1 && !in(ch, get_channel_bytype(obj1))
                warn_dialog("Incorrect list of channels.")
                ch = ch_init
                set_gtk_property!(entry_ch, :text, string(ch))
            elseif length(ch) == 1 && !in(ch, get_channel_bytype(obj2))
                warn_dialog("Incorrect list of channels.")
                ch = ch_init
                set_gtk_property!(entry_ch, :text, string(ch))
            elseif length(ch) > 1 && intersect(ch, get_channel_bytype(obj1)) != ch
                warn_dialog("Incorrect list of channels.")
                ch = ch_init
                set_gtk_property!(entry_ch, :text, string(ch))
            elseif length(ch) > 1 && intersect(ch, get_channel_bytype(obj2)) != ch
                warn_dialog("Incorrect list of channels.")
                ch = ch_init
                set_gtk_property!(entry_ch, :text, string(ch))
            end
            title = get_gtk_property(entry_title, :text, String)
            xlab = get_gtk_property(entry_xlab, :text, String)
            ylab = get_gtk_property(entry_ylab, :text, String)
            scale = get_gtk_property(cb_scale, :active, Bool)
            time1 = get_gtk_property(entry_time, :value, Float64)
            time2 = time1 + zoom
            time2 > obj1.time_pts[end] && (time2 = obj1.time_pts[end])
            p = NeuroAnalyzer.plot(obj1,
                                   obj2,
                                   ch=ch,
                                   seg=(time1, time2),
                                   scale=scale,
                                   title=title,
                                   xlabel=xlab,
                                   ylabel=ylab)
            if p.attr[:size][2] + 40 < 800
                Gtk.resize!(win, 1200, p.attr[:size][2] + 40)
                set_gtk_property!(win_view, :min_content_height, p.attr[:size][2])
            else
                Gtk.resize!(win, 1200, 800)
                set_gtk_property!(win_view, :min_content_height, 800)
            end
            ctx = getgc(can)
            show(io, MIME("image/png"), p)
            img = read_from_png(io)
            set_source_surface(ctx, img, 0, 0)
            paint(ctx)
        else
            warn_dialog("Incorrect channels!")
            set_gtk_property!(entry_ch, :text, string(ch_init))
        end
    end

    signal_connect(bt_refresh, "clicked") do widget
        draw(can)
    end
    signal_connect(entry_time, "value-changed") do widget
        draw(can)
    end
    signal_connect(cb_scale, "clicked") do widget
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
        if time_current >= obj1.time_pts[1] + 1
            time_current -= 1
            Gtk.@sigatom begin
                set_gtk_property!(entry_time, :value, time_current)
            end
        end
    end

    signal_connect(bt_prev5, "clicked") do widget
        time_current = get_gtk_property(entry_time, :value, Float64)
        if time_current >= obj1.time_pts[1] + zoom
            time_current = time_current - zoom
            Gtk.@sigatom begin
                set_gtk_property!(entry_time, :value, time_current)
            end
        end
    end

    signal_connect(bt_next, "clicked") do widget
        time_current = get_gtk_property(entry_time, :value, Float64)
        if time_current < obj1.time_pts[end] - zoom
            time_current += 1
            Gtk.@sigatom begin
                set_gtk_property!(entry_time, :value, time_current)
            end
        else
            time_current = obj1.time_pts[end] - zoom
            Gtk.@sigatom begin
                set_gtk_property!(entry_time, :value, time_current)
            end
        end
    end

    signal_connect(bt_next5, "clicked") do widget
        time_current = get_gtk_property(entry_time, :value, Float64)
        if time_current < obj1.time_pts[end] - zoom
            time_current += zoom
            Gtk.@sigatom begin
                set_gtk_property!(entry_time, :value, time_current)
            end
        else
            time_current = obj1.time_pts[end] - zoom
            Gtk.@sigatom begin
                set_gtk_property!(entry_time, :value, time_current)
            end
        end
    end

    signal_connect(bt_start, "clicked") do widget
        Gtk.@sigatom begin
            set_gtk_property!(entry_time, :value, obj1.time_pts[1])
        end
    end

    signal_connect(bt_end, "clicked") do widget
        time_current = obj1.time_pts[end] - zoom
        Gtk.@sigatom begin
            set_gtk_property!(entry_time, :value, time_current)
        end
    end

    signal_connect(bt_close, "clicked") do widget
        Gtk.destroy(win)
    end

    signal_connect(bt_help, "clicked") do widgete
        info_dialog("Keyboard shortcuts:\n\na\tgo to the signal beginning\ns\tgo to the signal end\nz\tgo back by 1 second\nx\tgo forward by 1 second\nc\tgo back by $zoom seconds\nv\tgo forward by $zoom seconds\n\nh\tthis info\nq\texit\n")
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        if k == 113 # q
            Gtk.destroy(win)
        elseif k == 104 # h
            info_dialog("Keyboard shortcuts:\n\na\tgo to the signal beginning\ns\tgo to the signal end\nz\tgo back by 1 second\nx\tgo forward by 1 second\nc\tgo back by $zoom seconds\nv\tgo forward by $zoom seconds\n\nh\tthis info\nq\texit\n")
        elseif k == 97 # a
            Gtk.@sigatom begin
                set_gtk_property!(entry_time, :value, obj1.time_pts[1])
            end
        elseif k == 115 # s
            time_current = obj1.time_pts[end] - zoom
            Gtk.@sigatom begin
                set_gtk_property!(entry_time, :value, time_current)
            end
        elseif k == 122 # z
            time_current = get_gtk_property(entry_time, :value, Float64)
            if time_current >= obj1.time_pts[1] + 1
                time_current -= 1
                Gtk.@sigatom begin
                    set_gtk_property!(entry_time, :value, time_current)
                end
            end
        elseif k == 99 # c
            time_current = get_gtk_property(entry_time, :value, Float64)
            if time_current >= obj1.time_pts[1] + zoom
                time_current = time_current - zoom
                Gtk.@sigatom begin
                    set_gtk_property!(entry_time, :value, time_current)
                end
            end
        elseif k == 120 # x
            time_current = get_gtk_property(entry_time, :value, Float64)
            if time_current < obj1.time_pts[end] - zoom
                time_current += 1
                Gtk.@sigatom begin
                    set_gtk_property!(entry_time, :value, time_current)
                end
            else
                time_current = obj1.time_pts[end] - zoom
                Gtk.@sigatom begin
                    set_gtk_property!(entry_time, :value, time_current)
                end
            end
        elseif k == 118 # v
            time_current = get_gtk_property(entry_time, :value, Float64)
            if time_current < obj1.time_pts[end] - zoom
                time_current += zoom
                Gtk.@sigatom begin
                    set_gtk_property!(entry_time, :value, time_current)
                end
            else
                time_current = obj1.time_pts[end] - zoom
                Gtk.@sigatom begin
                    set_gtk_property!(entry_time, :value, time_current)
                end
            end
        end
    end

    return nothing

end

"""
    iplot_ep(obj1, obj2, ch)

Interactive plot of two epoched signal.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `obj2::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: channel(s) to plot, default is all channels

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function iplot_ep(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=NeuroAnalyzer._c(nchannels(obj1)))

    @assert nepochs(obj1) > 1 "iplot_cont() should be used for continuous object."
    @assert nepochs(obj2) > 1 "iplot_cont() should be used for continuous object."
    _check_channels(obj1, ch)
    _check_channels(obj2, ch)

    p = NeuroAnalyzer.plot(obj1, obj2, ch=ch, ep=1)
    g = GtkGrid()
    g_opts = GtkGrid()
    win = GtkWindow("NeuroAnalyzer: iplot_cont()", 1200, (p.attr[:size][2] + 40))
    set_gtk_property!(win, :border_width, 20)
    set_gtk_property!(win, :resizable, true)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    can = GtkCanvas(Int32(p.attr[:size][1]), Int32(p.attr[:size][2]))
    set_gtk_property!(g, :column_homogeneous, false)
    set_gtk_property!(g_opts, :column_homogeneous, false)
    set_gtk_property!(g, :column_spacing, 10)
    set_gtk_property!(g, :row_spacing, 10)
    set_gtk_property!(g_opts, :row_spacing, 10)
    set_gtk_property!(g_opts, :column_spacing, 10)
    entry_epoch = GtkSpinButton(1, nepochs(obj1), 1)
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

    cb_scale = GtkCheckButton()
    set_gtk_property!(cb_scale, :tooltip_text, "Add scale to the plot")
    set_gtk_property!(cb_scale, :active, true)

    combo_save = GtkComboBoxText()
    file_types = ["PNG", "PDF"]
    for idx in file_types
        push!(combo_save, idx)
    end
    set_gtk_property!(combo_save, :active, 0)
    bt_save = GtkButton("Save as:")

    bt_refresh = GtkButton("Refresh")
    set_gtk_property!(bt_refresh, :tooltip_text, "Refresh the plot")

    lab_ch = GtkLabel("Channels:")
    set_gtk_property!(lab_ch, :halign, 2)
    lab_t = GtkLabel("Title:")
    set_gtk_property!(lab_t, :halign, 2)
    lab_x = GtkLabel("X lab:")
    set_gtk_property!(lab_x, :halign, 2)
    lab_y = GtkLabel("Y lab:")
    set_gtk_property!(lab_y, :halign, 2)
    lab_scale = GtkLabel("Draw scale:")
    set_gtk_property!(lab_scale, :halign, 2)
    g_opts[1, 1] = lab_ch
    g_opts[1, 2] = lab_t
    g_opts[1, 3] = lab_x
    g_opts[1, 4] = lab_y
    g_opts[1, 5] = lab_scale
    g_opts[1, 6] = bt_save
    g_opts[2, 1] = entry_ch
    g_opts[2, 2] = entry_title
    g_opts[2, 3] = entry_xlab
    g_opts[2, 4] = entry_ylab
    g_opts[2, 5] = cb_scale
    g_opts[2, 6] = combo_save
    g_opts[1:2, 7] = bt_refresh
    vbox = GtkBox(:v)
    push!(vbox, g_opts)

    g[1, 1] = vbox
    g[2:9, 1] = can
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
            if length(ch) == 1 && !in(ch, get_channel_bytype(obj1))
                warn_dialog("Incorrect list of channels.")
                ch = ch_init
                set_gtk_property!(entry_ch, :text, string(ch))
            elseif length(ch) == 1 && !in(ch, get_channel_bytype(obj2))
                warn_dialog("Incorrect list of channels.")
                ch = ch_init
                set_gtk_property!(entry_ch, :text, string(ch))
            elseif length(ch) > 1 && intersect(ch, get_channel_bytype(obj1)) != ch
                warn_dialog("Incorrect list of channels.")
                ch = ch_init
                set_gtk_property!(entry_ch, :text, string(ch))
            elseif length(ch) > 1 && intersect(ch, get_channel_bytype(obj2)) != ch
                warn_dialog("Incorrect list of channels.")
                ch = ch_init
                set_gtk_property!(entry_ch, :text, string(ch))
            end
            title = get_gtk_property(entry_title, :text, String)
            xlab = get_gtk_property(entry_xlab, :text, String)
            ylab = get_gtk_property(entry_ylab, :text, String)
            scale = get_gtk_property(cb_scale, :active, Bool)
            ep = get_gtk_property(entry_epoch, :value, Int64)
            p = NeuroAnalyzer.plot(obj1,
                                   obj2,
                                   ep=ep,
                                   ch=ch,
                                   scale=scale,
                                   title=title,
                                   xlabel=xlab,
                                   ylabel=ylab)
            Gtk.resize!(win, 1200, p.attr[:size][2] + 40)
            set_gtk_property!(can, :width_request, Int32(p.attr[:size][1]))
            set_gtk_property!(can, :height_request, Int32(p.attr[:size][2]))
            ctx = getgc(can)
            show(io, MIME("image/png"), NeuroAnalyzer.plot(obj1,
                                                           obj2,
                                                           ep=ep,
                                                           ch=ch,
                                                           scale=scale,
                                                           title=title,
                                                           xlabel=xlab,
                                                           ylabel=ylab))
            img = read_from_png(io)
            set_source_surface(ctx, img, 0, 0)
            paint(ctx)
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
    signal_connect(cb_scale, "clicked") do widget
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
        if ep < nepochs(obj1)
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
            set_gtk_property!(entry_epoch, :value, nepochs(obj1))
        end
    end

    signal_connect(bt_close, "clicked") do widget
        Gtk.destroy(win)
    end

    signal_connect(bt_help, "clicked") do widgete
        info_dialog("Keyboard shortcuts:\na\tgo to first epoch\ns\tgo to last epoch\nz\tprevious epoch\nx\tnext epoch\n\nh\tthis info\nq\texit\n")
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        if k == 113 # q
            Gtk.destroy(win)
        elseif k == 104 # h
            info_dialog("Keyboard shortcuts:\na\tgo to first epoch\ns\tgo to last epoch\nz\tprevious epoch\nx\tnext epoch\n\nh\tthis info\nq\texit\n")
        elseif k == 97 # a
            Gtk.@sigatom begin
                set_gtk_property!(entry_epoch, :value, 1)
            end
        elseif k == 115 # a
            Gtk.@sigatom begin
                set_gtk_property!(entry_epoch, :value, nepochs(obj1))
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
            if ep < nepochs(obj1)
                ep += 1
                Gtk.@sigatom begin
                    set_gtk_property!(entry_epoch, :value, ep)
                end
            end
        end
    end

    return nothing

end