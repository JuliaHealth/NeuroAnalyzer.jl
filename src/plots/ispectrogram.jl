export ispectrogram
export ispectrogram_cont
export ispectrogram_ep

"""
    ispectrogram(obj, ch, zoom)

Interactive spectrogram of continuous or epoched signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1`: index of channels, default is all signal channels
- `zoom::Int64=5`: how many seconds are displayed in one segment

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function ispectrogram(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1, zoom::Int64=5)

    if epoch_n(obj) == 1
        ispectrogram_cont(obj, ch=ch, zoom=zoom)
    else
        ispectrogram_ep(obj, ch=ch)
    end

end

"""
    ispectrogram_cont(obj, ch, zoom)

Interactive spectrogram of continuous signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1`: index of channels, default is all signal channels
- `zoom::Int64=5`: how many seconds are displayed in one segment

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function ispectrogram_cont(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1, zoom::Int64=5)

    @assert zoom >= 1 "zoom must be ≥ 1."
    @assert zoom <= signal_len(obj) / sr(obj) "zoom must be ≤ $(signal_len(obj) / sr(obj))."
    @assert epoch_n(obj) == 1 "ispectrogram_ep() should be used for epoched object."
    _check_channels(obj, ch)
    ch_init = ch

    p = NeuroAnalyzer.plot_spectrogram(obj, ch=ch)
    g = GtkGrid()
    g_opts = GtkGrid()
    win = GtkWindow("NeuroAnalyzer: ispectrogram_cont()", 1200, (p.attr[:size][2] + 40))
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

    cb_norm = GtkCheckButton()
    set_gtk_property!(cb_norm, :tooltip_text, "Normalize powers to dB")
    set_gtk_property!(cb_norm, :active, true)
    
    cb_hw = GtkCheckButton()
    set_gtk_property!(cb_hw, :tooltip_text, "Apply Hanning window")
    set_gtk_property!(cb_hw, :active, true)

    combo_method = GtkComboBoxText()
    psd_methods = ["standard periodogram", "short-time Fourier transform", "multi-taper", "Morlet wavelet"]
    for idx in psd_methods
        push!(combo_method, idx)
    end
    set_gtk_property!(combo_method, :active, 0)
    set_gtk_property!(combo_method, :tooltip_text, "Spectrogram method")

    entry_nt = GtkSpinButton(1, 128, 1)
    set_gtk_property!(entry_nt, :value, 8)
    set_gtk_property!(entry_nt, :tooltip_text, "Number of Slepian tapers")

    entry_ncyc = GtkSpinButton(1, 256, 1)
    set_gtk_property!(entry_ncyc, :value, 6)
    set_gtk_property!(entry_ncyc, :tooltip_text, "Number of Morlet wavelet cycles")

    entry_frq1 = GtkSpinButton(0.0, (sr(obj) / 2) - 0.5, 0.5)
    set_gtk_property!(entry_frq1, :value, 0.0)
    set_gtk_property!(entry_frq1, :tooltip_text, "Start frequency")

    entry_frq2 = GtkSpinButton(0.0, sr(obj) / 2, 0.5)
    set_gtk_property!(entry_frq2, :value, sr(obj) / 2)
    set_gtk_property!(entry_frq2, :tooltip_text, "End frequency")

    entry_wlen = GtkSpinButton(2, zoom * sr(obj) + 1, 1)
    set_gtk_property!(entry_wlen, :value, sr(obj))
    set_gtk_property!(entry_wlen, :tooltip_text, "Window length (samples)")

    entry_woverlap = GtkSpinButton(0, zoom * sr(obj), 1)
    set_gtk_property!(entry_woverlap, :value, round(Int64, sr(obj) * 0.97))
    set_gtk_property!(entry_woverlap, :tooltip_text, "Window overlap (samples)")

    combo_save = GtkComboBoxText()
    file_types = ["PNG", "PDF"]
    for idx in file_types
        push!(combo_save, idx)
    end
    set_gtk_property!(combo_save, :active, 0)
    bt_save = GtkButton("Save as:")

    bt_refresh = GtkButton("Refresh")
    set_gtk_property!(bt_refresh, :tooltip_text, "Refresh the plot")

    lab_method = GtkLabel("Spectrogram method:")
    set_gtk_property!(lab_method, :halign, 2)
    lab_ch = GtkLabel("Channels:")
    set_gtk_property!(lab_ch, :halign, 2)
    lab_t = GtkLabel("Title:")
    set_gtk_property!(lab_t, :halign, 2)
    lab_x = GtkLabel("X lab:")
    set_gtk_property!(lab_x, :halign, 2)
    lab_y = GtkLabel("Y lab:")
    set_gtk_property!(lab_y, :halign, 2)
    lab_mono = GtkLabel("Grayscale:")
    set_gtk_property!(lab_mono, :halign, 2)
    lab_norm = GtkLabel("Normalize:")
    set_gtk_property!(lab_norm, :halign, 2)
    lab_nt = GtkLabel("Slepians:")
    set_gtk_property!(lab_nt, :halign, 2)
    lab_wlen = GtkLabel("Window length:")
    set_gtk_property!(lab_wlen, :halign, 2)
    lab_woverlap = GtkLabel("Window overlap:")
    set_gtk_property!(lab_woverlap, :halign, 2)
    lab_frq1 = GtkLabel("Start frequency:")
    set_gtk_property!(lab_frq1, :halign, 2)
    lab_frq2 = GtkLabel("End frequency:")
    set_gtk_property!(lab_frq2, :halign, 2)
    lab_nc = GtkLabel("Cycles:")
    set_gtk_property!(lab_nc, :halign, 2)
    lab_hw = GtkLabel("Hanning:")
    set_gtk_property!(lab_hw, :halign, 2)
    g_opts[1, 1] = lab_method
    g_opts[1, 2] = lab_ch
    g_opts[1, 3] = lab_t
    g_opts[1, 4] = lab_x
    g_opts[1, 5] = lab_y
    g_opts[1, 6] = lab_frq1
    g_opts[1, 7] = lab_frq2
    g_opts[1, 8] = lab_nc
    g_opts[1, 9] = lab_nt
    g_opts[1, 10] = lab_wlen
    g_opts[1, 11] = lab_woverlap
    g_opts[1, 12] = lab_norm
    g_opts[1, 13] = lab_mono
    g_opts[1, 14] = lab_hw
    g_opts[1, 15] = bt_save
    g_opts[2, 1] = combo_method
    g_opts[2, 2] = entry_ch
    g_opts[2, 3] = entry_title
    g_opts[2, 4] = entry_xlab
    g_opts[2, 5] = entry_ylab
    g_opts[2, 6] = entry_frq1
    g_opts[2, 7] = entry_frq2
    g_opts[2, 8] = entry_ncyc
    g_opts[2, 9] = entry_nt
    g_opts[2, 10] = entry_wlen
    g_opts[2, 11] = entry_woverlap
    g_opts[2, 12] = cb_norm
    g_opts[2, 13] = cb_mono
    g_opts[2, 14] = cb_hw
    g_opts[2, 15] = combo_save
    g_opts[1:2, 16] = bt_refresh
    vbox = GtkBox(:v)
    push!(vbox, g_opts)

    g[1, 1] = vbox
    g[2:11, 1] = can
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
            hw = get_gtk_property(cb_hw, :active, Bool)
            method = get_gtk_property(combo_method, :active, String)
            method == "0" && (method = :standard)
            method == "1" && (method = :stft)
            method == "2" && (method = :mt)
            method == "3" && (method = :mw)
            frq1 = get_gtk_property(entry_frq1, :value, Float64)
            frq2 = get_gtk_property(entry_frq2, :value, Float64)
            nt = get_gtk_property(entry_nt, :value, Int64)
            ncyc = get_gtk_property(entry_ncyc, :value, Int64)
            wlen = get_gtk_property(entry_wlen, :value, Int64)
            woverlap = get_gtk_property(entry_woverlap, :value, Int64)
            if frq1 == frq2
                warn_dialog("Start and end frequencies must be different.")
            elseif frq1 > frq2
                warn_dialog("Start frequency must be < end frequency.")
            elseif woverlap >= wlen
                warn_dialog("Window overlap must be < window length.")
            else
                frq = (frq1, frq2)
                time1 = get_gtk_property(entry_time, :value, Float64)
                time2 = time1 + zoom
                time2 > obj.time_pts[end] && (time2 = obj.time_pts[end])
                p = NeuroAnalyzer.plot_spectrogram(obj,
                                                   ch=ch,
                                                   seg=(time1, time2),
                                                   mono=mono,
                                                   title=title,
                                                   xlabel=xlab,
                                                   ylabel=ylab,
                                                   norm=norm,
                                                   method=method,
                                                   frq_lim=frq,
                                                   ncyc=ncyc,
                                                   nt=nt,
                                                   wlen=wlen,
                                                   woverlap=woverlap,
                                                   w=hw)
                img = read_from_png(io)
                if typeof(p) == Plots.Plot{Plots.GRBackend}
                    Gtk.resize!(win, 1200, p.attr[:size][2] + 40)
                    set_gtk_property!(can, :width_request, Int32(p.attr[:size][1]))
                    set_gtk_property!(can, :height_request, Int32(p.attr[:size][2]))
                elseif typeof(p) == Makie.Figure
                    Gtk.resize!(win, 1000, 900 + 40)
                    set_gtk_property!(can, :width_request, Int32(900))
                    set_gtk_property!(can, :height_request, Int32(900))
                end
                ctx = getgc(can)
                show(io, MIME("image/png"), NeuroAnalyzer.plot_spectrogram(obj,
                                                                           ch=ch,
                                                                           seg=(time1, time2),
                                                                           mono=mono,
                                                                           title=title,
                                                                           xlabel=xlab,
                                                                           ylabel=ylab,
                                                                           norm=norm,
                                                                           method=method,
                                                                           frq_lim=frq,
                                                                           ncyc=ncyc,
                                                                           nt=nt,
                                                                           wlen=wlen,
                                                                           woverlap=woverlap,
                                                                           w=hw))
                img = read_from_png(io)
                set_source_surface(ctx, img, 0, 0)
                paint(ctx)
            end
        end
    end

    signal_connect(bt_refresh, "clicked") do widget
        draw(can)
    end
    signal_connect(entry_time, "value-changed") do widget
        draw(can)
    end
    signal_connect(combo_method, "changed") do widget
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
    signal_connect(cb_norm, "clicked") do widget
        draw(can)
    end
    signal_connect(cb_mono, "clicked") do widget
        draw(can)
    end
    signal_connect(cb_hw, "clicked") do widget
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
    ispectrogram_ep(obj, ch)

Interactive spectrogram of epoched signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1`: index of channels, default is all signal channels

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function ispectrogram_ep(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1)

    @assert epoch_n(obj) > 1 "ispectrogram_cont() should be used for continuous object."
    _check_channels(obj, ch)

    p = NeuroAnalyzer.plot_spectrogram(obj, ch=ch, ep=1)
    g = GtkGrid()
    g_opts = GtkGrid()
    win = GtkWindow("NeuroAnalyzer: ispectrogram_ep()", 1200, (p.attr[:size][2] + 40))
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
    entry_epoch = GtkSpinButton(1, epoch_n(obj), 1)
    set_gtk_property!(entry_epoch, :tooltip_text, "Epoch")
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

    cb_norm = GtkCheckButton()
    set_gtk_property!(cb_norm, :tooltip_text, "Normalize powers to dB")
    set_gtk_property!(cb_norm, :active, true)
    
    cb_hw = GtkCheckButton()
    set_gtk_property!(cb_hw, :tooltip_text, "Apply Hanning window")
    set_gtk_property!(cb_hw, :active, true)

    combo_method = GtkComboBoxText()
    psd_methods = ["standard periodogram", "short-time Fourier transform", "multi-taper", "Morlet wavelet"]
    for idx in psd_methods
        push!(combo_method, idx)
    end
    set_gtk_property!(combo_method, :active, 0)
    set_gtk_property!(combo_method, :tooltip_text, "Spectrogram method")

    entry_nt = GtkSpinButton(1, 128, 1)
    set_gtk_property!(entry_nt, :value, 8)
    set_gtk_property!(entry_nt, :tooltip_text, "Number of Slepian tapers")

    entry_ncyc = GtkSpinButton(1, 256, 1)
    set_gtk_property!(entry_ncyc, :value, 6)
    set_gtk_property!(entry_ncyc, :tooltip_text, "Number of Morlet wavelet cycles")

    entry_frq1 = GtkSpinButton(0.0, (sr(obj) / 2) - 0.5, 0.5)
    set_gtk_property!(entry_frq1, :value, 0.0)
    set_gtk_property!(entry_frq1, :tooltip_text, "Start frequency")

    entry_frq2 = GtkSpinButton(0.0, sr(obj) / 2, 0.5)
    set_gtk_property!(entry_frq2, :value, sr(obj) / 2)
    set_gtk_property!(entry_frq2, :tooltip_text, "End frequency")

    entry_wlen = GtkSpinButton(2, zoom * sr(obj) + 1, 1)
    set_gtk_property!(entry_wlen, :value, sr(obj))
    set_gtk_property!(entry_wlen, :tooltip_text, "Window length (samples)")

    entry_woverlap = GtkSpinButton(0, zoom * sr(obj), 1)
    set_gtk_property!(entry_woverlap, :value, round(Int64, sr(obj) * 0.97))
    set_gtk_property!(entry_woverlap, :tooltip_text, "Window overlap (samples)")

    combo_save = GtkComboBoxText()
    file_types = ["PNG", "PDF"]
    for idx in file_types
        push!(combo_save, idx)
    end
    set_gtk_property!(combo_save, :active, 0)
    bt_save = GtkButton("Save as:")

    bt_refresh = GtkButton("Refresh")
    set_gtk_property!(bt_refresh, :tooltip_text, "Refresh the plot")

    lab_method = GtkLabel("Spectrogram method:")
    set_gtk_property!(lab_method, :halign, 2)
    lab_ch = GtkLabel("Channels:")
    set_gtk_property!(lab_ch, :halign, 2)
    lab_t = GtkLabel("Title:")
    set_gtk_property!(lab_t, :halign, 2)
    lab_x = GtkLabel("X lab:")
    set_gtk_property!(lab_x, :halign, 2)
    lab_y = GtkLabel("Y lab:")
    set_gtk_property!(lab_y, :halign, 2)
    lab_mono = GtkLabel("Grayscale:")
    set_gtk_property!(lab_mono, :halign, 2)
    lab_norm = GtkLabel("Normalize:")
    set_gtk_property!(lab_norm, :halign, 2)
    lab_nt = GtkLabel("Slepians:")
    set_gtk_property!(lab_nt, :halign, 2)
    lab_wlen = GtkLabel("Window length:")
    set_gtk_property!(lab_wlen, :halign, 2)
    lab_woverlap = GtkLabel("Window overlap:")
    set_gtk_property!(lab_woverlap, :halign, 2)
    lab_frq1 = GtkLabel("Start frequency:")
    set_gtk_property!(lab_frq1, :halign, 2)
    lab_frq2 = GtkLabel("End frequency:")
    set_gtk_property!(lab_frq2, :halign, 2)
    lab_nc = GtkLabel("Cycles:")
    set_gtk_property!(lab_nc, :halign, 2)
    lab_hw = GtkLabel("Hanning:")
    set_gtk_property!(lab_hw, :halign, 2)
    g_opts[1, 1] = lab_method
    g_opts[1, 2] = lab_ch
    g_opts[1, 3] = lab_t
    g_opts[1, 4] = lab_x
    g_opts[1, 5] = lab_y
    g_opts[1, 6] = lab_frq1
    g_opts[1, 7] = lab_frq2
    g_opts[1, 8] = lab_nc
    g_opts[1, 9] = lab_nt
    g_opts[1, 10] = lab_wlen
    g_opts[1, 11] = lab_woverlap
    g_opts[1, 12] = lab_norm
    g_opts[1, 13] = lab_mono
    g_opts[1, 14] = lab_hw
    g_opts[1, 15] = bt_save
    g_opts[2, 1] = combo_method
    g_opts[2, 2] = entry_ch
    g_opts[2, 3] = entry_title
    g_opts[2, 4] = entry_xlab
    g_opts[2, 5] = entry_ylab
    g_opts[2, 6] = entry_frq1
    g_opts[2, 7] = entry_frq2
    g_opts[2, 8] = entry_ncyc
    g_opts[2, 9] = entry_nt
    g_opts[2, 10] = entry_wlen
    g_opts[2, 11] = entry_woverlap
    g_opts[2, 12] = cb_norm
    g_opts[2, 13] = cb_mono
    g_opts[2, 14] = cb_hw
    g_opts[2, 15] = combo_save
    g_opts[1:2, 16] = bt_refresh
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
            hw = get_gtk_property(cb_hw, :active, Bool)
            method = get_gtk_property(combo_method, :active, String)
            method == "0" && (method = :standard)
            method == "1" && (method = :stft)
            method == "2" && (method = :mt)
            method == "3" && (method = :mw)
            frq1 = get_gtk_property(entry_frq1, :value, Float64)
            frq2 = get_gtk_property(entry_frq2, :value, Float64)
            nt = get_gtk_property(entry_nt, :value, Int64)
            ncyc = get_gtk_property(entry_ncyc, :value, Int64)
            wlen = get_gtk_property(entry_wlen, :value, Int64)
            woverlap = get_gtk_property(entry_woverlap, :value, Int64)
            if frq1 == frq2
                warn_dialog("Start and end frequencies must be different.")
            elseif frq1 > frq2
                warn_dialog("Start frequency must be < end frequency.")
            elseif woverlap >= wlen
                warn_dialog("Window overlap must be < window length.")
            else
                frq = (frq1, frq2)
                ep = get_gtk_property(entry_epoch, :value, Int64)
                p = NeuroAnalyzer.plot_spectrogram(obj,
                                                   ch=ch,
                                                   ep=ep,
                                                   mono=mono,
                                                   title=title,
                                                   xlabel=xlab,
                                                   ylabel=ylab,
                                                   norm=norm,
                                                   method=method,
                                                   frq_lim=frq,
                                                   ncyc=ncyc,
                                                   nt=nt,
                                                   wlen=wlen,
                                                   woverlap=woverlap,
                                                   w=hw)
                img = read_from_png(io)
                if typeof(p) == Plots.Plot{Plots.GRBackend}
                    Gtk.resize!(win, 1200, p.attr[:size][2] + 40)
                    set_gtk_property!(can, :width_request, Int32(p.attr[:size][1]))
                    set_gtk_property!(can, :height_request, Int32(p.attr[:size][2]))
                elseif typeof(p) == Makie.Figure
                    Gtk.resize!(win, 1000, 900 + 40)
                    set_gtk_property!(can, :width_request, Int32(900))
                    set_gtk_property!(can, :height_request, Int32(900))
                end
                ctx = getgc(can)
                show(io, MIME("image/png"), NeuroAnalyzer.plot_spectrogram(obj,
                                                                           ch=ch,
                                                                           ep=ep,
                                                                           mono=mono,
                                                                           title=title,
                                                                           xlabel=xlab,
                                                                           ylabel=ylab,
                                                                           norm=norm,
                                                                           method=method,
                                                                           frq_lim=frq,
                                                                           ncyc=ncyc,
                                                                           nt=nt,
                                                                           wlen=wlen,
                                                                           woverlap=woverlap,
                                                                           w=hw))
                img = read_from_png(io)
                set_source_surface(ctx, img, 0, 0)
                paint(ctx)
            end
        end
    end

    signal_connect(bt_refresh, "clicked") do widget
        draw(can)
    end
    signal_connect(entry_epoch, "value-changed") do widget
        draw(can)
    end
    signal_connect(combo_method, "changed") do widget
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
    signal_connect(cb_norm, "clicked") do widget
        draw(can)
    end
    signal_connect(cb_mono, "clicked") do widget
        draw(can)
    end
    signal_connect(cb_hw, "clicked") do widget
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
        if ep < epoch_n(obj)
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
            set_gtk_property!(entry_epoch, :value, epoch_n(obj))
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
        elseif k == 115 # s
            Gtk.@sigatom begin
                set_gtk_property!(entry_epoch, :value, epoch_n(obj))
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
            if ep < epoch_n(obj)
                ep += 1
                Gtk.@sigatom begin
                    set_gtk_property!(entry_epoch, :value, ep)
                end
            end
        end
    end

    return nothing

end