export ipsd
export ipsd_cont
export ipsd_ep

"""
    ipsd(obj; <keyword arguments>)

Interactive PSD of continuous or epoched signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1`: index of channels, default is the first channel
- `zoom::Real=5`: how many seconds are displayed in one segment
"""
function ipsd(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1, zoom::Real=5)

    if nepochs(obj) == 1
        ipsd_cont(obj, ch=ch, zoom=zoom)
    else
        ipsd_ep(obj, ch=ch)
    end

    return nothing

end

"""
    ipsd_cont(obj, ch, zoom)

Interactive PSD of continuous signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1`: index of channels, default is the first channel
- `zoom::Real=5`: how many seconds are displayed in one segment
"""
function ipsd_cont(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1, zoom::Real=5)

    @assert zoom > 0 "zoom must be > 0."
    @assert zoom <= signal_len(obj) / sr(obj) "zoom must be ≤ $(signal_len(obj) / sr(obj))."
    @assert nepochs(obj) == 1 "ipsd_ep() should be used for epoched object."
    _check_channels(obj, ch)

    p = NeuroAnalyzer.plot_psd(obj, ch=ch)
    g = GtkGrid()
    g_opts = GtkGrid()
    win = GtkWindow("NeuroAnalyzer: plot_psd()", 1200, (p.attr[:size][2] + 40))
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
    psd_methods = ["Welch periodogram", "fast Fourier transform", "short-time Fourier transform", "multi-taper", "Morlet wavelet"]
    for idx in psd_methods
        push!(combo_method, idx)
    end
    set_gtk_property!(combo_method, :active, 0)
    set_gtk_property!(combo_method, :tooltip_text, "PSD method")

    combo_type = GtkComboBoxText()
    psd_types = ["normal", "butterfly", "mean", "w3d", "s3d", "topo"]
    for idx in psd_types
        push!(combo_type, idx)
    end
    set_gtk_property!(combo_type, :active, 0)
    set_gtk_property!(combo_type, :tooltip_text, "PSD type")

    combo_ref = GtkComboBoxText()
    ref_types = ["absolute", "total power", "delta", "theta", "alpha", "beta", "beta high", "gamma", "gamma 1", "gamma 2", "gamma lower", "gamma higher"]
    for idx in ref_types
        push!(combo_ref, idx)
    end
    set_gtk_property!(combo_ref, :active, 0)
    set_gtk_property!(combo_ref, :tooltip_text, "Powers referenced to")

    combo_ax = GtkComboBoxText()
    ax_types = ["linear-linear", "log10-linear", "linear-log10", "log10-log10"]
    for idx in ax_types
        push!(combo_ax, idx)
    end
    set_gtk_property!(combo_ax, :active, 0)
    set_gtk_property!(combo_ax, :tooltip_text, "Axes scaling (X-Y)")

    entry_nt = GtkSpinButton(1, 128, 1)
    set_gtk_property!(entry_nt, :value, 8)
    set_gtk_property!(entry_nt, :tooltip_text, "Number of Slepian tapers (for MT)")

    entry_ncyc = GtkSpinButton(1, 256, 1)
    set_gtk_property!(entry_ncyc, :value, 32)
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

    cb_frq = GtkCheckButton()
    set_gtk_property!(cb_frq, :tooltip_text, "Linear or logarithmic frequencies")
    set_gtk_property!(cb_frq, :active, true)
    
    entry_frq = GtkSpinButton(1, 2048, 1)
    set_gtk_property!(entry_frq, :value, _tlength((0.0, sr(obj) / 2)))
    set_gtk_property!(entry_frq, :tooltip_text, "Number of frequencies")

    lab_type = GtkLabel("PSD method:")
    set_gtk_property!(lab_type, :halign, 2)
    lab_method = GtkLabel("PSD method:")
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
    lab_ref = GtkLabel("PSD reference:")
    set_gtk_property!(lab_ref, :halign, 2)
    lab_ax = GtkLabel("Axes scaling:")
    set_gtk_property!(lab_ax, :halign, 2)
    lab_hw = GtkLabel("Hanning:")
    set_gtk_property!(lab_hw, :halign, 2)
    lab_frq = GtkLabel("Linear frequencies:")
    set_gtk_property!(lab_frq, :halign, 2)
    lab_frqn = GtkLabel("Frequencies:")
    set_gtk_property!(lab_frqn, :halign, 2)
    g_opts[1, 1] = lab_type
    g_opts[1, 2] = lab_method
    g_opts[1, 3] = lab_ref
    g_opts[1, 4] = lab_ax
    g_opts[1, 5] = lab_ch
    g_opts[1, 6] = lab_t
    g_opts[1, 7] = lab_x
    g_opts[1, 8] = lab_y
    g_opts[1, 9] = lab_frq1
    g_opts[1, 10] = lab_frq2
    g_opts[1, 11] = lab_nc
    g_opts[1, 12] = lab_nt
    g_opts[1, 13] = lab_wlen
    g_opts[1, 14] = lab_woverlap
    g_opts[1, 15] = lab_frqn
    g_opts[1, 16] = lab_norm
    g_opts[1, 17] = lab_mono
    g_opts[1, 18] = lab_hw
    g_opts[1, 19] = lab_frq
    g_opts[1, 20] = bt_save
    g_opts[2, 1] = combo_type
    g_opts[2, 2] = combo_method
    g_opts[2, 3] = combo_ref
    g_opts[2, 4] = combo_ax
    g_opts[2, 5] = entry_ch
    g_opts[2, 6] = entry_title
    g_opts[2, 7] = entry_xlab
    g_opts[2, 8] = entry_ylab
    g_opts[2, 9] = entry_frq1
    g_opts[2, 10] = entry_frq2
    g_opts[2, 11] = entry_ncyc
    g_opts[2, 12] = entry_nt
    g_opts[2, 13] = entry_wlen
    g_opts[2, 14] = entry_woverlap
    g_opts[2, 15] = entry_frq
    g_opts[2, 16] = cb_norm
    g_opts[2, 17] = cb_mono
    g_opts[2, 18] = cb_hw
    g_opts[2, 19] = cb_frq
    g_opts[2, 20] = combo_save
    g_opts[1:2, 21] = bt_refresh
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
            f = get_gtk_property(cb_frq, :active, Bool) ? :lin : :log
            method = get_gtk_property(combo_method, :active, String)
            method == "0" && (method = :welch)
            method == "1" && (method = :fft)
            method == "2" && (method = :stft)
            method == "3" && (method = :mt)
            method == "4" && (method = :mw)
            type = get_gtk_property(combo_type, :active, String)
            type == "0" && (type = :normal)
            type == "1" && (type = :butterfly)
            type == "2" && (type = :mean)
            type == "3" && (type = :w3d)
            type == "4" && (type = :s3d)
            type == "5" && (type = :topo)
            ref = get_gtk_property(combo_ref, :active, String)
            ref == "0" && (ref = :abs)
            ref == "1" && (ref = :total)
            ref == "2" && (ref = :delta)
            ref == "3" && (ref = :theta)
            ref == "4" && (ref = :alpha)
            ref == "5" && (ref = :beta)
            ref == "6" && (ref = :beta_high)
            ref == "7" && (ref = :gamma)
            ref == "8" && (ref = :gamma_1)
            ref == "9" && (ref = :gamma_2)
            ref == "10" && (ref = :gamma_lower)
            ref == "11" && (ref = :gamma_higher)
            ax = get_gtk_property(combo_ax, :active, String)
            ax == "0" && (ax = :linlin)
            ax == "1" && (ax = :loglin)
            ax == "2" && (ax = :linlog)
            ax == "3" && (ax = :loglog)
            frq1 = get_gtk_property(entry_frq1, :value, Float64)
            frq2 = get_gtk_property(entry_frq2, :value, Float64)
            nt = get_gtk_property(entry_nt, :value, Int64)
            ncyc = get_gtk_property(entry_ncyc, :value, Int64)
            wlen = get_gtk_property(entry_wlen, :value, Int64)
            woverlap = get_gtk_property(entry_woverlap, :value, Int64)
            frq_n = get_gtk_property(entry_frq, :value, Int64)
            if frq1 == frq2
                warn_dialog("Start and end frequencies must be different.")
            elseif frq1 > frq2
                warn_dialog("Start frequency must be < end frequency.")
            elseif woverlap >= wlen
                warn_dialog("Window overlap must be < window length.")
            elseif length(ch) < 2 && type === :butterfly
                warn_dialog("For butterfly plot, the signal must contain ≥ 2 channels.")
            elseif length(ch) < 2 && type === :mean
                warn_dialog("For mean plot, the signal must contain ≥ 2 channels.")
            elseif length(ch) < 2 && type === :w3d
                warn_dialog("For w3d plot, the signal must contain ≥ 2 channels.")
            elseif length(ch) < 2 && type === :s3d
                warn_dialog("For s3d plot, the signal must contain ≥ 2 channels.")
            elseif _has_locs(obj) == false && type === :topo
                warn_dialog("Electrode locations not available.")
            elseif length(unique(obj.header.recording[:channel_type][ch])) > 1 && type in [:butterfly, :mean, :w3d, :s3d, :topo]
                warn_dialog("For $(string(type)) plot\nall channels should be of the same type.")
            else
                frq = (frq1, frq2)
                time1 = get_gtk_property(entry_time, :value, Float64)
                time2 = time1 + zoom
                time2 > obj.time_pts[end] && (time2 = obj.time_pts[end])
                p = NeuroAnalyzer.plot_psd(obj,
                                           ch=ch,
                                           seg=(time1, time2),
                                           mono=mono,
                                           title=title,
                                           xlabel=xlab,
                                           ylabel=ylab,
                                           norm=norm,
                                           method=method,
                                           type=type,
                                           ax=ax,
                                           ref=ref,
                                           frq_lim=frq,
                                           frq=f,
                                           frq_n=frq_n,
                                           ncyc=ncyc,
                                           nt=nt,
                                           wlen=wlen,
                                           woverlap=woverlap,
                                           w=hw)
                img = read_from_png(io)
                Gtk.resize!(win, 1200, p.attr[:size][2] + 40)
                set_gtk_property!(can, :width_request, Int32(p.attr[:size][1]))
                set_gtk_property!(can, :height_request, Int32(p.attr[:size][2]))
                ctx = getgc(can)
                show(io, MIME("image/png"), p)
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
    signal_connect(combo_type, "changed") do widget
        draw(can)
    end
    signal_connect(combo_method, "changed") do widget
        draw(can)
    end
    signal_connect(combo_ref, "changed") do widget
        draw(can)
    end
    signal_connect(combo_ax, "changed") do widget
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
    signal_connect(entry_frq, "value-changed") do widget
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
    signal_connect(cb_frq, "clicked") do widget
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
    ipsd_ep(obj, ch)

Interactive PSD of epoched signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1`: index of channels, default is the first channel
"""
function ipsd_ep(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1)

    @assert nepochs(obj) > 1 "ipsd_cont() should be used for continuous object."
    _check_channels(obj, ch)

    p = NeuroAnalyzer.plot_psd(obj, ch=ch)
    g = GtkGrid()
    g_opts = GtkGrid()
    win = GtkWindow("NeuroAnalyzer: plot_psd()", 1200, (p.attr[:size][2] + 40))
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

    cb_hw = GtkCheckButton()
    set_gtk_property!(cb_hw, :tooltip_text, "Apply Hanning window")
    set_gtk_property!(cb_hw, :active, true)

    cb_norm = GtkCheckButton()
    set_gtk_property!(cb_norm, :tooltip_text, "Normalize powers to dB")
    set_gtk_property!(cb_norm, :active, true)

    combo_method = GtkComboBoxText()
    psd_methods = ["Welch periodogram", "fast Fourier transform", "short-time Fourier transform", "multi-taper", "Morlet wavelet"]
    for idx in psd_methods
        push!(combo_method, idx)
    end
    set_gtk_property!(combo_method, :active, 0)
    set_gtk_property!(combo_method, :tooltip_text, "PSD method")

    combo_type = GtkComboBoxText()
    psd_types = ["normal", "butterfly", "mean", "w3d", "s3d", "topo"]
    for idx in psd_types
        push!(combo_type, idx)
    end
    set_gtk_property!(combo_type, :active, 0)
    set_gtk_property!(combo_type, :tooltip_text, "PSD type")

    combo_ref = GtkComboBoxText()
    ref_types = ["absolute", "total power", "delta", "theta", "alpha", "beta", "beta high", "gamma", "gamma 1", "gamma 2", "gamma lower", "gamma higher"]
    for idx in ref_types
        push!(combo_ref, idx)
    end
    set_gtk_property!(combo_ref, :active, 0)
    set_gtk_property!(combo_ref, :tooltip_text, "Powers referenced to")

    combo_ax = GtkComboBoxText()
    ref_types = ["linear-linear", "log10-linear", "linear-log10", "log10-log10"]
    for idx in ref_types
        push!(combo_ax, idx)
    end
    set_gtk_property!(combo_ax, :active, 0)
    set_gtk_property!(combo_ax, :tooltip_text, "Axes scaling (X-Y)")

    entry_nt = GtkSpinButton(1, 128, 1)
    set_gtk_property!(entry_nt, :value, 8)
    set_gtk_property!(entry_nt, :tooltip_text, "Number of Slepian tapers")

    entry_ncyc = GtkSpinButton(1, 256, 1)
    set_gtk_property!(entry_ncyc, :value, 32)
    set_gtk_property!(entry_ncyc, :tooltip_text, "Number of Morlet wavelet cycles")

    entry_frq1 = GtkSpinButton(0.0, (sr(obj) / 2) - 0.5, 0.5)
    set_gtk_property!(entry_frq1, :value, 0.0)
    set_gtk_property!(entry_frq1, :tooltip_text, "Start frequency")

    entry_frq2 = GtkSpinButton(0.0, sr(obj) / 2, 0.5)
    set_gtk_property!(entry_frq2, :value, sr(obj) / 2)
    set_gtk_property!(entry_frq2, :tooltip_text, "End frequency")

    entry_wlen = GtkSpinButton(2, epoch_len(obj) * sr(obj) + 1, 1)
    set_gtk_property!(entry_wlen, :value, sr(obj))
    set_gtk_property!(entry_wlen, :tooltip_text, "Window length (samples)")

    entry_woverlap = GtkSpinButton(0, epoch_len(obj) * sr(obj), 1)
    set_gtk_property!(entry_woverlap, :value, round(Int64, sr(obj) * 0.97))
    set_gtk_property!(entry_woverlap, :tooltip_text, "Window overlap (samples)")

    cb_frq = GtkCheckButton()
    set_gtk_property!(cb_frq, :tooltip_text, "Linear or logarithmic frequencies")
    set_gtk_property!(cb_frq, :active, true)
    
    entry_frq = GtkSpinButton(1, 2048, 1)
    set_gtk_property!(entry_frq, :value, _tlength((0.0, sr(obj) / 2)))
    set_gtk_property!(entry_frq, :tooltip_text, "Number of frequencies")

    combo_save = GtkComboBoxText()
    file_types = ["PNG", "PDF"]
    for idx in file_types
        push!(combo_save, idx)
    end
    set_gtk_property!(combo_save, :active, 0)
    bt_save = GtkButton("Save as:")

    bt_refresh = GtkButton("Refresh")
    set_gtk_property!(bt_refresh, :tooltip_text, "Refresh the plot")

    lab_type = GtkLabel("PSD method:")
    set_gtk_property!(lab_type, :halign, 2)
    lab_method = GtkLabel("PSD method:")
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
    lab_ref = GtkLabel("PSD reference:")
    set_gtk_property!(lab_ref, :halign, 2)
    lab_ax = GtkLabel("Axes scaling:")
    set_gtk_property!(lab_ax, :halign, 2)
    lab_hw = GtkLabel("Hanning:")
    set_gtk_property!(lab_hw, :halign, 2)
    lab_frq = GtkLabel("Linear frequencies:")
    set_gtk_property!(lab_frq, :halign, 2)
    lab_frqn = GtkLabel("Frequencies:")
    set_gtk_property!(lab_frqn, :halign, 2)
    g_opts[1, 1] = lab_type
    g_opts[1, 2] = lab_method
    g_opts[1, 3] = lab_ref
    g_opts[1, 4] = lab_ax
    g_opts[1, 5] = lab_ch
    g_opts[1, 6] = lab_t
    g_opts[1, 7] = lab_x
    g_opts[1, 8] = lab_y
    g_opts[1, 9] = lab_frq1
    g_opts[1, 10] = lab_frq2
    g_opts[1, 11] = lab_nc
    g_opts[1, 12] = lab_nt
    g_opts[1, 13] = lab_wlen
    g_opts[1, 14] = lab_woverlap
    g_opts[1, 15] = lab_frqn
    g_opts[1, 16] = lab_norm
    g_opts[1, 17] = lab_mono
    g_opts[1, 18] = lab_hw
    g_opts[1, 19] = lab_frq
    g_opts[1, 20] = bt_save
    g_opts[2, 1] = combo_type
    g_opts[2, 2] = combo_method
    g_opts[2, 3] = combo_ref
    g_opts[2, 4] = combo_ax
    g_opts[2, 5] = entry_ch
    g_opts[2, 6] = entry_title
    g_opts[2, 7] = entry_xlab
    g_opts[2, 8] = entry_ylab
    g_opts[2, 9] = entry_frq1
    g_opts[2, 10] = entry_frq2
    g_opts[2, 11] = entry_ncyc
    g_opts[2, 12] = entry_nt
    g_opts[2, 13] = entry_wlen
    g_opts[2, 14] = entry_woverlap
    g_opts[2, 15] = entry_frq
    g_opts[2, 16] = cb_norm
    g_opts[2, 17] = cb_mono
    g_opts[2, 18] = cb_hw
    g_opts[2, 19] = cb_frq
    g_opts[2, 20] = combo_save
    g_opts[1:2, 21] = bt_refresh
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
            f = get_gtk_property(cb_frq, :active, Bool) ? :lin : :log
            method = get_gtk_property(combo_method, :active, String)
            method == "0" && (method = :welch)
            method == "1" && (method = :fft)
            method == "2" && (method = :stft)
            method == "3" && (method = :mt)
            method == "4" && (method = :mw)
            type = get_gtk_property(combo_type, :active, String)
            type == "0" && (type = :normal)
            type == "1" && (type = :butterfly)
            type == "2" && (type = :mean)
            type == "3" && (type = :w3d)
            type == "4" && (type = :s3d)
            type == "5" && (type = :topo)
            ref = get_gtk_property(combo_ref, :active, String)
            ref == "0" && (ref = :abs)
            ref == "1" && (ref = :total)
            ref == "2" && (ref = :delta)
            ref == "3" && (ref = :theta)
            ref == "4" && (ref = :alpha)
            ref == "5" && (ref = :beta)
            ref == "6" && (ref = :beta_high)
            ref == "7" && (ref = :gamma)
            ref == "8" && (ref = :gamma_1)
            ref == "9" && (ref = :gamma_2)
            ref == "10" && (ref = :gamma_lower)
            ref == "11" && (ref = :gamma_higher)
            ax = get_gtk_property(combo_ax, :active, String)
            ax == "0" && (ax = :linlin)
            ax == "1" && (ax = :loglin)
            ax == "2" && (ax = :linlog)
            ax == "3" && (ax = :loglog)
            frq1 = get_gtk_property(entry_frq1, :value, Float64)
            frq2 = get_gtk_property(entry_frq2, :value, Float64)
            nt = get_gtk_property(entry_nt, :value, Int64)
            ncyc = get_gtk_property(entry_ncyc, :value, Int64)
            wlen = get_gtk_property(entry_wlen, :value, Int64)
            woverlap = get_gtk_property(entry_woverlap, :value, Int64)
            frq_n = get_gtk_property(entry_frq, :value, Int64)
            if frq1 == frq2
                warn_dialog("Start and end frequencies must be different.")
            elseif frq1 > frq2
                warn_dialog("Start frequency must be < end frequency.")
            elseif woverlap >= wlen
                warn_dialog("Window overlap must be < window length.")
            elseif ref !== :abs && method === :mw
                warn_dialog("For relative plot, method must be Welch, STFT or multi-taper.")
            elseif length(ch) < 2 && type === :butterfly
                warn_dialog("For butterfly plot, the signal must contain ≥ 2 channels.")
            elseif length(ch) < 2 && type === :mean
                warn_dialog("For mean plot, the signal must contain ≥ 2 channels.")
            elseif length(ch) < 2 && type === :w3d
                warn_dialog("For w3d plot, the signal must contain ≥ 2 channels.")
            elseif length(ch) < 2 && type === :s3d
                warn_dialog("For s3d plot, the signal must contain ≥ 2 channels.")
            elseif _has_locs(obj) == false && type === :topo
                warn_dialog("Electrode locations not available.")
            elseif length(unique(obj.header.recording[:channel_type][ch])) > 1 && type in [:butterfly, :mean, :w3d, :s3d, :topo]
                warn_dialog("For $(string(type)) plot\nall channels should be of the same type.")
            else
                frq = (frq1, frq2)
                ep = get_gtk_property(entry_epoch, :value, Int64)
                p = NeuroAnalyzer.plot_psd(obj,
                                           ch=ch,
                                           ep=ep,
                                           mono=mono,
                                           title=title,
                                           xlabel=xlab,
                                           ylabel=ylab,
                                           norm=norm,
                                           method=method,
                                           type=type,
                                           ax=ax,
                                           ref=ref,
                                           frq_lim=frq,
                                           frq=f,
                                           frq_n=frq_n,
                                           ncyc=ncyc,
                                           nt=nt,
                                           wlen=wlen,
                                           woverlap=woverlap,
                                           w=hw)
                img = read_from_png(io)
                Gtk.resize!(win, 1200, p.attr[:size][2] + 40)
                set_gtk_property!(can, :width_request, Int32(p.attr[:size][1]))
                set_gtk_property!(can, :height_request, Int32(p.attr[:size][2]))
                ctx = getgc(can)
                show(io, MIME("image/png"), p)
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
    signal_connect(combo_type, "changed") do widget
        draw(can)
    end
    signal_connect(combo_method, "changed") do widget
        draw(can)
    end
    signal_connect(combo_ref, "changed") do widget
        draw(can)
    end
    signal_connect(combo_ax, "changed") do widget
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
    signal_connect(entry_frq, "value-changed") do widget
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
    signal_connect(cb_frq, "clicked") do widget
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