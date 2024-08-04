export ispectrogram
export ispectrogram_cont
export ispectrogram_ep

"""
    ispectrogram(obj; <keyword arguments>)

Interactive spectrogram of continuous or epoched signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::String`: channel name
- `zoom::Real=10`: how many seconds are displayed in one segment
"""
function ispectrogram(obj::NeuroAnalyzer.NEURO; ch::String, zoom::Real=10)

    if nepochs(obj) == 1
        ispectrogram_cont(obj, ch=ch, zoom=zoom)
    else
        ispectrogram_ep(obj, ch=ch)
    end

    return nothing

end

"""
    ispectrogram_cont(obj; <keyword arguments>)

Interactive spectrogram of continuous signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::String`: channel name
- `zoom::Real=10`: how many seconds are displayed in one segment
"""
function ispectrogram_cont(obj::NeuroAnalyzer.NEURO; ch::String, zoom::Real=10)

    @assert zoom > 0 "zoom must be > 0."
    @assert zoom <= signal_len(obj) / sr(obj) "zoom must be ≤ $(signal_len(obj) / sr(obj))."
    @assert nepochs(obj) == 1 "ispectrogram_ep() must be used for epoched object."

    ch_init = ch
    ch = get_channel(obj, ch=ch)
    clabels = labels(obj)

    p = NeuroAnalyzer.plot_spectrogram(obj, ch=clabels[ch])
    win = GtkWindow("NeuroAnalyzer: ispectrogram_cont()", 1200, 850)
    set_gtk_property!(win, :border_width, 5)
    set_gtk_property!(win, :resizable, true)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    set_gtk_property!(win, :startup_id, "org.neuroanalyzer")
    can = GtkCanvas(Int32(p.attr[:size][1]), Int32(p.attr[:size][2]))
    g = GtkGrid()
    g_opts = GtkGrid()
    set_gtk_property!(g, :column_homogeneous, false)
    set_gtk_property!(g_opts, :column_homogeneous, false)
    set_gtk_property!(g, :column_spacing, 5)
    set_gtk_property!(g, :row_spacing, 5)
    set_gtk_property!(g_opts, :row_spacing, 5)
    set_gtk_property!(g_opts, :column_spacing, 5)

    entry_time = GtkSpinButton(obj.time_pts[1], obj.time_pts[end] - zoom, 1)
    set_gtk_property!(entry_time, :digits, 2)
    set_gtk_property!(entry_time, :value, obj.time_pts[1])
    set_gtk_property!(entry_time, :tooltip_text, "Time position [s]")
    bt_start = GtkButton("⇤")
    set_gtk_property!(bt_start, :tooltip_text, "Go to the signal beginning")
    bt_prev5 = GtkButton("←")
    set_gtk_property!(bt_prev5, :tooltip_text, "Go back by $(round(zoom)) seconds")
    bt_next5 = GtkButton("→")
    set_gtk_property!(bt_next5, :tooltip_text, "Go forward by $(round(zoom)) seconds")
    bt_end = GtkButton("⇥")
    set_gtk_property!(bt_end, :tooltip_text, "Go to the signal end")
    bt_help = GtkButton("Help")
    set_gtk_property!(bt_help, :tooltip_text, "Show help")
    bt_close = GtkButton("Close")
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

    combo_ch = GtkComboBoxText()
    ctypes = uppercase.(unique(obj.header.recording[:channel_type]))
    ch_types = [ctypes; clabels]
    for idx in ch_types
        push!(combo_ch, idx)
    end
    if ch_init in lowercase.(ctypes)
        n = findfirst(isequal(ch_init), lowercase.(ctypes))
        set_gtk_property!(combo_ch, :active, n - 1)
    else
        set_gtk_property!(combo_ch, :active, (length(ctypes) + ch[1] - 1))
    end
    set_gtk_property!(combo_ch, :tooltip_text, "Channels")

    cb_mono = GtkCheckButton()
    set_gtk_property!(cb_mono, :tooltip_text, "Use color or gray palette")

    cb_db = GtkCheckButton()
    set_gtk_property!(cb_db, :tooltip_text, "Normalize powers to dB")
    set_gtk_property!(cb_db, :active, true)

    cb_hw = GtkCheckButton()
    set_gtk_property!(cb_hw, :tooltip_text, "Apply Hanning window")
    set_gtk_property!(cb_hw, :active, true)

    cb_frq = GtkCheckButton()
    set_gtk_property!(cb_frq, :tooltip_text, "Linear or logarithmic frequencies")
    set_gtk_property!(cb_frq, :active, true)

    cb_smooth = GtkCheckButton()
    set_gtk_property!(cb_smooth, :tooltip_text, "Smooth the spectrogram image")
    set_gtk_property!(cb_smooth, :active, true)

    combo_method = GtkComboBoxText()
    spectrogram_methods = ["short-time Fourier transform", "multi-taper", "Morlet wavelet", "Gaussian and Hilbert transform", "CWT"]
    for idx in spectrogram_methods
        push!(combo_method, idx)
    end
    set_gtk_property!(combo_method, :active, 0)
    set_gtk_property!(combo_method, :tooltip_text, "Spectrogram method")

    entry_wt = GtkEntry()
    set_gtk_property!(entry_wt, :text, "Morlet(2π), β=32, Q=128")
    set_gtk_property!(entry_wt, :tooltip_text, "Continuous wavelet formula")

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

    entry_wlen = GtkSpinButton(2, zoom * sr(obj) + 1, 1)
    set_gtk_property!(entry_wlen, :value, sr(obj))
    set_gtk_property!(entry_wlen, :tooltip_text, "Window length (samples)")

    entry_woverlap = GtkSpinButton(0, zoom * sr(obj), 1)
    set_gtk_property!(entry_woverlap, :value, round(Int64, sr(obj) * 0.97))
    set_gtk_property!(entry_woverlap, :tooltip_text, "Window overlap (samples)")

    entry_gw = GtkSpinButton(1, 128, 1)
    set_gtk_property!(entry_gw, :value, 6)
    set_gtk_property!(entry_gw, :tooltip_text, "Gaussian width in Hz")

    entry_n = GtkSpinButton(1, 64, 1)
    set_gtk_property!(entry_n, :value, 3)
    set_gtk_property!(entry_n, :tooltip_text, "Gaussian smoothing filter kernel size")

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
    lab_gw = GtkLabel("Gaussian width:")
    set_gtk_property!(lab_gw, :halign, 2)
    lab_wt = GtkLabel("Continuous wavelet:")
    set_gtk_property!(lab_wt, :halign, 2)
    lab_frq = GtkLabel("Linear frequencies:")
    set_gtk_property!(lab_frq, :halign, 2)
    lab_frqn = GtkLabel("Frequencies:")
    set_gtk_property!(lab_frqn, :halign, 2)
    lab_smooth = GtkLabel("Smooth:")
    set_gtk_property!(lab_smooth, :halign, 2)
    lab_n = GtkLabel("Kernel size:")
    set_gtk_property!(lab_n, :halign, 2)
    g_opts[1, 1] = lab_ch
    g_opts[1, 2] = lab_method
    g_opts[1, 3] = lab_t
    g_opts[1, 4] = lab_x
    g_opts[1, 5] = lab_y
    g_opts[1, 6] = lab_frq1
    g_opts[1, 7] = lab_frq2
    g_opts[1, 8] = lab_nc
    g_opts[1, 9] = lab_nt
    g_opts[1, 10] = lab_wlen
    g_opts[1, 11] = lab_woverlap
    g_opts[1, 12] = lab_gw
    g_opts[1, 13] = lab_wt
    g_opts[1, 14] = lab_norm
    g_opts[1, 15] = lab_hw
    g_opts[1, 16] = lab_frq
    g_opts[1, 17] = lab_mono
    g_opts[1, 18] = lab_smooth
    g_opts[1, 19] = lab_n
    g_opts[2, 1] = combo_ch
    g_opts[2, 2] = combo_method
    g_opts[2, 3] = entry_title
    g_opts[2, 4] = entry_xlab
    g_opts[2, 5] = entry_ylab
    g_opts[2, 6] = entry_frq1
    g_opts[2, 7] = entry_frq2
    g_opts[2, 8] = entry_ncyc
    g_opts[2, 9] = entry_nt
    g_opts[2, 10] = entry_wlen
    g_opts[2, 11] = entry_woverlap
    g_opts[2, 12] = entry_gw
    g_opts[2, 13] = entry_wt
    g_opts[2, 14] = cb_db
    g_opts[2, 15] = cb_hw
    g_opts[2, 16] = cb_frq
    g_opts[2, 17] = cb_mono
    g_opts[2, 18] = cb_smooth
    g_opts[2, 19] = entry_n
    g_opts[1:2, 20] = bt_refresh
    vbox = GtkBox(:v)
    push!(vbox, g_opts)

    g[1, 1] = vbox
    g[2:9, 1] = can
    g[2, 2] = bt_start
    g[3, 2] = bt_prev5
    g[4, 2] = entry_time
    g[5, 2] = bt_next5
    g[6, 2] = bt_end
    g[7, 2] = GtkLabel("")
    g[8, 2] = bt_help
    g[9, 2] = bt_close
    push!(win, g)

    showall(win)

    @guarded draw(can) do widget
        ch = ch_types[get_gtk_property(combo_ch, :active, Int64) + 1]
        ch in ctypes && (ch = lowercase(ch))
        title = get_gtk_property(entry_title, :text, String)
        xlab = get_gtk_property(entry_xlab, :text, String)
        ylab = get_gtk_property(entry_ylab, :text, String)
        mono = get_gtk_property(cb_mono, :active, Bool)
        db = get_gtk_property(cb_db, :active, Bool)
        frq = get_gtk_property(cb_frq, :active, Bool) ? :lin : :log
        hw = get_gtk_property(cb_hw, :active, Bool)
        smooth = get_gtk_property(cb_smooth, :active, Bool)
        n = get_gtk_property(entry_n, :value, Int64)
        method = get_gtk_property(combo_method, :active, String)
        method == "0" && (method = :stft)
        method == "1" && (method = :mt)
        method == "2" && (method = :mw)
        method == "3" && (method = :gh)
        method == "4" && (method = :cwt)
        wt = nothing
        try
            wt = eval(Meta.parse("wavelet(" * get_gtk_property(entry_wt, :text, String) * ")"))
        catch
        end
        frq1 = get_gtk_property(entry_frq1, :value, Float64)
        frq2 = get_gtk_property(entry_frq2, :value, Float64)
        nt = get_gtk_property(entry_nt, :value, Int64)
        ncyc = get_gtk_property(entry_ncyc, :value, Int64)
        wlen = get_gtk_property(entry_wlen, :value, Int64)
        woverlap = get_gtk_property(entry_woverlap, :value, Int64)
        gw = get_gtk_property(entry_gw, :value, Int64)

        no_error = true
        if frq1 == frq2
            warn_dialog("Start and end frequencies must be different.")
            no_error = false
        elseif wt === nothing
            warn_dialog("Incorrect wavelet formula.")
            no_error = false
        elseif frq1 > frq2
            warn_dialog("Start frequency must be < end frequency.")
            no_error = false
        elseif woverlap >= wlen
            warn_dialog("Window overlap must be < window length.")
            no_error = false
        elseif length(unique(obj.header.recording[:channel_type][get_channel(obj, ch=ch)])) > 1
            warn_dialog("For multi-channel spectrogram plot, all channels must be of the same type.")
            no_error = false
        end

        if no_error
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
                                               db=db,
                                               method=method,
                                               frq_lim=(frq1, frq2),
                                               ncyc=ncyc,
                                               nt=nt,
                                               wlen=wlen,
                                               woverlap=woverlap,
                                               w=hw,
                                               wt=wt,
                                               gw=gw,
                                               frq=frq,
                                               smooth=smooth,
                                               n=n)
            img = read_from_png(io)
            # Gtk.resize!(win, 1200, p.attr[:size][2] + 40)
            set_gtk_property!(can, :width_request, Int32(p.attr[:size][1]))
            set_gtk_property!(can, :height_request, Int32(p.attr[:size][2]))
            ctx = getgc(can)
            show(io, MIME("image/png"), p)
            img = read_from_png(io)
            set_source_surface(ctx, img, 0, 0)
            paint(ctx)
        end
    end

    can.mouse.scroll = @guarded (widget, event) -> begin
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
            end
        end
    end

    signal_connect(combo_ch, "changed") do widget
        draw(can)
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
    signal_connect(entry_gw, "value-changed") do widget
        draw(can)
    end
    signal_connect(entry_wlen, "value-changed") do widget
        draw(can)
    end
    signal_connect(entry_woverlap, "value-changed") do widget
        draw(can)
    end
    signal_connect(cb_db, "clicked") do widget
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
    signal_connect(cb_smooth, "clicked") do widget
        draw(can)
    end
    signal_connect(entry_n, "value-changed") do widget
        draw(can)
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

    help = "Keyboard shortcuts:\n\nHome\t\t\tGo to the signal beginning\nEnd\t\t\tGo to the signal end\nCtrl + ,\t\t\tGo back by 1 second\nCtrl + .\t\t\tGo forward by 1 second\nAlt + ,\t\t\tGo back by $(round(zoom)) seconds\nAlt + .\t\t\tGo forward by $(round(zoom)) seconds\n\n[\t\t\t\tZoom in\n]\t\t\t\tZoom out\n\nCtrl + s\t\t\tSave as PNG\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tExit\n"

    signal_connect(bt_help, "clicked") do widgete
        info_dialog(help)
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        s = event.state

        if k == 0x0000005b # [
            if zoom > 1
                zoom -= 1
                set_gtk_property!(bt_next5, :tooltip_text, "Go forward by $(round(zoom)) seconds")
                set_gtk_property!(bt_prev5, :tooltip_text, "Go back by $(round(zoom)) seconds")
                draw(can)
            end
            help = "Keyboard shortcuts:\n\nHome\t\t\tGo to the signal beginning\nEnd\t\t\tGo to the signal end\nCtrl + ,\t\t\tGo back by 1 second\nCtrl + .\t\t\tGo forward by 1 second\nAlt + ,\t\t\tGo back by $(round(zoom)) seconds\nAlt + .\t\t\tGo forward by $(round(zoom)) seconds\n\n[\t\t\t\tZoom in\n]\t\t\t\tZoom out\n\nCtrl + s\t\t\tSave as PNG\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tExit\n"
        elseif k == 0x0000005d # ]
            if zoom < 30 && zoom < obj.time_pts[end] - 1
                zoom += 1
                set_gtk_property!(bt_next5, :tooltip_text, "Go forward by $(round(zoom)) seconds")
                set_gtk_property!(bt_prev5, :tooltip_text, "Go back by $(round(zoom)) seconds")
                draw(can)
            else
                zoom = obj.time_pts[end]
                set_gtk_property!(bt_next5, :tooltip_text, "Go forward by $(round(zoom)) seconds")
                set_gtk_property!(bt_prev5, :tooltip_text, "Go back by $(round(zoom)) seconds")
                draw(can)
            end
            help = "Keyboard shortcuts:\n\nHome\t\t\tGo to the signal beginning\nEnd\t\t\tGo to the signal end\nCtrl + ,\t\t\tGo back by 1 second\nCtrl + .\t\t\tGo forward by 1 second\nAlt + ,\t\t\tGo back by $(round(zoom)) seconds\nAlt + .\t\t\tGo forward by $(round(zoom)) seconds\n\n[\t\t\t\tZoom in\n]\t\t\t\tZoom out\n\nCtrl + s\t\t\tSave as PNG\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tExit\n"
        elseif k == 0x0000ff50 # home
            Gtk.@sigatom begin
                set_gtk_property!(entry_time, :value, obj.time_pts[1])
            end
            draw(can)
        elseif k == 0x0000ff57 # end
            time_current = obj.time_pts[end] - zoom
            Gtk.@sigatom begin
                set_gtk_property!(entry_time, :value, time_current)
            end
            draw(can)
        end

        if s == 0x00000008 # alt
            if k == 0x0000002c # ,
                time_current = get_gtk_property(entry_time, :value, Float64)
                if time_current >= obj.time_pts[1] + zoom
                    time_current = time_current - zoom
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_time, :value, time_current)
                    end
                end
                draw(can)
            elseif k == 0x0000002e # .
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

        if s == 0x00000004 # ctrl
            if k == 0x00000071 # q
                Gtk.destroy(win)
            elseif k == 0x00000068 # h
                info_dialog(help)
            elseif k == 0x00000073 # s
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
            elseif k == 0x0000002c # ,
                time_current = get_gtk_property(entry_time, :value, Float64)
                if time_current >= obj.time_pts[1] + 1
                    time_current -= 1
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_time, :value, time_current)
                    end
                end
                draw(can)
            elseif k == 0x0000002e # .
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
    ispectrogram_ep(obj; <keyword arguments>)

Interactive spectrogram of epoched signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::String`: channel name
"""
function ispectrogram_ep(obj::NeuroAnalyzer.NEURO; ch::String)

    @assert nepochs(obj) > 1 "ispectrogram_cont() must be used for continuous object."

    ch_init = ch
    ch = get_channel(obj, ch=ch)
    clabels = labels(obj)

    p = NeuroAnalyzer.plot_spectrogram(obj, ch=clabels[ch], ep=1)
    g = GtkGrid()
    g_opts = GtkGrid()
    win = GtkWindow("NeuroAnalyzer: ispectrogram_ep()", 1200, 850)
    set_gtk_property!(win, :border_width, 5)
    set_gtk_property!(win, :resizable, true)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    set_gtk_property!(win, :startup_id, "org.neuroanalyzer")
    can = GtkCanvas(Int32(p.attr[:size][1]), Int32(p.attr[:size][2]))
    set_gtk_property!(g, :column_homogeneous, false)
    set_gtk_property!(g_opts, :column_homogeneous, false)
    set_gtk_property!(g, :column_spacing, 5)
    set_gtk_property!(g, :row_spacing, 5)
    set_gtk_property!(g_opts, :row_spacing, 5)
    set_gtk_property!(g_opts, :column_spacing, 5)
    entry_epoch = GtkSpinButton(1, nepochs(obj), 1)
    set_gtk_property!(entry_epoch, :tooltip_text, "Epoch")
    bt_start = GtkButton("⇤")
    set_gtk_property!(bt_start, :tooltip_text, "Go to the signal beginning")
    bt_end = GtkButton("⇥")
    set_gtk_property!(bt_end, :tooltip_text, "Go to the signal end")
    bt_help = GtkButton("Help")
    set_gtk_property!(bt_help, :tooltip_text, "Show help")
    bt_close = GtkButton("Close")
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

    combo_ch = GtkComboBoxText()
    ctypes = uppercase.(unique(obj.header.recording[:channel_type]))
    ch_types = [ctypes; clabels]
    for idx in ch_types
        push!(combo_ch, idx)
    end
    if ch_init in lowercase.(ctypes)
        n = findfirst(isequal(ch_init), lowercase.(ctypes))
        set_gtk_property!(combo_ch, :active, n - 1)
    else
        set_gtk_property!(combo_ch, :active, (length(ctypes) + ch[1] - 1))
    end
    set_gtk_property!(combo_ch, :tooltip_text, "Channels")

    cb_mono = GtkCheckButton()
    set_gtk_property!(cb_mono, :tooltip_text, "Use color or gray palette")

    cb_db = GtkCheckButton()
    set_gtk_property!(cb_db, :tooltip_text, "Normalize powers to dB")
    set_gtk_property!(cb_db, :active, true)

    cb_hw = GtkCheckButton()
    set_gtk_property!(cb_hw, :tooltip_text, "Apply Hanning window")
    set_gtk_property!(cb_hw, :active, true)

    cb_frq = GtkCheckButton()
    set_gtk_property!(cb_frq, :tooltip_text, "Linear or logarithmic frequencies")
    set_gtk_property!(cb_frq, :active, true)

    cb_smooth = GtkCheckButton()
    set_gtk_property!(cb_smooth, :tooltip_text, "Smooth the spectrogram image")
    set_gtk_property!(cb_smooth, :active, true)

    combo_method = GtkComboBoxText()
    spectrogram_methods = ["short-time Fourier transform", "multi-taper", "Morlet wavelet", "Gaussian and Hilbert transform", "CWT"]
    for idx in spectrogram_methods
        push!(combo_method, idx)
    end
    set_gtk_property!(combo_method, :active, 0)
    set_gtk_property!(combo_method, :tooltip_text, "Spectrogram method")

    entry_wt = GtkEntry()
    set_gtk_property!(entry_wt, :text, "Morlet(2π), β=32, Q=128")
    set_gtk_property!(entry_wt, :tooltip_text, "Continuous wavelet formula")

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

    entry_gw = GtkSpinButton(1, 128, 1)
    set_gtk_property!(entry_gw, :value, 6)
    set_gtk_property!(entry_gw, :tooltip_text, "Gaussian width in Hz")

    entry_n = GtkSpinButton(1, 64, 1)
    set_gtk_property!(entry_n, :value, 3)
    set_gtk_property!(entry_n, :tooltip_text, "Gaussian smoothing filter kernel size")

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
    lab_gw = GtkLabel("Gaussian width:")
    set_gtk_property!(lab_gw, :halign, 2)
    lab_wt = GtkLabel("Continuous wavelet:")
    set_gtk_property!(lab_wt, :halign, 2)
    lab_frq = GtkLabel("Linear frequencies:")
    set_gtk_property!(lab_frq, :halign, 2)
    lab_frqn = GtkLabel("Frequencies:")
    set_gtk_property!(lab_frqn, :halign, 2)
    lab_smooth = GtkLabel("Smooth:")
    set_gtk_property!(lab_smooth, :halign, 2)
    lab_n = GtkLabel("Kernel size:")
    set_gtk_property!(lab_n, :halign, 2)

    signal_slider = GtkScale(false, 1:nepochs(obj))
    set_gtk_property!(signal_slider, :draw_value, false)
    set_gtk_property!(signal_slider, :tooltip_text, "Current epoch")

    g_opts[1, 1] = lab_ch
    g_opts[1, 2] = lab_method
    g_opts[1, 3] = lab_t
    g_opts[1, 4] = lab_x
    g_opts[1, 5] = lab_y
    g_opts[1, 6] = lab_frq1
    g_opts[1, 7] = lab_frq2
    g_opts[1, 8] = lab_nc
    g_opts[1, 9] = lab_nt
    g_opts[1, 10] = lab_wlen
    g_opts[1, 11] = lab_woverlap
    g_opts[1, 12] = lab_gw
    g_opts[1, 13] = lab_wt
    g_opts[1, 14] = lab_norm
    g_opts[1, 15] = lab_hw
    g_opts[1, 16] = lab_frq
    g_opts[1, 17] = lab_mono
    g_opts[1, 18] = lab_smooth
    g_opts[1, 19] = lab_n
    g_opts[2, 1] = combo_ch
    g_opts[2, 2] = combo_method
    g_opts[2, 3] = entry_title
    g_opts[2, 4] = entry_xlab
    g_opts[2, 5] = entry_ylab
    g_opts[2, 6] = entry_frq1
    g_opts[2, 7] = entry_frq2
    g_opts[2, 8] = entry_ncyc
    g_opts[2, 9] = entry_nt
    g_opts[2, 10] = entry_wlen
    g_opts[2, 11] = entry_woverlap
    g_opts[2, 12] = entry_gw
    g_opts[2, 13] = entry_wt
    g_opts[2, 14] = cb_db
    g_opts[2, 15] = cb_hw
    g_opts[2, 16] = cb_frq
    g_opts[2, 17] = cb_mono
    g_opts[2, 18] = cb_smooth
    g_opts[2, 19] = entry_n
    g_opts[1:2, 20] = bt_refresh
    vbox = GtkBox(:v)
    push!(vbox, g_opts)

    g[1, 1] = vbox
    g[2:7, 1] = can
    g[2:7, 2] = signal_slider
    g[2, 3] = bt_start
    g[3, 3] = entry_epoch
    g[4, 3] = bt_end
    g[5, 3] = GtkLabel("")
    g[6, 3] = bt_help
    g[7, 3] = bt_close
    push!(win, g)

    showall(win)

    @guarded draw(can) do widget
        ch = ch_types[get_gtk_property(combo_ch, :active, Int64) + 1]
        ch in ctypes && (ch = lowercase(ch))
        title = get_gtk_property(entry_title, :text, String)
        xlab = get_gtk_property(entry_xlab, :text, String)
        ylab = get_gtk_property(entry_ylab, :text, String)
        mono = get_gtk_property(cb_mono, :active, Bool)
        db = get_gtk_property(cb_db, :active, Bool)
        frq = get_gtk_property(cb_frq, :active, Bool) ? :lin : :log
        hw = get_gtk_property(cb_hw, :active, Bool)
        smooth = get_gtk_property(cb_smooth, :active, Bool)
        n = get_gtk_property(entry_n, :value, Int64)
        method = get_gtk_property(combo_method, :active, String)
        method == "0" && (method = :stft)
        method == "1" && (method = :mt)
        method == "2" && (method = :mw)
        method == "3" && (method = :gh)
        method == "4" && (method = :cwt)
        wt = nothing
        try
            wt = eval(Meta.parse("wavelet(" * get_gtk_property(entry_wt, :text, String) * ")"))
        catch
        end
        frq1 = get_gtk_property(entry_frq1, :value, Float64)
        frq2 = get_gtk_property(entry_frq2, :value, Float64)
        nt = get_gtk_property(entry_nt, :value, Int64)
        ncyc = get_gtk_property(entry_ncyc, :value, Int64)
        wlen = get_gtk_property(entry_wlen, :value, Int64)
        woverlap = get_gtk_property(entry_woverlap, :value, Int64)
        gw = get_gtk_property(entry_gw, :value, Int64)

        no_error = true
        if frq1 == frq2
            warn_dialog("Start and end frequencies must be different.")
            no_error = false
        elseif wt === nothing
            warn_dialog("Incorrect wavelet formula.")
            no_error = false
        elseif frq1 > frq2
            warn_dialog("Start frequency must be < end frequency.")
            no_error = false
        elseif woverlap >= wlen
            warn_dialog("Window overlap must be < window length.")
            no_error = false
        elseif length(unique(obj.header.recording[:channel_type][get_channel(obj, ch=ch)])) > 1
            warn_dialog("For multi-channel spectrogram plot, all channels must be of the same type.")
            no_error = false
        end
        if no_error
            ep = get_gtk_property(entry_epoch, :value, Int64)
            p = NeuroAnalyzer.plot_spectrogram(obj,
                                               ch=ch,
                                               ep=ep,
                                               mono=mono,
                                               title=title,
                                               xlabel=xlab,
                                               ylabel=ylab,
                                               db=db,
                                               method=method,
                                               frq_lim=(frq1, frq2),
                                               ncyc=ncyc,
                                               nt=nt,
                                               wlen=wlen,
                                               woverlap=woverlap,
                                               w=hw,
                                               wt=wt,
                                               gw=gw,
                                               frq=frq,
                                               smooth=smooth,
                                               n=n)
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

    signal_connect(combo_ch, "changed") do widget
        draw(can)
    end
    signal_connect(bt_refresh, "clicked") do widget
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
    signal_connect(entry_gw, "value-changed") do widget
        draw(can)
    end
    signal_connect(entry_wlen, "value-changed") do widget
        draw(can)
    end
    signal_connect(entry_woverlap, "value-changed") do widget
        draw(can)
    end
    signal_connect(cb_db, "clicked") do widget
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
    signal_connect(cb_smooth, "clicked") do widget
        draw(can)
    end
    signal_connect(entry_n, "value-changed") do widget
        draw(can)
    end

    signal_connect(entry_epoch, "value-changed") do widget
         Gtk.@sigatom begin
            GAccessor.value(signal_slider, get_gtk_property(entry_epoch, :value, Int64))
        end
        draw(can)
    end

    signal_connect(signal_slider, "value-changed") do widget, others...
        Gtk.@sigatom begin
            set_gtk_property!(entry_epoch, :value, round(Int64, GAccessor.value(signal_slider)))
        end
        draw(can)
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

    help = "Keyboard shortcuts:\n\nHome\t\t\tGo to the first epoch\nEnd\t\t\tGo to last epoch\nCtrl + ,\t\t\tPrevious epoch\nCtrl + .\t\t\tNext epoch\n\nCtrl + s\t\t\tSave as PNG\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tExit\n"

    signal_connect(bt_help, "clicked") do widgete
        info_dialog(help)
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        s = event.state
        if k == 0x0000ff50 # home
            Gtk.@sigatom begin
                set_gtk_property!(entry_epoch, :value, 1)
            end
            draw(can)
        elseif k == 0x0000ff57 # end
            Gtk.@sigatom begin
                set_gtk_property!(entry_epoch, :value, nepochs(obj))
            end
            draw(can)
        end
        if s == 0x00000004
            if k == 0x00000071 # q
                Gtk.destroy(win)
            elseif k == 0x00000068 # h
                info_dialog(help)
            elseif k == 0x00000073 # s
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
            elseif k == 0x0000002c # ,
                ep = get_gtk_property(entry_epoch, :value, Int64)
                if ep >= 2
                    ep -= 1
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_epoch, :value, ep)
                    end
                end
            elseif k == 0x0000002e # .
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

    cnd = Condition()
    signal_connect(win, :destroy) do widget
        notify(cnd)
    end
    @async Gtk.gtk_main()
    wait(cnd)

    return nothing
end