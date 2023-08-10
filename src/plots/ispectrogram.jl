export ispectrogram
export ispectrogram_cont
export ispectrogram_ep

"""
    ispectrogram(obj; <keyword arguments>)

Interactive PSD.

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
    ispectrogram_cont(obj; <keyword arguments>)

Interactive PSD of continuous signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1`: index of channels, default is all signal channels
- `zoom::Int64=5`: how many seconds are displayed in one segment

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function ispectrogram_cont(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1, zoom::Int64=5)

    @assert epoch_n(obj) == 1 "ispectrogram_cont() should be used for epoched object."
    _check_channels(obj, ch)

    @assert zoom >= 1 "zoom must be ≥ 1."
    @assert zoom <= signal_len(obj) / sr(obj) "zoom must be ≤ $(signal_len(obj) / sr(obj))."

    p = NeuroAnalyzer.plot_spectrogram(obj, ch=ch)
    g = GtkGrid()
    g_opts = GtkGrid()
    win = GtkWindow("NeuroAnalyzer: ispectrogram_cont()", 1200, (p.attr[:size][2] + 40))
    set_gtk_property!(win, :border_width, 20)
    set_gtk_property!(win, :resizable, true)
    set_gtk_property!(win, :has_resize_grip, false)
    can = GtkCanvas(Int32(p.attr[:size][1]), Int32(p.attr[:size][2]))
    set_gtk_property!(g, :column_homogeneous, false)
    set_gtk_property!(g_opts, :column_homogeneous, false)
    set_gtk_property!(g, :column_spacing, 10)  # introduce a 10-pixel gap between columns
    set_gtk_property!(g, :row_spacing, 10)  # introduce a 10-pixel gap between columns
    set_gtk_property!(g_opts, :row_spacing, 10)  # introduce a 10-pixel gap between columns
    set_gtk_property!(g_opts, :column_spacing, 10)  # introduce a 10-pixel gap between columns
    entry_time = GtkButton(string(obj.time_pts[1]))
    set_gtk_property!(entry_time, :tooltip_text, "Time position [s]")
    bt_start = GtkButton("|<")
    set_gtk_property!(bt_start, :tooltip_text, "Go to the signal beginning")
    bt_prev5 = GtkButton("<<")
    set_gtk_property!(bt_prev5, :tooltip_text, "Go back by $zoom seconds")
    bt_prev = GtkButton("<")
    set_gtk_property!(bt_prev, :tooltip_text, "Go back by 1 second")
    bt_next = GtkButton(">")
    set_gtk_property!(bt_next, :tooltip_text, "Go forward by 1 second")
    bt_next5 = GtkButton(">>")
    set_gtk_property!(bt_next5, :tooltip_text, "Go forward by $zoom seconds")
    bt_end = GtkButton(">|")
    set_gtk_property!(bt_end, :tooltip_text, "Go to the signal end")
    bt_help = GtkButton("🛈")
    set_gtk_property!(bt_help, :tooltip_text, "Show keyboard shortcuts")
    bt_delete = GtkButton("DEL")
    set_gtk_property!(bt_delete, :tooltip_text, "Delete segment")
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
    set_gtk_property!(cb_mono, :tooltip_text, "Use color or grey palette")

    cb_norm = GtkCheckButton()
    set_gtk_property!(cb_norm, :tooltip_text, "Normalize powers to dB")
    set_gtk_property!(cb_norm, :active, true)

   combo_method = GtkComboBoxText()
    spectrogram_methods = ["Standard", "STFT", "MT", "MW"]
    for idx in spectrogram_methods
        push!(combo_method, idx)
    end
    set_gtk_property!(combo_method, :active, 0)
    set_gtk_property!(combo_method, :tooltip_text, "Spectrogram method")

    entry_nt = GtkEntry()
    set_gtk_property!(entry_nt, :text, "8")
    set_gtk_property!(entry_nt, :tooltip_text, "Number of Slepian tapers")

    entry_ncyc = GtkEntry()
    set_gtk_property!(entry_ncyc, :text, "6")
    set_gtk_property!(entry_ncyc, :tooltip_text, "Number of Morlet wavelet cycles")

    entry_frq = GtkEntry()
    set_gtk_property!(entry_frq, :text, "(0, $(sr(obj)/2))")
    set_gtk_property!(entry_frq, :tooltip_text, "Frequency limits")

    bt_png = GtkButton("PNG")
    set_gtk_property!(bt_png, :tooltip_text, "Save as PNG")

    bt_pdf = GtkButton("PDF")
    set_gtk_property!(bt_pdf, :tooltip_text, "Save as PDF")

    bt_refresh = GtkButton("Refresh")
    set_gtk_property!(bt_pdf, :tooltip_text, "Refresh the plot")

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
    lab_mono = GtkLabel("Greyscale:")
    set_gtk_property!(lab_mono, :halign, 2)
    lab_png = GtkLabel("Save as:")
    set_gtk_property!(lab_png, :halign, 2)
    lab_pdf = GtkLabel("Save as:")
    set_gtk_property!(lab_pdf, :halign, 2)
    lab_norm = GtkLabel("Normalize:")
    set_gtk_property!(lab_norm, :halign, 2)
    lab_nt = GtkLabel("Slepians:")
    set_gtk_property!(lab_nt, :halign, 2)
    lab_frq = GtkLabel("Frequencies:")
    set_gtk_property!(lab_frq, :halign, 2)
    lab_nc = GtkLabel("Cycles:")
    set_gtk_property!(lab_nc, :halign, 2)
    g_opts[1, 1] = lab_method
    g_opts[1, 2] = lab_ch
    g_opts[1, 3] = lab_t
    g_opts[1, 4] = lab_x
    g_opts[1, 5] = lab_y
    g_opts[1, 6] = lab_frq
    g_opts[1, 7] = lab_nc
    g_opts[1, 8] = lab_nt
    g_opts[1, 9] = lab_norm
    g_opts[1, 10] = lab_mono
    g_opts[1, 11] = lab_png
    g_opts[1, 12] = lab_pdf
    g_opts[2, 1] = combo_method
    g_opts[2, 2] = entry_ch
    g_opts[2, 3] = entry_title
    g_opts[2, 4] = entry_xlab
    g_opts[2, 5] = entry_ylab
    g_opts[2, 6] = entry_frq
    g_opts[2, 7] = entry_ncyc
    g_opts[2, 8] = entry_nt
    g_opts[2, 9] = cb_norm
    g_opts[2, 10] = cb_mono
    g_opts[2, 11] = bt_png
    g_opts[2, 12] = bt_pdf
    g_opts[1:2, 13] = bt_refresh
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
        if (occursin(", ", ch) && _check_svec(ch)) || (occursin(":", ch) && _check_srange(ch)) || _check_sint(ch)
            ch = _s2i(ch)
            _check_channels(obj, ch)
            title = get_gtk_property(entry_title, :text, String)
            xlab = get_gtk_property(entry_xlab, :text, String)
            ylab = get_gtk_property(entry_ylab, :text, String)
            mono = get_gtk_property(cb_mono, :active, Bool)
            norm = get_gtk_property(cb_norm, :active, Bool)
            method = get_gtk_property(combo_method, :active, String)
            method == "0" && (method = :standard)
            method == "1" && (method = :stft)
            method == "2" && (method = :mt)
            method == "3" && (method = :mw)
            frq = get_gtk_property(entry_frq, :text, String)
            if _check_stuplef(frq)
                frq = _s2tf(frq)
                fs = sr(obj)
                frq = tuple_order(frq)
                if frq[1] < 0 || frq[2] < 0 || frq[1] > sr(obj) / 2 || frq[2] > sr(obj) / 2
                    warn_dialog("frq must be in [0, $(fs / 2)].")
                    frq = nothing
                end
            end
            nt = get_gtk_property(entry_nt, :text, String)
            _check_sint(nt) && (nt = parse(Int64, nt))
            ncyc = get_gtk_property(entry_ncyc, :text, String)
            if _check_stuplei(ncyc)
                ncyc = _s2ti(ncyc)
                if ncyc[1] < 0 || ncyc[1] < 0
                    warn_dialog("ncyc must be in ≥ (1, 1).")
                    ncyc = nothing
                end
            elseif _check_sint(ncyc)
                ncyc = parse(Int64, ncyc)
                if ncyc < 1
                    warn_dialog("ncyc must be ≥ 1.")
                    ncyc = nothing
                end
            end
            if isa(ch, Union{Int64, Vector{Int64}}) && isa(frq, Tuple{Real, Real}) && isa(nt, Int64) && isa(ncyc, Union{Int64, Tuple{Int64, Int64}})
                time1 = parse(Float64, get_gtk_property(entry_time, :label, String))
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
                                                   nt=nt)
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
                                                                           nt=nt))
                img = read_from_png(io)
                set_source_surface(ctx, img, 0, 0)
                paint(ctx)
            end
        end
    end

    signal_connect(bt_refresh, "clicked") do widget
        draw(can)
    end

    signal_connect(bt_pdf, "clicked") do widget
        file_name = save_dialog("Save as PDF", GtkNullContainer(), (GtkFileFilter("*.pdf", name="All supported formats"), "*.pdf"))
        splitext(file_name)[2] == "" && (file_name *= ".pdf")
        if splitext(file_name)[2] == ".pdf"
            plot_save(p, file_name=file_name)
            _info("Plot saved as: $file_name")
        else
            warn_dialog("Incorrect file name!")
        end
    end

    signal_connect(bt_png, "clicked") do widget
        file_name = save_dialog("Save as PNG", GtkNullContainer(), (GtkFileFilter("*.png", name="All supported formats"), "*.png"))
        splitext(file_name)[2] == "" && (file_name *= ".png")
        if splitext(file_name)[2] == ".png"
            plot_save(p, file_name=file_name)
            _info("Plot saved as: $file_name")
        else
            warn_dialog("Incorrect file name!")
        end
    end

    signal_connect(bt_prev, "clicked") do widget
        time_current = parse(Float64, get_gtk_property(entry_time, :label, String))
        if time_current >= obj.time_pts[1] + 1
            time_current -= 1
            Gtk.@sigatom begin
                set_gtk_property!(entry_time, :label, string(time_current))
            end
        end
        draw(can)
    end

    signal_connect(bt_prev5, "clicked") do widget
        time_current = parse(Float64, get_gtk_property(entry_time, :label, String))
        if time_current >= obj.time_pts[1] + zoom
            time_current = time_current - zoom
            Gtk.@sigatom begin
                set_gtk_property!(entry_time, :label, string(time_current))
            end
        end
        draw(can)
    end

    signal_connect(bt_next, "clicked") do widget
        time_current = parse(Float64, get_gtk_property(entry_time, :label, String))
        if time_current < obj.time_pts[end] - zoom
            time_current += 1
            Gtk.@sigatom begin
                set_gtk_property!(entry_time, :label, string(time_current))
            end
        else
            time_current = obj.time_pts[end] - zoom
            Gtk.@sigatom begin
                set_gtk_property!(entry_time, :label, string(time_current))
            end
        end
        draw(can)
    end

    signal_connect(bt_next5, "clicked") do widget
        time_current = parse(Float64, get_gtk_property(entry_time, :label, String))
        if time_current < obj.time_pts[end] - zoom
            time_current += zoom
            Gtk.@sigatom begin
                set_gtk_property!(entry_time, :label, string(time_current))
            end
        else
            time_current = obj.time_pts[end] - zoom
            Gtk.@sigatom begin
                set_gtk_property!(entry_time, :label, string(time_current))
            end
        end
        draw(can)
    end

    signal_connect(bt_start, "clicked") do widget
        Gtk.@sigatom begin
            set_gtk_property!(entry_time, :label, string(obj.time_pts[1]))
        end
        draw(can)
    end

    signal_connect(bt_end, "clicked") do widget
        time_current = obj.time_pts[end] - zoom
        Gtk.@sigatom begin
            set_gtk_property!(entry_time, :label, string(time_current))
        end
        draw(can)
    end

    signal_connect(entry_time, "clicked") do widget
        value = parse(Float64, get_gtk_property(entry_time, :label, String))
        d_w = GtkWindow("Enter value", 200, 100)
        set_gtk_property!(d_w, :border_width, 20)
        set_gtk_property!(d_w, :resizable, true)
        d_g = GtkGrid()
        set_gtk_property!(d_g, :column_homogeneous, true)
        set_gtk_property!(d_g, :column_spacing, 10)  # introduce a 10-pixel gap between columns
        d_entry = GtkEntry()
        set_gtk_property!(d_entry, :text, string(value))
        d_bt_ok = GtkButton("Ok")
        d_bt_cancel = GtkButton("Cancel")
        d_g[1:2, 1] = d_entry
        d_g[1, 2] = d_bt_ok
        d_g[2, 2] = d_bt_cancel
        push!(d_w, d_g)
        showall(d_w)
        signal_connect(d_bt_ok, "clicked") do widget
            value_s = get_gtk_property(d_entry, :text, String)
            value_currect = true
            for idx in eachindex(value_s)
                string(value_s[idx]) in vcat(string.(0:9), ["."]) || (value_currect = false)
            end
            if value_currect
                v = parse(Float64, value_s)
                if v < obj.time_pts[1]
                    warn_dialog("Value must be ≥ $(obj.time_pts[1]).")
                elseif v > obj.time_pts[end] - zoom
                    warn_dialog("Value must be ≤ $(obj.time_pts[end] - zoom).")
                else
                    value_s = string(obj.time_pts[vsearch(parse(Float64, value_s), obj.time_pts)])
                    set_gtk_property!(entry_time, :label, value_s)
                    draw(can)
                    Gtk.destroy(d_w)
                end
            else
                warn_dialog("Incorrect value entered!")
            end
        end
        signal_connect(d_w, "key-press-event") do widget, event
            k = event.keyval
            if k == 65293 || k == 65421
                value_s = get_gtk_property(d_entry, :text, String)
                value_currect = true
                for idx in eachindex(value_s)
                    string(value_s[idx]) in vcat(string.(0:9), ["."]) || (value_currect = false)
                end
                if value_currect
                    v = parse(Float64, value_s)
                    if v < obj.time_pts[1]
                        warn_dialog("Value must be ≥ $(obj.time_pts[1]).")
                    elseif v > obj.time_pts[end] - zoom
                        warn_dialog("Value must be ≤ $(obj.time_pts[end] - zoom).")
                    else
                        value_s = string(obj.time_pts[vsearch(parse(Float64, value_s), obj.time_pts)])
                        set_gtk_property!(entry_time, :label, value_s)
                        draw(can)
                        Gtk.destroy(d_w)
                    end
                else
                    warn_dialog("Incorrect value entered!")
                end
            end
        end
        signal_connect(d_bt_cancel, "clicked") do widget
            Gtk.destroy(d_w)
        end
    end

    signal_connect(bt_close, "clicked") do widget
        Gtk.destroy(win)
    end

    signal_connect(bt_help, "clicked") do widgete
        info_dialog("Keyboard shortcuts:\n\ng\t\tgo to time point\nHOME\tgo to the signal beginning\nEND\t\tgo to the signal end\n,\t\tgo back by 1 second\n.\t\tgo forward by 1 second\n<\t\tgo back by $zoom seconds\n>\t\tgo forward by $zoom seconds\n\nh\t\tthis info\nq\t\texit\n")
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        if k == 113 # q
            Gtk.destroy(win)
        elseif k == 104 # h
            info_dialog("Keyboard shortcuts:\n\ng\t\tgo to time point\nHOME\tgo to the signal beginning\nEND\t\tgo to the signal end\n,\t\tgo back by 1 second\n.\t\tgo forward by 1 second\n<\t\tgo back by $zoom seconds\n>\t\tgo forward by $zoom seconds\n\nh\t\tthis info\nq\t\texit\n")
        elseif k == 103 # g
            value = parse(Float64, get_gtk_property(entry_time, :label, String))
            d_w = GtkWindow("Enter value", 200, 100)
            set_gtk_property!(d_w, :border_width, 20)
            set_gtk_property!(d_w, :resizable, true)
            d_g = GtkGrid()
            set_gtk_property!(d_g, :column_homogeneous, true)
            set_gtk_property!(d_g, :column_spacing, 10)  # introduce a 10-pixel gap between columns
            d_entry = GtkEntry()
            set_gtk_property!(d_entry, :text, string(value))
            d_bt_ok = GtkButton("Ok")
            d_bt_cancel = GtkButton("Cancel")
            d_g[1:2, 1] = d_entry
            d_g[1, 2] = d_bt_ok
            d_g[2, 2] = d_bt_cancel
            push!(d_w, d_g)
            showall(d_w)
            signal_connect(d_bt_ok, "clicked") do widget
                value_s = get_gtk_property(d_entry, :text, String)
                value_currect = true
                for idx in eachindex(value_s)
                    string(value_s[idx]) in vcat(string.(0:9), ["."]) || (value_currect = false)
                end
                if value_currect
                    v = parse(Float64, value_s)
                    if v < obj.time_pts[1]
                        warn_dialog("Value must be ≥ $(obj.time_pts[1]).")
                    elseif v > obj.time_pts[end] - zoom
                        warn_dialog("Value must be ≥ $(obj.time_pts[end] - zoom).")
                    else
                        value_s = string(obj.time_pts[vsearch(parse(Float64, value_s), obj.time_pts)])
                        set_gtk_property!(entry_time, :label, value_s)
                        draw(can)
                        Gtk.destroy(d_w)
                    end
                else
                    warn_dialog("Incorrect value entered!")
                end
            end
            signal_connect(d_w, "key-press-event") do widget, event
                k = event.keyval
                if k == 65293 || k == 65421
                    value_s = get_gtk_property(d_entry, :text, String)
                    value_currect = true
                    for idx in eachindex(value_s)
                        string(value_s[idx]) in vcat(string.(0:9), ["."]) || (value_currect = false)
                    end
                    if value_currect
                        v = parse(Float64, value_s)
                        if v < obj.time_pts[1]
                            warn_dialog("Value must be ≥ $(obj.time_pts[1]).")
                        elseif v > obj.time_pts[end] - zoom
                            warn_dialog("Value must be ≥ $(obj.time_pts[end] - zoom).")
                        else
                            value_s = string(obj.time_pts[vsearch(parse(Float64, value_s), obj.time_pts)])
                            set_gtk_property!(entry_time, :label, value_s)
                            draw(can)
                            Gtk.destroy(d_w)
                        end
                    else
                        warn_dialog("Incorrect value entered!")
                    end
                end
            end
            signal_connect(d_bt_cancel, "clicked") do widget
                Gtk.destroy(d_w)
            end
        elseif k == 65360 # HOME
            Gtk.@sigatom begin
                set_gtk_property!(entry_time, :label, string(obj.time_pts[1]))
            end
            draw(can)
        elseif k == 65367 # END
            time_current = obj.time_pts[end] - zoom
            Gtk.@sigatom begin
                set_gtk_property!(entry_time, :label, string(time_current))
            end
            draw(can)
        elseif k == 44 # ,
            time_current = parse(Float64, get_gtk_property(entry_time, :label, String))
            if time_current >= obj.time_pts[1] + 1
                time_current -= 1
                Gtk.@sigatom begin
                    set_gtk_property!(entry_time, :label, string(time_current))
                end
            end
            draw(can)
        elseif k == 60 # <
            time_current = parse(Float64, get_gtk_property(entry_time, :label, String))
            if time_current >= obj.time_pts[1] + zoom
                time_current = time_current - zoom
                Gtk.@sigatom begin
                    set_gtk_property!(entry_time, :label, string(time_current))
                end
            end
            draw(can)
        elseif k == 46 # .
            time_current = parse(Float64, get_gtk_property(entry_time, :label, String))
            if time_current < obj.time_pts[end] - zoom
                time_current += 1
                Gtk.@sigatom begin
                    set_gtk_property!(entry_time, :label, string(time_current))
                end
            else
                time_current = obj.time_pts[end] - zoom
                Gtk.@sigatom begin
                    set_gtk_property!(entry_time, :label, string(time_current))
                end
            end
            draw(can)
        elseif k == 62 # >
            time_current = parse(Float64, get_gtk_property(entry_time, :label, String))
            if time_current < obj.time_pts[end] - zoom
                time_current += zoom
                Gtk.@sigatom begin
                    set_gtk_property!(entry_time, :label, string(time_current))
                end
            else
                time_current = obj.time_pts[end] - zoom
                Gtk.@sigatom begin
                    set_gtk_property!(entry_time, :label, string(time_current))
                end
            end
            draw(can)
        end
    end

    return nothing

end

"""
    ispectrogram_ep(obj; <keyword arguments>)

Interactive plot of epoched signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1`: index of channels, default is all signal channels
- `zoom::Int64=5`: how many seconds are displayed in one segment

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function ispectrogram_ep(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1, zoom::Int64=5)

    @assert epoch_n(obj) > 1 "ispectrogram_cont() should be used for continuous object."
    _check_channels(obj, ch)

    @assert zoom >= 1 "zoom must be ≥ 1."
    @assert zoom <= signal_len(obj) / sr(obj) "zoom must be ≤ $(signal_len(obj) / sr(obj))."

    p = NeuroAnalyzer.plot_spectrogram(obj, ch=ch)
    g = GtkGrid()
    g_opts = GtkGrid()
    win = GtkWindow("NeuroAnalyzer: ispectrogram_ep()", 1200, (p.attr[:size][2] + 40))
    set_gtk_property!(win, :border_width, 20)
    set_gtk_property!(win, :resizable, true)
    set_gtk_property!(win, :has_resize_grip, false)
    can = GtkCanvas(Int32(p.attr[:size][1]), Int32(p.attr[:size][2]))
    set_gtk_property!(g, :column_homogeneous, false)
    set_gtk_property!(g_opts, :column_homogeneous, false)
    set_gtk_property!(g, :column_spacing, 10)  # introduce a 10-pixel gap between columns
    set_gtk_property!(g, :row_spacing, 10)  # introduce a 10-pixel gap between columns
    set_gtk_property!(g_opts, :row_spacing, 10)  # introduce a 10-pixel gap between columns
    set_gtk_property!(g_opts, :column_spacing, 10)  # introduce a 10-pixel gap between columns
    entry_epoch = GtkButton(string(1))
    set_gtk_property!(entry_epoch, :tooltip_text, "Epoch")
    bt_start = GtkButton("|<")
    set_gtk_property!(bt_start, :tooltip_text, "Go to the signal beginning")
    bt_prev = GtkButton("<")
    set_gtk_property!(bt_prev, :tooltip_text, "Go back by 1 epoch")
    bt_next = GtkButton(">")
    set_gtk_property!(bt_next, :tooltip_text, "Go forward by 1 epoch")
    bt_end = GtkButton(">|")
    set_gtk_property!(bt_end, :tooltip_text, "Go to the signal end")
    bt_help = GtkButton("🛈")
    set_gtk_property!(bt_help, :tooltip_text, "Show keyboard shortcuts")
    bt_delete = GtkButton("DEL")
    set_gtk_property!(bt_delete, :tooltip_text, "Delete epoch")
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
    set_gtk_property!(cb_mono, :tooltip_text, "Use color or grey palette")

    cb_norm = GtkCheckButton()
    set_gtk_property!(cb_norm, :tooltip_text, "Normalize powers to dB")
    set_gtk_property!(cb_norm, :active, true)

    combo_method = GtkComboBoxText()
    spectrogram_methods = ["Standard", "STFT", "MT", "MW"]
    for idx in spectrogram_methods
        push!(combo_method, idx)
    end
    set_gtk_property!(combo_method, :active, 0)
    set_gtk_property!(combo_method, :tooltip_text, "Spectrogram method")

    entry_nt = GtkEntry()
    set_gtk_property!(entry_nt, :text, "8")
    set_gtk_property!(entry_nt, :tooltip_text, "Number of Slepian tapers")

    entry_ncyc = GtkEntry()
    set_gtk_property!(entry_ncyc, :text, "6")
    set_gtk_property!(entry_ncyc, :tooltip_text, "Number of Morlet wavelet cycles")

    entry_frq = GtkEntry()
    set_gtk_property!(entry_frq, :text, "(0, $(sr(obj)/2))")
    set_gtk_property!(entry_frq, :tooltip_text, "Frequency limits")

    bt_png = GtkButton("PNG")
    set_gtk_property!(bt_png, :tooltip_text, "Save as PNG")

    bt_pdf = GtkButton("PDF")
    set_gtk_property!(bt_pdf, :tooltip_text, "Save as PDF")

    bt_refresh = GtkButton("Refresh")
    set_gtk_property!(bt_pdf, :tooltip_text, "Refresh the plot")

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
    lab_mono = GtkLabel("Greyscale:")
    set_gtk_property!(lab_mono, :halign, 2)
    lab_png = GtkLabel("Save as:")
    set_gtk_property!(lab_png, :halign, 2)
    lab_pdf = GtkLabel("Save as:")
    set_gtk_property!(lab_pdf, :halign, 2)
    lab_norm = GtkLabel("Normalize:")
    set_gtk_property!(lab_norm, :halign, 2)
    lab_nt = GtkLabel("Slepians:")
    set_gtk_property!(lab_nt, :halign, 2)
    lab_frq = GtkLabel("Frequencies:")
    set_gtk_property!(lab_frq, :halign, 2)
    lab_nc = GtkLabel("Cycles:")
    set_gtk_property!(lab_nc, :halign, 2)
    g_opts[1, 1] = lab_method
    g_opts[1, 2] = lab_ch
    g_opts[1, 3] = lab_t
    g_opts[1, 4] = lab_x
    g_opts[1, 5] = lab_y
    g_opts[1, 6] = lab_frq
    g_opts[1, 7] = lab_nc
    g_opts[1, 8] = lab_nt
    g_opts[1, 9] = lab_norm
    g_opts[1, 10] = lab_mono
    g_opts[1, 11] = lab_png
    g_opts[1, 12] = lab_pdf
    g_opts[2, 1] = combo_method
    g_opts[2, 2] = entry_ch
    g_opts[2, 3] = entry_title
    g_opts[2, 4] = entry_xlab
    g_opts[2, 5] = entry_ylab
    g_opts[2, 6] = entry_frq
    g_opts[2, 7] = entry_ncyc
    g_opts[2, 8] = entry_nt
    g_opts[2, 9] = cb_norm
    g_opts[2, 10] = cb_mono
    g_opts[2, 11] = bt_png
    g_opts[2, 12] = bt_pdf
    g_opts[1:2, 13] = bt_refresh
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
        if (occursin(", ", ch) && _check_svec(ch)) || (occursin(":", ch) && _check_srange(ch)) || _check_sint(ch)
            ch = _s2i(ch)
            title = get_gtk_property(entry_title, :text, String)
            xlab = get_gtk_property(entry_xlab, :text, String)
            ylab = get_gtk_property(entry_ylab, :text, String)
            mono = get_gtk_property(cb_mono, :active, Bool)
            norm = get_gtk_property(cb_norm, :active, Bool)
            method = get_gtk_property(combo_method, :active, String)
            method == "0" && (method = :standard)
            method == "1" && (method = :stft)
            method == "2" && (method = :mt)
            method == "3" && (method = :mw)
            frq = get_gtk_property(entry_frq, :text, String)
            if _check_stuplef(frq)
                frq = _s2tf(frq)
                fs = sr(obj)
                frq = tuple_order(frq)
                if frq[1] < 0 || frq[2] < 0 || frq[1] > sr(obj) / 2 || frq[2] > sr(obj) / 2
                    warn_dialog("frq must be in [0, $(fs / 2)].")
                    frq = nothing
                end
            end
            nt = get_gtk_property(entry_nt, :text, String)
            _check_sint(nt) && (nt = parse(Int64, nt))
            ncyc = get_gtk_property(entry_ncyc, :text, String)
            if _check_stuplei(ncyc)
                ncyc = _s2ti(ncyc)
                if ncyc[1] < 0 || ncyc[1] < 0
                    warn_dialog("ncyc must be in ≥ (1, 1).")
                    ncyc = nothing
                end
            elseif _check_sint(ncyc)
                ncyc = parse(Int64, ncyc)
                if ncyc < 1
                    warn_dialog("ncyc must be ≥ 1.")
                    ncyc = nothing
                end
            end
            if (isa(ch, Int64) || isa(ch, Vector{Int64})) && isa(frq, Tuple{Real, Real}) && isa(nt, Int64) && isa(ncyc, Union{Int64, Tuple{Int64, Int64}})
                ep = parse(Int64, get_gtk_property(entry_epoch, :label, String))
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
                                                   nt=nt)
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
                                                                           nt=nt))
                img = read_from_png(io)
                set_source_surface(ctx, img, 0, 0)
                paint(ctx)
            end
        end
    end

    signal_connect(bt_refresh, "clicked") do widget
        draw(can)
    end

    signal_connect(bt_pdf, "clicked") do widget
        file_name = save_dialog("Save as PDF", GtkNullContainer(), (GtkFileFilter("*.pdf", name="All supported formats"), "*.pdf"))
        splitext(file_name)[2] == "" && (file_name *= ".pdf")
        if splitext(file_name)[2] == ".pdf"
            plot_save(p, file_name=file_name)
            _info("Plot saved as: $file_name")
        elseif file_name != ".pdf"
            warn_dialog("Incorrect file name!")
        end
    end

    signal_connect(bt_png, "clicked") do widget
        file_name = save_dialog("Save as PNG", GtkNullContainer(), (GtkFileFilter("*.png", name="All supported formats"), "*.png"))
        splitext(file_name)[2] == "" && (file_name *= ".png")
        if splitext(file_name)[2] == ".png"
            plot_save(p, file_name=file_name)
            _info("Plot saved as: $file_name")
        elseif file_name != ".png"
            warn_dialog("Incorrect file name!")
        end
    end

    signal_connect(bt_prev, "clicked") do widget
        ep = parse(Int64, get_gtk_property(entry_epoch, :label, String))
        if ep >= 2
            ep -= 1
            Gtk.@sigatom begin
                set_gtk_property!(entry_epoch, :label, string(ep))
            end
        end
        draw(can)
    end

    signal_connect(bt_next, "clicked") do widget
        ep = parse(Int64, get_gtk_property(entry_epoch, :label, String))
        if ep < epoch_n(obj)
            ep += 1
            Gtk.@sigatom begin
                set_gtk_property!(entry_epoch, :label, string(ep))
            end
        end
        draw(can)
    end

    signal_connect(bt_start, "clicked") do widget
        Gtk.@sigatom begin
            set_gtk_property!(entry_epoch, :label, string(1))
        end
        draw(can)
    end

    signal_connect(bt_end, "clicked") do widget
        Gtk.@sigatom begin
            set_gtk_property!(entry_epoch, :label, string(epoch_n(obj)))
        end
        draw(can)
    end

    signal_connect(entry_epoch, "clicked") do widget
        value = parse(Int64, get_gtk_property(entry_epoch, :label, String))
        d_w = GtkWindow("Enter value", 200, 100)
        set_gtk_property!(d_w, :border_width, 20)
        set_gtk_property!(d_w, :resizable, true)
        d_g = GtkGrid()
        set_gtk_property!(d_g, :column_homogeneous, true)
        set_gtk_property!(g, :column_spacing, 10)
        set_gtk_property!(g, :row_spacing, 10)
        d_entry = GtkEntry()
        set_gtk_property!(d_entry, :text, string(value))
        d_bt_ok = GtkButton("Ok")
        d_bt_cancel = GtkButton("Cancel")
        d_g[1:2, 1] = d_entry
        d_g[1, 2] = d_bt_ok
        d_g[2, 2] = d_bt_cancel
        push!(d_w, d_g)
        showall(d_w)
        signal_connect(d_bt_ok, "clicked") do widget
            value_s = get_gtk_property(d_entry, :text, String)
            value_currect = true
            for idx in eachindex(value_s)
                string(value_s[idx]) in string.(0:9) || (value_currect = false)
            end
            if value_currect
                v = parse(Int64, value_s)
                if v < 1
                    warn_dialog("Value must be ≥ 1.")
                elseif v > epoch_n(obj)
                    warn_dialog("Value must be ≤ $(epoch_n(obj)).")
                else
                    set_gtk_property!(entry_epoch, :label, value_s)
                    draw(can)
                    Gtk.destroy(d_w)
                end
            else
                warn_dialog("Incorrect value entered!")
            end
        end
        signal_connect(d_w, "key-press-event") do widget, event
            k = event.keyval
            if k == 65293 || k == 65421
                value_s = get_gtk_property(d_entry, :text, String)
                value_currect = true
                for idx in eachindex(value_s)
                    string(value_s[idx]) in string.(0:9) || (value_currect = false)
                end
                if value_currect
                    v = parse(Int64, value_s)
                    if v < 1
                        warn_dialog("Value must be ≥ 1.")
                    elseif v > epoch_n(obj)
                        warn_dialog("Value must be ≤ $(epoch_n(obj)).")
                    else
                        set_gtk_property!(entry_epoch, :label, value_s)
                        draw(can)
                        Gtk.destroy(d_w)
                    end
                else
                    warn_dialog("Incorrect value entered!")
                end
            end
        end
        signal_connect(d_bt_cancel, "clicked") do widget
            Gtk.destroy(d_w)
        end
    end

    signal_connect(bt_close, "clicked") do widget
        Gtk.destroy(win)
    end

    signal_connect(bt_help, "clicked") do widgete
        info_dialog("Keyboard shortcuts:\n\nHOME\tgo to first epoch\nEND\t\tgo to last epoch\n,\t\tprevious epoch\n.\t\tnext epoch\n\nh\t\tthis info\nq\t\texit\n")
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        if k == 113 # q
            Gtk.destroy(win)
        elseif k == 104 # h
            info_dialog("Keyboard shortcuts:\n\nHOME\tgo to first epoch\nEND\t\tgo to last epoch\n,\t\tprevious epoch\n.\t\tnext epoch\n\nh\t\tthis info\nq\t\texit\n")
        elseif k == 103 # g
            value = parse(Int64, get_gtk_property(entry_epoch, :label, String))
            d_w = GtkWindow("Enter value", 200, 100)
            set_gtk_property!(d_w, :border_width, 20)
            set_gtk_property!(d_w, :resizable, true)
            d_g = GtkGrid()
            set_gtk_property!(d_g, :column_homogeneous, true)
            set_gtk_property!(g, :column_spacing, 10)
            set_gtk_property!(g, :row_spacing, 10)
            d_entry = GtkEntry()
            set_gtk_property!(d_entry, :text, string(value))
            d_bt_ok = GtkButton("Ok")
            d_bt_cancel = GtkButton("Cancel")
            d_g[1:2, 1] = d_entry
            d_g[1, 2] = d_bt_ok
            d_g[2, 2] = d_bt_cancel
            push!(d_w, d_g)
            showall(d_w)
            signal_connect(d_bt_ok, "clicked") do widget
                value_s = get_gtk_property(d_entry, :text, String)
                value_currect = true
                for idx in eachindex(value_s)
                    string(value_s[idx]) in string.(0:9) || (value_currect = false)
                end
                if value_currect
                    v = parse(Int64, value_s)
                    if v < 1
                        warn_dialog("Value must be ≥ 1.")
                    elseif v > epoch_n(obj)
                        warn_dialog("Value must be ≤ $(epoch_n(obj)).")
                    else
                        set_gtk_property!(entry_epoch, :label, value_s)
                        draw(can)
                        Gtk.destroy(d_w)
                    end
                else
                    warn_dialog("Incorrect value entered!")
                end
            end
            signal_connect(d_w, "key-press-event") do widget, event
                k = event.keyval
                if k == 65293 || k == 65421
                    value_s = get_gtk_property(d_entry, :text, String)
                    value_currect = true
                    for idx in eachindex(value_s)
                        string(value_s[idx]) in string.(0:9) || (value_currect = false)
                    end
                    if value_currect
                        v = parse(Int64, value_s)
                        if v < 1
                            warn_dialog("Value must be ≥ 1.")
                        elseif v > epoch_n(obj)
                            warn_dialog("Value must be ≤ $(epoch_n(obj)).")
                        else
                            set_gtk_property!(entry_epoch, :label, value_s)
                            draw(can)
                            Gtk.destroy(d_w)
                        end
                    else
                        warn_dialog("Incorrect value entered!")
                    end
                end
            end
            signal_connect(d_bt_cancel, "clicked") do widget
                Gtk.destroy(d_w)
            end
        elseif k == 65360 # HOME
            Gtk.@sigatom begin
                set_gtk_property!(entry_epoch, :label, string(1))
            end
            draw(can)
        elseif k == 65367 # END
            Gtk.@sigatom begin
                set_gtk_property!(entry_epoch, :label, string(epoch_n(obj)))
            end
            draw(can)
        elseif k == 44 # ,
            ep = parse(Int64, get_gtk_property(entry_epoch, :label, String))
            if ep >= 2
                ep -= 1
                Gtk.@sigatom begin
                    set_gtk_property!(entry_epoch, :label, string(ep))
                end
            end
            draw(can)
        elseif k == 46 # .
            ep = parse(Int64, get_gtk_property(entry_epoch, :label, String))
            if ep < epoch_n(obj)
                ep += 1
                Gtk.@sigatom begin
                    set_gtk_property!(entry_epoch, :label, string(ep))
                end
            end
            draw(can)
        end
    end

    return nothing

end