export iview
export iview_ep

"""
    iview(obj; <keyword arguments>)

Interactive view of continuous signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `mch::Bool=true`: draw multichannel signal (up to 20 channels in one plot)
- `zoom::Real=10`: how many seconds are displayed in one segment
- `bad::Bool=true`: list of bad channels; if not false -- plot bad channels using this list
- `snap::Bool=true`: snap region markers to grid at 0.0, 0.25, 0.5 and 0.75 time points

# Returns

- `seg::Union{Nothing, Tuple{Float64, Float64}}`
"""
function iview(obj::NeuroAnalyzer.NEURO; mch::Bool=true, zoom::Real=10, bad::Bool=true, snap::Bool=true)::Union{Nothing, Tuple{Float64, Float64}}

    obj.time_pts[end] < zoom && (zoom = round(obj.time_pts[end]) / 2)

    @assert zoom > 0 "zoom must be > 0."
    @assert zoom <= signal_len(obj) / sr(obj) "zoom must be ≤ $(signal_len(obj) / sr(obj))."

    ch_order = obj.header.recording[:channel_order]
    cl = labels(obj)[ch_order]
    ch = 1:nchannels(obj)

    mono = false
    quit = false
    scale = true
    plot_type = 0

    if mch
        if length(ch) > 20
            ch_first = 1
            ch_last = ch_first + 19
        else
            ch_first = 1
            ch_last = ch[end]
        end
    else
        ch_first = 1
        ch_last = ch[end]
    end

    if mch
        p = NeuroAnalyzer.plot(obj, ch=cl[ch_first:ch_last], title="", bad=bad)
    else
        p = NeuroAnalyzer.plot(obj, ch=cl[ch_first], title="Channel: $(cl[ch_first])")
    end
    win = GtkWindow("NeuroAnalyzer: iview()", Int32(p.attr[:size][1]) + 40, Int32(p.attr[:size][2]) + 40)
    set_gtk_property!(win, :border_width, 5)
    set_gtk_property!(win, :resizable, true)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    set_gtk_property!(win, :startup_id, "org.neuroanalyzer")
    can = GtkCanvas(Int32(p.attr[:size][1]), Int32(p.attr[:size][2]))
    g = GtkGrid()
    set_gtk_property!(g, :column_homogeneous, false)
    set_gtk_property!(g, :column_spacing, 5)
    set_gtk_property!(g, :row_spacing, 5)
    entry_time = GtkSpinButton(obj.time_pts[1], obj.time_pts[end] - zoom, 1)
    set_gtk_property!(entry_time, :digits, 2)
    set_gtk_property!(entry_time, :value, obj.time_pts[1])
    set_gtk_property!(entry_time, :tooltip_text, "Time position [s]")
    bt_start = GtkButton("⇤")
    set_gtk_property!(bt_start, :tooltip_text, "Go to the signal start")
    bt_prev = GtkButton("←")
    set_gtk_property!(bt_prev, :tooltip_text, "Go back by $(round(zoom)) seconds")
    bt_next = GtkButton("→")
    set_gtk_property!(bt_next, :tooltip_text, "Go forward by $(round(zoom)) seconds")
    bt_end = GtkButton("⇥")
    set_gtk_property!(bt_end, :tooltip_text, "Go to the signal end")
    entry_ts1 = GtkSpinButton(obj.time_pts[1], obj.time_pts[end], 0.5)
    bt_help = GtkButton("Help")
    set_gtk_property!(bt_help, :tooltip_text, "Show help")
    bt_close = GtkButton("Close")
    set_gtk_property!(bt_close, :tooltip_text, "Close this window")
    set_gtk_property!(entry_ts1, :tooltip_text, "Segment start [s]")
    set_gtk_property!(entry_ts1, :digits, 3)
    entry_ts2 = GtkSpinButton(obj.time_pts[1], obj.time_pts[end], 0.5)
    set_gtk_property!(entry_ts2, :digits, 3)
    set_gtk_property!(entry_ts2, :tooltip_text, "Segment end [s]")
    bt_ts = GtkButton("Return TS")
    set_gtk_property!(bt_ts, :tooltip_text, "Return selected time segment")
    bt_delete = GtkButton("Delete TS")
    set_gtk_property!(bt_delete, :tooltip_text, "Delete selected time segment")

    combo_ch = GtkComboBoxText()
    ch_types = uppercase.(unique(obj.header.recording[:channel_type]))
    for idx in ch_types
        length(get_channel(obj, type=lowercase(idx))) > 1 && push!(combo_ch, idx)
    end
    set_gtk_property!(combo_ch, :active, 0)
    set_gtk_property!(combo_ch, :sensitive, false)
    set_gtk_property!(combo_ch, :tooltip_text, "Channels")

    if !mch
        ch_slider = GtkScale(false, 1:nchannels(obj))
        set_gtk_property!(ch_slider, :draw_value, false)
    else
        if length(ch) > 20
            ch_slider = GtkScale(false, ch[ch_first]:(ch[end] - 19))
            set_gtk_property!(ch_slider, :draw_value, false)
        else
            ch_slider = GtkScale(false, ch[1]:ch[end])
            set_gtk_property!(ch_slider, :draw_value, false)
            set_gtk_property!(ch_slider, :sensitive, false)
        end
    end
    set_gtk_property!(ch_slider, :tooltip_text, "Scroll channels")
    set_gtk_property!(ch_slider, :vexpand, true)
    oc = GtkOrientable(ch_slider)
    set_gtk_property!(oc, :orientation, 1)

    signal_slider = GtkScale(false, obj.time_pts[1]:obj.time_pts[end] - zoom)
    set_gtk_property!(signal_slider, :draw_value, false)
    set_gtk_property!(signal_slider, :tooltip_text, "Time position")

    g[1:9, 1] = can
    g[10, 1] = ch_slider
    g[1:9, 2] = signal_slider
    g[1, 3] = bt_start
    g[2, 3] = bt_prev
    g[3, 3] = entry_time
    g[4, 3] = bt_next
    g[5, 3] = bt_end
    g[6, 3] = combo_ch
    g[7, 3] = entry_ts1
    g[8, 3] = bt_ts
    g[9, 3] = bt_help
    g[7, 4] = entry_ts2
    g[8, 4] = bt_delete
    g[9, 4] = bt_close
    push!(win, g)

    showall(win)

    @guarded draw(can) do widget
        time1 = get_gtk_property(entry_time, :value, Float64)
        time2 = time1 + zoom
        time2 > obj.time_pts[end] && (time2 = obj.time_pts[end])
        ts1 = get_gtk_property(entry_ts1, :value, Float64)
        ts2 = get_gtk_property(entry_ts2, :value, Float64)
        ctx = getgc(can)
        channel_type = lowercase.(ch_types)[get_gtk_property(combo_ch, :active, Int64) + 1]
        if mch
            if plot_type == 0
                p = NeuroAnalyzer.plot(obj,
                                       ch=cl[ch_first:ch_last],
                                       seg=(time1, time2),
                                       s_pos=(ts1, ts2),
                                       mono=mono,
                                       title="",
                                       scale=scale)
            elseif plot_type == 1
                p =  NeuroAnalyzer.plot(obj,
                                        ch=get_channel(obj, type=channel_type),
                                        seg=(time1, time2),
                                        s_pos=(ts1, ts2),
                                        mono=mono,
                                        type=:butterfly,
                                        avg=false)
            elseif plot_type == 2
                p = NeuroAnalyzer.plot(obj,
                                       ch=get_channel(obj, type=channel_type),
                                       seg=(time1, time2),
                                       s_pos=(ts1, ts2),
                                       mono=mono,
                                       type=:mean)
            end
        else
            p = NeuroAnalyzer.plot(obj,
                                   ch=cl[ch_first],
                                   seg=(time1, time2),
                                   s_pos=(ts1, ts2),
                                   mono=mono,
                                   title="Channel: $(cl[ch_first])",
                                   scale=scale,
                                   bad=bad)
        end
        show(io, MIME("image/png"), p)
        Gtk.resize!(win, p.attr[:size][2] + 40, p.attr[:size][2] + 40)
        set_gtk_property!(can, :width_request, Int32(p.attr[:size][1]))
        set_gtk_property!(can, :height_request, Int32(p.attr[:size][2]))
        img = read_from_png(io)
        set_source_surface(ctx, img, 0, 0)
        paint(ctx)
    end

    can.mouse.button1press = @guarded (widget, event) -> begin
        x_pos = event.x
        if x_pos < 82
            if mch
                y_pos = event.y
                ch_y = collect(50:39:793)[1:length(ch_first:ch_last)]
                ch_idx = nothing
                for idx in eachindex(ch_y)
                    if y_pos > ch_y[idx] && y_pos < ch_y[idx] + 15 && x_pos > 0 && x_pos < 82
                        ch_idx = idx + ch_first - 1
                    end
                end
                !isnothing(ch_idx) && channel_info(obj, ch=obj.header.recording[:channel_order][ch_idx])
            end
        else
            time_current = get_gtk_property(entry_time, :value, Float64)
            x_pos > 1172 && (x_pos = 1172)
            if time_current + zoom < obj.time_pts[end]
                ts1 = time_current + round((x_pos - 82) / (1090 / zoom), digits=3)
            else
                ts1 = time_current + round((x_pos - 82) / (1090 / (obj.time_pts[end] - time_current)), digits=3)
            end
            snap && (ts1 = round(ts1 * 4) / 4)
            round(ts1, digits=3)
            Gtk.@sigatom begin
                set_gtk_property!(entry_ts1, :value, ts1)
            end
        end
    end

    can.mouse.button3press = @guarded (widget, event) -> begin
        x_pos = event.x
        if x_pos < 82
            if mch
                y_pos = event.y
                ch_y = collect(50:39:793)[1:length(ch_first:ch_last)]
                ch_idx = nothing
                for idx in eachindex(ch_y)
                    if y_pos > ch_y[idx] && y_pos < ch_y[idx] + 15 && x_pos > 0 && x_pos < 82
                        ch_idx = idx + ch_first - 1
                    end
                end
                !isnothing(ch_idx) && (obj.header.recording[:bad_channel][obj.header.recording[:channel_order][ch_idx, 1]] = !obj.header.recording[:bad_channel][obj.header.recording[:channel_order][ch_idx, 1]])
                draw(can)
            end
        else
            time_current = get_gtk_property(entry_time, :value, Float64)
            x_pos = event.x
            x_pos > 1172 && (x_pos = 1172)
            if time_current + zoom < obj.time_pts[end]
                ts2 = time_current + ((x_pos - 82) / (1090 / zoom))
            else
                ts2 = time_current + ((x_pos - 82) / (1090 / (obj.time_pts[end] - time_current)))
            end
            snap && (ts2 = round(ts2 * 4) / 4)
            round(ts2, digits=3)
            Gtk.@sigatom begin
                set_gtk_property!(entry_ts2, :value, round(ts2, digits=3))
            end
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
            else
                if mch
                    if ch_last < length(ch)
                        ch_first += 1
                        ch_last += 1
                        Gtk.@sigatom begin
                            GAccessor.value(ch_slider, ch_first)
                        end
                        draw(can)
                    end
                else
                    if ch_first < length(ch)
                        ch_first += 1
                        Gtk.@sigatom begin
                            GAccessor.value(ch_slider, ch_first)
                        end
                        draw(can)
                    end
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
            else
                if ch_first > 1
                    ch_first -= 1
                    ch_last -= 1
                    Gtk.@sigatom begin
                        GAccessor.value(ch_slider, ch_first)
                    end
                    draw(can)
                end
            end
        end
    end

    @guarded signal_connect(bt_delete, "clicked") do widget
        time_current = get_gtk_property(entry_time, :value, Float64)
        time1 = obj.time_pts[vsearch(get_gtk_property(entry_ts1, :value, Float64), obj.time_pts)]
        time2 = obj.time_pts[vsearch(get_gtk_property(entry_ts2, :value, Float64), obj.time_pts)]
        if time1 > time2
            warn_dialog("Cannot delete!\nSegment start is larger than segment end.")
        elseif time1 == time2
            warn_dialog("Cannot delete!\nSegment start must be different from segment end.")
        elseif ask_dialog("Delete segment $time1:$time2 ?", "No", "Yes")
            trim!(obj, seg=(time1, time2), remove_epochs=false)
            _info("Deleted segment: $time1:$time2")
            if time1 == time_current && time2 > obj.time_pts[end]
                time_current = obj.time_pts[end] - zoom
                time_current < obj.time_pts[1] && (time_current = obj.time_pts[1])
            else
                if obj.time_pts[end] % zoom == 0
                    time_current >= (obj.time_pts[end] - zoom) && (time_current = obj.time_pts[end] - zoom)
                else
                    time_current >= (obj.time_pts[end] - (obj.time_pts[end] % zoom)) && (time_current = obj.time_pts[end] - (obj.time_pts[end] % zoom))
                end
                time_current < obj.time_pts[1] && (time_current = obj.time_pts[1])
            end
            Gtk.@sigatom begin
                set_gtk_property!(entry_time, :value, time_current)
                set_gtk_property!(entry_ts1, :value, time_current)
                set_gtk_property!(entry_ts2, :value, time_current)
                GAccessor.range(signal_slider, obj.time_pts[1], obj.time_pts[end] - zoom)
                GAccessor.range(entry_time, obj.time_pts[1], obj.time_pts[end] - zoom)
                GAccessor.range(entry_ts1, obj.time_pts[1], obj.time_pts[end])
                GAccessor.range(entry_ts2, obj.time_pts[1], obj.time_pts[end])
            end
        end
    end

    signal_connect(signal_slider, "value-changed") do widget, others...
        Gtk.@sigatom begin
            set_gtk_property!(entry_time, :value, round(GAccessor.value(signal_slider)))
        end
        draw(can)
    end

    signal_connect(ch_slider, "value-changed") do widget, others...
        ch_first = round(Int64, GAccessor.value(ch_slider))
        length(ch) > 20 && (ch_last = ch_first + 19)
        draw(can)
    end

    signal_connect(entry_ts1, "value-changed") do widget
        Gtk.@sigatom begin
            set_gtk_property!(entry_ts1, :value, obj.time_pts[vsearch(get_gtk_property(entry_ts1, :value, Float64), obj.time_pts)])
        end
        draw(can)
    end

    signal_connect(entry_ts2, "value-changed") do widget
        Gtk.@sigatom begin
            set_gtk_property!(entry_ts2, :value, obj.time_pts[vsearch(get_gtk_property(entry_ts2, :value, Float64), obj.time_pts)])
        end
        draw(can)
    end

    signal_connect(entry_time, "value-changed") do widget
        Gtk.@sigatom begin
            GAccessor.value(signal_slider, get_gtk_property(entry_time, :value, Float64))
        end
        draw(can)
    end

    signal_connect(bt_prev, "clicked") do widget
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
            time_current += zoom
        else
            time_current = obj.time_pts[end] - zoom
        end
        Gtk.@sigatom begin
            set_gtk_property!(entry_time, :value, time_current)
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

    signal_connect(bt_ts, "clicked") do widget
        Gtk.destroy(win)
    end

    signal_connect(bt_close, "clicked") do widget
        quit = true
        Gtk.destroy(win)
    end

    if mch
        help = "Keyboard shortcuts:\n\nCtrl + b\t\t\tToggle butterfly plot\nCtrl + m\t\tToggle mean plot\n\nPage Up\t\tScroll channels up\nPage Down\tScroll channels down\n\nHome\t\t\tGo to the signal start\nEnd\t\t\tGo to the signal end\nCtrl + ,\t\t\tGo back by 1 second\nCtrl + .\t\t\tGo forward by 1 second\nAlt + ,\t\t\tGo back by $(round(zoom)) seconds\nAlt + .\t\t\tGo forward by $(round(zoom)) seconds\n\n[\t\t\t\tZoom in\n]\t\t\t\tZoom out\n\nCtrl + Enter\tReturn selected time segment\nCtrl + d\t\t\tDelete selected time segment\n\nCtrl + s\t\t\tToggle snapping\nAlt + s\t\t\tToggle scales\nAlt + m\t\t\tToggle monochromatic mode\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tClose\n"
    else
        help = "Keyboard shortcuts:\n\nPage Up\t\tScroll channels up\nPage Down\tScroll channels down\n\nHome\t\t\tGo to the signal start\nEnd\t\t\tGo to the signal end\nCtrl + ,\t\t\tGo back by 1 second\nCtrl + .\t\t\tGo forward by 1 second\nAlt + ,\t\t\tGo back by $(round(zoom)) seconds\nAlt + .\t\t\tGo forward by $(round(zoom)) seconds\n\n[\t\t\t\tZoom in\n]\t\t\t\tZoom out\n\nCtrl + Enter\tReturn selected time segment\nCtrl + d\t\t\tDelete selected time segment\n\nCtrl + s\t\t\tToggle snapping\nAlt + s\t\t\tToggle scales\nAlt + m\t\t\tToggle monochromatic mode\n\nCtrl + p\t\t\tPlay current segment as audio\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tClose\n"
    end

    signal_connect(combo_ch, "changed") do widgete
        draw(can)
    end

    signal_connect(bt_help, "clicked") do widgete
        info_dialog(help)
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        s = event.state
        if k == 0x0000ff55 # Page Up
            if ch_first > 1
                ch_first -= 1
                ch_last -= 1
                Gtk.@sigatom begin
                    GAccessor.value(ch_slider, ch_first)
                end
                draw(can)
            end
        elseif k == 0x0000ff56 # Page Down
            if mch
                if ch_last < length(ch)
                    ch_first += 1
                    ch_last += 1
                    Gtk.@sigatom begin
                        GAccessor.value(ch_slider, ch_first)
                    end
                    draw(can)
                end
            else
                if ch_first < length(ch)
                    ch_first += 1
                    Gtk.@sigatom begin
                        GAccessor.value(ch_slider, ch_first)
                    end
                    draw(can)
                end
            end
        elseif k == 0x0000005b # [
            if zoom > 1
                zoom -= 1
                Gtk.@sigatom begin
                    set_gtk_property!(bt_next, :tooltip_text, "Go forward by $(round(zoom)) seconds")
                    set_gtk_property!(bt_prev, :tooltip_text, "Go back by $(round(zoom)) seconds")
                    GAccessor.range(signal_slider, obj.time_pts[1], obj.time_pts[end] - zoom)
                    GAccessor.range(entry_time, obj.time_pts[1], obj.time_pts[end] - zoom)
                end
                draw(can)
            end
            if mch
                help = "Keyboard shortcuts:\n\nPage Up\t\tScroll channels up\nPage Down\tScroll channels down\n\nHome\t\t\tGo to the signal start\nEnd\t\t\tGo to the signal end\nCtrl + ,\t\t\tGo back by 1 second\nCtrl + .\t\t\tGo forward by 1 second\nAlt + ,\t\t\tGo back by $(round(zoom)) seconds\nAlt + .\t\t\tGo forward by $(round(zoom)) seconds\n\n[\t\t\t\tZoom in\n]\t\t\t\tZoom out\n\nCtrl + Enter\tReturn selected time segment\nCtrl + d\t\t\tDelete selected time segment\n\nCtrl + s\t\t\tToggle snapping\nAlt + s\t\t\tToggle scales\nAlt + m\t\t\tToggle monochromatic mode\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tClose\n"
            else
                help = "Keyboard shortcuts:\n\nPage Up\t\tScroll channels up\nPage Down\tScroll channels down\n\nHome\t\t\tGo to the signal start\nEnd\t\t\tGo to the signal end\nCtrl + ,\t\t\tGo back by 1 second\nCtrl + .\t\t\tGo forward by 1 second\nAlt + ,\t\t\tGo back by $(round(zoom)) seconds\nAlt + .\t\t\tGo forward by $(round(zoom)) seconds\n\n[\t\t\t\tZoom in\n]\t\t\t\tZoom out\n\nCtrl + Enter\tReturn selected time segment\nCtrl + d\t\t\tDelete selected time segment\n\nCtrl + s\t\t\tToggle snapping\nAlt + s\t\t\tToggle scales\nAlt + m\t\t\tToggle monochromatic mode\n\nCtrl + p\t\t\tPlay current segment as audio\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tClose\n"

            end
        elseif k == 0x0000005d # ]
            if zoom < 30 && zoom < obj.time_pts[end] - 1
                zoom += 1
                Gtk.@sigatom begin
                    set_gtk_property!(bt_next, :tooltip_text, "Go forward by $(round(zoom)) seconds")
                    set_gtk_property!(bt_prev, :tooltip_text, "Go back by $(round(zoom)) seconds")
                    GAccessor.range(signal_slider, obj.time_pts[1], obj.time_pts[end] - zoom)
                    GAccessor.range(entry_time, obj.time_pts[1], obj.time_pts[end] - zoom)
                end
                draw(can)
            else
                zoom = obj.time_pts[end]
                Gtk.@sigatom begin
                    set_gtk_property!(bt_next, :tooltip_text, "Go forward by $(round(zoom)) seconds")
                    set_gtk_property!(bt_prev, :tooltip_text, "Go back by $(round(zoom)) seconds")
                    GAccessor.range(signal_slider, obj.time_pts[1], obj.time_pts[end] - zoom)
                    GAccessor.range(entry_time, obj.time_pts[1], obj.time_pts[end] - zoom)
                end
                draw(can)
            end
            if mch
                help = "Keyboard shortcuts:\n\nPage Up\t\tScroll channels up\nPage Down\tScroll channels down\n\nHome\t\t\tGo to the signal start\nEnd\t\t\tGo to the signal end\nCtrl + ,\t\t\tGo back by 1 second\nCtrl + .\t\t\tGo forward by 1 second\nAlt + ,\t\t\tGo back by $(round(zoom)) seconds\nAlt + .\t\t\tGo forward by $(round(zoom)) seconds\n\n[\t\t\t\tZoom in\n]\t\t\t\tZoom out\n\nCtrl + Enter\tReturn selected time segment\nCtrl + d\t\t\tDelete selected time segment\n\nCtrl + s\t\t\tToggle snapping\nAlt + s\t\t\tToggle scales\nAlt + m\t\t\tToggle monochromatic mode\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tClose\n"
            else
                help = "Keyboard shortcuts:\n\nPage Up\t\tScroll channels up\nPage Down\tScroll channels down\n\nHome\t\t\tGo to the signal start\nEnd\t\t\tGo to the signal end\nCtrl + ,\t\t\tGo back by 1 second\nCtrl + .\t\t\tGo forward by 1 second\nAlt + ,\t\t\tGo back by $(round(zoom)) seconds\nAlt + .\t\t\tGo forward by $(round(zoom)) seconds\n\n[\t\t\t\tZoom in\n]\t\t\t\tZoom out\n\nCtrl + Enter\tReturn selected time segment\nCtrl + d\t\t\tDelete selected time segment\n\nCtrl + s\t\t\tToggle snapping\nAlt + s\t\t\tToggle scales\nAlt + m\t\t\tToggle monochromatic mode\n\nCtrl + p\t\t\tPlay current segment as audio\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tClose\n"

            end
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

        if s == 0x00000008 || s == 0x00000010 # alt
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
            elseif k == 0x0000006d # m
                mono = !mono
                draw(can)
            elseif k == 0x00000073 # s
                scale = !scale
                draw(can)
            end
        end

        if s == 0x00000004 || s == 0x00000014 # ctrl
            if k == 0x00000071 # q
                quit = true
                Gtk.destroy(win)
            elseif k == 104 # h
                info_dialog(help)
            elseif k == 0x00000062 # b
                if mch
                    if plot_type != 1
                        Gtk.@sigatom begin
                            set_gtk_property!(ch_slider, :sensitive, false)
                            set_gtk_property!(combo_ch, :sensitive, true)
                        end
                        plot_type = 1
                    else
                        set_gtk_property!(ch_slider, :sensitive, true)
                        set_gtk_property!(combo_ch, :sensitive, false)
                        plot_type = 0
                    end
                    draw(can)
                end
            elseif k == 0x0000006d # m
                if mch
                    if plot_type != 2
                        Gtk.@sigatom begin
                            set_gtk_property!(ch_slider, :sensitive, false)
                            set_gtk_property!(combo_ch, :sensitive, true)
                        end
                        plot_type = 2
                    else
                        set_gtk_property!(ch_slider, :sensitive, true)
                        set_gtk_property!(combo_ch, :sensitive, false)
                        plot_type = 0
                    end
                    draw(can)
                end
            elseif k == 0x00000070 # p
                if !mch
                    _info("Playing current segment as audio")
                    time_current = get_gtk_property(entry_time, :value, Float64)
                    play(obj, ch=cl[ch_first], seg=(time_current, time_current+zoom), ep=1)
                end
            elseif k == 0x0000ff0d # Enter
                Gtk.destroy(win)
            elseif k == 104 # h
                info_dialog(help)
            elseif k == 0x00000073 # s
                snap = !snap
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
            elseif k == 100 # d
                time_current = get_gtk_property(entry_time, :value, Float64)
                time1 = obj.time_pts[vsearch(get_gtk_property(entry_ts1, :value, Float64), obj.time_pts)]
                time2 = obj.time_pts[vsearch(get_gtk_property(entry_ts2, :value, Float64), obj.time_pts)]
                if time1 > time2
                    warn_dialog("Cannot delete!\nSegment start is larger than segment end.")
                elseif time1 == time2
                    warn_dialog("Cannot delete!\nSegment start must be different from segment end.")
                elseif ask_dialog("Delete segment $time1:$time2 ?", "No", "Yes")
                    trim!(obj, seg=(time1, time2), remove_epochs=false)
                    _info("Deleted segment: $time1:$time2")
                    if time1 == time_current && time2 > obj.time_pts[end]
                        time_current = obj.time_pts[end] - zoom
                        time_current < obj.time_pts[1] && (time_current = obj.time_pts[1])
                    else
                        if obj.time_pts[end] % zoom == 0
                            time_current >= (obj.time_pts[end] - zoom) && (time_current = obj.time_pts[end] - zoom)
                        else
                            time_current >= obj.time_pts[end] - (obj.time_pts[end] % zoom) && (time_current = obj.time_pts[end] - (obj.time_pts[end] % zoom))
                        end
                        time_current < obj.time_pts[1] && (time_current = obj.time_pts[1])
                    end
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_time, :value, time_current)
                        set_gtk_property!(entry_ts1, :value, time_current)
                        set_gtk_property!(entry_ts2, :value, time_current)
                        GAccessor.range(entry_time, obj.time_pts[1], obj.time_pts[end] - zoom)
                        GAccessor.range(entry_ts1, obj.time_pts[1], obj.time_pts[end])
                        GAccessor.range(entry_ts2, obj.time_pts[1], obj.time_pts[end])
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

    seg = nothing
    time1 = obj.time_pts[vsearch(get_gtk_property(entry_ts1, :value, Float64), obj.time_pts)]
    time2 = obj.time_pts[vsearch(get_gtk_property(entry_ts2, :value, Float64), obj.time_pts)]
    seg = (time1, time2)
    seg == (0.0, 0.0) && (seg = nothing)

    return seg

end

"""
    iview_ep(obj; <keyword arguments>)

Interactive view of epoched signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `mch::Bool=true`: draw multichannel signal (up to 20 channels in one plot)
- `ep::Int64=1`: initial epoch to display
- `bad::Bool=true`: list of bad channels; if not false -- plot bad channels using this list
- `snap::Bool=true`: snap region markers to grid at 0.0, 0.25, 0.5 and 0.75 time points

# Returns

- `seg::Union{Nothing, Tuple{Float64, Float64}}`
"""
function iview_ep(obj::NeuroAnalyzer.NEURO; mch::Bool=true, ep::Int64=1, bad::Bool=true, snap::Bool=true)::Union{Nothing, Tuple{Float64, Float64}}

    @assert nepochs(obj) > 1 "iview() must be used for continuous object."
    _check_epochs(obj, ep)

    ch_order = obj.header.recording[:channel_order]
    cl = labels(obj)[ch_order]
    ch = 1:nchannels(obj)

    mono = false
    quit = false
    scale = true
    plot_type = 0
    zoom = epoch_len(obj)

    if mch
        if length(ch) > 20
            ch_first = 1
            ch_last = ch_first + 19
        else
            ch_first = 1
            ch_last = ch[end]
        end
    else
        ch_first = 1
        ch_last = ch[end]
    end

    if mch
        p = NeuroAnalyzer.plot(obj, ch=cl[ch_first:ch_last], ep=ep, title="", bad=bad)
    else
        p = NeuroAnalyzer.plot(obj, ch=cl[ch_first], ep=ep, title="Channel: $(cl[ch_first])")
    end
    win = GtkWindow("NeuroAnalyzer: iview_ep()", Int32(p.attr[:size][1]) + 40, Int32(p.attr[:size][2]) + 40)
    set_gtk_property!(win, :border_width, 5)
    set_gtk_property!(win, :resizable, true)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    set_gtk_property!(win, :startup_id, "org.neuroanalyzer")
    can = GtkCanvas(Int32(p.attr[:size][1]), Int32(p.attr[:size][2]))
    g = GtkGrid()
    set_gtk_property!(g, :column_homogeneous, false)
    set_gtk_property!(g, :column_spacing, 5)
    set_gtk_property!(g, :row_spacing, 5)
    entry_epoch = GtkSpinButton(1, nepochs(obj), 1)
    set_gtk_property!(entry_epoch, :digits, 0)
    set_gtk_property!(entry_epoch, :value, ep)
    set_gtk_property!(entry_epoch, :tooltip_text, "Epoch")
    bt_start = GtkButton("⇤")
    set_gtk_property!(bt_start, :tooltip_text, "Go to the signal start")
    bt_end = GtkButton("⇥")
    set_gtk_property!(bt_end, :tooltip_text, "Go to the signal end")
    bt_help = GtkButton("Help")
    set_gtk_property!(bt_help, :tooltip_text, "Show help")
    bt_close = GtkButton("Close")
    set_gtk_property!(bt_close, :tooltip_text, "Close this window")
    entry_ts1 = GtkSpinButton(obj.time_pts[1], obj.time_pts[end], 0.5)
    set_gtk_property!(entry_ts1, :tooltip_text, "Segment start [s]")
    set_gtk_property!(entry_ts1, :digits, 3)
    entry_ts2 = GtkSpinButton(obj.time_pts[1], obj.time_pts[end], 0.5)
    set_gtk_property!(entry_ts2, :digits, 3)
    set_gtk_property!(entry_ts2, :tooltip_text, "Segment end [s]")
    bt_ts = GtkButton("Return TS")
    set_gtk_property!(bt_ts, :tooltip_text, "Return selected time segment")
    bt_delete = GtkButton("Delete epoch")
    set_gtk_property!(bt_delete, :tooltip_text, "Delete current epoch")
    if !mch
        ch_slider = GtkScale(false, 1:nchannels(obj))
        set_gtk_property!(ch_slider, :draw_value, false)
    else
        if length(ch) > 20
            ch_slider = GtkScale(false, ch[ch_first]:(ch[end] - 19))
            set_gtk_property!(ch_slider, :draw_value, false)
        else
            ch_slider = GtkScale(false, ch[1]:ch[end])
            set_gtk_property!(ch_slider, :draw_value, false)
            set_gtk_property!(ch_slider, :sensitive, false)
        end
    end

    combo_ch = GtkComboBoxText()
    ch_types = uppercase.(unique(obj.header.recording[:channel_type]))
    for idx in ch_types
        length(get_channel(obj, type=lowercase(idx))) > 1 && push!(combo_ch, idx)
    end
    set_gtk_property!(combo_ch, :active, 0)
    set_gtk_property!(combo_ch, :sensitive, false)
    set_gtk_property!(combo_ch, :tooltip_text, "Channels")

    set_gtk_property!(ch_slider, :tooltip_text, "Scroll channels")
    set_gtk_property!(ch_slider, :vexpand, true)
    oc = GtkOrientable(ch_slider)
    set_gtk_property!(oc, :orientation, 1)

    signal_slider = GtkScale(false, 1:1:nepochs(obj))
    set_gtk_property!(signal_slider, :draw_value, false)
    set_gtk_property!(signal_slider, :tooltip_text, "Current epoch")

    g[1:7, 1] = can
    g[8, 1] = ch_slider
    g[1:7, 2] = signal_slider
    g[1, 3] = bt_start
    g[2, 3] = entry_epoch
    g[3, 3] = bt_end
    g[4, 3] = combo_ch
    g[5, 3] = entry_ts1
    g[6, 3] = bt_ts
    g[7, 3] = bt_help
    g[5, 4] = entry_ts2
    g[6, 4] = bt_delete
    g[7, 4] = bt_close
    push!(win, g)

    showall(win)

    @guarded draw(can) do widget
        ep = get_gtk_property(entry_epoch, :value, Int64)
        ts1 = get_gtk_property(entry_ts1, :value, Float64)
        ts2 = get_gtk_property(entry_ts2, :value, Float64)
        ctx = getgc(can)
        channel_type = lowercase.(ch_types)[get_gtk_property(combo_ch, :active, Int64) + 1]
        if mch
            if plot_type == 0
                p = NeuroAnalyzer.plot(obj,
                                       ch=cl[ch_first:ch_last],
                                       ep=ep,
                                       s_pos=(ts1, ts2),
                                       mono=mono,
                                       title="",
                                       scale=scale)
            elseif plot_type == 1
                p =  NeuroAnalyzer.plot(obj,
                                        ch=get_channel(obj, type=channel_type),
                                        ep=ep,
                                        s_pos=(ts1, ts2),
                                        mono=mono,
                                        type=:butterfly,
                                        avg=false)
            elseif plot_type == 2
                p = NeuroAnalyzer.plot(obj,
                                       ch=get_channel(obj, type=channel_type),
                                       ep=ep,
                                       s_pos=(ts1, ts2),
                                       mono=mono,
                                       type=:mean)
            end
        else
            p = NeuroAnalyzer.plot(obj,
                                   ch=cl[ch_first],
                                   ep=ep,
                                   s_pos=(ts1, ts2),
                                   mono=mono,
                                   title="Channel: $(cl[ch_first])",
                                   scale=scale,
                                   bad=bad)
        end
        show(io, MIME("image/png"), p)
        Gtk.resize!(win, p.attr[:size][2] + 40, p.attr[:size][2] + 40)
        set_gtk_property!(can, :width_request, Int32(p.attr[:size][1]))
        set_gtk_property!(can, :height_request, Int32(p.attr[:size][2]))
        img = read_from_png(io)
        set_source_surface(ctx, img, 0, 0)
        paint(ctx)
    end

    can.mouse.button1press = @guarded (widget, event) -> begin
        x_pos = event.x
        if x_pos < 82
            if mch
                y_pos = event.y
                ch_y = collect(50:39:793)[1:length(ch_first:ch_last)]
                ch_idx = nothing
                for idx in eachindex(ch_y)
                    if y_pos > ch_y[idx] && y_pos < ch_y[idx] + 15 && x_pos > 0 && x_pos < 82
                        ch_idx = idx + ch_first - 1
                    end
                end
                !isnothing(ch_idx) && channel_info(obj, ch=obj.header.recording[:channel_order][ch_idx])
            end
        else
            x_pos > 1172 && (x_pos = 1172)
            ep = get_gtk_property(entry_epoch, :value, Int64)
            ts1 = ((epoch_len(obj) / sr(obj)) * (ep - 1)) + obj.epoch_time[1] + (x_pos - 82) / (1090 / (obj.epoch_time[end] - obj.epoch_time[1]))
            snap && (ts1 = round(ts1 * 4) / 4)
            round(ts1, digits=3)
            Gtk.@sigatom begin
                set_gtk_property!(entry_ts1, :value, ts1)
            end
        end
    end

    can.mouse.button3press = @guarded (widget, event) -> begin
        x_pos = event.x
        if x_pos < 82
            if mch
                ep = get_gtk_property(entry_epoch, :value, Int64)
                y_pos = event.y
                ch_y = collect(50:39:793)[1:length(ch_first:ch_last)]
                ch_idx = nothing
                for idx in eachindex(ch_y)
                    if y_pos > ch_y[idx] && y_pos < ch_y[idx] + 15 && x_pos > 0 && x_pos < 82
                        ch_idx = idx + ch_first - 1
                    end
                end
                !isnothing(ch_idx) && (obj.header.recording[:bad_channel][obj.header.recording[:channel_order][ch_idx], ep] = !obj.header.recording[:bad_channel][obj.header.recording[:channel_order][ch_idx], ep])
                draw(can)
            end
        else
            x_pos > 1172 && (x_pos = 1172)
            ep = get_gtk_property(entry_epoch, :value, Int64)
            ts2 = ((epoch_len(obj) / sr(obj)) * (ep - 1)) + obj.epoch_time[1] + ((x_pos - 82) / (1090 / (obj.epoch_time[end] - obj.epoch_time[1])))
            snap && (ts2 = round(ts2 * 4) / 4)
            round(ts2, digits=3)
            Gtk.@sigatom begin
                set_gtk_property!(entry_ts2, :value, round(ts2, digits=3))
            end
        end
    end

    can.mouse.scroll = @guarded (widget, event) -> begin
        s = event.state
        if event.direction == 1 # down
            if s == 0x00000001
                ep = get_gtk_property(entry_epoch, :value, Int64)
                if ep < nepochs(obj)
                    ep += 1
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_epoch, :value, ep)
                    end
                end
            else
                if mch
                    if ch_last < length(ch)
                        ch_first += 1
                        ch_last += 1
                        Gtk.@sigatom begin
                            GAccessor.value(ch_slider, ch_first)
                        end
                        draw(can)
                    end
                else
                    if ch_first < length(ch)
                        ch_first += 1
                        Gtk.@sigatom begin
                            GAccessor.value(ch_slider, ch_first)
                        end
                        draw(can)
                    end
                end
            end
        elseif event.direction == 0 # up
            if s == 0x00000001
                ep = get_gtk_property(entry_epoch, :value, Int64)
                if ep > 1
                    ep -= 1
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_epoch, :value, ep)
                    end
                end
            else
                if ch_first > 1
                    ch_first -= 1
                    ch_last -= 1
                    Gtk.@sigatom begin
                        GAccessor.value(ch_slider, ch_first)
                    end
                    draw(can)
                end
            end
        end
    end

    signal_connect(signal_slider, "value-changed") do widget, others...
        Gtk.@sigatom begin
            set_gtk_property!(entry_epoch, :value, round(Int64, GAccessor.value(signal_slider)))
        end
        draw(can)
    end

    signal_connect(ch_slider, "value-changed") do widget, others...
        ch_first = round(Int64, GAccessor.value(ch_slider))
        ch_last = ch_first + 19
        draw(can)
    end

    signal_connect(entry_epoch, "value-changed") do widget
        Gtk.@sigatom begin
            GAccessor.value(signal_slider, get_gtk_property(entry_epoch, :value, Int64))
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

    signal_connect(entry_ts1, "value-changed") do widget
        Gtk.@sigatom begin
            set_gtk_property!(entry_ts1, :value, obj.time_pts[vsearch(get_gtk_property(entry_ts1, :value, Float64), obj.time_pts)])
        end
        draw(can)
    end

    signal_connect(entry_ts2, "value-changed") do widget
        Gtk.@sigatom begin
            set_gtk_property!(entry_ts2, :value, obj.time_pts[vsearch(get_gtk_property(entry_ts2, :value, Float64), obj.time_pts)])
        end
        draw(can)
    end

    signal_connect(bt_delete, "clicked") do widget
        if nepochs(obj) > 1
            ep = get_gtk_property(entry_epoch, :value, Int64)
            if ask_dialog("Delete epoch $ep ?", "No", "Yes")
                delete_epoch!(obj, ep=ep)
                _info("Deleted epoch: $ep")
                Gtk.@sigatom begin
                    set_gtk_property!(entry_epoch, :value, ep)
                    GAccessor.range(entry_epoch, 1, nepochs(obj))
                    GAccessor.range(signal_slider, 1, nepochs(obj))
                end
            end
            draw(can)
        else
            error_dialog("You cannot delete the last epoch.")
        end
    end

    signal_connect(bt_ts, "clicked") do widget
        Gtk.destroy(win)
    end

    signal_connect(bt_close, "clicked") do widget
        quit = true
        Gtk.destroy(win)
    end

    signal_connect(combo_ch, "changed") do widget
        draw(can)
    end

    if mch
        help = "Keyboard shortcuts:\n\nCtrl + b\t\t\tToggle butterfly plot\nCtrl + m\t\tToggle mean plot\n\nPage Up\t\tScroll channels up\nPage Down\tScroll channels down\n\nHome\t\t\tGo to first epoch\nEnd\t\t\tGo to last epoch\nCtrl + ,\t\t\tPrevious epoch\nCtrl + .\t\t\tNext epoch\n\nCtrl + Enter\tReturn selected time segment\nCtrl + d\t\t\tDelete selected time segment\n\nCtrl + s\t\t\tToggle snapping\nAlt + s\t\t\tToggle scales\nAlt + m\t\t\tToggle monochromatic mode\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tClose\n"
    else
        help = "Keyboard shortcuts:\n\nCtrl + b\t\t\tToggle butterfly plot\nCtrl + m\t\tToggle mean plot\n\nPage Up\t\tScroll channels up\nPage Down\tScroll channels down\n\nHome\t\t\tGo to first epoch\nEnd\t\t\tGo to last epoch\nCtrl + ,\t\t\tPrevious epoch\nCtrl + .\t\t\tNext epoch\n\nCtrl + Enter\tReturn selected time segment\nCtrl + d\t\t\tDelete selected time segment\n\nCtrl + s\t\t\tToggle snapping\nAlt + s\t\t\tToggle scales\nAlt + m\t\t\tToggle monochromatic mode\n\nplay current epoch as audio\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tClose\n"
    end

    signal_connect(bt_help, "clicked") do widgete
        info_dialog(help)
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        s = event.state

        if k == 0x0000ff55 # Page Up
            if ch_first > 1
                ch_first -= 1
                ch_last -= 1
                Gtk.@sigatom begin
                    GAccessor.value(ch_slider, ch_first)
                end
                draw(can)
            end
        elseif k == 0x0000ff56 # Page Down
            if mch
                if ch_last < length(ch)
                    ch_first += 1
                    ch_last += 1
                    Gtk.@sigatom begin
                        GAccessor.value(ch_slider, ch_first)
                    end
                    draw(can)
                end
            else
                if ch_first < length(ch)
                    ch_first += 1
                    Gtk.@sigatom begin
                        GAccessor.value(ch_slider, ch_first)
                    end
                    draw(can)
                end
            end
        elseif k == 0x0000ff50 # home
            Gtk.@sigatom begin
                set_gtk_property!(entry_epoch, :value, 1)
            end
        elseif k == 0x0000ff57 # end
            Gtk.@sigatom begin
                set_gtk_property!(entry_epoch, :value, nepochs(obj))
            end
        end

        if s == 0x00000008 || s == 0x00000010 # alt
            if k == 0x00000073 # s
                scale = !scale
                draw(can)
            elseif k == 0x0000006d # m
                mono = !mono
                draw(can)
            end
        end

        if s == 0x00000004 || s == 0x00000014 # ctrl
            if k == 0x00000071 # q
                quit = true
                Gtk.destroy(win)
            elseif k == 104 # h
                info_dialog(help)
            elseif k == 0x00000062 # b
                if mch
                    if plot_type != 1
                        Gtk.@sigatom begin
                            set_gtk_property!(ch_slider, :sensitive, false)
                            set_gtk_property!(combo_ch, :sensitive, true)
                        end
                        plot_type = 1
                    else
                        set_gtk_property!(ch_slider, :sensitive, true)
                        set_gtk_property!(combo_ch, :sensitive, false)
                        plot_type = 0
                    end
                    draw(can)
                end
            elseif k == 0x0000006d # m
                if mch
                    if plot_type != 2
                        Gtk.@sigatom begin
                            set_gtk_property!(ch_slider, :sensitive, false)
                            set_gtk_property!(combo_ch, :sensitive, true)
                        end
                        plot_type = 2
                    else
                        set_gtk_property!(ch_slider, :sensitive, true)
                        set_gtk_property!(combo_ch, :sensitive, false)
                        plot_type = 0
                    end
                    draw(can)
                end
            elseif k == 0x00000070 # p
                if !mch
                    _info("Playing current epoch as audio")
                    ep = get_gtk_property(entry_epoch, :value, Int64)
                    play(obj, ch=cl[ch_first], seg=(0, epoch_len(obj) / sr(obj)), ep=ep)
                end
            elseif k == 0x0000ff0d # Enter
                quit = false
                Gtk.destroy(win)
            elseif k == 0x00000073 # s
                snap = !snap
            elseif k == 0x0000002c # ,
                ep = get_gtk_property(entry_epoch, :value, Int64)
                if ep > 1
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
            elseif k == 100 # d
                ep = get_gtk_property(entry_epoch, :value, Int64)
                if ask_dialog("Delete epoch $ep ?", "No", "Yes")
                    delete_epoch!(obj, ep=ep)
                    _info("Deleted epoch: $ep")
                    ep = ep > 1 ? ep -= 1 : ep = 1
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_epoch, :value, ep)
                        GAccessor.range(entry_epoch, 1, nepochs(obj))
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
    if !quit
        time1 = obj.time_pts[vsearch(get_gtk_property(entry_ts1, :value, Float64), obj.epoch_time)]
        time2 = obj.time_pts[vsearch(get_gtk_property(entry_ts2, :value, Float64), obj.epoch_time)]
        seg = (time1, time2)
        return seg
    else
        return nothing
    end

end

"""
    iview(obj1, obj2; <keyword arguments>)

Interactive view of two continuous signals.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `obj2::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `zoom::Real=10`: how many seconds are displayed in one segment

# Returns

Nothing
"""
function iview(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; zoom::Real=10)::Nothing

    obj1.time_pts[end] < zoom && (zoom = round(obj1.time_pts[end]) / 2)

    @assert size(obj1) == size(obj2) "Both signals must have the same size."
    @assert obj1.header.recording[:channel_order] == obj2.header.recording[:channel_order] "Both signals must have the same order."
    @assert sr(obj1) == sr(obj2) "Both signals must have the same sampling rate."
    @assert zoom > 0 "zoom must be > 0."
    @assert zoom <= signal_len(obj1) / sr(obj1) "zoom must be ≤ $(signal_len(obj1) / sr(obj1))."

    ch_order = obj1.header.recording[:channel_order]
    cl = labels(obj1)[ch_order]
    ch = 1:nchannels(obj1)

    scale = true

    if length(ch) > 20
        ch_first = 1
        ch_last = ch_first + 19
    else
        ch_first = 1
        ch_last = length(ch)
    end

    p = NeuroAnalyzer.plot(obj1, obj2, ch=cl[ch_first:ch_last], title="")
    win = GtkWindow("NeuroAnalyzer: iview()", Int32(p.attr[:size][1]) + 40, Int32(p.attr[:size][2]) + 40)
    set_gtk_property!(win, :border_width, 5)
    set_gtk_property!(win, :resizable, true)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    set_gtk_property!(win, :startup_id, "org.neuroanalyzer")
    can = GtkCanvas(Int32(p.attr[:size][1]), Int32(p.attr[:size][2]))
    g = GtkGrid()
    set_gtk_property!(g, :column_homogeneous, false)
    set_gtk_property!(g, :column_spacing, 5)
    set_gtk_property!(g, :row_spacing, 5)
    entry_time = GtkSpinButton(obj1.time_pts[1], obj1.time_pts[end] - zoom, 1)
    set_gtk_property!(entry_time, :digits, 2)
    set_gtk_property!(entry_time, :value, obj1.time_pts[1])
    set_gtk_property!(entry_time, :tooltip_text, "Time position [s]")
    bt_start = GtkButton("⇤")
    set_gtk_property!(bt_start, :tooltip_text, "Go to the signal start")
    bt_prev = GtkButton("←")
    set_gtk_property!(bt_prev, :tooltip_text, "Go back by $(round(zoom)) seconds")
    bt_next = GtkButton("→")
    set_gtk_property!(bt_next, :tooltip_text, "Go forward by $(round(zoom)) seconds")
    bt_end = GtkButton("⇥")
    set_gtk_property!(bt_end, :tooltip_text, "Go to the signal end")
    bt_help = GtkButton("Help")
    set_gtk_property!(bt_help, :tooltip_text, "Show help")
    bt_close = GtkButton("Close")
    set_gtk_property!(bt_close, :tooltip_text, "Close this window")

    if length(ch) > 20
        ch_slider = GtkScale(false, ch[ch_first]:(ch[end] - 19))
        set_gtk_property!(ch_slider, :draw_value, false)
    else
        ch_slider = GtkScale(false, ch[ch_first]:ch[end])
        set_gtk_property!(ch_slider, :draw_value, false)
        set_gtk_property!(ch_slider, :sensitive, false)
    end
    set_gtk_property!(ch_slider, :tooltip_text, "Scroll channels")
    set_gtk_property!(ch_slider, :vexpand, true)
    oc = GtkOrientable(ch_slider)
    set_gtk_property!(oc, :orientation, 1)

    signal_slider = GtkScale(false, 1:obj1.time_pts[end] - zoom)
    set_gtk_property!(signal_slider, :draw_value, false)
    set_gtk_property!(signal_slider, :tooltip_text, "Time position")

    g[1:7, 1] = can
    g[8, 1] = ch_slider
    g[1:7, 2] = signal_slider
    g[1, 3] = bt_start
    g[2, 3] = bt_prev
    g[3, 3] = entry_time
    g[4, 3] = bt_next
    g[5, 3] = bt_end
    g[6, 3] = bt_help
    g[7, 3] = bt_close
    push!(win, g)

    showall(win)

    @guarded draw(can) do widget
        time1 = get_gtk_property(entry_time, :value, Float64)
        time2 = time1 + zoom
        time2 > obj1.time_pts[end] && (time2 = obj1.time_pts[end])
        ctx = getgc(can)
        p = NeuroAnalyzer.plot(obj1, obj2,
                               ch=cl[ch_first:ch_last],
                               seg=(time1, time2),
                               title="",
                               scale=scale)
        show(io, MIME("image/png"), p)
        img = read_from_png(io)
        set_source_surface(ctx, img, 0, 0)
        paint(ctx)
    end

    can.mouse.scroll = @guarded (widget, event) -> begin
        s = event.state
        if event.direction == 1 # down
            if s == 0x00000001
                time_current = get_gtk_property(entry_time, :value, Float64)
                if time_current < obj1.time_pts[end] - zoom
                    time_current += 1
                else
                    time_current = obj1.time_pts[end] - zoom
                end
                Gtk.@sigatom begin
                    set_gtk_property!(entry_time, :value, time_current)
                end
            elseif s == 0x00000004
                time_current = get_gtk_property(entry_time, :value, Float64)
                if time_current < obj1.time_pts[end] - zoom
                    time_current += zoom
                else
                    time_current = obj1.time_pts[end] - zoom
                end
                Gtk.@sigatom begin
                    set_gtk_property!(entry_time, :value, time_current)
                end
            else
                if ch_last < length(ch)
                    ch_first += 1
                    ch_last += 1
                    Gtk.@sigatom begin
                        GAccessor.value(ch_slider, ch_first)
                    end
                    draw(can)
                end
            end
        elseif event.direction == 0 # up
            if s == 0x00000001
                time_current = get_gtk_property(entry_time, :value, Float64)
                if time_current >= obj1.time_pts[1] + 1
                    time_current -= 1
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_time, :value, time_current)
                    end
                end
            elseif s == 0x00000004
                time_current = get_gtk_property(entry_time, :value, Float64)
                if time_current >= obj1.time_pts[1] + zoom
                    time_current = time_current - zoom
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_time, :value, time_current)
                    end
                end
            else
                if ch_first > 1
                    ch_first -= 1
                    ch_last -= 1
                    Gtk.@sigatom begin
                        GAccessor.value(ch_slider, ch_first)
                    end
                    draw(can)
                end
            end
        end
    end

    signal_connect(signal_slider, "value-changed") do widget, others...
        Gtk.@sigatom begin
            set_gtk_property!(entry_time, :value, round(GAccessor.value(signal_slider)))
        end
        draw(can)
    end

    signal_connect(ch_slider, "value-changed") do widget, others...
        ch_first = round(Int64, GAccessor.value(ch_slider))
        ch_last = ch_first + 19
        draw(can)
    end

    signal_connect(entry_time, "value-changed") do widget
        Gtk.@sigatom begin
            GAccessor.value(signal_slider, get_gtk_property(entry_time, :value, Float64))
        end
        draw(can)
    end

    signal_connect(bt_prev, "clicked") do widget
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
            time_current += zoom
        else
            time_current = obj1.time_pts[end] - zoom
        end
        Gtk.@sigatom begin
            set_gtk_property!(entry_time, :value, time_current)
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

    help = "Keyboard shortcuts:\n\nPage Up\t\tScroll channels up\nPage Down\tScroll channels down\n\nHome\t\t\tGo to the signal start\nEnd\t\t\tGo to the signal end\nCtrl + ,\t\t\tGo back by 1 second\nCtrl + .\t\t\tGo forward by 1 second\nAlt + ,\t\t\tGo back by $(round(zoom)) seconds\nAlt + .\t\t\tGo forward by $(round(zoom)) seconds\n\n[\t\t\t\tZoom in\n]\t\t\t\tZoom out\n\nAlt + s\t\t\tToggle scales\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tClose\n"

    signal_connect(bt_help, "clicked") do widgete
        info_dialog(help)
    end

    signal_connect(bt_close, "clicked") do widget
        Gtk.destroy(win)
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        s = event.state

        if k == 0x0000ff55 # Page Up
            if ch_first > 1
                ch_first -= 1
                ch_last -= 1
                Gtk.@sigatom begin
                    GAccessor.value(ch_slider, ch_first)
                end
                draw(can)
            end
        elseif k == 0x0000ff56 # Page Down
            if ch_last < length(ch)
                ch_first += 1
                ch_last += 1
                Gtk.@sigatom begin
                    GAccessor.value(ch_slider, ch_first)
                end
                draw(can)
            end
        elseif k == 0x0000005b # [
            if zoom > 1
                zoom -= 1
                set_gtk_property!(bt_next, :tooltip_text, "Go forward by $(round(zoom)) seconds")
                set_gtk_property!(bt_prev, :tooltip_text, "Go back by $(round(zoom)) seconds")
                GAccessor.range(entry_time, obj1.time_pts[1], obj1.time_pts[end] - zoom)
                draw(can)
            end
            help = "Keyboard shortcuts:\n\nPage Up\t\tScroll channels up\nPage Down\tScroll channels down\n\nHome\t\t\tGo to the signal start\nEnd\t\t\tGo to the signal end\nCtrl + ,\t\t\tGo back by 1 second\nCtrl + .\t\t\tGo forward by 1 second\nAlt + ,\t\t\tGo back by $(round(zoom)) seconds\nAlt + .\t\t\tGo forward by $(round(zoom)) seconds\n\n[\t\t\t\tZoom in\n]\t\t\t\tZoom out\n\nAlt + s\t\t\tToggle scales\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tClose\n"
        elseif k == 0x0000005d # ]
            if zoom < 30 && zoom < obj1.time_pts[end] - 1
                zoom += 1
                set_gtk_property!(bt_next, :tooltip_text, "Go forward by $(round(zoom)) seconds")
                set_gtk_property!(bt_prev, :tooltip_text, "Go back by $(round(zoom)) seconds")
                GAccessor.range(entry_time, obj1.time_pts[1], obj1.time_pts[end] - zoom)
                draw(can)
            else
                zoom = obj1.time_pts[end]
                set_gtk_property!(bt_next, :tooltip_text, "Go forward by $(round(zoom)) seconds")
                set_gtk_property!(bt_prev, :tooltip_text, "Go back by $(round(zoom)) seconds")
                GAccessor.range(entry_time, obj1.time_pts[1], obj1.time_pts[end] - zoom)
                draw(can)
            end
            help = "Keyboard shortcuts:\n\nPage Up\t\tScroll channels up\nPage Down\tScroll channels down\n\nHome\t\t\tGo to the signal start\nEnd\t\t\tGo to the signal end\nCtrl + ,\t\t\tGo back by 1 second\nCtrl + .\t\t\tGo forward by 1 second\nAlt + ,\t\t\tGo back by $(round(zoom)) seconds\nAlt + .\t\t\tGo forward by $(round(zoom)) seconds\n\n[\t\t\t\tZoom in\n]\t\t\t\tZoom out\n\nAlt + s\t\t\tToggle scales\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tClose\n"
        elseif k == 0x0000ff50 # home
            Gtk.@sigatom begin
                set_gtk_property!(entry_time, :value, obj1.time_pts[1])
            end
            draw(can)
        elseif k == 0x0000ff57 # end
            time_current = obj1.time_pts[end] - zoom
            Gtk.@sigatom begin
                set_gtk_property!(entry_time, :value, time_current)
            end
            draw(can)
        end

        if s == 0x00000008 || s == 0x00000010 # alt
            if k == 0x0000002c # ,
                time_current = get_gtk_property(entry_time, :value, Float64)
                if time_current >= obj1.time_pts[1] + zoom
                    time_current = time_current - zoom
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_time, :value, time_current)
                    end
                end
                draw(can)
            elseif k == 0x0000002e # .
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
            elseif k == 0x00000073 # s
                scale = !scale
                draw(can)
            end
        end

        if s == 0x00000004 || s == 0x00000014 # ctrl
            if k == 0x00000071 # q
                Gtk.destroy(win)
            elseif k == 104 # h
                info_dialog(help)
            elseif k == 0x0000002c # ,
                time_current = get_gtk_property(entry_time, :value, Float64)
                if time_current >= obj1.time_pts[1] + 1
                    time_current -= 1
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_time, :value, time_current)
                    end
                end
                draw(can)
            elseif k == 0x0000002e # .
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
    iview_ep(obj1, obj2; <keyword arguments>)

Interactive view of two epoched signals.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `obj2::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ep::Int64=1`: initial epoch to display

# Returns

Nothing
"""
function iview_ep(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ep::Int64=1)::Nothing

    @assert obj1.header.recording[:channel_order] == obj2.header.recording[:channel_order] "Both signals must have the same order."
    @assert size(obj1) == size(obj2) "Both signals must have the same size."
    @assert sr(obj1) == sr(obj2) "Both signals must have the same sampling rate."
    @assert nepochs(obj1) > 1 "iview() must be used for continuous object."

    ch_order = obj1.header.recording[:channel_order]
    cl = labels(obj1)[ch_order]
    ch = 1:nchannels(obj1)
    _check_epochs(obj1, ep)

    scale = true

    if length(ch) > 20
        ch_first = 1
        ch_last = ch_first + 19
    else
        ch_first = 1
        ch_last = length(ch)
    end

    p = NeuroAnalyzer.plot(obj1, obj2, ch=cl[ch_first:ch_last], ep=ep, title="")
    win = GtkWindow("NeuroAnalyzer: iview_ep()", Int32(p.attr[:size][1]) + 40, Int32(p.attr[:size][2]) + 40)
    set_gtk_property!(win, :border_width, 5)
    set_gtk_property!(win, :resizable, true)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    set_gtk_property!(win, :startup_id, "org.neuroanalyzer")
    can = GtkCanvas(Int32(p.attr[:size][1]), Int32(p.attr[:size][2]))
    g = GtkGrid()
    set_gtk_property!(g, :column_homogeneous, false)
    set_gtk_property!(g, :column_spacing, 5)
    set_gtk_property!(g, :row_spacing, 5)
    entry_epoch = GtkSpinButton(1, nepochs(obj1), 1)
    set_gtk_property!(entry_epoch, :value, ep)
    set_gtk_property!(entry_epoch, :tooltip_text, "Epoch")
    bt_start = GtkButton("⇤")
    set_gtk_property!(bt_start, :tooltip_text, "Go to the signal start")
    bt_end = GtkButton("⇥")
    set_gtk_property!(bt_end, :tooltip_text, "Go to the signal end")
    bt_help = GtkButton("Help")
    set_gtk_property!(bt_help, :tooltip_text, "Show help")
    bt_close = GtkButton("Close")
    set_gtk_property!(bt_close, :tooltip_text, "Close this window")

    if length(ch) > 20
        ch_slider = GtkScale(false, ch[ch_first]:(ch[end] - 19))
        set_gtk_property!(ch_slider, :draw_value, false)
    else
        ch_slider = GtkScale(false, ch[ch_first]:ch[end])
        set_gtk_property!(ch_slider, :draw_value, false)
        set_gtk_property!(ch_slider, :sensitive, false)
    end
    set_gtk_property!(ch_slider, :tooltip_text, "Scroll channels")
    set_gtk_property!(ch_slider, :vexpand, true)
    oc = GtkOrientable(ch_slider)
    set_gtk_property!(oc, :orientation, 1)

    signal_slider = GtkScale(false, 1:nepochs(obj1))
    set_gtk_property!(signal_slider, :draw_value, false)
    set_gtk_property!(signal_slider, :tooltip_text, "Current epoch")

    g[1:5, 1] = can
    g[6, 1] = ch_slider
    g[1:5, 2] = signal_slider
    g[1, 3] = bt_start
    g[2, 3] = entry_epoch
    g[3, 3] = bt_end
    g[4, 3] = bt_help
    g[5, 3] = bt_close
    push!(win, g)

    showall(win)

    @guarded draw(can) do widget
        ep = get_gtk_property(entry_epoch, :value, Int64)
        ctx = getgc(can)
        p = NeuroAnalyzer.plot(obj1, obj2,
                               ch=cl[ch_first:ch_last],
                               ep=ep,
                               title="",
                               scale=scale)
        show(io, MIME("image/png"), p)
        img = read_from_png(io)
        set_source_surface(ctx, img, 0, 0)
        paint(ctx)
    end

    can.mouse.scroll = @guarded (widget, event) -> begin
        s = event.state
        if event.direction == 1 # down
            if s == 1
                ep = get_gtk_property(entry_epoch, :value, Int64)
                if ep < nepochs(obj1)
                    ep += 1
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_epoch, :value, ep)
                    end
                end
            else
                if ch_last < length(ch)
                    ch_first += 1
                    ch_last += 1
                    Gtk.@sigatom begin
                        GAccessor.value(ch_slider, ch_first)
                    end
                    draw(can)
                end
            end
        elseif event.direction == 0 # up
            if s == 1
                ep = get_gtk_property(entry_epoch, :value, Int64)
                if ep > 1
                    ep -= 1
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_epoch, :value, ep)
                    end
                end
            else
                if ch_first > 1
                    ch_first -= 1
                    ch_last -= 1
                    Gtk.@sigatom begin
                        GAccessor.value(ch_slider, ch_first)
                    end
                    draw(can)
                end
            end
        end
    end

    signal_connect(signal_slider, "value-changed") do widget, others...
        Gtk.@sigatom begin
            set_gtk_property!(entry_epoch, :value, round(Int64, GAccessor.value(signal_slider)))
        end
        draw(can)
    end

    signal_connect(ch_slider, "value-changed") do widget, others...
        ch_first = round(Int64, GAccessor.value(ch_slider))
        ch_last = ch_first + 19
        draw(can)
    end

    signal_connect(entry_epoch, "value-changed") do widget
        Gtk.@sigatom begin
            GAccessor.value(signal_slider, get_gtk_property(entry_epoch, :value, Int64))
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
            set_gtk_property!(entry_epoch, :value, nepochs(obj1))
        end
    end

    signal_connect(bt_close, "clicked") do widget
        Gtk.destroy(win)
    end

    help = "Keyboard shortcuts:\n\nPage Up\t\tScroll channels up\nPage Down\tScroll channels down\n\nHome\t\t\tGo to first epoch\nEnd\t\t\tGo to last epoch\nCtrl + ,\t\t\tPrevious epoch\nCtrl + .\t\t\tNext epoch\n\nAlt + s\t\t\tToggle scales\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tClose\n"

    signal_connect(bt_help, "clicked") do widgete
        info_dialog(help)
    end

    signal_connect(bt_close, "clicked") do widget
        Gtk.destroy(win)
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        s = event.state
        if k == 0x0000ff55 # Page Up
            if ch_first > 1
                ch_first -= 1
                ch_last -= 1
                Gtk.@sigatom begin
                    GAccessor.value(ch_slider, ch_first)
                end
                draw(can)
            end
        elseif k == 0x0000ff56 # Page Down
            if ch_last < length(ch)
                ch_first += 1
                ch_last += 1
                Gtk.@sigatom begin
                    GAccessor.value(ch_slider, ch_first)
                end
                draw(can)
            end
        elseif k == 0x0000ff50 # home
            Gtk.@sigatom begin
                set_gtk_property!(entry_epoch, :value, 1)
            end
        elseif k == 0x0000ff57 # end
            Gtk.@sigatom begin
                set_gtk_property!(entry_epoch, :value, nepochs(obj1))
            end
        end
        if s == 0x00000008 || s == 0x00000010 # alt
            if k == 0x00000073 # s
                scale = !scale
                draw(can)
            end
        end
        if s == 0x00000004 || s == 0x00000014 # ctrl
            if k == 0x00000071 # q
                Gtk.destroy(win)
            elseif k == 104 # h
                info_dialog(help)
            elseif k == 97 # a
                Gtk.@sigatom begin
                    set_gtk_property!(entry_epoch, :value, 1)
                end
            elseif k == 0x0000002c # ,
                ep = get_gtk_property(entry_epoch, :value, Int64)
                if ep > 1
                    ep -= 1
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_epoch, :value, ep)
                    end
                end
            elseif k == 0x0000002e # .
                ep = get_gtk_property(entry_epoch, :value, Int64)
                if ep < nepochs(obj1)
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

"""
    iview(obj, c; <keyword arguments>)

Interactive view of embedded or external component of continuous signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `c::Union{Symbol, AbstractArray}`: component to plot
- `zoom::Real=10`: how many seconds are displayed in one segment

# Returns

Nothing
"""
function iview(obj::NeuroAnalyzer.NEURO, c::Union{Symbol, AbstractArray}; zoom::Real=10)::Nothing

    @assert zoom > 0 "zoom must be > 0."
    @assert zoom <= signal_len(obj) / sr(obj) "zoom must be ≤ $(signal_len(obj) / sr(obj))."

    c isa Symbol && (c = _get_component(obj, c))
    ch = axes(c, 1)

    mono = false
    quit = false
    scale = true

    if length(ch) > 20
        ch_first = 1
        ch_last = 20
    else
        ch_first = 1
        ch_last = length(ch)
    end

    length(ch) == 1 && (ch = [ch])

    p = NeuroAnalyzer.plot(obj, c, c_idx=ch[ch_first:ch_last], mono=mono, scale=scale, title="")
    win = GtkWindow("NeuroAnalyzer: iview()", Int32(p.attr[:size][1]) + 40, Int32(p.attr[:size][2]) + 40)
    set_gtk_property!(win, :border_width, 5)
    set_gtk_property!(win, :resizable, true)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    set_gtk_property!(win, :startup_id, "org.neuroanalyzer")
    can = GtkCanvas(Int32(p.attr[:size][1]), Int32(p.attr[:size][2]))
    g = GtkGrid()
    set_gtk_property!(g, :column_homogeneous, false)
    set_gtk_property!(g, :column_spacing, 5)
    set_gtk_property!(g, :row_spacing, 5)
    entry_time = GtkSpinButton(obj.time_pts[1], obj.time_pts[end] - zoom, 1)
    set_gtk_property!(entry_time, :digits, 2)
    set_gtk_property!(entry_time, :value, obj.time_pts[1])
    set_gtk_property!(entry_time, :tooltip_text, "Time position [s]")
    bt_start = GtkButton("⇤")
    set_gtk_property!(bt_start, :tooltip_text, "Go to the signal start")
    bt_prev = GtkButton("←")
    set_gtk_property!(bt_prev, :tooltip_text, "Go back by $(round(zoom)) seconds")
    bt_next = GtkButton("→")
    set_gtk_property!(bt_next, :tooltip_text, "Go forward by $(round(zoom)) seconds")
    bt_end = GtkButton("⇥")
    set_gtk_property!(bt_end, :tooltip_text, "Go to the signal end")
    bt_help = GtkButton("Help")
    set_gtk_property!(bt_help, :tooltip_text, "Show help")
    bt_close = GtkButton("Close")
    set_gtk_property!(bt_close, :tooltip_text, "Close this window")
    if length(ch) > 20
        ch_slider = GtkScale(false, ch[ch_first]:(ch[end] - 19))
        set_gtk_property!(ch_slider, :draw_value, false)
    else
        ch_slider = GtkScale(false, ch[ch_first]:ch[end])
        set_gtk_property!(ch_slider, :draw_value, false)
        set_gtk_property!(ch_slider, :sensitive, false)
    end
    set_gtk_property!(ch_slider, :tooltip_text, "Scroll components")
    set_gtk_property!(ch_slider, :vexpand, true)
    oc = GtkOrientable(ch_slider)
    set_gtk_property!(oc, :orientation, 1)

    signal_slider = GtkScale(false, obj.time_pts[1]:obj.time_pts[end] - zoom)
    set_gtk_property!(signal_slider, :draw_value, false)
    set_gtk_property!(signal_slider, :tooltip_text, "Time position")

    g[1:7, 1] = can
    g[8, 1] = ch_slider
    g[1:7, 2] = signal_slider
    g[1, 3] = bt_start
    g[2, 3] = bt_prev
    g[3, 3] = entry_time
    g[4, 3] = bt_next
    g[5, 3] = bt_end
    g[6, 3] = bt_help
    g[7, 3] = bt_close
    push!(win, g)

    showall(win)

    @guarded draw(can) do widget
        time1 = get_gtk_property(entry_time, :value, Float64)
        time2 = time1 + zoom
        time2 > obj.time_pts[end] && (time2 = obj.time_pts[end])
        ctx = getgc(can)
        p = NeuroAnalyzer.plot(obj,
                               c,
                               c_idx=ch[ch_first]:ch[ch_last],
                               seg=(time1, time2),
                               mono=mono,
                               scale=scale,
                               title="")
        show(io, MIME("image/png"), p)
        img = read_from_png(io)
        set_source_surface(ctx, img, 0, 0)
        paint(ctx)
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
            else
                if ch_last < ch[end]
                    ch_first += 1
                    ch_last += 1
                    Gtk.@sigatom begin
                        GAccessor.value(ch_slider, ch_first)
                    end
                    draw(can)
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
            else
                if ch_first > 1
                    ch_first -= 1
                    ch_last -= 1
                    Gtk.@sigatom begin
                        GAccessor.value(ch_slider, ch_first)
                    end
                    draw(can)
                end
            end
        end
    end

    signal_connect(signal_slider, "value-changed") do widget, others...
        Gtk.@sigatom begin
            set_gtk_property!(entry_time, :value, round(GAccessor.value(signal_slider)))
        end
        draw(can)
    end

    signal_connect(ch_slider, "value-changed") do widget
        ch_first = round(Int64, GAccessor.value(ch_slider))
        ch_last = ch_first + 19
        draw(can)
    end

    signal_connect(entry_time, "value-changed") do widget
        Gtk.@sigatom begin
            GAccessor.value(signal_slider, get_gtk_property(entry_time, :value, Float64))
        end
        draw(can)
    end

    signal_connect(bt_prev, "clicked") do widget
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
            time_current += zoom
        else
            time_current = obj.time_pts[end] - zoom
        end
        Gtk.@sigatom begin
            set_gtk_property!(entry_time, :value, time_current)
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

    help = "Keyboard shortcuts:\n\nPage Up\t\tScroll components up\nPage Down\tScroll components down\n\nHome\t\t\tGo to the signal start\nEnd\t\t\tGo to the signal end\nCtrl + ,\t\t\tGo back by 1 second\nCtrl + .\t\t\tGo forward by 1 second\nAlt + ,\t\t\tGo back by $(round(zoom)) seconds\nAlt + .\t\t\tGo forward by $(round(zoom)) seconds\n\n[\t\t\t\tZoom in\n]\t\t\t\tZoom out\n\nAlt + s\t\t\tToggle scales\nAlt + m\t\t\tToggle monochromatic mode\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tClose\n"

    signal_connect(bt_help, "clicked") do widgete
        info_dialog(help)
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        s = event.state

        if k == 0x0000ff55 # Page Up
            if ch_first > 1
                ch_first -= 1
                ch_last -= 1
                Gtk.@sigatom begin
                    GAccessor.value(ch_slider, ch_first)
                end
                draw(can)
            end
        elseif k == 0x0000ff56 # Page Down
            if ch_last < ch[end]
                ch_first += 1
                ch_last += 1
                Gtk.@sigatom begin
                    GAccessor.value(ch_slider, ch_first)
                end
                draw(can)
            end
        elseif k == 0x0000005b # [
            if zoom > 1
                zoom -= 1
                set_gtk_property!(bt_next, :tooltip_text, "Go forward by $(round(zoom)) seconds")
                set_gtk_property!(bt_prev, :tooltip_text, "Go back by $(round(zoom)) seconds")
                GAccessor.range(signal_slider, obj.time_pts[1], obj.time_pts[end] - zoom)
                GAccessor.range(entry_time, obj.time_pts[1], obj.time_pts[end] - zoom)
                draw(can)
            end
            help = "Keyboard shortcuts:\n\nPage Up\t\tScroll channels up\nPage Down\tScroll channels down\n\nHome\t\t\tGo to the signal start\nEnd\t\t\tGo to the signal end\nCtrl + ,\t\t\tGo back by 1 second\nCtrl + .\t\t\tGo forward by 1 second\nAlt + ,\t\t\tGo back by $(round(zoom)) seconds\nAlt + .\t\t\tGo forward by $(round(zoom)) seconds\n\n[\t\t\t\tZoom in\n]\t\t\t\tZoom out\n\nAlt + s\t\t\tToggle scales\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tClose\n"
        elseif k == 0x0000005d # ]
            if zoom < 30 && zoom < obj.time_pts[end] - 1
                zoom += 1
                set_gtk_property!(bt_next, :tooltip_text, "Go forward by $(round(zoom)) seconds")
                set_gtk_property!(bt_prev, :tooltip_text, "Go back by $(round(zoom)) seconds")
                GAccessor.range(signal_slider, obj.time_pts[1], obj.time_pts[end] - zoom)
                GAccessor.range(entry_time, obj.time_pts[1], obj.time_pts[end] - zoom)
                draw(can)
            else
                zoom = obj.time_pts[end]
                set_gtk_property!(bt_next, :tooltip_text, "Go forward by $(round(zoom)) seconds")
                set_gtk_property!(bt_prev, :tooltip_text, "Go back by $(round(zoom)) seconds")
                GAccessor.range(signal_slider, obj.time_pts[1], obj.time_pts[end] - zoom)
                GAccessor.range(entry_time, obj.time_pts[1], obj.time_pts[end] - zoom)
                draw(can)
            end
            help = "Keyboard shortcuts:\n\nPage Up\t\tScroll channels up\nPage Down\tScroll channels down\n\nHome\t\t\tGo to the signal start\nEnd\t\t\tGo to the signal end\nCtrl + ,\t\t\tGo back by 1 second\nCtrl + .\t\t\tGo forward by 1 second\nAlt + ,\t\t\tGo back by $(round(zoom)) seconds\nAlt + .\t\t\tGo forward by $(round(zoom)) seconds\n\n[\t\t\t\tZoom in\n]\t\t\t\tZoom out\n\nAlt + s\t\t\tToggle scales\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tClose\n"
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

        if s == 0x00000008 || s == 0x00000010 # alt
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
            elseif k == 0x0000006d # m
                mono = !mono
                draw(can)
            elseif k == 0x00000073 # s
                scale = !scale
                draw(can)
            end
        end

        if s == 0x00000004 || s == 0x00000014 # ctrl
            if k == 0x00000071 # q
                Gtk.destroy(win)
            elseif k == 104 # h
                info_dialog(help)
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
    iview_ep(obj, c; <keyword arguments>)

Interactive view of embedded or external component of epoched signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `c::Union{Symbol, AbstractArray}`: component to plot
- `ep::Int64=1`: initial epoch to display

# Returns

Nothing
"""
function iview_ep(obj::NeuroAnalyzer.NEURO, c::Union{Symbol, AbstractArray}; ep::Int64=1)::Nothing

    @assert nepochs(obj) > 1 "iview() must be used for continuous object."

    c isa Symbol && (c = _get_component(obj, c))
    ch = axes(c, 1)

    _check_epochs(obj, ep)

    mono = false
    scale = true

    if length(ch) > 20
        ch_first = 1
        ch_last = 20
    else
        ch_first = 1
        ch_last = length(ch)
    end

    p = NeuroAnalyzer.plot(obj, c, c_idx=ch[ch_first:ch_last], ep=ep, mono=mono, scale=scale, title="")
    win = GtkWindow("NeuroAnalyzer: iview_ep()", Int32(p.attr[:size][1]) + 40, Int32(p.attr[:size][2]) + 40)
    set_gtk_property!(win, :border_width, 5)
    set_gtk_property!(win, :resizable, true)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    set_gtk_property!(win, :startup_id, "org.neuroanalyzer")
    can = GtkCanvas(Int32(p.attr[:size][1]), Int32(p.attr[:size][2]))
    g = GtkGrid()
    set_gtk_property!(g, :column_homogeneous, false)
    set_gtk_property!(g, :column_spacing, 5)
    set_gtk_property!(g, :row_spacing, 5)
    entry_epoch = GtkSpinButton(1, nepochs(obj), 1)
    set_gtk_property!(entry_epoch, :value, ep)
    set_gtk_property!(entry_epoch, :tooltip_text, "Epoch")
    bt_start = GtkButton("⇤")
    set_gtk_property!(bt_start, :tooltip_text, "Go to the signal start")
    bt_end = GtkButton("⇥")
    set_gtk_property!(bt_end, :tooltip_text, "Go to the signal end")
    bt_help = GtkButton("Help")
    set_gtk_property!(bt_help, :tooltip_text, "Show help")
    bt_close = GtkButton("Close")
    set_gtk_property!(bt_close, :tooltip_text, "Close this window")

    if length(ch) > 20
        ch_slider = GtkScale(false, ch[ch_first]:(ch[end] - 19))
        set_gtk_property!(ch_slider, :draw_value, false)
    else
        ch_slider = GtkScale(false, ch[ch_first]:ch[end])
        set_gtk_property!(ch_slider, :draw_value, false)
        set_gtk_property!(ch_slider, :sensitive, false)
    end
    set_gtk_property!(ch_slider, :tooltip_text, "Scroll components")
    set_gtk_property!(ch_slider, :vexpand, true)
    oc = GtkOrientable(ch_slider)
    set_gtk_property!(oc, :orientation, 1)

    signal_slider = GtkScale(false, 1:nepochs(obj))
    set_gtk_property!(signal_slider, :draw_value, false)
    set_gtk_property!(signal_slider, :tooltip_text, "Current epoch")

    g[1:5, 1] = can
    g[6, 1] = ch_slider
    g[1:5, 2] = signal_slider
    g[1, 3] = bt_start
    g[2, 3] = entry_epoch
    g[3, 3] = bt_end
    g[4, 3] = bt_help
    g[5, 3] = bt_close
    push!(win, g)

    showall(win)

    @guarded draw(can) do widget
        ep = get_gtk_property(entry_epoch, :value, Int64)
        ctx = getgc(can)
        p = NeuroAnalyzer.plot(obj,
                               c,
                               c_idx=ch[ch_first]:ch[ch_last],
                               ep=ep,
                               mono=mono,
                               scale=scale,
                               title="")
        show(io, MIME("image/png"), p)
        img = read_from_png(io)
        set_source_surface(ctx, img, 0, 0)
        paint(ctx)
    end

    can.mouse.scroll = @guarded (widget, event) -> begin
        s = event.state
        if event.direction == 1 # down
            if s == 0x00000001
                ep = get_gtk_property(entry_epoch, :value, Int64)
                if ep < nepochs(obj)
                    ep += 1
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_epoch, :value, ep)
                    end
                end
            else
                if ch_last < ch[end]
                    ch_first += 1
                    ch_last += 1
                    Gtk.@sigatom begin
                        GAccessor.value(ch_slider, ch_first)
                    end
                    draw(can)
                end
            end
        elseif event.direction == 0 # up
            if s == 0x00000001
                ep = get_gtk_property(entry_epoch, :value, Int64)
                if ep > 1
                    ep -= 1
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_epoch, :value, ep)
                    end
                end
            else
                if ch_first > 1
                    ch_first -= 1
                    ch_last -= 1
                    Gtk.@sigatom begin
                        GAccessor.value(ch_slider, ch_first)
                    end
                    draw(can)
                end
            end
        end
    end

    signal_connect(signal_slider, "value-changed") do widget, others...
        Gtk.@sigatom begin
            set_gtk_property!(entry_epoch, :value, round(Int64, GAccessor.value(signal_slider)))
        end
        draw(can)
    end

    signal_connect(ch_slider, "value-changed") do widget, others...
        ch_first = round(Int64, GAccessor.value(ch_slider))
        ch_last = ch_first + 19
        draw(can)
    end

    signal_connect(entry_epoch, "value-changed") do widget
        Gtk.@sigatom begin
            GAccessor.value(signal_slider, get_gtk_property(entry_epoch, :value, Int64))
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

    help = "Keyboard shortcuts:\n\nPage Up\t\tScroll components up\nPage Down\tScroll components down\n\nHome\t\t\tGo to first epoch\nEnd\t\t\tGo to last epoch\nCtrl + ,\t\t\tPrevious epoch\nCtrl + .\t\t\tNext epoch\n\nAlt + s\t\t\tToggle scales\nAlt + m\t\t\tToggle monochromatic mode\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tClose\n"

    signal_connect(bt_help, "clicked") do widgete
        info_dialog(help)
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        s = event.state

        if k == 0x0000ff55 # Page Up
            if ch_first > 1
                ch_first -= 1
                ch_last -= 1
                Gtk.@sigatom begin
                    GAccessor.value(ch_slider, ch_first)
                end
                draw(can)
            end
        elseif k == 0x0000ff56 # Page Down
            if ch_last < ch[end]
                ch_first += 1
                ch_last += 1
                Gtk.@sigatom begin
                    GAccessor.value(ch_slider, ch_first)
                end
                draw(can)
            end
        elseif k == 0x0000ff50 # home
            Gtk.@sigatom begin
                set_gtk_property!(entry_epoch, :value, 1)
            end
        elseif k == 0x0000ff57 # end
            Gtk.@sigatom begin
                set_gtk_property!(entry_epoch, :value, nepochs(obj))
            end
        end
        if s == 0x00000008 || s == 0x00000010 # alt
            if k == 0x00000073 # s
                scale = !scale
                draw(can)
            elseif k == 0x0000006d # m
                mono = !mono
                draw(can)
            end
        end

        if s == 0x00000004 || s == 0x00000014 # ctrl
            if k == 0x00000071 # q
                Gtk.destroy(win)
            elseif k == 104 # h
                info_dialog(help)
            elseif k == 0x0000002c # ,
                ep = get_gtk_property(entry_epoch, :value, Int64)
                if ep > 1
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


"""
    iview(p)

View plot object.

# Arguments

- `p::Plots.Plot{Plots.GRBackend}`

# Returns

Nothing
"""
function iview(p::Plots.Plot{Plots.GRBackend})::Nothing

    win = GtkWindow("NeuroAnalyzer: iview()", p.attr[:size][1] + 2, p.attr[:size][2] + 2)
    set_gtk_property!(win, :border_width, 0)
    set_gtk_property!(win, :resizable, false)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    set_gtk_property!(win, :startup_id, "org.neuroanalyzer")
    can = GtkCanvas(p.attr[:size][1] + 2, p.attr[:size][2] + 2)
    push!(win, can)
    showall(win)

    @guarded draw(can) do widget
        show(io, MIME("image/png"), p)
        img = read_from_png(io)
        w = img.width
        h = img.height
        ctx = getgc(can)
        # Cairo.scale(ctx, 0.75, 0.75)
        Cairo.set_source_surface(ctx, img, 1, 2)
        Cairo.paint(ctx)
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        s = event.state
        if s == 0x00000004 || s == 0x00000014 # ctrl
            if k == 115 # s
                file_name = save_dialog("Pick image file", GtkNullContainer(), (GtkFileFilter("*.png", name="All supported formats"), "*.png"))
                    if file_name != ""
                        if splitext(file_name)[2] in [".png"]
                            try
                                surface_buf = Gtk.cairo_surface(can)
                                Cairo.write_to_png(surface_buf, file_name)
                                _info("Plot saved as: $file_name")
                            catch
                                warn_dialog("File $file_name cannot be written!")
                            end
                        else
                            warn_dialog("Incorrect filename!")
                        end
                    end
            elseif k == 0x00000071 # q
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

"""
    iview(file_name)

View PNG image.

# Arguments

- `file_name::String`

# Returns

Nothing
"""
function iview(file_name::String)::Nothing

    @assert isfile(file_name) "File $file_name cannot be opened."
    if splitext(file_name)[2] != ".png"
        @error "Incorrect filename!"
        return nothing
    end

    img = read_from_png(file_name)

    win = GtkWindow("NeuroAnalyzer: iview()", img.width, img.height)
    set_gtk_property!(win, :border_width, 0)
    set_gtk_property!(win, :resizable, false)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    set_gtk_property!(win, :startup_id, "org.neuroanalyzer")
    can = GtkCanvas(img.width, img.height)
    push!(win, can)
    showall(win)

    @guarded draw(can) do widget
        ctx = getgc(can)
        Cairo.set_source_surface(ctx, img, 0, 0)
        Cairo.paint(ctx)
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        s = event.state
        if s == 0x00000004 || s == 0x00000014 # ctrl
            if k == 0x00000071 # q
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

"""
    iview(c)

View Cairo surface object.

# Arguments

- `c::Cairo.CairoSurfaceBase{UInt32}`

# Returns

Nothing
"""
function iview(c::Cairo.CairoSurfaceBase{UInt32})::Nothing

    win = GtkWindow("NeuroAnalyzer: iview()", c.width, c.height)
    set_gtk_property!(win, :border_width, 0)
    set_gtk_property!(win, :resizable, false)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    set_gtk_property!(win, :startup_id, "org.neuroanalyzer")
    can = GtkCanvas(c.width, c.height)
    push!(win, can)
    showall(win)

    @guarded draw(can) do widget
        ctx = getgc(can)
        Cairo.set_source_rgb(ctx, 255, 255, 255)
        Cairo.set_source_surface(ctx, c, 1, 1)
        Cairo.paint(ctx)
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        s = event.state
        if s == 0x00000004
            if k == 115 # s
                file_name = save_dialog("Pick image file", GtkNullContainer(), (GtkFileFilter("*.png", name="All supported formats"), "*.png"))
                    if file_name != ""
                        if splitext(file_name)[2] in [".png"]
                            try
                                surface_buf = Gtk.cairo_surface(can)
                                Cairo.write_to_png(surface_buf, file_name)
                                _info("Plot saved as: $file_name")
                            catch
                                warn_dialog("File $file_name cannot be written!")
                            end
                        else
                            warn_dialog("Incorrect filename!")
                        end
                    end
            elseif k == 0x00000071 # q
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
