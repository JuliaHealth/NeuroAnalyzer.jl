export iselect_ts

"""
    iselect_ts(obj; <keyword arguments>)

Select time segment.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=get_channel_bytype(obj, type=datatype(obj))`: channel(s) to plot, default is EEG/MEG/ERP channels
- `mono::Bool=true`: use color or gray palette
- `zoom::Real=5`: how many seconds are displayed in one segment
- `snap::Bool=true`: snap region markers to grid at 0.0, 0.25, 0.5 and 0.75 time points

# Returns

- `seg::Tuple{Float64, Float64}`
"""
function iselect_ts(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=get_channel_bytype(obj, type=datatype(obj)), mono::Bool=true, zoom::Real=5, snap::Bool=true)

    quit = false

    _check_datatype(obj, ["eeg", "meg", "erp"])

    nepochs(obj) > 1 && (zoom = epoch_len(obj) / sr(obj))

    (signal_len(obj) / sr(obj)) < zoom && (zoom = obj.time_pts[end])

    @assert zoom > 0 "zoom must be > 0."
    @assert zoom <= signal_len(obj) / sr(obj) "zoom must be ≤ $(signal_len(obj) / sr(obj))."
    _check_channels(obj, ch)

    p = NeuroAnalyzer.plot(obj, ch=ch, mono=mono, title="")
    win = GtkWindow("NeuroAnalyzer: iselect_ts()", Int32(p.attr[:size][1]), Int32(p.attr[:size][2]) + 40)
    win_view = GtkScrolledWindow()
    set_gtk_property!(win_view, :min_content_width, Int32(p.attr[:size][1]))
    set_gtk_property!(win_view, :min_content_height, Int32(p.attr[:size][2]))
    set_gtk_property!(win, :border_width, 20)
    set_gtk_property!(win, :resizable, false)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    set_gtk_property!(win, :startup_id, "org.neuroanalyzer")
    can = GtkCanvas(Int32(p.attr[:size][1]), Int32(p.attr[:size][2]))
    push!(win_view, can)
    g = GtkGrid()
    set_gtk_property!(g, :column_homogeneous, false)
    set_gtk_property!(g, :column_spacing, 10)
    set_gtk_property!(g, :row_spacing, 10)
    entry_time = GtkSpinButton(obj.time_pts[1], obj.time_pts[end] - zoom, zoom)
    set_gtk_property!(entry_time, :digits, 2)
    set_gtk_property!(entry_time, :value, obj.time_pts[1])
    set_gtk_property!(entry_time, :tooltip_text, "Time position [s]")
    entry_ts1 = GtkSpinButton(obj.time_pts[1], obj.time_pts[end], 0.5)
    set_gtk_property!(entry_ts1, :tooltip_text, "Segment start [s]")
    set_gtk_property!(entry_ts1, :digits, 3)
    entry_ts2 = GtkSpinButton(obj.time_pts[1], obj.time_pts[end], 0.5)
    set_gtk_property!(entry_ts2, :digits, 3)
    set_gtk_property!(entry_ts2, :tooltip_text, "Segment end [s]")
    bt_start = GtkButton("⇤")
    set_gtk_property!(bt_start, :tooltip_text, "Go to the signal beginning")
    bt_prev5 = GtkButton("↞")
    nepochs(obj) == 1 && set_gtk_property!(bt_prev5, :tooltip_text, "Go back by $zoom seconds")
    nepochs(obj) > 1 && set_gtk_property!(bt_prev5, :tooltip_text, "Go back by 1 epoch")
    bt_prev = GtkButton("←")
    set_gtk_property!(bt_prev, :tooltip_text, "Go back by 1 second")
    bt_next = GtkButton("→")
    set_gtk_property!(bt_next, :tooltip_text, "Go forward by 1 second")
    bt_next5 = GtkButton("↠")
    nepochs(obj) == 1 && set_gtk_property!(bt_next5, :tooltip_text, "Go forward by $zoom seconds")
    nepochs(obj) > 1 && set_gtk_property!(bt_next5, :tooltip_text, "Go forward by 1 epoch")
    bt_end = GtkButton("⇥")
    set_gtk_property!(bt_end, :tooltip_text, "Go to the signal end")
    bt_help = GtkButton("🛈")
    set_gtk_property!(bt_help, :tooltip_text, "Show keyboard shortcuts")
    bt_ts = GtkButton("Return TS")
    set_gtk_property!(bt_ts, :tooltip_text, "Return selected time segment")
    bt_close = GtkButton("✖")
    set_gtk_property!(bt_close, :tooltip_text, "Close this window")
    g[1:16, 1] = win_view
    g[1, 2] = bt_start
    g[2, 2] = bt_prev5
    g[3, 2] = bt_prev
    g[4, 2] = entry_time
    g[5, 2] = bt_next
    g[6, 2] = bt_next5
    g[7, 2] = bt_end
    g[8, 2] = GtkLabel("")
    g[9, 2] = entry_ts1
    g[10, 2] = GtkLabel("|")
    g[11, 2] = entry_ts2
    g[12, 2] = GtkLabel("")
    g[13, 2] = bt_ts
    g[14, 2] = GtkLabel("")
    g[15, 2] = bt_help
    g[16, 2] = bt_close
    push!(win, g)

    showall(win)

    @guarded draw(can) do widget
        time1 = get_gtk_property(entry_time, :value, Float64)
        time2 = time1 + zoom
        time2 > obj.time_pts[end] && (time2 = obj.time_pts[end])
        ts1 = get_gtk_property(entry_ts1, :value, Float64)
        ts2 = get_gtk_property(entry_ts2, :value, Float64)
        ctx = getgc(can)
        show(io, MIME("image/png"), NeuroAnalyzer.plot(obj,
                                                       ch=ch,
                                                       seg=(time1, time2),
                                                       s_pos=(ts1, ts2),
                                                       mono=mono,
                                                       title=""))
        img = read_from_png(io)
        set_source_surface(ctx, img, 0, 0)
        paint(ctx)
    end

    signal_connect(entry_time, "value-changed") do widget
        time_current = get_gtk_property(entry_time, :value, Float64)
        Gtk.@sigatom begin
            set_gtk_property!(entry_ts1, :value, time_current)
            set_gtk_property!(entry_ts2, :value, time_current)
        end
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

    can.mouse.button1press = @guarded (widget, event) -> begin
        time_current = get_gtk_property(entry_time, :value, Float64)
        x_pos = event.x
        x_pos < 52 && (x_pos = 52)
        x_pos > 1182 && (x_pos = 1182)
        if time_current + zoom < obj.time_pts[end]
            ts1 = time_current + round((x_pos - 52) / (1130 / zoom), digits=3)
        else
            ts1 = time_current + round((x_pos - 52) / (1130 / (obj.time_pts[end] - time_current)), digits=3)
        end
        snap && (ts1 = round(ts1 * 4) / 4)
        round(ts1, digits=3)
        Gtk.@sigatom begin
            set_gtk_property!(entry_ts1, :value, ts1)
        end
    end

    can.mouse.button3press = @guarded (widget, event) -> begin
        time_current = get_gtk_property(entry_time, :value, Float64)
        x_pos = event.x
        x_pos < 52 && (x_pos = 52)
        x_pos > 1182 && (x_pos = 1182)
        if time_current + zoom < obj.time_pts[end]
            ts2 = time_current + ((x_pos - 52) / (1130 / zoom))
        else
            ts2 = time_current + ((x_pos - 52) / (1130 / (obj.time_pts[end] - time_current)))
        end
        snap && (ts2 = round(ts2 * 4) / 4)
        round(ts2, digits=3)
        Gtk.@sigatom begin
            set_gtk_property!(entry_ts2, :value, round(ts2, digits=3))
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
        else
            time_current = obj.time_pts[end] - zoom
        end
        Gtk.@sigatom begin
            set_gtk_property!(entry_time, :value, time_current)
        end
    end

    signal_connect(bt_next5, "clicked") do widget
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

    signal_connect(bt_help, "clicked") do widgete
        nepochs(obj) == 1 && info_dialog("Keyboard shortcuts:\n\nctrl-a\tgo to the signal beginning\nctrl-s\tgo to the signal end\nctrl-z\tgo back by 1 second\nctrl-x\tgo forward by 1 second\nctrl-c\tgo back by $zoom seconds\nctrl-v\tgo forward by $zoom seconds\n\nctrl-p\treturn selected time segment\n\nctrl-\\\tswitch snapping\n\nctrl-h\tthis info\nctrl-q\texit\n")
        nepochs(obj) > 1 && info_dialog("Keyboard shortcuts:\n\nctrl-a\tgo to the signal beginning\nctrl-s\tgo to the signal end\nctrl-z\tgo back by 1 second\nctrl-x\tgo forward by 1 second\nctrl-c\tgo back by 1 epoch\nctrl-v\tgo forward by 1 epoch\n\nctrl-p\treturn selected time segment\n\nctrl-\\\tswitch snapping\n\nctrl-h\tthis info\nctrl-q\texit\n")
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        s = event.state
        if s == 4
            if k == 113 # q
                quit = true
                Gtk.destroy(win)
            elseif k == 104 # h
                nepochs(obj) == 1 && info_dialog("Keyboard shortcuts:\n\nctrl-a\tgo to the signal beginning\nctrl-s\tgo to the signal end\nctrl-z\tgo back by 1 second\nctrl-x\tgo forward by 1 second\nctrl-c\tgo back by $zoom seconds\nctrl-v\tgo forward by $zoom seconds\n\nctrl-p\treturn selected time segment\n\nctrl-\\\tswitch snapping\n\nctrl-h\tthis info\nctrl-q\texit\n")
                nepochs(obj) > 1 && info_dialog("Keyboard shortcuts:\n\nctrl-a\tgo to the signal beginning\nctrl-s\tgo to the signal end\nctrl-z\tgo back by 1 second\nctrl-x\tgo forward by 1 second\nctrl-c\tgo back by 1 epoch\nctrl-v\tgo forward by 1 epoch\n\nctrl-p\treturn selected time segment\n\nctrl-\\\tswitch snapping\n\nctrl-h\tthis info\nctrl-q\texit\n")
            elseif k == 92 # \
                snap = !snap
            elseif k == 97 # a
                Gtk.@sigatom begin
                    set_gtk_property!(entry_time, :value, obj.time_pts[1])
                    set_gtk_property!(entry_ts1, :value, obj.time_pts[1])
                    set_gtk_property!(entry_ts1, :value, obj.time_pts[1])
                end
                draw(can)
            elseif k == 115 # s
                time_current = obj.time_pts[end] - zoom
                Gtk.@sigatom begin
                    set_gtk_property!(entry_time, :value, time_current)
                    set_gtk_property!(entry_ts1, :value, time_current)
                    set_gtk_property!(entry_ts2, :value, time_current)
                end
                draw(can)
            elseif k == 122 # z
                time_current = get_gtk_property(entry_time, :value, Float64)
                if time_current >= obj.time_pts[1] + 1
                    time_current -= 1
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_time, :value, time_current)
                        set_gtk_property!(entry_ts1, :value, time_current)
                        set_gtk_property!(entry_ts2, :value, time_current)
                    end
                end
                draw(can)
            elseif k == 99 # c
                time_current = get_gtk_property(entry_time, :value, Float64)
                if time_current >= obj.time_pts[1] + zoom
                    time_current = time_current - zoom
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_time, :value, time_current)
                        set_gtk_property!(entry_ts1, :value, time_current)
                        set_gtk_property!(entry_ts2, :value, time_current)
                    end
                end
                draw(can)
            elseif k == 120 # x
                time_current = get_gtk_property(entry_time, :value, Float64)
                if time_current < obj.time_pts[end] - zoom
                    time_current += 1
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_time, :value, time_current)
                        set_gtk_property!(entry_ts1, :value, time_current)
                        set_gtk_property!(entry_ts2, :value, time_current)
                    end
                else
                    time_current = obj.time_pts[end] - zoom
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_time, :value, time_current)
                        set_gtk_property!(entry_ts1, :value, time_current)
                        set_gtk_property!(entry_ts2, :value, time_current)
                    end
                end
            elseif k == 118 # v
                time_current = get_gtk_property(entry_time, :value, Float64)
                if time_current < obj.time_pts[end] - zoom
                    time_current += zoom
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_time, :value, time_current)
                        set_gtk_property!(entry_ts1, :value, time_current)
                        set_gtk_property!(entry_ts2, :value, time_current)
                    end
                else
                    time_current = obj.time_pts[end] - zoom
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_time, :value, time_current)
                        set_gtk_property!(entry_ts1, :value, time_current)
                        set_gtk_property!(entry_ts2, :value, time_current)
                    end
                end
            elseif k == 112 # p
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
    if !quit
        time_current = get_gtk_property(entry_time, :value, Float64)
        time1 = obj.time_pts[vsearch(get_gtk_property(entry_ts1, :value, Float64), obj.time_pts)]
        time2 = obj.time_pts[vsearch(get_gtk_property(entry_ts2, :value, Float64), obj.time_pts)]
        seg = (time1, time2)
        return seg
    else
        return nothing
    end

end
