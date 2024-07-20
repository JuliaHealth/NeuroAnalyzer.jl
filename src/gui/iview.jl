export iview
export iview_ep
export iview_cont

"""
    iview(obj; <keyword arguments>)

Interactive view of continuous or epoched signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Vector{Int64}, AbstractRange}=_c(nchannels(obj))`: channels to plot, default is all channels
- `ep::Int64=1`: initial epoch to display
- `zoom::Real=10`: how many seconds are displayed in one segment
- `bad::Bool=true`: show bad channels
- `snap::Bool=true`: snap region markers to grid at 0.0, 0.25, 0.5 and 0.75 time points

# Returns

- `seg::Union{Nothing, Tuple{Float64, Float64}}`
"""
function iview(obj::NeuroAnalyzer.NEURO; ch::Union{Vector{Int64}, AbstractRange}=_c(nchannels(obj)), ep::Int64=1, zoom::Real=10, bad::Bool=true, snap::Bool=true)

    @assert length(ch) > 1 "iview() requires > 1 channels."
    seg = nothing
    if nepochs(obj) == 1
        seg = iview_cont(obj, ch=ch, zoom=zoom, bad=bad, snap=snap)
    else
        seg = iview_ep(obj, ch=ch, ep=ep, bad=bad)
    end

    return seg

end

"""
    iview_cont(obj; <keyword arguments>)

Interactive view of continuous signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Vector{Int64}, AbstractRange}=_c(nchannels(obj))`: channels to plot, default is all channels
- `zoom::Real=10`: how many seconds are displayed in one segment
- `bad::Bool=true`: list of bad channels; if not false -- plot bad channels using this list
- `snap::Bool=true`: snap region markers to grid at 0.0, 0.25, 0.5 and 0.75 time points

# Returns

- `seg::Union{Nothing, Tuple{Float64, Float64}}`
"""
function iview_cont(obj::NeuroAnalyzer.NEURO; ch::Union{Vector{Int64}, AbstractRange}=_c(nchannels(obj)), zoom::Real=10, bad::Bool=true, snap::Bool=true)

    obj.time_pts[end] < zoom && (zoom = obj.time_pts[end])

    @assert length(ch) > 1 "iview() requires > 1 channels."
    @assert zoom > 0 "zoom must be > 0."
    @assert zoom <= signal_len(obj) / sr(obj) "zoom must be ≤ $(signal_len(obj) / sr(obj))."
    @assert nepochs(obj) == 1 "iview_ep() should be used for epoched object."
    _check_channels(obj, ch)

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

    if length(ch) > 1
        p = NeuroAnalyzer.plot(obj, ch=ch[ch_first:ch_last], title="", bad=bad)
    else
        p = NeuroAnalyzer.plot(obj, ch=ch, title="", bad=bad)
    end

    win = GtkWindow("NeuroAnalyzer: iview_cont()", Int32(p.attr[:size][1]) + 40, Int32(p.attr[:size][2]) + 40)
    set_gtk_property!(win, :border_width, 20)
    set_gtk_property!(win, :resizable, true)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    set_gtk_property!(win, :startup_id, "org.neuroanalyzer")
    can = GtkCanvas(Int32(p.attr[:size][1]), Int32(p.attr[:size][2]))
    g = GtkGrid()
    set_gtk_property!(g, :column_homogeneous, false)
    set_gtk_property!(g, :column_spacing, 10)
    set_gtk_property!(g, :row_spacing, 10)
    entry_time = GtkSpinButton(obj.time_pts[1], obj.time_pts[end] - zoom, 1)
    set_gtk_property!(entry_time, :digits, 2)
    set_gtk_property!(entry_time, :value, obj.time_pts[1])
    set_gtk_property!(entry_time, :tooltip_text, "Time position [s]")
    bt_start = GtkButton("⇤")
    set_gtk_property!(bt_start, :tooltip_text, "Go to the signal beginning")
    bt_prev5 = GtkButton("↞")
    set_gtk_property!(bt_prev5, :tooltip_text, "Go back by $zoom seconds")
    bt_next5 = GtkButton("↠")
    set_gtk_property!(bt_next5, :tooltip_text, "Go forward by $zoom seconds")
    bt_end = GtkButton("⇥")
    set_gtk_property!(bt_end, :tooltip_text, "Go to the signal end")
    entry_ts1 = GtkSpinButton(obj.time_pts[1], obj.time_pts[end], 0.5)
    bt_help = GtkButton("🛈")
    set_gtk_property!(bt_help, :tooltip_text, "Show keyboard shortcuts")
    bt_close = GtkButton("✖")
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
    set_gtk_property!(entry_ts1, :tooltip_text, "Segment start [s]")
    set_gtk_property!(entry_ts1, :digits, 3)
    entry_ts2 = GtkSpinButton(obj.time_pts[1], obj.time_pts[end], 0.5)
    set_gtk_property!(entry_ts2, :digits, 3)
    set_gtk_property!(entry_ts2, :tooltip_text, "Segment end [s]")
    bt_ts = GtkButton("Return TS")
    set_gtk_property!(bt_ts, :tooltip_text, "Return selected time segment")
    bt_delete = GtkButton("Delete TS")
    set_gtk_property!(bt_delete, :tooltip_text, "Delete selected time segment")
    set_gtk_property!(ch_slider, :vexpand, true)
    oc = GtkOrientable(ch_slider)
    set_gtk_property!(oc, :orientation, 1)
    g[1:10, 1] = can
    g[11, 1] = ch_slider
    g[1, 2] = bt_start
    g[2, 2] = bt_prev5
    g[4, 2] = entry_time
    g[6, 2] = bt_next5
    g[7, 2] = bt_end
    g[8, 2] = entry_ts1
    g[9, 2] = bt_ts
    g[10, 2] = bt_help
    g[8, 3] = entry_ts2
    g[9, 3] = bt_delete
    g[10, 3] = bt_close
    push!(win, g)

    showall(win)

    @guarded draw(can) do widget
        time1 = get_gtk_property(entry_time, :value, Float64)
        time2 = time1 + zoom
        time2 > obj.time_pts[end] && (time2 = obj.time_pts[end])
        ts1 = get_gtk_property(entry_ts1, :value, Float64)
        ts2 = get_gtk_property(entry_ts2, :value, Float64)
        ctx = getgc(can)

        if length(ch) > 1
            show(io, MIME("image/png"), NeuroAnalyzer.plot(obj,
                                                           ch=ch[ch_first]:ch[ch_last],
                                                           seg=(time1, time2),
                                                           s_pos=(ts1, ts2),
                                                           mono=mono,
                                                           title="",
                                                           scale=scale,
                                                           bad=bad))
        else
            show(io, MIME("image/png"), NeuroAnalyzer.plot(obj,
                                                           ch=ch,
                                                           seg=(time1, time2),
                                                           s_pos=(ts1, ts2),
                                                           mono=mono,
                                                           title="",
                                                           scale=scale,
                                                           bad=bad))
        end

        img = read_from_png(io)
        set_source_surface(ctx, img, 0, 0)
        paint(ctx)
    end

    can.mouse.button2press = @guarded (widget, event) -> begin
        x_pos = event.x
        y_pos = event.y
        ch_n = length(ch)
        ch_x1 = 0
        ch_x2 = 82
        ch_y = collect(50:39:793)[1:length(ch_first:ch_last)]
        ch_idx = nothing
        for idx in eachindex(ch_y)
            if y_pos > ch_y[idx] && y_pos < ch_y[idx] + 15 && x_pos > 0 && x_pos < 82
                ch_idx = idx
            end
        end
    end

    can.mouse.button1press = @guarded (widget, event) -> begin
        x_pos = event.x
        if x_pos < 82
            y_pos = event.y
            ch_y = collect(50:39:793)[1:length(ch_first:ch_last)]
            ch_idx = nothing
            for idx in eachindex(ch_y)
                if y_pos > ch_y[idx] && y_pos < ch_y[idx] + 15 && x_pos > 0 && x_pos < 82
                    ch_idx = idx + ch_first - 1
                end
            end
            !isnothing(ch_idx) && channel_info(obj, ch=obj.header.recording[:channel_order][ch_idx])
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
            y_pos = event.y
            ch_y = collect(50:39:793)[1:length(ch_first:ch_last)]
            ch_idx = nothing
            for idx in eachindex(ch_y)
                if y_pos > ch_y[idx] && y_pos < ch_y[idx] + 15 && x_pos > 0 && x_pos < 82
                    ch_idx = idx + ch_first - 1
                end
            end
            !isnothing(ch_idx) && (obj.header.recording[:bad_channels][obj.header.recording[:channel_order][ch_idx, 1]] = !obj.header.recording[:bad_channels][obj.header.recording[:channel_order][ch_idx, 1]])
            draw(can)
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
            if s == 0x00000011
                time_current = get_gtk_property(entry_time, :value, Float64)
                if time_current < obj.time_pts[end] - zoom
                    time_current += 1
                else
                    time_current = obj.time_pts[end] - zoom
                end
                Gtk.@sigatom begin
                    set_gtk_property!(entry_time, :value, time_current)
                end
            elseif s == 0x00000014
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
            if s == 0x00000011
                time_current = get_gtk_property(entry_time, :value, Float64)
                if time_current >= obj.time_pts[1] + 1
                    time_current -= 1
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_time, :value, time_current)
                    end
                end
            elseif s == 0x00000014
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

    signal_connect(ch_slider, "value-changed") do widget, others...
        ch_first = round(Int64, GAccessor.value(ch_slider))
        ch_last = ch_first + 19
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

    help = "Keyboard shortcuts:\n\nPage Up\t\tscroll channels up\nPage Down\tscroll channels down\n\nHome\tgo to the signal beginning\nEnd\tgo to the signal end\nctrl-,\tgo back by 1 second\nctrl-.\tgo forward by 1 second\nalt-,\tgo back by $zoom seconds\nalt-.\tgo forward by $zoom seconds\n\n[\t zoom in\n]\tzoom out\n\nctrl-Enter\treturn selected time segment\nctrl-d\t\tdelete selected time segment\n\nalt-s\ttoggle snapping\nctrl-s\ttoggle scales\nctrl-m\ttoggle monochromatic mode\n\nctrl-h\tthis info\nctrl-q\texit\n"

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
                set_gtk_property!(bt_next5, :tooltip_text, "Go forward by $zoom seconds")
                set_gtk_property!(bt_prev5, :tooltip_text, "Go back by $zoom seconds")
                draw(can)
            end
            help = "Keyboard shortcuts:\n\nPage Up\t\tscroll channels up\nPage Down\tscroll channels down\n\nHome\tgo to the signal beginning\nEnd\tgo to the signal end\nctrl-,\tgo back by 1 second\nctrl-.\tgo forward by 1 second\nalt-,\tgo back by $zoom seconds\nalt-.\tgo forward by $zoom seconds\n\n[\t zoom in\n]\tzoom out\n\nctrl-s\ttoggle scales\n\nctrl-h\tthis info\nctrl-q\texit\n"
        elseif k == 0x0000005d # ]
            if zoom < 30 && zoom < obj.time_pts[end] - 1
                zoom += 1
                set_gtk_property!(bt_next5, :tooltip_text, "Go forward by $zoom seconds")
                set_gtk_property!(bt_prev5, :tooltip_text, "Go back by $zoom seconds")
                draw(can)
            else
                zoom = obj.time_pts[end]
                set_gtk_property!(bt_next5, :tooltip_text, "Go forward by $zoom seconds")
                set_gtk_property!(bt_prev5, :tooltip_text, "Go back by $zoom seconds")
                draw(can)
            end
            help = "Keyboard shortcuts:\n\nPage Up\t\tscroll channels up\nPage Down\tscroll channels down\n\nHome\tgo to the signal beginning\nEnd\tgo to the signal end\nctrl-,\tgo back by 1 second\nctrl-.\tgo forward by 1 second\nalt-,\tgo back by $zoom seconds\nalt-.\tgo forward by $zoom seconds\n\n[\t zoom in\n]\tzoom out\n\nctrl-s\ttoggle scales\n\nctrl-h\tthis info\nctrl-q\texit\n"
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
        if s == 0x00000018 # alt
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
            elseif k == 0x00000073 # s
                snap = !snap
            end
        end
        if s == 0x00000014 # ctrl
            if k == 113 # q
                quit = true
                Gtk.destroy(win)
            elseif k == 0x0000ff0d # Enter
                Gtk.destroy(win)
            elseif k == 104 # h
                info_dialog(help)
            elseif k == 0x0000006d # m
                mono = !mono
                draw(can)
            elseif k == 0x00000073 # s
                scale = !scale
                draw(can)
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
    if !quit
        time1 = obj.time_pts[vsearch(get_gtk_property(entry_ts1, :value, Float64), obj.time_pts)]
        time2 = obj.time_pts[vsearch(get_gtk_property(entry_ts2, :value, Float64), obj.time_pts)]
        seg = (time1, time2)
        return seg
    else
        return nothing
    end

end

"""
    iview_ep(obj; <keyword arguments>)

Interactive view of epoched signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Vector{Int64}, AbstractRange}=_c(nchannels(obj))`: channels to plot, default is all channels
- `ep::Int64=1`: initial epoch to display
- `bad::Bool=true`: list of bad channels; if not false -- plot bad channels using this list
- `snap::Bool=true`: snap region markers to grid at 0.0, 0.25, 0.5 and 0.75 time points

# Returns

- `seg::Union{Nothing, Tuple{Float64, Float64}}`
"""
function iview_ep(obj::NeuroAnalyzer.NEURO; ch::Union{Vector{Int64}, AbstractRange}=_c(nchannels(obj)), ep::Int64=1, bad::Bool=true, snap::Bool=true)

    @assert length(ch) > 1 "iview() requires > 1 channels."
    @assert nepochs(obj) > 1 "iview_cont() should be used for continuous object."
    _check_channels(obj, ch)
    _check_epochs(obj, ep)

    mono = false
    quit = false
    scale = true
    zoom = epoch_len(obj)

    if length(ch) > 20
        ch_first = 1
        ch_last = 20
    else
        ch_first = 1
        ch_last = length(ch)
    end

    if length(ch) > 1
        p = NeuroAnalyzer.plot(obj, ch=ch[ch_first:ch_last], ep=ep, title="", bad=bad)
    else
        p = NeuroAnalyzer.plot(obj, ch=ch, ep=ep, mono=mono, title="", bad=bad)
    end
    win = GtkWindow("NeuroAnalyzer: iview_ep()", Int32(p.attr[:size][1]) + 40, Int32(p.attr[:size][2]) + 40)
    set_gtk_property!(win, :border_width, 20)
    set_gtk_property!(win, :resizable, true)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    set_gtk_property!(win, :startup_id, "org.neuroanalyzer")
    can = GtkCanvas(Int32(p.attr[:size][1]), Int32(p.attr[:size][2]))
    g = GtkGrid()
    set_gtk_property!(g, :column_homogeneous, false)
    set_gtk_property!(g, :column_spacing, 10)
    set_gtk_property!(g, :row_spacing, 10)
    entry_epoch = GtkSpinButton(1, nepochs(obj), 1)
    set_gtk_property!(entry_epoch, :value, ep)
    set_gtk_property!(entry_epoch, :tooltip_text, "Epoch")
    bt_start = GtkButton("⇤")
    set_gtk_property!(bt_start, :tooltip_text, "Go to the signal beginning")
    bt_end = GtkButton("⇥")
    set_gtk_property!(bt_end, :tooltip_text, "Go to the signal end")
    bt_help = GtkButton("🛈")
    set_gtk_property!(bt_help, :tooltip_text, "Show keyboard shortcuts")
    bt_close = GtkButton("✖")
    set_gtk_property!(bt_close, :tooltip_text, "Close this window")
    entry_ts1 = GtkSpinButton(obj.epoch_time[1], obj.epoch_time[end], 0.5)
    set_gtk_property!(entry_ts1, :tooltip_text, "Segment start [s]")
    set_gtk_property!(entry_ts1, :digits, 3)
    entry_ts2 = GtkSpinButton(obj.epoch_time[1], obj.epoch_time[end], 0.5)
    set_gtk_property!(entry_ts2, :digits, 3)
    set_gtk_property!(entry_ts2, :tooltip_text, "Segment end [s]")
    bt_ts = GtkButton("Return TS")
    set_gtk_property!(bt_ts, :tooltip_text, "Return selected time segment")
    bt_delete = GtkButton("Delete epoch")
    set_gtk_property!(bt_delete, :tooltip_text, "Delete current epoch")
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
    g[1:6, 1] = can
    g[7, 1] = ch_slider
    g[1, 2] = bt_start
    g[2, 2] = entry_epoch
    g[3, 2] = bt_end
    g[4, 2] = entry_ts1
    g[5, 2] = bt_ts
    g[6, 2] = bt_help
    g[4, 3] = entry_ts2
    g[5, 3] = bt_delete
    g[6, 3] = bt_close
    push!(win, g)

    showall(win)

    @guarded draw(can) do widget
        ep = get_gtk_property(entry_epoch, :value, Int64)
        ts1 = get_gtk_property(entry_ts1, :value, Float64)
        ts2 = get_gtk_property(entry_ts2, :value, Float64)
        ctx = getgc(can)
        if length(ch) > 1
            show(io, MIME("image/png"), NeuroAnalyzer.plot(obj,
                                                           ch=ch[ch_first]:ch[ch_last],
                                                           ep=ep,
                                                           s_pos=(ts1, ts2),
                                                           mono=mono,
                                                           title="",
                                                           scale=scale,
                                                           bad=bad))
        else
            show(io, MIME("image/png"), NeuroAnalyzer.plot(obj,
                                                           ch=ch,
                                                           ep=ep,
                                                           s_pos=(ts1, ts2),
                                                           mono=mono,
                                                           title="",
                                                           scale=scale,
                                                           bad=bad))
        end

        img = read_from_png(io)
        set_source_surface(ctx, img, 0, 0)
        paint(ctx)
    end

    can.mouse.button1press = @guarded (widget, event) -> begin
        x_pos = event.x
        if x_pos < 82
            y_pos = event.y
            ch_y = collect(50:39:793)[1:length(ch_first:ch_last)]
            ch_idx = nothing
            for idx in eachindex(ch_y)
                if y_pos > ch_y[idx] && y_pos < ch_y[idx] + 15 && x_pos > 0 && x_pos < 82
                    ch_idx = idx + ch_first - 1
                end
            end
            !isnothing(ch_idx) && channel_info(obj, ch=obj.header.recording[:channel_order][ch_idx])
        else
            x_pos > 1172 && (x_pos = 1172)
            ts1 = obj.epoch_time[1] + (x_pos - 82) / (1090 / (obj.epoch_time[end] - obj.epoch_time[1]))
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
            ep = get_gtk_property(entry_epoch, :value, Int64)
            y_pos = event.y
            ch_y = collect(50:39:793)[1:length(ch_first:ch_last)]
            ch_idx = nothing
            for idx in eachindex(ch_y)
                if y_pos > ch_y[idx] && y_pos < ch_y[idx] + 15 && x_pos > 0 && x_pos < 82
                    ch_idx = idx + ch_first - 1
                end
            end
            !isnothing(ch_idx) && (obj.header.recording[:bad_channels][obj.header.recording[:channel_order][ch_idx], ep] = !obj.header.recording[:bad_channels][obj.header.recording[:channel_order][ch_idx], ep])
            draw(can)
        else
            x_pos > 1172 && (x_pos = 1172)
            ts2 = obj.epoch_time[1] + ((x_pos - 82) / (1090 / (obj.epoch_time[end] - obj.epoch_time[1])))
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
            if s == 0x00000011
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
            if s == 0x00000011
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

    signal_connect(ch_slider, "value-changed") do widget, others...
        ch_first = round(Int64, GAccessor.value(ch_slider))
        ch_last = ch_first + 19
        draw(can)
    end

    signal_connect(entry_epoch, "value-changed") do widget
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
        quit = true
        Gtk.destroy(win)
    end

    signal_connect(entry_ts1, "value-changed") do widget
        Gtk.@sigatom begin
            set_gtk_property!(entry_ts1, :value, obj.epoch_time[vsearch(get_gtk_property(entry_ts1, :value, Float64), obj.epoch_time)])
        end
        draw(can)
    end

    signal_connect(entry_ts2, "value-changed") do widget
        Gtk.@sigatom begin
            set_gtk_property!(entry_ts2, :value, obj.epoch_time[vsearch(get_gtk_property(entry_ts2, :value, Float64), obj.epoch_time)])
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
                end
            end
            draw(can)
        else
            error_dialog("You cannot delete the last epoch.")
        end
    end

    signal_connect(bt_ts, "clicked") do widget
        quit = true
        Gtk.destroy(win)
    end

    help = "Keyboard shortcuts:\n\nPage Up\t\tscroll channels up\nPage Down\tscroll channels down\n\nHome\tgo to first epoch\nEnd\tgo to last epoch\nctrl-,\tprevious epoch\nctrl-.\tnext epoch\n\nctrl-Enter\treturn selected time segment\nctrl-d\t\tdelete current epoch\n\nalt-s\ttoggle snapping\nctrl-s\ttoggle scales\nctrl-m\ttoggle monochromatic mode\n\nctrl-h\tthis info\nctrl-q\texit\n"

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
        if s == 0x00000018 # alt
            if k == 0x00000073 # s
                snap = !snap
            end
        end
        if s == 0x00000014 # ctrl
            if k == 113 # q
                quit = true
                Gtk.destroy(win)
            elseif k == 0x0000ff0d # Enter
                quit = false
                Gtk.destroy(win)
            elseif k == 104 # h
                info_dialog(help)
            elseif k == 0x0000006d # m
                mono = !mono
                draw(can)
            elseif k == 0x00000073 # s
                scale = !scale
                draw(can)
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

Interactive view of two continuous or epoched signals.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `obj2::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Vector{Int64}, AbstractRange}=_c(nchannels(obj1))`: channels to plot, default is all channels
- `ep::Int64=1`: initial epoch to display
- `zoom::Real=10`: how many seconds are displayed in one segment
"""
function iview(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch::Union{Vector{Int64}, AbstractRange}=_c(nchannels(obj1)), ep::Int64=1, zoom::Real=10)

    @assert length(ch) > 1 "iview() requires > 1 channels."
    @assert size(obj1) == size(obj2) "Both signals must have the same size."
    @assert sr(obj1) == sr(obj2) "Both signals must have the same sampling rate."

    if nepochs(obj1) == 1
        iview_cont(obj1, obj2, ch=ch, zoom=zoom)
    else
        iview_ep(obj1, obj2, ch=ch, ep=ep)
    end

    return nothing

end

"""
    iview_cont(obj1, obj2; <keyword arguments>)

Interactive view of two continuous signals.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `obj2::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Vector{Int64}, AbstractRange}=_c(nchannels(obj1))`: channels to plot, default is all channels
- `zoom::Real=10`: how many seconds are displayed in one segment
"""
function iview_cont(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch::Union{Vector{Int64}, AbstractRange}=_c(nchannels(obj1)), zoom::Real=10)

    (signal_len(obj1) / sr(obj1)) < zoom && (zoom = obj1.time_pts[end])

    @assert length(ch) > 1 "iview() requires > 1 channels."
    @assert size(obj1) == size(obj2) "Both signals must have the same size."
    @assert sr(obj1) == sr(obj2) "Both signals must have the same sampling rate."
    @assert zoom > 0 "zoom must be > 0."
    @assert zoom <= signal_len(obj1) / sr(obj1) "zoom must be ≤ $(signal_len(obj1) / sr(obj1))."
    @assert nepochs(obj1) == 1 "iview_ep() should be used for epoched object."
    _check_channels(obj1, ch)

    scale = true

    if length(ch) > 20
        ch_first = 1
        ch_last = 20
    else
        ch_first = 1
        ch_last = length(ch)
    end

    if length(ch) > 1
        p = NeuroAnalyzer.plot(obj1, obj2, ch=ch[ch_first]:ch[ch_last], title="")
    else
        p = NeuroAnalyzer.plot(obj1, obj2, ch=ch, title="")
    end

    win = GtkWindow("NeuroAnalyzer: iview_cont()", Int32(p.attr[:size][1]) + 40, Int32(p.attr[:size][2]) + 40)
    set_gtk_property!(win, :border_width, 20)
    set_gtk_property!(win, :resizable, true)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    set_gtk_property!(win, :startup_id, "org.neuroanalyzer")
    can = GtkCanvas(Int32(p.attr[:size][1]), Int32(p.attr[:size][2]))
    g = GtkGrid()
    set_gtk_property!(g, :column_homogeneous, false)
    set_gtk_property!(g, :column_spacing, 10)
    set_gtk_property!(g, :row_spacing, 10)
    entry_time = GtkSpinButton(obj1.time_pts[1], obj1.time_pts[end] - zoom, 1)
    set_gtk_property!(entry_time, :digits, 2)
    set_gtk_property!(entry_time, :value, obj1.time_pts[1])
    set_gtk_property!(entry_time, :tooltip_text, "Time position [s]")
    bt_start = GtkButton("⇤")
    set_gtk_property!(bt_start, :tooltip_text, "Go to the signal beginning")
    bt_prev5 = GtkButton("↞")
    set_gtk_property!(bt_prev5, :tooltip_text, "Go back by $zoom seconds")
    bt_next5 = GtkButton("↠")
    set_gtk_property!(bt_next5, :tooltip_text, "Go forward by $zoom seconds")
    bt_end = GtkButton("⇥")
    set_gtk_property!(bt_end, :tooltip_text, "Go to the signal end")
    bt_help = GtkButton("🛈")
    set_gtk_property!(bt_help, :tooltip_text, "Show keyboard shortcuts")
    bt_close = GtkButton("✖")
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
    g[1:7, 1] = can
    g[8, 1] = ch_slider
    g[1, 2] = bt_start
    g[2, 2] = bt_prev5
    g[3, 2] = entry_time
    g[4, 2] = bt_next5
    g[5, 2] = bt_end
    g[6, 2] = bt_help
    g[7, 2] = bt_close
    push!(win, g)

    showall(win)

    @guarded draw(can) do widget
        time1 = get_gtk_property(entry_time, :value, Float64)
        time2 = time1 + zoom
        time2 > obj1.time_pts[end] && (time2 = obj1.time_pts[end])
        ctx = getgc(can)

        if length(ch) > 1
            show(io, MIME("image/png"), NeuroAnalyzer.plot(obj1, obj2,
                                                           ch=ch[ch_first]:ch[ch_last],
                                                           seg=(time1, time2),
                                                           title="",
                                                           scale=scale))
        else
            show(io, MIME("image/png"), NeuroAnalyzer.plot(obj1, obj2,
                                                           ch=ch,
                                                           seg=(time1, time2),
                                                           title="",
                                                           scale=scale))
        end

        img = read_from_png(io)
        set_source_surface(ctx, img, 0, 0)
        paint(ctx)
    end

    can.mouse.scroll = @guarded (widget, event) -> begin
        s = event.state
        if event.direction == 1 # down
            if s == 0x00000011
                time_current = get_gtk_property(entry_time, :value, Float64)
                if time_current < obj1.time_pts[end] - zoom
                    time_current += 1
                else
                    time_current = obj1.time_pts[end] - zoom
                end
                Gtk.@sigatom begin
                    set_gtk_property!(entry_time, :value, time_current)
                end
            elseif s == 0x00000014
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
            if s == 0x00000011
                time_current = get_gtk_property(entry_time, :value, Float64)
                if time_current >= obj1.time_pts[1] + 1
                    time_current -= 1
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_time, :value, time_current)
                    end
                end
            elseif s == 0x00000014
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

    signal_connect(ch_slider, "value-changed") do widget, others...
        ch_first = round(Int64, GAccessor.value(ch_slider))
        ch_last = ch_first + 19
        draw(can)
    end

    signal_connect(entry_time, "value-changed") do widget
        draw(can)
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

    signal_connect(bt_next5, "clicked") do widget
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

    signal_connect(bt_close, "clicked") do widget
        Gtk.destroy(win)
    end

    help = "Keyboard shortcuts:\n\nPage Up\t\tscroll channels up\nPage Down\tscroll channels down\n\nHome\tgo to the signal beginning\nEnd\tgo to the signal end\nctrl-,\tgo back by 1 second\nctrl-.\tgo forward by 1 second\nalt-,\tgo back by $zoom seconds\nalt-.\tgo forward by $zoom seconds\n\n[\t zoom in\n]\tzoom out\n\nctrl-s\ttoggle scales\n\nctrl-h\tthis info\nctrl-q\texit\n"

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
                set_gtk_property!(bt_next5, :tooltip_text, "Go forward by $zoom seconds")
                set_gtk_property!(bt_prev5, :tooltip_text, "Go back by $zoom seconds")
                draw(can)
            end
            help = "Keyboard shortcuts:\n\nPage Up\t\tscroll channels up\nPage Down\tscroll channels down\n\nHome\tgo to the signal beginning\nEnd\tgo to the signal end\nctrl-,\tgo back by 1 second\nctrl-.\tgo forward by 1 second\nalt-,\tgo back by $zoom seconds\nalt-.\tgo forward by $zoom seconds\n\n[\t zoom in\n]\tzoom out\n\nctrl-s\ttoggle scales\n\nctrl-h\tthis info\nctrl-q\texit\n"
        elseif k == 0x0000005d # ]
            if zoom < 30 && zoom < obj1.time_pts[end] - 1
                zoom += 1
                set_gtk_property!(bt_next5, :tooltip_text, "Go forward by $zoom seconds")
                set_gtk_property!(bt_prev5, :tooltip_text, "Go back by $zoom seconds")
                draw(can)
            else
                zoom = obj1.time_pts[end]
                set_gtk_property!(bt_next5, :tooltip_text, "Go forward by $zoom seconds")
                set_gtk_property!(bt_prev5, :tooltip_text, "Go back by $zoom seconds")
                draw(can)
            end
            help = "Keyboard shortcuts:\n\nPage Up\t\tscroll channels up\nPage Down\tscroll channels down\n\nHome\tgo to the signal beginning\nEnd\tgo to the signal end\nctrl-,\tgo back by 1 second\nctrl-.\tgo forward by 1 second\nalt-,\tgo back by $zoom seconds\nalt-.\tgo forward by $zoom seconds\n\n[\t zoom in\n]\tzoom out\n\nctrl-s\ttoggle scales\n\nctrl-h\tthis info\nctrl-q\texit\n"
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
        if s == 0x00000018 # alt
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
            end
        end
        if s == 0x00000014 # ctrl
            if k == 113 # q
                Gtk.destroy(win)
            elseif k == 104 # h
                info_dialog(help)
            elseif k == 0x00000073 # s
                scale = !scale
                draw(can)
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
- `ch::Union{Vector{Int64}, AbstractRange}=_c(nchannels(obj1))`: channels to plot, default is all channels
- `ep::Int64=1`: initial epoch to display
"""
function iview_ep(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch::Union{Vector{Int64}, AbstractRange}=_c(nchannels(obj1)), ep::Int64=1)

    @assert length(ch) > 1 "iview() requires > 1 channels."
    @assert size(obj1) == size(obj2) "Both signals must have the same size."
    @assert sr(obj1) == sr(obj2) "Both signals must have the same sampling rate."
    @assert nepochs(obj1) > 1 "iview_cont() should be used for continuous object."
    _check_channels(obj1, ch)
    _check_epochs(obj1, ep)

    if length(ch) > 20
        ch_first = 1
        ch_last = 20
    else
        ch_first = 1
        ch_last = length(ch)
    end

    if length(ch) > 1
        p = NeuroAnalyzer.plot(obj1, obj2, ch=ch[ch_first]:ch[ch_last], ep=ep, title="")
    else
        p = NeuroAnalyzer.plot(obj1, obj2, ch=ch, ep=ep, title="")
    end

    win = GtkWindow("NeuroAnalyzer: iview_ep()", Int32(p.attr[:size][1]) + 40, Int32(p.attr[:size][2]) + 40)
    set_gtk_property!(win, :border_width, 20)
    set_gtk_property!(win, :resizable, true)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    set_gtk_property!(win, :startup_id, "org.neuroanalyzer")
    can = GtkCanvas(Int32(p.attr[:size][1]), Int32(p.attr[:size][2]))
    g = GtkGrid()
    set_gtk_property!(g, :column_homogeneous, false)
    set_gtk_property!(g, :column_spacing, 10)
    set_gtk_property!(g, :row_spacing, 10)
    entry_epoch = GtkSpinButton(1, nepochs(obj1), 1)
    set_gtk_property!(entry_epoch, :value, ep)
    set_gtk_property!(entry_epoch, :tooltip_text, "Epoch")
    bt_start = GtkButton("⇤")
    set_gtk_property!(bt_start, :tooltip_text, "Go to the signal beginning")
    bt_end = GtkButton("⇥")
    set_gtk_property!(bt_end, :tooltip_text, "Go to the signal end")
    bt_help = GtkButton("🛈")
    set_gtk_property!(bt_help, :tooltip_text, "Show keyboard shortcuts")
    bt_close = GtkButton("✖")
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
    g[1:5, 1] = can
    g[6, 1] = ch_slider
    g[1, 2] = bt_start
    g[2, 2] = entry_epoch
    g[3, 2] = bt_end
    g[4, 2] = bt_help
    g[5, 2] = bt_close
    push!(win, g)

    showall(win)

    @guarded draw(can) do widget
        ep = get_gtk_property(entry_epoch, :value, Int64)
        ctx = getgc(can)

        if length(ch) > 1
            show(io, MIME("image/png"), NeuroAnalyzer.plot(obj1, obj2,
                                                           ch=ch[ch_first]:ch[ch_last],
                                                           ep=ep,
                                                           title=""))
        else
            show(io, MIME("image/png"), NeuroAnalyzer.plot(obj1, obj2,
                                                           ch=ch,
                                                           ep=ep,
                                                           title=""))
        end

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

    signal_connect(ch_slider, "value-changed") do widget, others...
        ch_first = round(Int64, GAccessor.value(ch_slider))
        ch_last = ch_first + 19
        draw(can)
    end

    signal_connect(entry_epoch, "value-changed") do widget
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

    help = "Keyboard shortcuts:\n\nPage Up\t\tscroll channels up\nPage Down\tscroll channels down\n\nctrl-a\tgo to first epoch\nctrl-s\tgo to last epoch\nctrl-z\tprevious epoch\nctrl-x\tnext epoch\n\nctrl-h\tthis info\nctrl-q\texit\n"

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
                set_gtk_property!(entry_epoch, :value, nepochs(obj1))
            end
        end
        if s == 0x00000014
            if k == 113 # q
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

export iview
export iview_ep
export iview_cont

"""
    iview(obj, c; <keyword arguments>)

Interactive view of embedded or external component of continuous or epoched signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `c::Union{Symbol, AbstractArray}`: component to plot
- `c_idx::Union{Vector{Int64}, AbstractRange}=_c(nchannels(obj))`: component channel to display, default is all component channels
- `ep::Int64=1`: initial epoch to display
- `zoom::Real=10`: how many seconds are displayed in one segment
"""
function iview(obj::NeuroAnalyzer.NEURO, c::Union{Symbol, AbstractArray}; c_idx::Union{Vector{Int64}, AbstractRange}=_c(nchannels(obj)), ep::Int64=1, zoom::Real=10)

    @assert length(c_idx) > 1 "iview() requires > 1 channels."

    if nepochs(obj) == 1
        iview_cont(obj, c, c_idx=c_idx, zoom=zoom)
    else
        iview_ep(obj, c, c_idx=c_idx, ep=ep)
    end

    return nothing

end

"""
    iview_cont(obj, c; <keyword arguments>)

Interactive view of embedded or external component of continuous signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `c::Union{Symbol, AbstractArray}`: component to plot
- `c_idx::Union{Vector{Int64}, AbstractRange}=_c(nchannels(obj))`: component channel to display, default is all component channels
- `zoom::Real=10`: how many seconds are displayed in one segment
"""
function iview_cont(obj::NeuroAnalyzer.NEURO, c::Union{Symbol, AbstractArray}; c_idx::Union{Vector{Int64}, AbstractRange}=_c(nchannels(obj)), zoom::Real=10)

    @assert length(c_idx) > 1 "iview() requires > 1 component."
    @assert zoom > 0 "zoom must be > 0."
    @assert zoom <= signal_len(obj) / sr(obj) "zoom must be ≤ $(signal_len(obj) / sr(obj))."
    @assert nepochs(obj) == 1 "iview_ep() should be used for epoched object."

    mono = false
    quit = false
    scale = true

    if length(c_idx) > 20
        c_idx_first = 1
        c_idx_last = 20
    else
        c_idx_first = 1
        c_idx_last = length(c_idx)
    end

    length(c_idx) == 1 && (c_idx = [c_idx])

    p = NeuroAnalyzer.plot(obj, c, c_idx=c_idx[c_idx_first:c_idx_last], mono=mono, scale=scale, title="")
    win = GtkWindow("NeuroAnalyzer: iview_cont()", Int32(p.attr[:size][1]) + 40, Int32(p.attr[:size][2]) + 40)
    set_gtk_property!(win, :border_width, 20)
    set_gtk_property!(win, :resizable, true)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    set_gtk_property!(win, :startup_id, "org.neuroanalyzer")
    can = GtkCanvas(Int32(p.attr[:size][1]), Int32(p.attr[:size][2]))
    g = GtkGrid()
    set_gtk_property!(g, :column_homogeneous, false)
    set_gtk_property!(g, :column_spacing, 10)
    set_gtk_property!(g, :row_spacing, 10)
    entry_time = GtkSpinButton(obj.time_pts[1], obj.time_pts[end] - zoom, zoom)
    set_gtk_property!(entry_time, :digits, 2)
    set_gtk_property!(entry_time, :value, obj.time_pts[1])
    set_gtk_property!(entry_time, :tooltip_text, "Time position [s]")
    bt_start = GtkButton("⇤")
    set_gtk_property!(bt_start, :tooltip_text, "Go to the signal beginning")
    bt_prev5 = GtkButton("↞")
    set_gtk_property!(bt_prev5, :tooltip_text, "Go back by $zoom seconds")
    bt_next5 = GtkButton("↠")
    set_gtk_property!(bt_next5, :tooltip_text, "Go forward by $zoom seconds")
    bt_end = GtkButton("⇥")
    set_gtk_property!(bt_end, :tooltip_text, "Go to the signal end")
    bt_help = GtkButton("🛈")
    set_gtk_property!(bt_help, :tooltip_text, "Show keyboard shortcuts")
    bt_close = GtkButton("✖")
    set_gtk_property!(bt_close, :tooltip_text, "Close this window")
    if length(c_idx) > 20
        ch_slider = GtkScale(false, c_idx[c_idx_first]:(c_idx[end] - 19))
        set_gtk_property!(ch_slider, :draw_value, false)
    else
        ch_slider = GtkScale(false, c_idx[c_idx_first]:c_idx[end])
        set_gtk_property!(ch_slider, :draw_value, false)
        set_gtk_property!(ch_slider, :sensitive, false)
    end
    set_gtk_property!(ch_slider, :tooltip_text, "Scroll components")
    set_gtk_property!(ch_slider, :vexpand, true)
    oc = GtkOrientable(ch_slider)
    set_gtk_property!(oc, :orientation, 1)
    g[1:7, 1] = can
    g[8, 1] = ch_slider
    g[1, 2] = bt_start
    g[2, 2] = bt_prev5
    g[3, 2] = entry_time
    g[4, 2] = bt_next5
    g[5, 2] = bt_end
    g[6, 2] = bt_help
    g[7, 2] = bt_close
    push!(win, g)

    showall(win)

    @guarded draw(can) do widget
        time1 = get_gtk_property(entry_time, :value, Float64)
        time2 = time1 + zoom
        time2 > obj.time_pts[end] && (time2 = obj.time_pts[end])
        ctx = getgc(can)
        show(io, MIME("image/png"), NeuroAnalyzer.plot(obj,
                                                       c,
                                                       c_idx=c_idx[c_idx_first]:c_idx[c_idx_last],
                                                       seg=(time1, time2),
                                                       mono=mono,
                                                       scale=scale,
                                                       title=""))
        img = read_from_png(io)
        set_source_surface(ctx, img, 0, 0)
        paint(ctx)
    end

    can.mouse.scroll = @guarded (widget, event) -> begin
        s = event.state
        if event.direction == 1 # down
            if s == 0x00000011
                time_current = get_gtk_property(entry_time, :value, Float64)
                if time_current < obj.time_pts[end] - zoom
                    time_current += 1
                else
                    time_current = obj.time_pts[end] - zoom
                end
                Gtk.@sigatom begin
                    set_gtk_property!(entry_time, :value, time_current)
                end
            elseif s == 0x00000014
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
                if c_idx_last < c_idx[end]
                    c_idx_first += 1
                    c_idx_last += 1
                    Gtk.@sigatom begin
                        GAccessor.value(ch_slider, c_idx_first)
                    end
                    draw(can)
                end
            end
        elseif event.direction == 0 # up
            if s == 0x00000011
                time_current = get_gtk_property(entry_time, :value, Float64)
                if time_current >= obj.time_pts[1] + 1
                    time_current -= 1
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_time, :value, time_current)
                    end
                end
            elseif s == 0x00000014
                time_current = get_gtk_property(entry_time, :value, Float64)
                if time_current >= obj.time_pts[1] + zoom
                    time_current = time_current - zoom
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_time, :value, time_current)
                    end
                end
            else
                if c_idx_first > 1
                    c_idx_first -= 1
                    c_idx_last -= 1
                    Gtk.@sigatom begin
                        GAccessor.value(ch_slider, c_idx_first)
                    end
                    draw(can)
                end
            end
        end
    end

    signal_connect(entry_time, "value-changed") do widget
        c_idx_first = round(Int64, GAccessor.value(c_idx_slider))
        c_idx_last = c_idx_first + 19
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

    help = "Keyboard shortcuts:\n\nPage Up\t\tscroll channels up\nPage Down\tscroll channels down\n\nHome\tgo to the signal beginning\nEnd\tgo to the signal end\nctrl-,\tgo back by 1 second\nctrl-.\tgo forward by 1 second\nalt-,\tgo back by $zoom seconds\nalt-.\tgo forward by $zoom seconds\n\n[\t zoom in\n]\tzoom out\n\nctrl-s\ttoggle scales\nctrl-m\ttoggle monochromatic mode\n\nctrl-h\tthis info\nctrl-q\texit\n"

    signal_connect(bt_help, "clicked") do widgete
        info_dialog()
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        s = event.state
        if k == 0x0000ff55 # Page Up
            if c_idx_first > 1
                c_idx_first -= 1
                c_idx_last -= 1
                Gtk.@sigatom begin
                    GAccessor.value(ch_slider, c_idx_first)
                end
                draw(can)
            end
        elseif k == 0x0000ff56 # Page Down
            if c_idx_last < ch[end]
                c_idx_first += 1
                c_idx_last += 1
                Gtk.@sigatom begin
                    GAccessor.value(ch_slider, c_idx_first)
                end
                draw(can)
            end
        elseif k == 0x0000005b # [
            if zoom > 1
                zoom -= 1
                set_gtk_property!(bt_next5, :tooltip_text, "Go forward by $zoom seconds")
                set_gtk_property!(bt_prev5, :tooltip_text, "Go back by $zoom seconds")
                draw(can)
            end
            help = "Keyboard shortcuts:\n\nPage Up\t\tscroll channels up\nPage Down\tscroll channels down\n\nHome\tgo to the signal beginning\nEnd\tgo to the signal end\nctrl-,\tgo back by 1 second\nctrl-.\tgo forward by 1 second\nalt-,\tgo back by $zoom seconds\nalt-.\tgo forward by $zoom seconds\n\n[\t zoom in\n]\tzoom out\n\nctrl-s\ttoggle scales\n\nctrl-h\tthis info\nctrl-q\texit\n"
        elseif k == 0x0000005d # ]
            if zoom < 30 && zoom < obj.time_pts[end] - 1
                zoom += 1
                set_gtk_property!(bt_next5, :tooltip_text, "Go forward by $zoom seconds")
                set_gtk_property!(bt_prev5, :tooltip_text, "Go back by $zoom seconds")
                draw(can)
            else
                zoom = obj.time_pts[end]
                set_gtk_property!(bt_next5, :tooltip_text, "Go forward by $zoom seconds")
                set_gtk_property!(bt_prev5, :tooltip_text, "Go back by $zoom seconds")
                draw(can)
            end
            help = "Keyboard shortcuts:\n\nPage Up\t\tscroll channels up\nPage Down\tscroll channels down\n\nHome\tgo to the signal beginning\nEnd\tgo to the signal end\nctrl-,\tgo back by 1 second\nctrl-.\tgo forward by 1 second\nalt-,\tgo back by $zoom seconds\nalt-.\tgo forward by $zoom seconds\n\n[\t zoom in\n]\tzoom out\n\nctrl-s\ttoggle scales\n\nctrl-h\tthis info\nctrl-q\texit\n"
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
        if s == 0x00000018 # alt
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
        if s == 0x00000014 # ctrl
            if k == 113 # q
                Gtk.destroy(win)
            elseif k == 104 # h
                info_dialog(help)
            elseif k == 0x0000006d # m
                mono = !mono
                draw(can)
            elseif k == 0x00000073 # s
                scale = !scale
                draw(can)
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
- `c_idx::Union{Vector{Int64}, AbstractRange}=_c(nchannels(obj))`: component channel to display, default is all component channels
- `ep::Int64=1`: initial epoch to display
"""
function iview_ep(obj::NeuroAnalyzer.NEURO, c::Union{Symbol, AbstractArray}; c_idx::Union{Vector{Int64}, AbstractRange}=_c(nchannels(obj)), ep::Int64=1)

    @assert length(c_idx) > 1 "iview() requires > 1 channels."
    @assert nepochs(obj) > 1 "iview_cont() should be used for continuous object."

    _check_epochs(obj, ep)

    mono = false
    scale = true

    if length(c_idx) > 20
        c_idx_first = 1
        c_idx_last = 20
    else
        c_idx_first = 1
        c_idx_last = length(c_idx)
    end

    p = NeuroAnalyzer.plot(obj, c, c_idx=c_idx[c_idx_first:c_idx_last], ep=ep, mono=mono, scale=scale, title="")
    win = GtkWindow("NeuroAnalyzer: iview_ep()", Int32(p.attr[:size][1]) + 40, Int32(p.attr[:size][2]) + 40)
    set_gtk_property!(win, :border_width, 20)
    set_gtk_property!(win, :resizable, true)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    set_gtk_property!(win, :startup_id, "org.neuroanalyzer")
    can = GtkCanvas(Int32(p.attr[:size][1]), Int32(p.attr[:size][2]))
    g = GtkGrid()
    set_gtk_property!(g, :column_homogeneous, false)
    set_gtk_property!(g, :column_spacing, 10)
    set_gtk_property!(g, :row_spacing, 10)
    entry_epoch = GtkSpinButton(1, nepochs(obj), 1)
    set_gtk_property!(entry_epoch, :value, ep)
    set_gtk_property!(entry_epoch, :tooltip_text, "Epoch")
    bt_start = GtkButton("⇤")
    set_gtk_property!(bt_start, :tooltip_text, "Go to the signal beginning")
    bt_end = GtkButton("⇥")
    set_gtk_property!(bt_end, :tooltip_text, "Go to the signal end")
    bt_help = GtkButton("🛈")
    set_gtk_property!(bt_help, :tooltip_text, "Show keyboard shortcuts")
    bt_close = GtkButton("✖")
    set_gtk_property!(bt_close, :tooltip_text, "Close this window")
    if length(c_idx) > 20
        ch_slider = GtkScale(false, c_idx[c_idx_first]:(c_idx[end] - 19))
        set_gtk_property!(ch_slider, :draw_value, false)
    else
        ch_slider = GtkScale(false, c_idx[c_idx_first]:c_idx[end])
        set_gtk_property!(ch_slider, :draw_value, false)
        set_gtk_property!(ch_slider, :sensitive, false)
    end
    set_gtk_property!(ch_slider, :tooltip_text, "Scroll components")
    set_gtk_property!(ch_slider, :vexpand, true)
    oc = GtkOrientable(ch_slider)
    set_gtk_property!(oc, :orientation, 1)
    g[1:5, 1] = can
    g[6, 1] = ch_slider
    g[1, 2] = bt_start
    g[2, 2] = entry_epoch
    g[3, 2] = bt_end
    g[4, 2] = bt_help
    g[5, 2] = bt_close
    push!(win, g)

    showall(win)

    @guarded draw(can) do widget
        ep = get_gtk_property(entry_epoch, :value, Int64)
        ctx = getgc(can)
        show(io, MIME("image/png"), NeuroAnalyzer.plot(obj,
                                                       c,
                                                       c_idx=c_idx[c_idx_first]:c_idx[c_idx_last],
                                                       ep=ep,
                                                       mono=mono,
                                                       scale=scale,
                                                       title=""))
        img = read_from_png(io)
        set_source_surface(ctx, img, 0, 0)
        paint(ctx)
    end

    signal_connect(ch_slider, "value-changed") do widget, others...
        c_idx_first = round(Int64, GAccessor.value(ch_slider))
        c_idx_last = c_idx_first + 19
        draw(can)
    end

    can.mouse.scroll = @guarded (widget, event) -> begin
        s = event.state
        if event.direction == 1 # down
            if s == 0x00000011
                ep = get_gtk_property(entry_epoch, :value, Int64)
                if ep < nepochs(obj)
                    ep += 1
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_epoch, :value, ep)
                    end
                end
            else
                if c_idx_last < c_idx[end]
                    c_idx_first += 1
                    c_idx_last += 1
                    Gtk.@sigatom begin
                        GAccessor.value(ch_slider, c_idx_first)
                    end
                    draw(can)
                end
            end
        elseif event.direction == 0 # up
            if s == 0x00000011
                ep = get_gtk_property(entry_epoch, :value, Int64)
                if ep > 1
                    ep -= 1
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_epoch, :value, ep)
                    end
                end
            else
                if c_idx_first > 1
                    c_idx_first -= 1
                    c_idx_last -= 1
                    Gtk.@sigatom begin
                        GAccessor.value(ch_slider, c_idx_first)
                    end
                    draw(can)
                end
            end
        end
    end

    signal_connect(entry_epoch, "value-changed") do widget
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

    help = "Keyboard shortcuts:\n\nPage Up\t\tscroll channels up\nPage Down\tscroll channels down\n\nHome\tgo to first epoch\nEnd\tgo to last epoch\nctrl-,\tprevious epoch\nctrl-.\tnext epoch\n\nctrl-s\ttoggle scales\nctrl-m\ttoggle monochromatic mode\n\nctrl-h\tthis info\nctrl-q\texit\n"

    signal_connect(bt_help, "clicked") do widgete
        info_dialog(help)
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        s = event.state
        if k == 0x0000ff55 # Page Up
            if c_idx_first > 1
                c_idx_first -= 1
                ch_last -= 1
                Gtk.@sigatom begin
                    GAccessor.value(ch_slider, c_idx_first)
                end
                draw(can)
            end
        elseif k == 0x0000ff56 # Page Down
            if ch_last < c_idx[end]
                c_idx_first += 1
                ch_last += 1
                Gtk.@sigatom begin
                    GAccessor.value(ch_slider, c_idx_first)
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
        if s == 0x00000014
            if k == 113 # q
                Gtk.destroy(win)
            elseif k == 104 # h
                info_dialog(help)
            elseif k == 0x0000006d # m
                mono = !mono
                draw(can)
            elseif k == 0x00000073 # s
                scale = !scale
                draw(can)
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
"""
function iview(p::Plots.Plot{Plots.GRBackend})

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
        # Cairo.scale(ctx, 242/w, 241/h)
        Cairo.set_source_surface(ctx, img, 1, 2)
        Cairo.paint(ctx)
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        s = event.state
        if s == 0x00000014
            if k == 115 # s
                file_name = save_dialog("Pick image file", GtkNullContainer(), (GtkFileFilter("*.png", name="All supported formats"), "*.png"))
                    if file_name != ""
                        if splitext(file_name)[2] in [".png"]
                            surface_buf = Gtk.cairo_surface(can)
                            Cairo.write_to_png(surface_buf, file_name)
                        else
                            warn_dialog("Incorrect filename!")
                        end
                    end
            elseif k == 113 # q
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
"""
function iview(file_name::String)

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
        if s == 0x00000014
            if k == 115 # s
                file_name = save_dialog("Pick image file", GtkNullContainer(), (GtkFileFilter("*.png", name="All supported formats"), "*.png"))
                    if file_name != ""
                        if splitext(file_name)[2] in [".png"]
                            surface_buf = Gtk.cairo_surface(can)
                            Cairo.write_to_png(surface_buf, file_name)
                        else
                            warn_dialog("Incorrect filename!")
                        end
                    end
            elseif k == 111 # o
                info_dialog("This feature has not been implemented yet.")
            elseif k == 113 # q
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
"""
function iview(c::Cairo.CairoSurfaceBase{UInt32})

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
        if s == 0x00000014
            if k == 115 # s
                file_name = save_dialog("Pick image file", GtkNullContainer(), (GtkFileFilter("*.png", name="All supported formats"), "*.png"))
                    if file_name != ""
                        if splitext(file_name)[2] in [".png"]
                            surface_buf = Gtk.cairo_surface(can)
                            Cairo.write_to_png(surface_buf, file_name)
                        else
                            warn_dialog("Incorrect filename!")
                        end
                    end
            elseif k == 113 # q
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