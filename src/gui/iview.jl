export iview
export iview_ep
export iview_cont

"""
    iview(obj; ch, zoom)

Interactive view of continuous or epoched signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: channel(s) to plot, default is all channels
- `ep::Int64=1`: initial epoch to display
- `zoom::Real=5`: how many seconds are displayed in one segment
"""
function iview(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj)), ep::Int64=1, zoom::Real=5)

    if nepochs(obj) == 1
        iview_cont(obj, ch=ch, zoom=zoom)
    else
        iview_ep(obj, ch=ch, ep=ep)
    end

    return nothing

end

"""
    iview_cont(obj; ch, zoom)

Interactive view of continuous signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: channel(s) to plot, default is all channels
- `zoom::Real=5`: how many seconds are displayed in one segment
"""
function iview_cont(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj)), zoom::Real=5)

    (signal_len(obj) / sr(obj)) < zoom && (zoom = obj.time_pts[end])

    @assert zoom > 0 "zoom must be > 0."
    @assert zoom <= signal_len(obj) / sr(obj) "zoom must be ≤ $(signal_len(obj) / sr(obj))."
    @assert nepochs(obj) == 1 "iview_ep() should be used for epoched object."
    _check_channels(obj, ch)

    if length(ch) > 10
        ch_first = 1
        ch_last = 10
    else
        ch_first = 1
        ch_last = length(ch)
    end

    if length(ch) > 1
        p = NeuroAnalyzer.plot(obj, ch=ch[ch_first:ch_last], mono=true, title="")
    else
        p = NeuroAnalyzer.plot(obj, ch=ch, mono=true, title="")
    end
    
    win = GtkWindow("NeuroAnalyzer: iview_cont()", Int32(p.attr[:size][1]), Int32(p.attr[:size][2]) + 40)
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
    bt_chup = GtkButton("△")
    set_gtk_property!(bt_chup, :tooltip_text, "Slide channels up")
    length(ch) < 11 && set_gtk_property!(bt_chup, :sensitive, false)
    ch[ch_first] == ch[1] && set_gtk_property!(bt_chup, :sensitive, false)
    bt_chdown = GtkButton("▽")
    set_gtk_property!(bt_chdown, :tooltip_text, "Slide channels down")
    length(ch) < 11 && set_gtk_property!(bt_chdown, :sensitive, false)
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
    lab_ch = GtkLabel("$(lpad(string(ch[ch_first]), 2, '0')):$(lpad(string(ch[ch_last]), 2, '0'))")
    set_gtk_property!(lab_ch, :halign, 0)
    g[1:14, 1] = can
    g[1, 2] = bt_chup
    g[2, 2] = lab_ch
    g[3, 2] = bt_chdown
    g[4, 2] = GtkLabel("")
    g[5, 2] = bt_start
    g[6, 2] = bt_prev5
    g[7, 2] = bt_prev
    g[8, 2] = entry_time
    g[9, 2] = bt_next
    g[10, 2] = bt_next5
    g[11, 2] = bt_end
    g[12, 2] = GtkLabel("")
    g[13, 2] = bt_help
    g[14, 2] = bt_close
    push!(win, g)

    showall(win)

    @guarded draw(can) do widget
        time1 = get_gtk_property(entry_time, :value, Float64)
        time2 = time1 + zoom
        time2 > obj.time_pts[end] && (time2 = obj.time_pts[end])
        ctx = getgc(can)

        if length(ch) > 1
            show(io, MIME("image/png"), NeuroAnalyzer.plot(obj,
                                                           ch=ch[ch_first]:ch[ch_last],
                                                           seg=(time1, time2),
                                                           mono=true,
                                                           title=""))
        else
            show(io, MIME("image/png"), NeuroAnalyzer.plot(obj,
                                                           ch=ch,
                                                           seg=(time1, time2),
                                                           mono=true,
                                                           title=""))
        end

        img = read_from_png(io)
        set_source_surface(ctx, img, 0, 0)
        paint(ctx)
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

    can.mouse.scroll = @guarded (widget, event) -> begin
        if event.direction == 1 # down
            if ch_last < ch[end]
                ch_first += 1
                ch_last += 1
                Gtk.@sigatom begin
                    set_gtk_property!(lab_ch, :label, "$(lpad(string(ch[ch_first]), 2, '0')):$(lpad(string(ch[ch_last]), 2, '0'))")
                    ch[ch_first] != ch[1] && set_gtk_property!(bt_chup, :sensitive, true)
                    ch[ch_last] == ch[end] && set_gtk_property!(bt_chdown, :sensitive, false)
                end
                draw(can)
            end
        elseif event.direction == 0 # up
            if ch_first > 1
                ch_first -= 1
                ch_last -= 1
                Gtk.@sigatom begin
                    set_gtk_property!(lab_ch, :label, "$(lpad(string(ch[ch_first]), 2, '0')):$(lpad(string(ch[ch_last]), 2, '0'))")
                    ch[ch_first] == ch[1] && set_gtk_property!(bt_chup, :sensitive, false)
                    ch[ch_last] != ch[end] && set_gtk_property!(bt_chdown, :sensitive, true)
                end
                draw(can)
            end
        end
    end

    signal_connect(bt_chdown, "clicked") do widget
        if ch_last < ch[end]
            ch_first += 1
            ch_last += 1
            Gtk.@sigatom begin
                set_gtk_property!(lab_ch, :label, "$(lpad(string(ch[ch_first]), 2, '0')):$(lpad(string(ch[ch_last]), 2, '0'))")
                ch[ch_first] != ch[1] && set_gtk_property!(bt_chup, :sensitive, true)
                ch[ch_last] == ch[end] && set_gtk_property!(bt_chdown, :sensitive, false)
            end
            draw(can)
        end
    end

    signal_connect(bt_chup, "clicked") do widget
        if ch_first > 1
            ch_first -= 1
            ch_last -= 1
            Gtk.@sigatom begin
                set_gtk_property!(lab_ch, :label, "$(lpad(string(ch[ch_first]), 2, '0')):$(lpad(string(ch[ch_last]), 2, '0'))")
                ch[ch_first] == ch[1] && set_gtk_property!(bt_chup, :sensitive, false)
                ch[ch_last] != ch[end] && set_gtk_property!(bt_chdown, :sensitive, true)
            end
            draw(can)
        end
    end

    signal_connect(entry_time, "value-changed") do widget
        draw(can)
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

    signal_connect(bt_close, "clicked") do widget
        Gtk.destroy(win)
    end

    signal_connect(bt_help, "clicked") do widgete
        info_dialog("Keyboard shortcuts:\n\nctrl-b\tslide channels up\nctrl-n\tslide channels down\n\nctrl-a\tgo to the signal beginning\nctrl-s\tgo to the signal end\nctrl-z\tgo back by 1 second\nctrl-x\tgo forward by 1 second\nctrl-c\tgo back by $zoom seconds\nctrl-v\tgo forward by $zoom seconds\n\nctrl-h\tthis info\nctrl-q\texit\n")
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        s = event.state
        if s == 4
            if k == 113 # q
                Gtk.destroy(win)
            elseif k == 104 # h
                info_dialog("Keyboard shortcuts:\n\nctrl-b\tslide channels up\nctrl-n\tslide channels down\n\nctrl-a\tgo to the signal beginning\nctrl-s\tgo to the signal end\nctrl-z\tgo back by 1 second\nctrl-x\tgo forward by 1 second\nctrl-c\tgo back by $zoom seconds\nctrl-v\tgo forward by $zoom seconds\n\nctrl-h\tthis info\nctrl-q\texit\n")
            elseif k == 98 # b
                if ch_first > 1
                    ch_first -= 1
                    ch_last -= 1
                    Gtk.@sigatom begin
                        set_gtk_property!(lab_ch, :label, "$(lpad(string(ch[ch_first]), 2, '0')):$(lpad(string(ch[ch_last]), 2, '0'))")
                        ch[ch_first] == ch[1] && set_gtk_property!(bt_chup, :sensitive, false)
                        ch[ch_last] != ch[end] && set_gtk_property!(bt_chdown, :sensitive, true)
                    end
                    draw(can)
                end
            elseif k == 110 # n
                if ch_last < ch[end]
                    ch_first += 1
                    ch_last += 1
                    Gtk.@sigatom begin
                        set_gtk_property!(lab_ch, :label, "$(lpad(string(ch[ch_first]), 2, '0')):$(lpad(string(ch[ch_last]), 2, '0'))")
                        ch[ch_first] != ch[1] && set_gtk_property!(bt_chup, :sensitive, true)
                        ch[ch_last] == ch[end] && set_gtk_property!(bt_chdown, :sensitive, false)
                    end
                    draw(can)
                end
            elseif k == 97 # a
                Gtk.@sigatom begin
                    set_gtk_property!(entry_time, :value, obj.time_pts[1])
                end
                draw(can)
            elseif k == 115 # s
                time_current = obj.time_pts[end] - zoom
                Gtk.@sigatom begin
                    set_gtk_property!(entry_time, :value, time_current)
                end
                draw(can)
            elseif k == 122 # z
                time_current = get_gtk_property(entry_time, :value, Float64)
                if time_current >= obj.time_pts[1] + 1
                    time_current -= 1
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_time, :value, time_current)
                    end
                end
                draw(can)
            elseif k == 99 # c
                time_current = get_gtk_property(entry_time, :value, Float64)
                if time_current >= obj.time_pts[1] + zoom
                    time_current = time_current - zoom
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_time, :value, time_current)
                    end
                end
                draw(can)
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

    cnd = Condition()
    signal_connect(win, :destroy) do widget
        notify(cnd)
    end
    @async Gtk.gtk_main()
    wait(cnd)

    return nothing

end

"""
    iview_ep(obj, ch, ep)

Interactive view of epoched signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: channel(s) to plot, default is all channels
- `ep::Int64=1`: initial epoch to display
"""
function iview_ep(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj)), ep::Int64=1)

    @assert nepochs(obj) > 1 "iview_cont() should be used for continuous object."
    _check_channels(obj, ch)
    _check_epochs(obj, ep)

    if length(ch) > 10
        ch_first = 1
        ch_last = 10
    else
        ch_first = 1
        ch_last = length(ch)
    end

    if length(ch) > 1
        p = NeuroAnalyzer.plot(obj, ch=ch[ch_first:ch_last], ep=ep, mono=true, title="")
    else
        p = NeuroAnalyzer.plot(obj, ch=ch, ep=ep, mono=true, title="")
    end
    win = GtkWindow("NeuroAnalyzer: iview_ep()", Int32(p.attr[:size][1]), Int32(p.attr[:size][2]) + 40)
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
    bt_chup = GtkButton("△")
    set_gtk_property!(bt_chup, :tooltip_text, "Slide channels up")
    length(ch) < 11 && set_gtk_property!(bt_chup, :sensitive, false)
    ch[ch_first] == ch[1] && set_gtk_property!(bt_chup, :sensitive, false)
    bt_chdown = GtkButton("▽")
    set_gtk_property!(bt_chdown, :tooltip_text, "Slide channels down")
    length(ch) < 11 && set_gtk_property!(bt_chdown, :sensitive, false)
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
    lab_ch = GtkLabel("$(lpad(string(ch[ch_first]), 2, '0')):$(lpad(string(ch[ch_last]), 2, '0'))")
    set_gtk_property!(lab_ch, :halign, 0)
    g[1:14, 1] = can
    g[1, 2] = bt_chup
    g[2, 2] = lab_ch
    g[3, 2] = bt_chdown
    g[4, 2] = GtkLabel("")
    g[5, 2] = bt_start
    g[6, 2] = bt_prev
    g[7, 2] = entry_epoch
    g[8, 2] = bt_next
    g[9, 2] = bt_end
    g[10, 2] = GtkLabel("")
    g[11, 2] = bt_help
    g[12, 2] = bt_close
    push!(win, g)

    showall(win)

    @guarded draw(can) do widget
        ep = get_gtk_property(entry_epoch, :value, Int64)
        ctx = getgc(can)

        if length(ch) > 1
            show(io, MIME("image/png"), NeuroAnalyzer.plot(obj,
                                                           ch=ch[ch_first]:ch[ch_last],
                                                           ep=ep,
                                                           mono=true,
                                                           title=""))
        else
            show(io, MIME("image/png"), NeuroAnalyzer.plot(obj,
                                                           ch=ch,
                                                           ep=ep,
                                                           mono=true,
                                                           title=""))
        end

        img = read_from_png(io)
        set_source_surface(ctx, img, 0, 0)
        paint(ctx)
    end

    can.mouse.scroll = @guarded (widget, event) -> begin
        if event.direction == 1 # down
            if ch_last < ch[end]
                ch_first += 1
                ch_last += 1
                Gtk.@sigatom begin
                    set_gtk_property!(lab_ch, :label, "$(lpad(string(ch[ch_first]), 2, '0')):$(lpad(string(ch[ch_last]), 2, '0'))")
                    ch[ch_first] != ch[1] && set_gtk_property!(bt_chup, :sensitive, true)
                    ch[ch_last] == ch[end] && set_gtk_property!(bt_chdown, :sensitive, false)
                end
                draw(can)
            end
        elseif event.direction == 0 # up
            if ch_first > 1
                ch_first -= 1
                ch_last -= 1
                Gtk.@sigatom begin
                    set_gtk_property!(lab_ch, :label, "$(lpad(string(ch[ch_first]), 2, '0')):$(lpad(string(ch[ch_last]), 2, '0'))")
                    ch[ch_first] == ch[1] && set_gtk_property!(bt_chup, :sensitive, false)
                    ch[ch_last] != ch[end] && set_gtk_property!(bt_chdown, :sensitive, true)
                end
                draw(can)
            end
        end
    end

    signal_connect(bt_chdown, "clicked") do widget
        if ch_last < ch[end]
            ch_first += 1
            ch_last += 1
            Gtk.@sigatom begin
                set_gtk_property!(lab_ch, :label, "$(lpad(string(ch[ch_first]), 2, '0')):$(lpad(string(ch[ch_last]), 2, '0'))")
                ch[ch_first] != ch[1] && set_gtk_property!(bt_chup, :sensitive, true)
                ch[ch_last] == ch[end] && set_gtk_property!(bt_chdown, :sensitive, false)
            end
            draw(can)
        end
    end

    signal_connect(bt_chup, "clicked") do widget
        if ch_first > 1
            ch_first -= 1
            ch_last -= 1
            Gtk.@sigatom begin
                set_gtk_property!(lab_ch, :label, "$(lpad(string(ch[ch_first]), 2, '0')):$(lpad(string(ch[ch_last]), 2, '0'))")
                ch[ch_first] == ch[1] && set_gtk_property!(bt_chup, :sensitive, false)
                ch[ch_last] != ch[end] && set_gtk_property!(bt_chdown, :sensitive, true)
            end
            draw(can)
        end
    end

    signal_connect(entry_epoch, "value-changed") do widget
        draw(can)
    end

    signal_connect(bt_prev, "clicked") do widget
        ep = get_gtk_property(entry_epoch, :value, Int64)
        if ep > 1
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
        info_dialog("Keyboard shortcuts:\n\nctrl-b\tslide channels up\nctrl-n\tslide channels down\n\nctrl-a\tgo to first epoch\nctrl-s\tgo to last epoch\nctrl-z\tprevious epoch\nctrl-x\tnext epoch\n\nctrl-h\tthis info\nctrl-q\texit\n")
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        s = event.state
        if s == 4
            if k == 113 # q
                Gtk.destroy(win)
            elseif k == 104 # h
                info_dialog("Keyboard shortcuts:\n\nctrl-b\tslide channels up\nctrl-n\tslide channels down\n\nctrl-a\tgo to first epoch\nctrl-s\tgo to last epoch\nctrl-z\tprevious epoch\nctrl-x\tnext epoch\n\nctrl-h\tthis info\nctrl-q\texit\n")
            elseif k == 98 # b
                if ch_first > 1
                    ch_first -= 1
                    ch_last -= 1
                    Gtk.@sigatom begin
                        set_gtk_property!(lab_ch, :label, "$(lpad(string(ch[ch_first]), 2, '0')):$(lpad(string(ch[ch_last]), 2, '0'))")
                        ch[ch_first] == ch[1] && set_gtk_property!(bt_chup, :sensitive, false)
                        ch[ch_last] != ch[end] && set_gtk_property!(bt_chdown, :sensitive, true)
                    end
                    draw(can)
                end
            elseif k == 110 # n
                if ch_last < ch[end]
                    ch_first += 1
                    ch_last += 1
                    Gtk.@sigatom begin
                        set_gtk_property!(lab_ch, :label, "$(lpad(string(ch[ch_first]), 2, '0')):$(lpad(string(ch[ch_last]), 2, '0'))")
                        ch[ch_first] != ch[1] && set_gtk_property!(bt_chup, :sensitive, true)
                        ch[ch_last] == ch[end] && set_gtk_property!(bt_chdown, :sensitive, false)
                    end
                    draw(can)
                end
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
                if ep > 1
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

    cnd = Condition()
    signal_connect(win, :destroy) do widget
        notify(cnd)
    end
    @async Gtk.gtk_main()
    wait(cnd)

    return nothing

end

"""
    iview(obj1, obj2; ch, zoom)

Interactive view of continuous or epoched signal.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `obj2::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj1))`: channel(s) to plot, default is all channels
- `ep::Int64=1`: initial epoch to display
- `zoom::Real=5`: how many seconds are displayed in one segment
"""
function iview(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj1)), ep::Int64=1, zoom::Real=5)

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
    iview_cont(obj1, obj2; ch, zoom)

Interactive view of continuous signal.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `obj2::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj1))`: channel(s) to plot, default is all channels
- `zoom::Real=5`: how many seconds are displayed in one segment
"""
function iview_cont(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj1)), zoom::Real=5)

    (signal_len(obj1) / sr(obj1)) < zoom && (zoom = obj1.time_pts[end])

    @assert size(obj1) == size(obj2) "Both signals must have the same size."
    @assert sr(obj1) == sr(obj2) "Both signals must have the same sampling rate."
    @assert zoom > 0 "zoom must be > 0."
    @assert zoom <= signal_len(obj1) / sr(obj1) "zoom must be ≤ $(signal_len(obj1) / sr(obj1))."
    @assert nepochs(obj1) == 1 "iview_ep() should be used for epoched object."
    _check_channels(obj1, ch)

    if length(ch) > 10
        ch_first = 1
        ch_last = 10
    else
        ch_first = 1
        ch_last = length(ch)
    end

    if length(ch) > 1
        p = NeuroAnalyzer.plot(obj1, obj2, ch=ch[ch_first]:ch[ch_last], title="")
    else
        p = NeuroAnalyzer.plot(obj1, obj2, ch=ch, title="")
    end

    win = GtkWindow("NeuroAnalyzer: iview_cont()", Int32(p.attr[:size][1]), Int32(p.attr[:size][2]) + 40)
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
    entry_time = GtkSpinButton(obj1.time_pts[1], obj1.time_pts[end] - zoom, zoom)
    set_gtk_property!(entry_time, :digits, 2)
    set_gtk_property!(entry_time, :value, obj1.time_pts[1])
    set_gtk_property!(entry_time, :tooltip_text, "Time position [s]") 
    bt_chup = GtkButton("△")
    set_gtk_property!(bt_chup, :tooltip_text, "Slide channels up")
    length(ch) < 11 && set_gtk_property!(bt_chup, :sensitive, false)
    ch[ch_first] == ch[1] && set_gtk_property!(bt_chup, :sensitive, false)
    bt_chdown = GtkButton("▽")
    set_gtk_property!(bt_chdown, :tooltip_text, "Slide channels down")
    length(ch) < 11 && set_gtk_property!(bt_chdown, :sensitive, false)
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
    lab_ch = GtkLabel("$(lpad(string(ch[ch_first]), 2, '0')):$(lpad(string(ch[ch_last]), 2, '0'))")
    set_gtk_property!(lab_ch, :halign, 0)
    g[1:14, 1] = can
    g[1, 2] = bt_chup
    g[2, 2] = lab_ch
    g[3, 2] = bt_chdown
    g[4, 2] = GtkLabel("")
    g[5, 2] = bt_start
    g[6, 2] = bt_prev5
    g[7, 2] = bt_prev
    g[8, 2] = entry_time
    g[9, 2] = bt_next
    g[10, 2] = bt_next5
    g[11, 2] = bt_end
    g[12, 2] = GtkLabel("")
    g[13, 2] = bt_help
    g[14, 2] = bt_close
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
                                                           title=""))
        else
            show(io, MIME("image/png"), NeuroAnalyzer.plot(obj1, obj2,
                                                           ch=ch,
                                                           seg=(time1, time2),
                                                           title=""))
        end

        img = read_from_png(io)
        set_source_surface(ctx, img, 0, 0)
        paint(ctx)
    end

    can.mouse.scroll = @guarded (widget, event) -> begin
        if event.direction == 1 # down
            if ch_last < ch[end]
                ch_first += 1
                ch_last += 1
                Gtk.@sigatom begin
                    set_gtk_property!(lab_ch, :label, "$(lpad(string(ch[ch_first]), 2, '0')):$(lpad(string(ch[ch_last]), 2, '0'))")
                    ch[ch_first] != ch[1] && set_gtk_property!(bt_chup, :sensitive, true)
                    ch[ch_last] == ch[end] && set_gtk_property!(bt_chdown, :sensitive, false)
                end
                draw(can)
            end
        elseif event.direction == 0 # up
            if ch_first > 1
                ch_first -= 1
                ch_last -= 1
                Gtk.@sigatom begin
                    set_gtk_property!(lab_ch, :label, "$(lpad(string(ch[ch_first]), 2, '0')):$(lpad(string(ch[ch_last]), 2, '0'))")
                    ch[ch_first] == ch[1] && set_gtk_property!(bt_chup, :sensitive, false)
                    ch[ch_last] != ch[end] && set_gtk_property!(bt_chdown, :sensitive, true)
                end
                draw(can)
            end
        end
    end

    signal_connect(bt_chdown, "clicked") do widget
        if ch_last < ch[end]
            ch_first += 1
            ch_last += 1
            Gtk.@sigatom begin
                set_gtk_property!(lab_ch, :label, "$(lpad(string(ch[ch_first]), 2, '0')):$(lpad(string(ch[ch_last]), 2, '0'))")
                ch[ch_first] != ch[1] && set_gtk_property!(bt_chup, :sensitive, true)
                ch[ch_last] == ch[end] && set_gtk_property!(bt_chdown, :sensitive, false)
            end
            draw(can)
        end
    end

    signal_connect(bt_chup, "clicked") do widget
        if ch_first > 1
            ch_first -= 1
            ch_last -= 1
            Gtk.@sigatom begin
                set_gtk_property!(lab_ch, :label, "$(lpad(string(ch[ch_first]), 2, '0')):$(lpad(string(ch[ch_last]), 2, '0'))")
                ch[ch_first] == ch[1] && set_gtk_property!(bt_chup, :sensitive, false)
                ch[ch_last] != ch[end] && set_gtk_property!(bt_chdown, :sensitive, true)
            end
            draw(can)
        end
    end

    signal_connect(entry_time, "value-changed") do widget
        draw(can)
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
        else
            time_current = obj1.time_pts[end] - zoom
        end
        Gtk.@sigatom begin
            set_gtk_property!(entry_time, :value, time_current)
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

    signal_connect(bt_help, "clicked") do widgete
        info_dialog("Keyboard shortcuts:\n\nctrl-b\tslide channels up\nctrl-n\tslide channels down\n\nctrl-a\tgo to the signal beginning\nctrl-s\tgo to the signal end\nctrl-z\tgo back by 1 second\nctrl-x\tgo forward by 1 second\nctrl-c\tgo back by $zoom seconds\nctrl-v\tgo forward by $zoom seconds\n\nctrl-h\tthis info\nctrl-q\texit\n")
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        s = event.state
        if s == 4
            if k == 113 # q
                Gtk.destroy(win)
            elseif k == 104 # h
                info_dialog("Keyboard shortcuts:\n\nctrl-b\tslide channels up\nctrl-n\tslide channels down\n\nctrl-a\tgo to the signal beginning\nctrl-s\tgo to the signal end\nctrl-z\tgo back by 1 second\nctrl-x\tgo forward by 1 second\nctrl-c\tgo back by $zoom seconds\nctrl-v\tgo forward by $zoom seconds\n\nctrl-h\tthis info\nctrl-q\texit\n")
            elseif k == 98 # b
                if ch_first > 1
                    ch_first -= 1
                    ch_last -= 1
                    Gtk.@sigatom begin
                        set_gtk_property!(lab_ch, :label, "$(lpad(string(ch[ch_first]), 2, '0')):$(lpad(string(ch[ch_last]), 2, '0'))")
                        ch[ch_first] == ch[1] && set_gtk_property!(bt_chup, :sensitive, false)
                        ch[ch_last] != ch[end] && set_gtk_property!(bt_chdown, :sensitive, true)
                    end
                    draw(can)
                end
            elseif k == 110 # n
                if ch_last < ch[end]
                    ch_first += 1
                    ch_last += 1
                    Gtk.@sigatom begin
                        set_gtk_property!(lab_ch, :label, "$(lpad(string(ch[ch_first]), 2, '0')):$(lpad(string(ch[ch_last]), 2, '0'))")
                        ch[ch_first] != ch[1] && set_gtk_property!(bt_chup, :sensitive, true)
                        ch[ch_last] == ch[end] && set_gtk_property!(bt_chdown, :sensitive, false)
                    end
                    draw(can)
                end
            elseif k == 97 # a
                Gtk.@sigatom begin
                    set_gtk_property!(entry_time, :value, obj1.time_pts[1])
                end
                draw(can)
            elseif k == 115 # s
                time_current = obj1.time_pts[end] - zoom
                Gtk.@sigatom begin
                    set_gtk_property!(entry_time, :value, time_current)
                end
                draw(can)
            elseif k == 122 # z
                time_current = get_gtk_property(entry_time, :value, Float64)
                if time_current >= obj1.time_pts[1] + 1
                    time_current -= 1
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_time, :value, time_current)
                    end
                end
                draw(can)
            elseif k == 99 # c
                time_current = get_gtk_property(entry_time, :value, Float64)
                if time_current >= obj1.time_pts[1] + zoom
                    time_current = time_current - zoom
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_time, :value, time_current)
                    end
                end
                draw(can)
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
    iview_ep(obj1, obj2, ch)

Interactive view of epoched signal.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `obj2::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj1))`: channel(s) to plot, default is all channels
- `ep::Int64=1`: initial epoch to display
"""
function iview_ep(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj1)), ep::Int64=1)

    @assert size(obj1) == size(obj2) "Both signals must have the same size."
    @assert sr(obj1) == sr(obj2) "Both signals must have the same sampling rate."
    @assert nepochs(obj1) > 1 "iview_cont() should be used for continuous object."
    _check_channels(obj1, ch)
    _check_epochs(obj1, ep)

    if length(ch) > 10
        ch_first = 1
        ch_last = 10
    else
        ch_first = 1
        ch_last = length(ch)
    end

    if length(ch) > 1
        p = NeuroAnalyzer.plot(obj1, obj2, ch=ch[ch_first]:ch[ch_last], ep=ep, title="")
    else
        p = NeuroAnalyzer.plot(obj1, obj2, ch=ch, ep=ep, title="")
    end

    win = GtkWindow("NeuroAnalyzer: iview_ep()", Int32(p.attr[:size][1]), Int32(p.attr[:size][2]) + 40)
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
    bt_chup = GtkButton("△")
    set_gtk_property!(bt_chup, :tooltip_text, "Slide channels up")
    length(ch) < 11 && set_gtk_property!(bt_chup, :sensitive, false)
    ch[ch_first] == ch[1] && set_gtk_property!(bt_chup, :sensitive, false)
    bt_chdown = GtkButton("▽")
    set_gtk_property!(bt_chdown, :tooltip_text, "Slide channels down")
    length(ch) < 11 && set_gtk_property!(bt_chdown, :sensitive, false)
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
    lab_ch = GtkLabel("$(lpad(string(ch[ch_first]), 2, '0')):$(lpad(string(ch[ch_last]), 2, '0'))")
    set_gtk_property!(lab_ch, :halign, 0)
    g[1:14, 1] = can
    g[1, 2] = bt_chup
    g[2, 2] = lab_ch
    g[3, 2] = bt_chdown
    g[4, 2] = GtkLabel("")
    g[5, 2] = bt_start
    g[6, 2] = bt_prev
    g[7, 2] = entry_epoch
    g[8, 2] = bt_next
    g[9, 2] = bt_end
    g[10, 2] = GtkLabel("")
    g[11, 2] = bt_help
    g[12, 2] = bt_close
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
        if event.direction == 1 # down
            if ch_last < ch[end]
                ch_first += 1
                ch_last += 1
                Gtk.@sigatom begin
                    set_gtk_property!(lab_ch, :label, "$(lpad(string(ch[ch_first]), 2, '0')):$(lpad(string(ch[ch_last]), 2, '0'))")
                    ch[ch_first] != ch[1] && set_gtk_property!(bt_chup, :sensitive, true)
                    ch[ch_last] == ch[end] && set_gtk_property!(bt_chdown, :sensitive, false)
                end
                draw(can)
            end
        elseif event.direction == 0 # up
            if ch_first > 1
                ch_first -= 1
                ch_last -= 1
                Gtk.@sigatom begin
                    set_gtk_property!(lab_ch, :label, "$(lpad(string(ch[ch_first]), 2, '0')):$(lpad(string(ch[ch_last]), 2, '0'))")
                    ch[ch_first] == ch[1] && set_gtk_property!(bt_chup, :sensitive, false)
                    ch[ch_last] != ch[end] && set_gtk_property!(bt_chdown, :sensitive, true)
                end
                draw(can)
            end
        end
    end

    signal_connect(bt_chdown, "clicked") do widget
        if ch_last < ch[end]
            ch_first += 1
            ch_last += 1
            Gtk.@sigatom begin
                set_gtk_property!(lab_ch, :label, "$(lpad(string(ch[ch_first]), 2, '0')):$(lpad(string(ch[ch_last]), 2, '0'))")
                ch[ch_first] != ch[1] && set_gtk_property!(bt_chup, :sensitive, true)
                ch[ch_last] == ch[end] && set_gtk_property!(bt_chdown, :sensitive, false)
            end
            draw(can)
        end
    end

    signal_connect(bt_chup, "clicked") do widget
        if ch_first > 1
            ch_first -= 1
            ch_last -= 1
            Gtk.@sigatom begin
                set_gtk_property!(lab_ch, :label, "$(lpad(string(ch[ch_first]), 2, '0')):$(lpad(string(ch[ch_last]), 2, '0'))")
                ch[ch_first] == ch[1] && set_gtk_property!(bt_chup, :sensitive, false)
                ch[ch_last] != ch[end] && set_gtk_property!(bt_chdown, :sensitive, true)
            end
            draw(can)
        end
    end

    signal_connect(entry_epoch, "value-changed") do widget
        draw(can)
    end

    signal_connect(bt_prev, "clicked") do widget
        ep = get_gtk_property(entry_epoch, :value, Int64)
        if ep > 1
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
        info_dialog("Keyboard shortcuts:\n\nctrl-b\tslide channels up\nctrl-n\tslide channels down\n\nctrl-a\tgo to first epoch\nctrl-s\tgo to last epoch\nctrl-z\tprevious epoch\nctrl-x\tnext epoch\n\nctrl-h\tthis info\nctrl-q\texit\n")
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        s = event.state
        if s == 4
            if k == 113 # q
                Gtk.destroy(win)
            elseif k == 104 # h
                info_dialog("Keyboard shortcuts:\n\nctrl-b\tslide channels up\nctrl-n\tslide channels down\n\nctrl-a\tgo to first epoch\nctrl-s\tgo to last epoch\nctrl-z\tprevious epoch\nctrl-x\tnext epoch\n\nctrl-h\tthis info\nctrl-q\texit\n")
            elseif k == 98 # b
                if ch_first > 1
                    ch_first -= 1
                    ch_last -= 1
                    Gtk.@sigatom begin
                        set_gtk_property!(lab_ch, :label, "$(lpad(string(ch[ch_first]), 2, '0')):$(lpad(string(ch[ch_last]), 2, '0'))")
                        ch[ch_first] == ch[1] && set_gtk_property!(bt_chup, :sensitive, false)
                        ch[ch_last] != ch[end] && set_gtk_property!(bt_chdown, :sensitive, true)
                    end
                    draw(can)
                end
            elseif k == 110 # n
                if ch_last < ch[end]
                    ch_first += 1
                    ch_last += 1
                    Gtk.@sigatom begin
                        set_gtk_property!(lab_ch, :label, "$(lpad(string(ch[ch_first]), 2, '0')):$(lpad(string(ch[ch_last]), 2, '0'))")
                        ch[ch_first] != ch[1] && set_gtk_property!(bt_chup, :sensitive, true)
                        ch[ch_last] == ch[end] && set_gtk_property!(bt_chdown, :sensitive, false)
                    end
                    draw(can)
                end
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
                if ep > 1
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
    iview(obj, c; c_idx, zoom)

Interactive view of embedded or external component of continuous or epoched signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `c::Union{Symbol, AbstractArray}`: component to plot
- `c_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0`: component channel to display, default is all component channels
- `ep::Int64=1`: initial epoch to display
- `zoom::Real=5`: how many seconds are displayed in one segment
"""
function iview(obj::NeuroAnalyzer.NEURO, c::Union{Symbol, AbstractArray}; c_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0, ep::Int64=1, zoom::Real=5)

    if nepochs(obj) == 1
        iview_cont(obj, c, c_idx=c_idx, zoom=zoom)
    else
        iview_ep(obj, c, c_idx=c_idx, ep=ep)
    end

    return nothing

end

"""
    iview_cont(obj, c; c_idx, zoom)

Interactive view of embedded or external component of continuous signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `c::Union{Symbol, AbstractArray}`: component to plot
- `c_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0`: component channel to display, default is all component channels
- `zoom::Real=5`: how many seconds are displayed in one segment
"""
function iview_cont(obj::NeuroAnalyzer.NEURO, c::Union{Symbol, AbstractArray}; c_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0, zoom::Real=5)

    @assert zoom > 0 "zoom must be > 0."
    @assert zoom <= signal_len(obj) / sr(obj) "zoom must be ≤ $(signal_len(obj) / sr(obj))."
    @assert nepochs(obj) == 1 "iview_ep() should be used for epoched object."

    if length(c_idx) > 10
        c_idx_first = 1
        c_idx_last = 10
    else
        c_idx_first = 1
        c_idx_last = length(c_idx)
    end

    length(c_idx) == 1 && (c_idx = [c_idx])

    p = NeuroAnalyzer.plot(obj, c, c_idx=c_idx[c_idx_first:c_idx_last], mono=true, title="")
    win = GtkWindow("NeuroAnalyzer: iview_cont()", Int32(p.attr[:size][1]), Int32(p.attr[:size][2]) + 40)
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
    bt_chup = GtkButton("△")
    set_gtk_property!(bt_chup, :tooltip_text, "Slide channels up")
    length(c_idx) < 11 && set_gtk_property!(bt_chup, :sensitive, false)
    c_idx[c_idx_first] == c_idx[1] && set_gtk_property!(bt_chup, :sensitive, false)
    bt_chdown = GtkButton("▽")
    set_gtk_property!(bt_chdown, :tooltip_text, "Slide channels down")
    length(c_idx) < 11 && set_gtk_property!(bt_chdown, :sensitive, false)
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
    lab_ch = GtkLabel("$(lpad(string(c_idx[c_idx_first]), 2, '0')):$(lpad(string(c_idx[c_idx_last]), 2, '0'))")
    set_gtk_property!(lab_ch, :halign, 0)
    g[1:14, 1] = can
    g[1, 2] = bt_chup
    g[2, 2] = lab_ch
    g[3, 2] = bt_chdown
    g[4, 2] = GtkLabel("")
    g[5, 2] = bt_start
    g[6, 2] = bt_prev5
    g[7, 2] = bt_prev
    g[8, 2] = entry_time
    g[9, 2] = bt_next
    g[10, 2] = bt_next5
    g[11, 2] = bt_end
    g[12, 2] = GtkLabel("")
    g[13, 2] = bt_help
    g[14, 2] = bt_close
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
                                                       mono=true,
                                                       title=""))
        img = read_from_png(io)
        set_source_surface(ctx, img, 0, 0)
        paint(ctx)
    end

    can.mouse.scroll = @guarded (widget, event) -> begin
        if event.direction == 1 # down
            if ch_last < ch[end]
                ch_first += 1
                ch_last += 1
                Gtk.@sigatom begin
                    set_gtk_property!(lab_ch, :label, "$(lpad(string(ch[ch_first]), 2, '0')):$(lpad(string(ch[ch_last]), 2, '0'))")
                    ch[ch_first] != ch[1] && set_gtk_property!(bt_chup, :sensitive, true)
                    ch[ch_last] == ch[end] && set_gtk_property!(bt_chdown, :sensitive, false)
                end
                draw(can)
            end
        elseif event.direction == 0 # up
            if ch_first > 1
                ch_first -= 1
                ch_last -= 1
                Gtk.@sigatom begin
                    set_gtk_property!(lab_ch, :label, "$(lpad(string(ch[ch_first]), 2, '0')):$(lpad(string(ch[ch_last]), 2, '0'))")
                    ch[ch_first] == ch[1] && set_gtk_property!(bt_chup, :sensitive, false)
                    ch[ch_last] != ch[end] && set_gtk_property!(bt_chdown, :sensitive, true)
                end
                draw(can)
            end
        end
    end

    signal_connect(bt_chdown, "clicked") do widget
        if c_idx_last < c_idx[end]
            c_idx_first += 1
            c_idx_last += 1
            Gtk.@sigatom begin
                set_gtk_property!(lab_ch, :label, "$(lpad(string(c_idx[c_idx_first]), 2, '0')):$(lpad(string(c_idx[c_idx_last]), 2, '0'))")
                c_idx[c_idx_first] != c_idx[1] && set_gtk_property!(bt_chup, :sensitive, true)
                c_idx[c_idx_last] == c_idx[end] && set_gtk_property!(bt_chdown, :sensitive, false)
            end
            draw(can)
        end
    end

    signal_connect(bt_chup, "clicked") do widget
        if c_idx_first > 1
            c_idx_first -= 1
            c_idx_last -= 1
            Gtk.@sigatom begin
                set_gtk_property!(lab_ch, :label, "$(lpad(string(c_idx[c_idx_first]), 2, '0')):$(lpad(string(c_idx[c_idx_last]), 2, '0'))")
                c_idx[c_idx_first] == c_idx[1] && set_gtk_property!(bt_chup, :sensitive, false)
                c_idx[c_idx_last] != c_idx[end] && set_gtk_property!(bt_chdown, :sensitive, true)
            end
            draw(can)
        end
    end

    signal_connect(entry_time, "value-changed") do widget
        draw(can)
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

    signal_connect(bt_close, "clicked") do widget
        Gtk.destroy(win)
    end

    signal_connect(bt_help, "clicked") do widgete
        info_dialog("Keyboard shortcuts:\n\nctrl-b\tslide channels up\nctrl-n\tslide channels down\n\nctrl-a\tgo to the signal beginning\nctrl-s\tgo to the signal end\nctrl-z\tgo back by 1 second\nctrl-x\tgo forward by 1 second\nctrl-c\tgo back by $zoom seconds\nctrl-v\tgo forward by $zoom seconds\n\nctrl-h\tthis info\nctrl-q\texit\n")
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        s = event.state
        if s == 4
            if k == 113 # q
                Gtk.destroy(win)
            elseif k == 104 # h
                info_dialog("Keyboard shortcuts:\n\nctrl-b\tslide channels up\nctrl-n\tslide channels down\n\nctrl-a\tgo to the signal beginning\nctrl-s\tgo to the signal end\nctrl-z\tgo back by 1 second\nctrl-x\tgo forward by 1 second\nctrl-c\tgo back by $zoom seconds\nctrl-v\tgo forward by $zoom seconds\n\nctrl-h\tthis info\nctrl-q\texit\n")
            elseif k == 98 # b
                if c_idx_first > 1
                    c_idx_first -= 1
                    c_idx_last -= 1
                    Gtk.@sigatom begin
                        set_gtk_property!(lab_ch, :label, "$(lpad(string(c_idx[c_idx_first]), 2, '0')):$(lpad(string(c_idx[c_idx_last]), 2, '0'))")
                        c_idx[c_idx_first] == c_idx[1] && set_gtk_property!(bt_chup, :sensitive, false)
                        c_idx[c_idx_last] != c_idx[end] && set_gtk_property!(bt_chdown, :sensitive, true)
                    end
                    draw(can)
                end
            elseif k == 110 # n
                if c_idx_last < c_idx[end]
                    c_idx_first += 1
                    c_idx_last += 1
                    Gtk.@sigatom begin
                        set_gtk_property!(lab_ch, :label, "$(lpad(string(c_idx[c_idx_first]), 2, '0')):$(lpad(string(c_idx[c_idx_last]), 2, '0'))")
                        c_idx[c_idx_first] != c_idx[1] && set_gtk_property!(bt_chup, :sensitive, true)
                        c_idx[c_idx_last] == c_idx[end] && set_gtk_property!(bt_chdown, :sensitive, false)
                    end
                    draw(can)
                end
            elseif k == 97 # a
                Gtk.@sigatom begin
                    set_gtk_property!(entry_time, :value, obj.time_pts[1])
                end
                draw(can)
            elseif k == 115 # s
                time_current = obj.time_pts[end] - zoom
                Gtk.@sigatom begin
                    set_gtk_property!(entry_time, :value, time_current)
                end
                draw(can)
            elseif k == 122 # z
                time_current = get_gtk_property(entry_time, :value, Float64)
                if time_current >= obj.time_pts[1] + 1
                    time_current -= 1
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_time, :value, time_current)
                    end
                end
                draw(can)
            elseif k == 99 # c
                time_current = get_gtk_property(entry_time, :value, Float64)
                if time_current >= obj.time_pts[1] + zoom
                    time_current = time_current - zoom
                    Gtk.@sigatom begin
                        set_gtk_property!(entry_time, :value, time_current)
                    end
                end
                draw(can)
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

    cnd = Condition()
    signal_connect(win, :destroy) do widget
        notify(cnd)
    end
    @async Gtk.gtk_main()
    wait(cnd)

    return nothing

end

"""
    iview_ep(obj, c; c_idx, mono)

Interactive view of embedded or external component of epoched signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `c::Union{Symbol, AbstractArray}`: component to plot
- `c_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0`: component channel to display, default is all component channels
- `ep::Int64=1`: initial epoch to display
"""
function iview_ep(obj::NeuroAnalyzer.NEURO, c::Union{Symbol, AbstractArray}; c_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0, ep::Int64=1)

    @assert nepochs(obj) > 1 "iview_cont() should be used for continuous object."

    _check_epochs(obj, ep)

    if length(c_idx) > 10
        c_idx_first = 1
        c_idx_last = 10
    else
        c_idx_first = 1
        c_idx_last = length(c_idx)
    end

    p = NeuroAnalyzer.plot(obj, c, c_idx=c_idx[c_idx_first:c_idx_last], ep=ep, mono=true, title="")
    win = GtkWindow("NeuroAnalyzer: iview_ep()", Int32(p.attr[:size][1]), Int32(p.attr[:size][2]) + 40)
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
    bt_chup = GtkButton("△")
    set_gtk_property!(bt_chup, :tooltip_text, "Slide channels up")
    length(c_idx) < 11 && set_gtk_property!(bt_chup, :sensitive, false)
    c_idx[c_idx_first] == c_idx[1] && set_gtk_property!(bt_chup, :sensitive, false)
    bt_chdown = GtkButton("▽")
    set_gtk_property!(bt_chdown, :tooltip_text, "Slide channels down")
    length(c_idx) < 11 && set_gtk_property!(bt_chdown, :sensitive, false)
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
    lab_ch = GtkLabel("$(lpad(string(c_idx[c_idx_first]), 2, '0')):$(lpad(string(c_idx[c_idx_last]), 2, '0'))")
    set_gtk_property!(lab_ch, :halign, 0)
    g[1:14, 1] = can
    g[1, 2] = bt_chup
    g[2, 2] = lab_ch
    g[3, 2] = bt_chdown
    g[4, 2] = GtkLabel("")
    g[5, 2] = bt_start
    g[6, 2] = bt_prev
    g[7, 2] = entry_epoch
    g[8, 2] = bt_next
    g[9, 2] = bt_end
    g[10, 2] = GtkLabel("")
    g[11, 2] = bt_help
    g[12, 2] = bt_close
    push!(win, g)

    showall(win)

    @guarded draw(can) do widget
        ep = get_gtk_property(entry_epoch, :value, Int64)
        ctx = getgc(can)
        show(io, MIME("image/png"), NeuroAnalyzer.plot(obj,
                                                       c,
                                                       c_idx=c_idx[c_idx_first]:c_idx[c_idx_last],
                                                       ep=ep,
                                                       mono=true,
                                                       title=""))
        img = read_from_png(io)
        set_source_surface(ctx, img, 0, 0)
        paint(ctx)
    end

    can.mouse.scroll = @guarded (widget, event) -> begin
        if event.direction == 1 # down
            if ch_last < ch[end]
                ch_first += 1
                ch_last += 1
                Gtk.@sigatom begin
                    set_gtk_property!(lab_ch, :label, "$(lpad(string(ch[ch_first]), 2, '0')):$(lpad(string(ch[ch_last]), 2, '0'))")
                    ch[ch_first] != ch[1] && set_gtk_property!(bt_chup, :sensitive, true)
                    ch[ch_last] == ch[end] && set_gtk_property!(bt_chdown, :sensitive, false)
                end
                draw(can)
            end
        elseif event.direction == 0 # up
            if ch_first > 1
                ch_first -= 1
                ch_last -= 1
                Gtk.@sigatom begin
                    set_gtk_property!(lab_ch, :label, "$(lpad(string(ch[ch_first]), 2, '0')):$(lpad(string(ch[ch_last]), 2, '0'))")
                    ch[ch_first] == ch[1] && set_gtk_property!(bt_chup, :sensitive, false)
                    ch[ch_last] != ch[end] && set_gtk_property!(bt_chdown, :sensitive, true)
                end
                draw(can)
            end
        end
    end

    signal_connect(bt_chdown, "clicked") do widget
        if c_idx_last < c_idx[end]
            c_idx_first += 1
            c_idx_last += 1
            Gtk.@sigatom begin
                set_gtk_property!(lab_ch, :label, "$(lpad(string(c_idx[c_idx_first]), 2, '0')):$(lpad(string(c_idx[c_idx_last]), 2, '0'))")
                c_idx[c_idx_first] != c_idx[1] && set_gtk_property!(bt_chup, :sensitive, true)
                c_idx[c_idx_last] == c_idx[end] && set_gtk_property!(bt_chdown, :sensitive, false)
            end
            draw(can)
        end
    end

    signal_connect(bt_chup, "clicked") do widget
        if c_idx_first > 1
            c_idx_first -= 1
            c_idx_last -= 1
            Gtk.@sigatom begin
                set_gtk_property!(lab_ch, :label, "$(lpad(string(c_idx[c_idx_first]), 2, '0')):$(lpad(string(c_idx[c_idx_last]), 2, '0'))")
                c_idx[c_idx_first] == c_idx[1] && set_gtk_property!(bt_chup, :sensitive, false)
                c_idx[c_idx_last] != c_idx[end] && set_gtk_property!(bt_chdown, :sensitive, true)
            end
            draw(can)
        end
    end

    signal_connect(entry_epoch, "value-changed") do widget
        draw(can)
    end

    signal_connect(bt_prev, "clicked") do widget
        ep = get_gtk_property(entry_epoch, :value, Int64)
        if ep > 1
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
        info_dialog("Keyboard shortcuts:\n\nctrl-b\tslide channels up\nctrl-n\tslide channels down\n\nctrl-a\tgo to first epoch\nctrl-s\tgo to last epoch\nctrl-z\tprevious epoch\nctrl-x\tnext epoch\n\nctrl-h\tthis info\nctrl-q\texit\n")
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        s = event.state
        if s == 4
            if k == 113 # q
                Gtk.destroy(win)
            elseif k == 104 # h
                info_dialog("Keyboard shortcuts:\n\nctrl-b\tslide channels up\nctrl-n\tslide channels down\n\nctrl-a\tgo to first epoch\nctrl-s\tgo to last epoch\nctrl-z\tprevious epoch\nctrl-x\tnext epoch\n\nctrl-h\tthis info\nctrl-q\texit\n")
            elseif k == 98 # b
                if c_idx_first > 1
                    c_idx_first -= 1
                    c_idx_last -= 1
                    Gtk.@sigatom begin
                        set_gtk_property!(lab_ch, :label, "$(lpad(string(c_idx[c_idx_first]), 2, '0')):$(lpad(string(c_idx[c_idx_last]), 2, '0'))")
                        c_idx[c_idx_first] == c_idx[1] && set_gtk_property!(bt_chup, :sensitive, false)
                        c_idx[c_idx_last] != c_idx[end] && set_gtk_property!(bt_chdown, :sensitive, true)
                    end
                    draw(can)
                end
            elseif k == 110 # n
                if c_idx_last < c_idx[end]
                    c_idx_first += 1
                    c_idx_last += 1
                    Gtk.@sigatom begin
                        set_gtk_property!(lab_ch, :label, "$(lpad(string(c_idx[c_idx_first]), 2, '0')):$(lpad(string(c_idx[c_idx_last]), 2, '0'))")
                        c_idx[c_idx_first] != c_idx[1] && set_gtk_property!(bt_chup, :sensitive, true)
                        c_idx[c_idx_last] == c_idx[end] && set_gtk_property!(bt_chdown, :sensitive, false)
                    end
                    draw(can)
                end
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
                if ep > 1
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

    cnd = Condition()
    signal_connect(win, :destroy) do widget
        notify(cnd)
    end
    @async Gtk.gtk_main()
    wait(cnd)

    return nothing

end
