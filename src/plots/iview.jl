export iview
export iview_ep
export iview_cont

"""
    iview(obj, ch, zoom)

Interactive view of continuous or epoched signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: channel(s) to plot, default is all channels
- `zoom::Int64=5`: how many seconds are displayed in one segment

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function iview(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=NeuroAnalyzer._c(nchannels(obj)), zoom::Int64=5)

    if nepochs(obj) == 1
        iview_cont(obj, ch=ch, zoom=zoom)
    else
        iview_ep(obj, ch=ch)
    end

end

"""
    iview_cont(obj, ch, zoom)

Interactive view of continuous signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: channel(s) to plot, default is all channels
- `zoom::Int64=5`: how many seconds are displayed in one segment

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function iview_cont(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=NeuroAnalyzer._c(nchannels(obj)), zoom::Int64=5)

    @assert zoom >= 1 "zoom must be ≥ 1."
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

    p = NeuroAnalyzer.plot(obj, ch=ch[ch_first:ch_last], mono=true, title="")
    win = GtkWindow("NeuroAnalyzer: iview_cont()", Int32(p.attr[:size][1]), Int32(p.attr[:size][2]) + 40)
    set_gtk_property!(win, :border_width, 20)
    set_gtk_property!(win, :resizable, true)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
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
        show(io, MIME("image/png"), NeuroAnalyzer.plot(obj,
                                                       ch=ch[ch_first]:ch[ch_last],
                                                       seg=(time1, time2),
                                                       mono=true,
                                                       title=""))
        img = read_from_png(io)
        set_source_surface(ctx, img, 0, 0)
        paint(ctx)
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
        info_dialog("Keyboard shortcuts:\n\nb\tslide channels up\nn\tslide channels down\n\na\tgo to the signal beginning\ns\tgo to the signal end\nz\tgo back by 1 second\nx\tgo forward by 1 second\nc\tgo back by $zoom seconds\nv\tgo forward by $zoom seconds\n\nh\tthis info\nq\texit\n")
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        if k == 113 # q
            Gtk.destroy(win)
        elseif k == 104 # h
            info_dialog("Keyboard shortcuts:\n\nb\tslide channels up\nn\tslide channels down\n\na\tgo to the signal beginning\ns\tgo to the signal end\nz\tgo back by 1 second\nx\tgo forward by 1 second\nc\tgo back by $zoom seconds\nv\tgo forward by $zoom seconds\n\nh\tthis info\nq\texit\n")
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

    return nothing

end

"""
    iview_ep(obj, ch, mono)

Interactive view of epoched signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: channel(s) to plot, default is all channels

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function iview_ep(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=NeuroAnalyzer._c(nchannels(obj)))

    @assert nepochs(obj) > 1 "iview_cont() should be used for continuous object."
    _check_channels(obj, ch)


    if length(ch) > 10
        ch_first = 1
        ch_last = 10
    else
        ch_first = 1
        ch_last = length(ch)
    end

    p = NeuroAnalyzer.plot(obj, ch=ch[ch_first:ch_last], ep=1, mono=true, title="")
    win = GtkWindow("NeuroAnalyzer: iview_ep()", Int32(p.attr[:size][1]), Int32(p.attr[:size][2]) + 40)
    set_gtk_property!(win, :border_width, 20)
    set_gtk_property!(win, :resizable, true)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    can = GtkCanvas(Int32(p.attr[:size][1]), Int32(p.attr[:size][2]))
    g = GtkGrid()
    set_gtk_property!(g, :column_homogeneous, false)
    set_gtk_property!(g, :column_spacing, 10)
    set_gtk_property!(g, :row_spacing, 10)
    entry_epoch = GtkSpinButton(1, nepochs(obj), 1)
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
    bt_delete = GtkButton("DEL")
    set_gtk_property!(bt_delete, :tooltip_text, "Delete epoch")
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
    g[11, 2] = bt_delete
    g[12, 2] = GtkLabel("")
    g[13, 2] = bt_help
    g[14, 2] = bt_close
    push!(win, g)

    showall(win)

    @guarded draw(can) do widget
        ep = get_gtk_property(entry_epoch, :value, Int64)
        ctx = getgc(can)
        show(io, MIME("image/png"), NeuroAnalyzer.plot(obj,
                                                       ch=ch[ch_first]:ch[ch_last],
                                                       ep=ep,
                                                       mono=true,
                                                       title=""))
        img = read_from_png(io)
        set_source_surface(ctx, img, 0, 0)
        paint(ctx)
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

    signal_connect(bt_delete, "clicked") do widget
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

    signal_connect(bt_close, "clicked") do widget
        Gtk.destroy(win)
    end

    signal_connect(bt_help, "clicked") do widgete
        info_dialog("Keyboard shortcuts:\n\nb\tslide channels up\nn\tslide channels down\n\na\tgo to first epoch\ns\tgo to last epoch\nz\tprevious epoch\nx\tnext epoch\n\nh\tthis info\nq\texit\n")
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        if k == 113 # q
            Gtk.destroy(win)
        elseif k == 104 # h
            info_dialog("Keyboard shortcuts:\n\nb\tslide channels up\nn\tslide channels down\n\na\tgo to first epoch\ns\tgo to last epoch\nz\tprevious epoch\nx\tnext epoch\n\nh\tthis info\nq\texit\n")
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

    return nothing

end
