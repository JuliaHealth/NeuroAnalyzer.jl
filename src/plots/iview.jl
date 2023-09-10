export iview
export iview_ep
export iview_cont

"""
    iview(obj, ch, mono, zoom)

Interactive view of continuous or epoched signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: channel(s) to plot, default is all channels
- `mono::Bool=false`: Use color or gray palette
- `zoom::Int64=5`: how many seconds are displayed in one segment

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function iview(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=NeuroAnalyzer._c(nchannels(obj)), mono::Bool=false, zoom::Int64=5)

    if nepochs(obj) == 1
        iview_cont(obj, ch=ch, mono=mono, zoom=zoom)
    else
        iview_ep(obj, ch=ch, mono=mono)
    end

end

"""
    iview_cont(obj, ch, mono, zoom)

Interactive view of continuous signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: channel(s) to plot, default is all channels
- `mono::Bool=false`: Use color or gray palette
- `zoom::Int64=5`: how many seconds are displayed in one segment

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function iview_cont(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=NeuroAnalyzer._c(nchannels(obj)), mono::Bool=false, zoom::Int64=5)

    @assert zoom >= 1 "zoom must be ≥ 1."
    @assert zoom <= signal_len(obj) / sr(obj) "zoom must be ≤ $(signal_len(obj) / sr(obj))."
    @assert nepochs(obj) == 1 "iview_ep() should be used for epoched object."
    _check_channels(obj, ch)

    p = NeuroAnalyzer.plot(obj, ch=ch, mono=mono, title="")
    win = GtkWindow("NeuroAnalyzer: iview_cont()", 1200, 800)
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
    g[1:10, 1] = win_view
    g[1, 2] = bt_start
    g[2, 2] = bt_prev5
    g[3, 2] = bt_prev
    g[4, 2] = entry_time
    g[5, 2] = bt_next
    g[6, 2] = bt_next5
    g[7, 2] = bt_end
    g[8, 2] = GtkLabel("")
    g[9, 2] = bt_help
    g[10, 2] = bt_close
    push!(win, g)

    showall(win)

    @guarded draw(can) do widget
        time1 = get_gtk_property(entry_time, :value, Float64)
        time2 = time1 + zoom
        time2 > obj.time_pts[end] && (time2 = obj.time_pts[end])
        ctx = getgc(can)
        show(io, MIME("image/png"), NeuroAnalyzer.plot(obj,
                                                       ch=ch,
                                                       seg=(time1, time2),
                                                       mono=mono,
                                                       title=""))
        img = read_from_png(io)
        set_source_surface(ctx, img, 0, 0)
        paint(ctx)
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
- `mono::Bool=false`: Use color or gray palette

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function iview_ep(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=NeuroAnalyzer._c(nchannels(obj)), mono::Bool=false)

    @assert nepochs(obj) > 1 "iview_cont() should be used for continuous object."
    _check_channels(obj, ch)

    p = NeuroAnalyzer.plot(obj, ch=ch, ep=1, mono=mono, title="")
    win = GtkWindow("NeuroAnalyzer: iview_ep()", 1200, 800)
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
    set_gtk_property!(g, :column_homogeneous, false)
    set_gtk_property!(g, :column_spacing, 10)
    set_gtk_property!(g, :row_spacing, 10)
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
    bt_delete = GtkButton("DEL")
    set_gtk_property!(bt_delete, :tooltip_text, "Delete epoch")
    bt_close = GtkButton("✖")
    set_gtk_property!(bt_close, :tooltip_text, "Close this window")
    g[1:10, 1] = win_view
    g[1, 2] = bt_start
    g[2, 2] = bt_prev
    g[3, 2] = entry_epoch
    g[4, 2] = bt_next
    g[5, 2] = bt_end
    g[6, 2] = GtkLabel("")
    g[7, 2] = bt_delete
    g[8, 2] = GtkLabel("")
    g[9, 2] = bt_help
    g[10, 2] = bt_close
    push!(win, g)

    showall(win)

    @guarded draw(can) do widget
        ep = get_gtk_property(entry_epoch, :value, Int64)
        ctx = getgc(can)
        show(io, MIME("image/png"), NeuroAnalyzer.plot(obj,
                                                       ch=ch,
                                                       ep=ep,
                                                       mono=mono,
                                                       title=""))
        img = read_from_png(io)
        set_source_surface(ctx, img, 0, 0)
        paint(ctx)
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
        info_dialog("Keyboard shortcuts:\n\nHOME\tgo to first epoch\nEND\t\tgo to last epoch\n,\t\tprevious epoch\n.\t\tnext epoch\n\nDEL\t\tdelete current epoch\n\nh\t\tthis info\nq\t\texit\n")
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        if k == 113 # q
            Gtk.destroy(win)
        elseif k == 104 # h
            info_dialog("Keyboard shortcuts:\n\nHOME\tgo to first epoch\nEND\t\tgo to last epoch\n,\t\tprevious epoch\n.\t\tnext epoch\n\nDEL\t\tdelete current epoch\n\nh\t\tthis info\nq\t\texit\n")
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
