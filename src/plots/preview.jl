export preview
export preview_ep
export preview_cont

"""
    preview(obj; <keyword arguments>)

Interactive preview of continuous or epoched signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: channel(s) to plot, default is all channels
- `mono::Bool=false`: use color or grey palette

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function preview(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=NeuroAnalyzer._c(channel_n(obj)), mono::Bool=false)

    if epoch_n(obj) == 1
        preview_cont(obj, ch=ch, mono=mono)
    else
        preview_ep(obj, ch=ch, mono=mono)
    end

end

"""
    preview_ep(obj; <keyword arguments>)

Interactive preview of epoched signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: channel(s) to plot, default is all channels
- `mono::Bool=false`: use color or grey palette

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function preview_ep(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=NeuroAnalyzer._c(channel_n(obj)), mono::Bool=false)

    epoch_n(obj) < 2 && @error "preview_cont() should be used for continuous object."
    _check_channels(obj, ch)

    p = NeuroAnalyzer.plot(obj, ch=ch)
    win = GtkWindow("NeuroAnalyzer", 1200, (p.attr[:size][2] + 90))
    can = GtkCanvas(Int32(p.attr[:size][1]), Int32(p.attr[:size][2]))

    g = GtkGrid()
    slider_ep = GtkScale(false, 1:epoch_n(obj))
    label_ep = GtkLabel("Epoch")
    bt_start = GtkButton("⇤")
    bt_prev = GtkButton("←")
    bt_next = GtkButton("→")
    bt_end = GtkButton("⇥")
    bt_help = GtkButton("🛈")
    bt_delete = GtkButton("DEL")
    bt_close = GtkButton("✖")

    @guarded draw(can) do widget
        ctx = getgc(can)
        ep = Int(GAccessor.value(slider_ep))
        show(io, MIME("image/png"), NeuroAnalyzer.plot(obj, ch=ch, ep=ep, mono=mono))
        img = read_from_png(io)
        set_source_surface(ctx, img, 0, 0)
        paint(ctx)
    end
    signal_connect((w) -> draw(can), slider_ep, "value-changed")

    signal_connect(bt_prev, "clicked") do widget
        Gtk.@sigatom begin
            ep = Int(GAccessor.value(slider_ep))
            if ep > 1
                GAccessor.value(slider_ep, ep - 1)
            end
        end
    end

    signal_connect(bt_next, "clicked") do widget
        Gtk.@sigatom begin
            ep = Int(GAccessor.value(slider_ep))
            if ep < epoch_n(obj)
                GAccessor.value(slider_ep, ep + 1)
            end
        end
    end

    signal_connect(bt_start, "clicked") do widget
        Gtk.@sigatom begin
            GAccessor.value(slider_ep, 1)
        end
    end

    signal_connect(bt_end, "clicked") do widget
        Gtk.@sigatom begin
            GAccessor.value(slider_ep, epoch_n(obj))
        end
    end

    signal_connect(bt_close, "clicked") do widget
        Gtk.destroy(win)
    end

    signal_connect(bt_help, "clicked") do widget
        info_dialog("Keyboard shortcuts:\n\nHOME\tgo to first epoch\nEND\t\tgo to last epoch\n,\t\tprevious epoch\n.\t\tnext epoch\n\nDEL\t\tdelete current epoch\n\nh\t\tthis info\nq\t\texit\n")
    end

    signal_connect(bt_delete, "clicked") do widget
        Gtk.@sigatom begin
            ep = Int(GAccessor.value(slider_ep))
            if ask_dialog("Delete epoch $ep ?", "No", "Yes")
                delete_epoch!(obj, ep=ep)
                _info("Deleted epoch: $ep")
                GAccessor.range(slider_ep, 1, epoch_n(obj))
                if ep != 1
                    GAccessor.value(slider_ep, ep - 1)
                else
                    GAccessor.value(slider_ep, 1)
                end
                draw(can)
            end
        end
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        if k == 113 # q
            Gtk.destroy(win)
        elseif k == 65360 # HOME
            Gtk.@sigatom begin
                GAccessor.value(slider_ep, 1)
            end
        elseif k == 65367 # END
            Gtk.@sigatom begin
                GAccessor.value(slider_ep, epoch_n(obj))
            end
        elseif k == 104 # h
            info_dialog("Keyboard shortcuts:\n\nHOME\tgo to first epoch\nEND\t\tgo to last epoch\n,\t\tprevious epoch\n.\t\tnext epoch\n\nDEL\t\tdelete current epoch\n\nh\t\tthis info\nq\t\texit\n")
        elseif k == 44 # ,
            Gtk.@sigatom begin
                ep = Int(GAccessor.value(slider_ep))
                if ep > 1
                    GAccessor.value(slider_ep, ep - 1)
                else
                    GAccessor.value(slider_ep, 1)
                end
            end
        elseif k == 46 # .
            Gtk.@sigatom begin
                ep = Int(GAccessor.value(slider_ep))
                if ep < epoch_n(obj)
                    GAccessor.value(slider_ep, ep + 1)
                else
                    GAccessor.value(slider_ep, epoch_n(obj))
                end
            end
        elseif k == 65535 # DEL
            Gtk.@sigatom begin
                ep = Int(GAccessor.value(slider_ep))
                if ask_dialog("Delete epoch $ep ?", "No", "Yes")
                    delete_epoch!(obj, ep=ep)
                    _info("Deleted epoch: $ep")
                    GAccessor.range(slider_ep, 1, epoch_n(obj))
                    if ep != 1
                        GAccessor.value(slider_ep, ep - 1)
                    else
                        GAccessor.value(slider_ep, 1)
                    end
                    draw(can)
                end
            end
        end
    end

    g[1:7, 1] = can         # spans all columns
    g[1, 2] = label_ep
    g[2:7, 2] = slider_ep   # spans all columns
    g[1, 3] = bt_start      # Cartesian coordinates, g[x,y]
    g[2, 3] = bt_prev
    g[3, 3] = bt_next
    g[4, 3] = bt_end
    g[5, 3] = bt_delete
    g[6, 3] = bt_help
    g[7, 3] = bt_close
    set_gtk_property!(g, :column_homogeneous, true)
    set_gtk_property!(g, :column_spacing, 10)  # introduce a 5-pixel gap between columns
    set_gtk_property!(win, :border_width, 20)
    push!(win, g)

    showall(win)

end

"""
    preview_cont(obj; <keyword arguments>)

Interactive preview of continuous signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: channel(s) to plot, default is all channels
- `mono::Bool=false`: use color or grey palette

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function preview_cont(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=NeuroAnalyzer._c(channel_n(obj)), mono::Bool=false)

    epoch_n(obj) > 1 && @error "preview_ep() should be used for epoched object."
    _check_channels(obj, ch)

    zoom_level = 10

    p = NeuroAnalyzer.plot(obj, ch=ch)
    win = GtkWindow("NeuroAnalyzer", 1200, (p.attr[:size][2] + 180))
    can = GtkCanvas(Int32(p.attr[:size][1]), Int32(p.attr[:size][2]))
    g = GtkGrid()
    slider_time = GtkScale(false, obj.time_pts[1]:obj.time_pts[end])
    slider_ts1 = GtkScale(false, obj.time_pts[1]:obj.time_pts[1] + 10)
    slider_ts2 = GtkScale(false, obj.time_pts[1]:obj.time_pts[1] + 10)
    label_time = GtkLabel("Signal time")
    label_ts1 = GtkLabel("Segment start")
    label_ts2 = GtkLabel("Segment end")
    bt_start = GtkButton("⇤")
    bt_prev10 = GtkButton("⇐")
    bt_prev = GtkButton("←")
    bt_next = GtkButton("→")
    bt_next10 = GtkButton("⇒")
    bt_end = GtkButton("⇥")
    bt_zoomin = GtkButton("▼")
    label_zoomlevel = GtkLabel("10")
    bt_zoomout = GtkButton("▲")
    bt_help = GtkButton("🛈")
    bt_delete = GtkButton("DEL")
    bt_close = GtkButton("✖")
    separator1 = GtkLabel("")
    separator2 = GtkLabel("")
    separator3 = GtkLabel("")
    separator4 = GtkLabel("")

    @guarded draw(can) do widget

        time1 = GAccessor.value(slider_time)
        time2 = time1 + zoom_level
        time2 > obj.time_pts[end] && (time2 = obj.time_pts[end])
        ctx = getgc(can)
        show(io, MIME("image/png"), NeuroAnalyzer.plot(obj, ch=ch, seg=(time1, time2), mono=mono))
        img = read_from_png(io)
        set_source_surface(ctx, img, 0, 0)
        paint(ctx)
    end

    # signal_connect((w) -> draw(can), slider_time, "value-changed")
    signal_connect(slider_time, "value-changed") do widget
        time1 = GAccessor.value(slider_time)
        time2 = time1 + zoom_level
        time2 > obj.time_pts[end] && (time2 = obj.time_pts[end])
        Gtk.@sigatom begin
            GAccessor.range(slider_ts1, time1, time2)
            GAccessor.range(slider_ts2, time1, time2)
            GAccessor.value(slider_ts1, time1)
            GAccessor.value(slider_ts2, time1)
        end
        draw(can)
    end

    signal_connect(bt_prev, "clicked") do widget
        if GAccessor.value(slider_time) > obj.time_pts[1] + 1
            time_current = GAccessor.value(slider_time) - 1
            Gtk.@sigatom begin
                GAccessor.value(slider_time, time_current)
                GAccessor.value(slider_ts1, time_current)
                GAccessor.value(slider_ts2, time_current)
            end
        else
            Gtk.@sigatom begin
                GAccessor.value(slider_time, obj.time_pts[1])
                GAccessor.value(slider_ts1, obj.time_pts[1])
                GAccessor.value(slider_ts2, obj.time_pts[1])
            end
        end
    end

    signal_connect(bt_next, "clicked") do widget
        if GAccessor.value(slider_time) < obj.time_pts[end] - 1
            time_current = GAccessor.value(slider_time) + 1
            Gtk.@sigatom begin
                GAccessor.value(slider_time, time_current)
                GAccessor.value(slider_ts1, time_current)
                GAccessor.value(slider_ts2, time_current)
            end
        else
            Gtk.@sigatom begin
                GAccessor.value(slider_time, round(obj.time_pts[end] - 1))
                GAccessor.value(slider_ts1, round(obj.time_pts[end] - 1))
                GAccessor.value(slider_ts2, round(obj.time_pts[end] - 1))
            end
        end
    end

    signal_connect(bt_prev10, "clicked") do widget
        time_current = GAccessor.value(slider_time)
        if time_current > obj.time_pts[1]
            time_current = time_current - 10
            Gtk.@sigatom begin
                GAccessor.value(slider_time, time_current)
                GAccessor.value(slider_ts1, time_current)
                GAccessor.value(slider_ts2, time_current)
            end
        end
    end

    signal_connect(bt_next10, "clicked") do widget
        if obj.time_pts[end] % zoom_level == 0
            if GAccessor.value(slider_time) < obj.time_pts[end] - 10
                time_current = GAccessor.value(slider_time) + 10
                Gtk.@sigatom begin
                    GAccessor.value(slider_time, time_current)
                    GAccessor.value(slider_ts1, time_current)
                    GAccessor.value(slider_ts2, time_current)
                end
            end
        else
            if GAccessor.value(slider_time) < obj.time_pts[end] - (obj.time_pts[end] % 10)
                time_current = GAccessor.value(slider_time) + 10
                Gtk.@sigatom begin
                    GAccessor.value(slider_time, time_current)
                    GAccessor.value(slider_ts1, time_current)
                    GAccessor.value(slider_ts2, time_current)
                end
            else
                time_current = obj.time_pts[end] - (obj.time_pts[end] % 10)
                Gtk.@sigatom begin
                    GAccessor.value(slider_time, time_current)
                    GAccessor.value(slider_ts1, time_current)
                    GAccessor.value(slider_ts2, time_current)
                end
            end
        end
    end

    signal_connect(bt_start, "clicked") do widget
        Gtk.@sigatom begin
            GAccessor.value(slider_time, obj.time_pts[1])
            GAccessor.value(slider_ts1, obj.time_pts[1])
            GAccessor.value(slider_ts2, obj.time_pts[1])
        end
    end

    signal_connect(bt_end, "clicked") do widget
        if obj.time_pts[end] % zoom_level == 0
            Gtk.@sigatom begin
                GAccessor.value(slider_time, obj.time_pts[end] - zoom_level)
                GAccessor.value(slider_ts1, obj.time_pts[end] - zoom_level)
                GAccessor.value(slider_ts2, obj.time_pts[end] - zoom_level)
            end
        else
            Gtk.@sigatom begin
                GAccessor.value(slider_time, obj.time_pts[end] - obj.time_pts[end] % zoom_level)
                GAccessor.value(slider_ts1, obj.time_pts[end] - obj.time_pts[end] % zoom_level)
                GAccessor.value(slider_ts2, obj.time_pts[end] - obj.time_pts[end] % zoom_level)
            end
        end
    end

    signal_connect(bt_zoomin, "clicked") do widget
        time_current = GAccessor.value(slider_time)
        if zoom_level == 10
            zoom_level = 5
        elseif zoom_level == 5
            zoom_level = 1
        end
        GAccessor.text(label_zoomlevel, lpad(string(zoom_level), 2, '0'))
        if time_current + zoom_level <= obj.time_pts[end]
            Gtk.@sigatom begin
                GAccessor.range(slider_ts1, time_current, time_current + zoom_level)
                GAccessor.range(slider_ts2, time_current, time_current + zoom_level)
                GAccessor.value(slider_ts1, time_current)
                GAccessor.value(slider_ts2, time_current)
            end
            draw(can)
        else
            Gtk.@sigatom begin
                GAccessor.range(slider_ts1, obj.time_pts[end])
                GAccessor.range(slider_ts2, obj.time_pts[end])
                GAccessor.value(slider_ts1, time_current)
                GAccessor.value(slider_ts2, time_current)
            end
            draw(can)
        end
    end

    signal_connect(bt_zoomout, "clicked") do widget
        time_current = GAccessor.value(slider_time)
        if zoom_level == 1
            zoom_level = 5
        elseif zoom_level == 5
            zoom_level = 10
        end
        GAccessor.text(label_zoomlevel, lpad(string(zoom_level), 2, '0'))
        if time_current + zoom_level < obj.time_pts[end]
            Gtk.@sigatom begin
                GAccessor.range(slider_ts1, time_current, time_current + zoom_level)
                GAccessor.range(slider_ts2, time_current, time_current + zoom_level)
                GAccessor.value(slider_ts1, time_current)
                GAccessor.value(slider_ts2, time_current)
            end
            draw(can)
        else
            Gtk.@sigatom begin
                GAccessor.range(slider_ts1, obj.time_pts[end])
                GAccessor.range(slider_ts2, obj.time_pts[end])
                GAccessor.value(slider_ts1, time_current)
                GAccessor.value(slider_ts2, time_current)
            end
            draw(can)
        end
    end

    signal_connect(bt_delete, "clicked") do widget
        time_current = GAccessor.value(slider_time)
        time1 = GAccessor.value(slider_ts1)
        time2 = GAccessor.value(slider_ts2)
        if time1 < time2
            if ask_dialog("Delete segment $time1:$time2 ?", "No", "Yes")
                trim!(obj, seg=(time1, time2), remove_epochs=false)
                _info("Deleted segment: $time1:$time2")
                Gtk.@sigatom begin
                    GAccessor.range(slider_time, obj.time_pts[1], obj.time_pts[end])
                    if obj.time_pts[end] % zoom_level == 0
                        time_current >= (obj.time_pts[end] - zoom_level) && (time_current = obj.time_pts[end] - zoom_level)
                    else
                        time_current >= obj.time_pts[end] - (obj.time_pts[end] % zoom_level) && (time_current = obj.time_pts[end] - (obj.time_pts[end] % zoom_level))
                    end
                    time_current < obj.time_pts[1] && (time_current = obj.time_pts[1])
                    GAccessor.value(slider_time, time_current)
                    GAccessor.value(slider_ts1, time_current)
                    GAccessor.value(slider_ts2, time_current)
                end
                draw(can)
            end
        else
            info_dialog("Time segment start must be < time segment end.")
        end
    end

    signal_connect(bt_close, "clicked") do widget
        Gtk.destroy(win)
    end

    signal_connect(bt_help, "clicked") do widget
        info_dialog("Keyboard shortcuts:\n\nHOME\tgo to the signal beginning\nEND\t\tgo to the signal end\n,\t\tgo back by 1 second\n.\t\tgo forward by 1 second\n<\t\tgo back by 10 seconds\n>\t\tgo forward by 10 seconds\n\n+\t\tzoom in\n−\t\tzoom out\n\nDEL\t\tdelete current segment\n\nh\t\tthis info\nq\t\texit\n")
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        if k == 113 # q
            Gtk.destroy(win)
        elseif k == 65360 # HOME
            Gtk.@sigatom begin
                GAccessor.value(slider_time, obj.time_pts[1])
                GAccessor.value(slider_ts1, obj.time_pts[1])
                GAccessor.value(slider_ts2, obj.time_pts[1])
            end
        elseif k == 65367 # END
            if obj.time_pts[end] % 10 == 0
                Gtk.@sigatom begin
                    GAccessor.value(slider_time, obj.time_pts[end] - 10)
                    GAccessor.value(slider_ts1, obj.time_pts[end] - 10)
                    GAccessor.value(slider_ts2, obj.time_pts[end] - 10)
                end
            else
                Gtk.@sigatom begin
                    GAccessor.value(slider_time, obj.time_pts[end] - obj.time_pts[end] % 10)
                    GAccessor.value(slider_ts1, obj.time_pts[end] - obj.time_pts[end] % 10)
                    GAccessor.value(slider_ts2, obj.time_pts[end] - obj.time_pts[end] % 10)
                end
            end
        elseif k == 104 # h
            info_dialog("Keyboard shortcuts:\n\nHOME\tgo to the signal beginning\nEND\t\tgo to the signal end\n,\t\tgo back by 1 second\n.\t\tgo forward by 1 second\n<\t\tgo back by 10 seconds\n>\t\tgo forward by 10 seconds\n\n+\t\tzoom in\n−\t\tzoom out\n\nDEL\t\tdelete current segment\n\nh\t\tthis info\nq\t\texit\n")
        elseif k == 65453 || k == 45 # -
            time_current = GAccessor.value(slider_time)
            if zoom_level == 10
                zoom_level = 5
            elseif zoom_level == 5
                zoom_level = 1
            end
            GAccessor.text(label_zoomlevel, lpad(string(zoom_level), 2, '0'))
            if time_current + zoom_level < obj.time_pts[end]
                Gtk.@sigatom begin
                    GAccessor.range(slider_ts1, time_current, time_current + zoom_level)
                    GAccessor.range(slider_ts2, time_current, time_current + zoom_level)
                    GAccessor.value(slider_ts1, time_current)
                    GAccessor.value(slider_ts2, time_current)
                end
                draw(can)
            else
                Gtk.@sigatom begin
                    GAccessor.range(slider_ts1, obj.time_pts[end])
                    GAccessor.range(slider_ts2, obj.time_pts[end])
                    GAccessor.value(slider_ts1, time_current)
                    GAccessor.value(slider_ts2, time_current)
                end
                draw(can)
            end
        elseif k == 65451 || k == 43 # +
            time_current = GAccessor.value(slider_time)
            if zoom_level == 1
                zoom_level = 5
            elseif zoom_level == 5
                zoom_level = 10
            end
            GAccessor.text(label_zoomlevel, lpad(string(zoom_level), 2, '0'))
            if time_current + zoom_level < obj.time_pts[end]
                Gtk.@sigatom begin
                    GAccessor.range(slider_ts1, time_current, time_current + zoom_level)
                    GAccessor.range(slider_ts2, time_current, time_current + zoom_level)
                    GAccessor.value(slider_ts1, time_current)
                    GAccessor.value(slider_ts2, time_current)
                end
                draw(can)
            else
                Gtk.@sigatom begin
                    GAccessor.range(slider_ts1, obj.time_pts[end])
                    GAccessor.range(slider_ts2, obj.time_pts[end])
                    GAccessor.value(slider_ts1, time_current)
                    GAccessor.value(slider_ts2, time_current)
                end
                draw(can)
            end
        elseif k == 44 # ,
            if GAccessor.value(slider_time) > obj.time_pts[1] + 1
                time_current = GAccessor.value(slider_time) - 1
                Gtk.@sigatom begin
                    GAccessor.value(slider_time, time_current)
                    GAccessor.value(slider_ts1, time_current)
                    GAccessor.value(slider_ts2, time_current)
                end
            else
                Gtk.@sigatom begin
                    GAccessor.value(slider_time, obj.time_pts[1])
                    GAccessor.value(slider_ts1, obj.time_pts[1])
                    GAccessor.value(slider_ts2, obj.time_pts[1])
                end
            end
        elseif k == 46 # .
            if GAccessor.value(slider_time) < obj.time_pts[end] - 1
                time_current = GAccessor.value(slider_time) + 1
                Gtk.@sigatom begin
                    GAccessor.value(slider_time, time_current)
                    GAccessor.value(slider_ts1, time_current)
                    GAccessor.value(slider_ts2, time_current)
                end
            end
        elseif k == 60 # <
            time_current = GAccessor.value(slider_time)
            if time_current > obj.time_pts[1]
                time_current = time_current - 10
                Gtk.@sigatom begin
                    GAccessor.value(slider_time, time_current)
                    GAccessor.value(slider_ts1, time_current)
                    GAccessor.value(slider_ts2, time_current)
                end
            end
        elseif k == 62 # >
            if obj.time_pts[end] % 10 == 0
                if GAccessor.value(slider_time) < obj.time_pts[end] - 10
                    time_current = GAccessor.value(slider_time) + 10
                    Gtk.@sigatom begin
                        GAccessor.value(slider_time, time_current)
                        GAccessor.value(slider_ts1, time_current)
                        GAccessor.value(slider_ts2, time_current)
                    end
                end
            else
                if GAccessor.value(slider_time) < obj.time_pts[end] - (obj.time_pts[end] % 10)
                    time_current = GAccessor.value(slider_time) + 10
                    Gtk.@sigatom begin
                        GAccessor.value(slider_time, time_current)
                        GAccessor.value(slider_ts1, time_current)
                        GAccessor.value(slider_ts2, time_current)
                    end
                else
                    time_current = obj.time_pts[end] - (obj.time_pts[end] % 10)
                    Gtk.@sigatom begin
                        GAccessor.value(slider_time, time_current)
                        GAccessor.value(slider_ts1, time_current)
                        GAccessor.value(slider_ts2, time_current)
                    end
                end
            end
        elseif k == 65535 # DEL
            time_current = GAccessor.value(slider_time)
            time1 = GAccessor.value(slider_ts1)
            time2 = GAccessor.value(slider_ts2)
            if time1 < time2
                if ask_dialog("Delete segment $time1:$time2 ?", "No", "Yes")
                    trim!(obj, seg=(time1, time2), remove_epochs=false)
                    _info("Deleted segment: $time1:$time2")
                    Gtk.@sigatom begin
                        GAccessor.range(slider_time, obj.time_pts[1], obj.time_pts[end])
                        if obj.time_pts[end] % zoom_level == 0
                            time_current >= (obj.time_pts[end] - zoom_level) && (time_current = obj.time_pts[end] - zoom_level)
                        else
                            time_current >= obj.time_pts[end] - (obj.time_pts[end] % zoom_level) && (time_current = obj.time_pts[end] - (obj.time_pts[end] % zoom_level))
                        end
                        time_current < obj.time_pts[1] && (time_current = obj.time_pts[1])
                        GAccessor.value(slider_time, time_current)
                        GAccessor.value(slider_ts1, time_current)
                        GAccessor.value(slider_ts2, time_current)
                    end
                    draw(can)
                end
            else
                info_dialog("Time segment start must be < time segment end.")
            end
        end
    end

    g[1:16, 1] = can
    g[1, 2] = label_time
    g[2:16, 2] = slider_time
    g[1, 3] = label_ts1
    g[2:16, 3] = slider_ts1
    g[1, 4] = label_ts2
    g[2:16, 4] = slider_ts2
    g[1, 5] = separator1
    g[2, 5] = bt_start
    g[3, 5] = bt_prev10
    g[4, 5] = bt_prev
    g[5, 5] = bt_next
    g[6, 5] = bt_next10
    g[7, 5] = bt_end
    g[8, 5] = separator2
    g[9, 5] = bt_zoomin
    g[10, 5] = label_zoomlevel
    g[11, 5] = bt_zoomout
    g[12, 5] = separator3
    g[13, 5] = bt_delete
    g[14, 5] = separator4
    g[15, 5] = bt_help
    g[16, 5] = bt_close
    set_gtk_property!(g, :column_homogeneous, false)
    set_gtk_property!(g, :column_spacing, 10)  # introduce a 10-pixel gap between columns
    set_gtk_property!(win, :border_width, 20)
    set_gtk_property!(win, :resizable, false)
    set_gtk_property!(slider_time, :round_digits, 1)
    set_gtk_property!(slider_ts1, :round_digits, 1)
    set_gtk_property!(slider_ts2, :round_digits, 1)
    set_gtk_property!(slider_time, :digits, 1)
    set_gtk_property!(slider_ts1, :digits, 1)
    set_gtk_property!(slider_ts2, :digits, 1)
    push!(win, g)

    showall(win)

    return nothing

end

"""
    preview(obj1, obj2; <keyword arguments>)

Interactive preview of two continuous or epoched signals.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object (before) - drawn in black
- `obj2::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object (after) - drawn in red
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj1))`: channel(s) to plot, default is all channels

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function preview(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=NeuroAnalyzer._c(channel_n(obj1)))

    if epoch_n(obj1) == 1
        preview_cont(obj1, obj2, ch=ch)
    else
        preview_ep(obj1, obj2, ch=ch)
    end

end

"""
    preview_ep(obj1, obj2; <keyword arguments>)

Interactive preview of epoched signal.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object (before) - drawn in black
- `obj2::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object (after) - drawn in red
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj1))`: channel(s) to plot, default is all channels

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function preview_ep(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=NeuroAnalyzer._c(channel_n(obj1)))

    epoch_n(obj1) < 2 && @error "preview_cont() should be used for continuous object."
    _check_channels(obj1, ch1)

    p = NeuroAnalyzer.plot(obj1, obj2, ch=ch)
    win = GtkWindow("NeuroAnalyzer", 1200, (p.attr[:size][2] + 90))
    can = GtkCanvas(Int32(p.attr[:size][1]), Int32(p.attr[:size][2]))

    g = GtkGrid()
    slider_ep = GtkScale(false, 1:epoch_n(obj1))
    label_ep = GtkLabel("Epoch")
    bt_start = GtkButton("⇤")
    bt_prev = GtkButton("←")
    bt_next = GtkButton("→")
    bt_end = GtkButton("⇥")
    bt_help = GtkButton("🛈")
    bt_delete = GtkButton("DEL")
    bt_close = GtkButton("✖")

    @guarded draw(can) do widget
        ctx = getgc(can)
        ep = Int(GAccessor.value(slider_ep))
        show(io, MIME("image/png"), NeuroAnalyzer.plot(obj1, obj2, ch=ch, ep=ep))
        img = read_from_png(io)
        set_source_surface(ctx, img, 0, 0)
        paint(ctx)
    end
    signal_connect((w) -> draw(can), slider_ep, "value-changed")

    signal_connect(bt_prev, "clicked") do widget
        Gtk.@sigatom begin
            ep = Int(GAccessor.value(slider_ep))
            if ep > 1
                GAccessor.value(slider_ep, ep - 1)
            end
        end
    end

    signal_connect(bt_next, "clicked") do widget
        Gtk.@sigatom begin
            ep = Int(GAccessor.value(slider_ep))
            if ep < epoch_n(obj1)
                GAccessor.value(slider_ep, ep + 1)
            end
        end
    end

    signal_connect(bt_start, "clicked") do widget
        Gtk.@sigatom begin
            GAccessor.value(slider_ep, 1)
        end
    end

    signal_connect(bt_end, "clicked") do widget
        Gtk.@sigatom begin
            GAccessor.value(slider_ep, epoch_n(obj1))
        end
    end

    signal_connect(bt_close, "clicked") do widget
        Gtk.destroy(win)
    end

    signal_connect(bt_help, "clicked") do widget
        info_dialog("Keyboard shortcuts:\n\nHOME\tgo to first epoch\nEND\t\tgo to last epoch\n,\t\tprevious epoch\n.\t\tnext epoch\n\nDEL\t\tdelete current epoch\n\nh\t\tthis info\nq\t\texit\n")
    end

    signal_connect(bt_delete, "clicked") do widget
        Gtk.@sigatom begin
            ep = Int(GAccessor.value(slider_ep))
            if ask_dialog("Delete epoch $ep ?", "No", "Yes")
                delete_epoch!(obj1, ep=ep)
                delete_epoch!(obj2, ep=ep)
                _info("Deleted epoch: $ep")
                GAccessor.range(slider_ep, 1, epoch_n(obj1))
                if ep != 1
                    GAccessor.value(slider_ep, ep - 1)
                else
                    GAccessor.value(slider_ep, 1)
                end
                draw(can)
            end
        end
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        if k == 113 # q
            Gtk.destroy(win)
        elseif k == 65360 # HOME
            Gtk.@sigatom begin
                GAccessor.value(slider_ep, 1)
            end
        elseif k == 65367 # END
            Gtk.@sigatom begin
                GAccessor.value(slider_ep, epoch_n(obj1))
            end
        elseif k == 104 # h
            info_dialog("Keyboard shortcuts:\n\nHOME\tgo to first epoch\nEND\t\tgo to last epoch\n,\t\tprevious epoch\n.\t\tnext epoch\n\nDEL\t\tdelete current epoch\n\nh\t\tthis info\nq\t\texit\n")
        elseif k == 44 # ,
            Gtk.@sigatom begin
                ep = Int(GAccessor.value(slider_ep))
                if ep > 1
                    GAccessor.value(slider_ep, ep - 1)
                else
                    GAccessor.value(slider_ep, 1)
                end
            end
        elseif k == 46 # .
            Gtk.@sigatom begin
                ep = Int(GAccessor.value(slider_ep))
                if ep < epoch_n(obj1)
                    GAccessor.value(slider_ep, ep + 1)
                else
                    GAccessor.value(slider_ep, epoch_n(obj1))
                end
            end
        elseif k == 65535 # DEL
            Gtk.@sigatom begin
                ep = Int(GAccessor.value(slider_ep))
                if ask_dialog("Delete epoch $ep ?", "No", "Yes")
                    delete_epoch!(obj1, ep=ep)
                    delete_epoch!(obj2, ep=ep)
                    _info("Deleted epoch: $ep")
                    GAccessor.range(slider_ep, 1, epoch_n(obj1))
                    if ep != 1
                        GAccessor.value(slider_ep, ep - 1)
                    else
                        GAccessor.value(slider_ep, 1)
                    end
                    draw(can)
                end
            end
        end
    end

    g[1:7, 1] = can         # spans all columns
    g[1, 2] = label_ep
    g[2:7, 2] = slider_ep   # spans all columns
    g[1, 3] = bt_start      # Cartesian coordinates, g[x,y]
    g[2, 3] = bt_prev
    g[3, 3] = bt_next
    g[4, 3] = bt_end
    g[5, 3] = bt_delete
    g[6, 3] = bt_help
    g[7, 3] = bt_close
    set_gtk_property!(g, :column_homogeneous, true)
    set_gtk_property!(g, :column_spacing, 10)  # introduce a 5-pixel gap between columns
    set_gtk_property!(win, :border_width, 20)
    push!(win, g)

    showall(win)

end

"""
    preview_cont(obj1, obj2; <keyword arguments>)

Interactive preview of continuous signal.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object (before) - drawn in black
- `obj2::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object (after) - drawn in red
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj1))`: channel(s) to plot, default is all channels

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function preview_cont(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=NeuroAnalyzer._c(channel_n(obj1)))

    epoch_n(obj1) > 1 && @error "preview_ep() should be used for epoched object."
    _check_channels(obj1, ch)

    zoom_level = 10

    p = NeuroAnalyzer.plot(obj1, obj2, ch=ch)
    win = GtkWindow("NeuroAnalyzer", 1200, (p.attr[:size][2] + 180))
    can = GtkCanvas(Int32(p.attr[:size][1]), Int32(p.attr[:size][2]))
    g = GtkGrid()
    slider_time = GtkScale(false, obj1.time_pts[1]:obj1.time_pts[end])
    slider_ts1 = GtkScale(false, obj1.time_pts[1]:obj1.time_pts[1] + 10)
    slider_ts2 = GtkScale(false, obj1.time_pts[1]:obj1.time_pts[1] + 10)
    label_time = GtkLabel("Signal time")
    label_ts1 = GtkLabel("Segment start")
    label_ts2 = GtkLabel("Segment end")
    bt_start = GtkButton("⇤")
    bt_prev10 = GtkButton("⇐")
    bt_prev = GtkButton("←")
    bt_next = GtkButton("→")
    bt_next10 = GtkButton("⇒")
    bt_end = GtkButton("⇥")
    bt_zoomin = GtkButton("▼")
    label_zoomlevel = GtkLabel("10")
    bt_zoomout = GtkButton("▲")
    bt_help = GtkButton("🛈")
    bt_delete = GtkButton("DEL")
    bt_close = GtkButton("✖")
    separator1 = GtkLabel("")
    separator2 = GtkLabel("")
    separator3 = GtkLabel("")
    separator4 = GtkLabel("")

    @guarded draw(can) do widget
        time1 = GAccessor.value(slider_time)
        time2 = time1 + zoom_level
        time2 > obj1.time_pts[end] && (time2 = obj1.time_pts[end])
        ctx = getgc(can)
        show(io, MIME("image/png"), NeuroAnalyzer.plot(obj1, obj2, ch=ch, seg=(time1, time2)))
        img = read_from_png(io)
        set_source_surface(ctx, img, 0, 0)
        paint(ctx)
    end

    # signal_connect((w) -> draw(can), slider_time, "value-changed")
    signal_connect(slider_time, "value-changed") do widget
        time1 = GAccessor.value(slider_time)
        time2 = time1 + zoom_level
        time2 > obj1.time_pts[end] && (time2 = obj1.time_pts[end])
        Gtk.@sigatom begin
            GAccessor.range(slider_ts1, time1, time2)
            GAccessor.range(slider_ts2, time1, time2)
            GAccessor.value(slider_ts1, time1)
            GAccessor.value(slider_ts2, time1)
        end
        draw(can)
    end

    signal_connect(bt_prev, "clicked") do widget
        if GAccessor.value(slider_time) > obj1.time_pts[1] + 1
            time_current = GAccessor.value(slider_time) - 1
            Gtk.@sigatom begin
                GAccessor.value(slider_time, time_current)
                GAccessor.value(slider_ts1, time_current)
                GAccessor.value(slider_ts2, time_current)
            end
        else
            Gtk.@sigatom begin
                GAccessor.value(slider_time, obj1.time_pts[1])
                GAccessor.value(slider_ts1, obj1.time_pts[1])
                GAccessor.value(slider_ts2, obj1.time_pts[1])
            end
        end
    end

    signal_connect(bt_next, "clicked") do widget
        if GAccessor.value(slider_time) < obj1.time_pts[end] - 1
            time_current = GAccessor.value(slider_time) + 1
            Gtk.@sigatom begin
                GAccessor.value(slider_time, time_current)
                GAccessor.value(slider_ts1, time_current)
                GAccessor.value(slider_ts2, time_current)
            end
        else
            Gtk.@sigatom begin
                GAccessor.value(slider_time, round(obj1.time_pts[end] - 1))
                GAccessor.value(slider_ts1, round(obj1.time_pts[end] - 1))
                GAccessor.value(slider_ts2, round(obj1.time_pts[end] - 1))
            end
        end
    end

    signal_connect(bt_prev10, "clicked") do widget
        time_current = GAccessor.value(slider_time)
        if time_current > obj1.time_pts[1]
            time_current = time_current - 10
            Gtk.@sigatom begin
                GAccessor.value(slider_time, time_current)
                GAccessor.value(slider_ts1, time_current)
                GAccessor.value(slider_ts2, time_current)
            end
        end
    end

    signal_connect(bt_next10, "clicked") do widget
        if obj1.time_pts[end] % zoom_level == 0
            if GAccessor.value(slider_time) < obj1.time_pts[end] - 10
                time_current = GAccessor.value(slider_time) + 10
                Gtk.@sigatom begin
                    GAccessor.value(slider_time, time_current)
                    GAccessor.value(slider_ts1, time_current)
                    GAccessor.value(slider_ts2, time_current)
                end
            end
        else
            if GAccessor.value(slider_time) < obj1.time_pts[end] - (obj1.time_pts[end] % 10)
                time_current = GAccessor.value(slider_time) + 10
                Gtk.@sigatom begin
                    GAccessor.value(slider_time, time_current)
                    GAccessor.value(slider_ts1, time_current)
                    GAccessor.value(slider_ts2, time_current)
                end
            else
                time_current = obj1.time_pts[end] - (obj1.time_pts[end] % 10)
                Gtk.@sigatom begin
                    GAccessor.value(slider_time, time_current)
                    GAccessor.value(slider_ts1, time_current)
                    GAccessor.value(slider_ts2, time_current)
                end
            end
        end
    end

    signal_connect(bt_start, "clicked") do widget
        Gtk.@sigatom begin
            GAccessor.value(slider_time, obj1.time_pts[1])
            GAccessor.value(slider_ts1, obj1.time_pts[1])
            GAccessor.value(slider_ts2, obj1.time_pts[1])
        end
    end

    signal_connect(bt_end, "clicked") do widget
        if obj1.time_pts[end] % zoom_level == 0
            Gtk.@sigatom begin
                GAccessor.value(slider_time, obj1.time_pts[end] - zoom_level)
                GAccessor.value(slider_ts1, obj1.time_pts[end] - zoom_level)
                GAccessor.value(slider_ts2, obj1.time_pts[end] - zoom_level)
            end
        else
            Gtk.@sigatom begin
                GAccessor.value(slider_time, obj1.time_pts[end] - obj1.time_pts[end] % zoom_level)
                GAccessor.value(slider_ts1, obj1.time_pts[end] - obj1.time_pts[end] % zoom_level)
                GAccessor.value(slider_ts2, obj1.time_pts[end] - obj1.time_pts[end] % zoom_level)
            end
        end
    end

    signal_connect(bt_zoomin, "clicked") do widget
        time_current = GAccessor.value(slider_time)
        if zoom_level == 10
            zoom_level = 5
        elseif zoom_level == 5
            zoom_level = 1
        end
        GAccessor.text(label_zoomlevel, lpad(string(zoom_level), 2, '0'))
        if time_current + zoom_level < obj1.time_pts[end]
            Gtk.@sigatom begin
                GAccessor.range(slider_ts1, time_current, time_current + zoom_level)
                GAccessor.range(slider_ts2, time_current, time_current + zoom_level)
                GAccessor.value(slider_ts1, time_current)
                GAccessor.value(slider_ts2, time_current)
            end
            draw(can)
        else
            Gtk.@sigatom begin
                GAccessor.range(slider_ts1, obj1.time_pts[end])
                GAccessor.range(slider_ts2, obj1.time_pts[end])
                GAccessor.value(slider_ts1, time_current)
                GAccessor.value(slider_ts2, time_current)
            end
            draw(can)
        end
    end

    signal_connect(bt_zoomout, "clicked") do widget
        time_current = GAccessor.value(slider_time)
        if zoom_level == 1
            zoom_level = 5
        elseif zoom_level == 5
            zoom_level = 10
        end
        GAccessor.text(label_zoomlevel, lpad(string(zoom_level), 2, '0'))
        if time_current + zoom_level < obj1.time_pts[end]
            Gtk.@sigatom begin
                GAccessor.range(slider_ts1, time_current, time_current + zoom_level)
                GAccessor.range(slider_ts2, time_current, time_current + zoom_level)
                GAccessor.value(slider_ts1, time_current)
                GAccessor.value(slider_ts2, time_current)
            end
            draw(can)
        else
            Gtk.@sigatom begin
                GAccessor.range(slider_ts1, obj1.time_pts[end])
                GAccessor.range(slider_ts2, obj1.time_pts[end])
                GAccessor.value(slider_ts1, time_current)
                GAccessor.value(slider_ts2, time_current)
            end
            draw(can)
        end
    end

    signal_connect(bt_delete, "clicked") do widget
        time_current = GAccessor.value(slider_time)
        time1 = GAccessor.value(slider_ts1)
        time2 = GAccessor.value(slider_ts2)
        if time1 < time2
            if ask_dialog("Delete segment $time1:$time2 ?", "No", "Yes")
                trim!(obj1, seg=(time1, time2), remove_epochs=false)
                trim!(obj2, seg=(time1, time2), remove_epochs=false)
                _info("Deleted segment: $time1:$time2")
                Gtk.@sigatom begin
                    GAccessor.range(slider_time, obj1.time_pts[1], obj1.time_pts[end])
                    if obj1.time_pts[end] % zoom_level == 0
                        time_current >= (obj1.time_pts[end] - zoom_level) && (time_current = obj1.time_pts[end] - zoom_level)
                    else
                        time_current >= obj1.time_pts[end] - (obj1.time_pts[end] % zoom_level) && (time_current = obj1.time_pts[end] - (obj1.time_pts[end] % zoom_level))
                    end
                    time_current < obj1.time_pts[1] && (time_current = obj1.time_pts[1])
                    GAccessor.value(slider_time, time_current)
                    GAccessor.value(slider_ts1, time_current)
                    GAccessor.value(slider_ts2, time_current)
                end
                draw(can)
            end
        else
            info_dialog("Time segment start must be < time segment end.")
        end
    end

    signal_connect(bt_close, "clicked") do widget
        Gtk.destroy(win)
    end

    signal_connect(bt_help, "clicked") do widget
        info_dialog("Keyboard shortcuts:\n\nHOME\tgo to the signal beginning\nEND\t\tgo to the signal end\n,\t\tgo back by 1 second\n.\t\tgo forward by 1 second\n<\t\tgo back by 10 seconds\n>\t\tgo forward by 10 seconds\n\n+\t\tzoom in\n−\t\tzoom out\n\nDEL\t\tdelete current segment\n\nh\t\tthis info\nq\t\texit\n")
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        if k == 113 # q
            Gtk.destroy(win)
        elseif k == 65360 # HOME
            Gtk.@sigatom begin
                GAccessor.value(slider_time, obj1.time_pts[1])
                GAccessor.value(slider_ts1, obj1.time_pts[1])
                GAccessor.value(slider_ts2, obj1.time_pts[1])
            end
        elseif k == 65367 # END
            if obj1.time_pts[end] % 10 == 0
                Gtk.@sigatom begin
                    GAccessor.value(slider_time, obj1.time_pts[end] - 10)
                    GAccessor.value(slider_ts1, obj1.time_pts[end] - 10)
                    GAccessor.value(slider_ts2, obj1.time_pts[end] - 10)
                end
            else
                Gtk.@sigatom begin
                    GAccessor.value(slider_time, obj1.time_pts[end] - obj1.time_pts[end] % 10)
                    GAccessor.value(slider_ts1, obj1.time_pts[end] - obj1.time_pts[end] % 10)
                    GAccessor.value(slider_ts2, obj1.time_pts[end] - obj1.time_pts[end] % 10)
                end
            end
        elseif k == 104 # h
            info_dialog("Keyboard shortcuts:\n\nHOME\tgo to the signal beginning\nEND\t\tgo to the signal end\n,\t\tgo back by 1 second\n.\t\tgo forward by 1 second\n<\t\tgo back by 10 seconds\n>\t\tgo forward by 10 seconds\n\n+\t\tzoom in\n−\t\tzoom out\n\nDEL\t\tdelete current segment\n\nh\t\tthis info\nq\t\texit\n")
        elseif k == 65453 || k == 45 # -
            time_current = GAccessor.value(slider_time)
            if zoom_level == 10
                zoom_level = 5
            elseif zoom_level == 5
                zoom_level = 1
            end
            GAccessor.text(label_zoomlevel, lpad(string(zoom_level), 2, '0'))
            if time_current + zoom_level < obj1.time_pts[end]
                Gtk.@sigatom begin
                    GAccessor.range(slider_ts1, time_current, time_current + zoom_level)
                    GAccessor.range(slider_ts2, time_current, time_current + zoom_level)
                    GAccessor.value(slider_ts1, time_current)
                    GAccessor.value(slider_ts2, time_current)
                end
                draw(can)
            else
                Gtk.@sigatom begin
                    GAccessor.range(slider_ts1, obj1.time_pts[end])
                    GAccessor.range(slider_ts2, obj1.time_pts[end])
                    GAccessor.value(slider_ts1, time_current)
                    GAccessor.value(slider_ts2, time_current)
                end
                draw(can)
            end
        elseif k == 65451 || k == 43 # +
            time_current = GAccessor.value(slider_time)
            if zoom_level == 1
                zoom_level = 5
            elseif zoom_level == 5
                zoom_level = 10
            end
            GAccessor.text(label_zoomlevel, lpad(string(zoom_level), 2, '0'))
            if time_current + zoom_level < obj1.time_pts[end]
                Gtk.@sigatom begin
                    GAccessor.range(slider_ts1, time_current, time_current + zoom_level)
                    GAccessor.range(slider_ts2, time_current, time_current + zoom_level)
                    GAccessor.value(slider_ts1, time_current)
                    GAccessor.value(slider_ts2, time_current)
                end
                draw(can)
            else
                Gtk.@sigatom begin
                    GAccessor.range(slider_ts1, obj1.time_pts[end])
                    GAccessor.range(slider_ts2, obj1.time_pts[end])
                    GAccessor.value(slider_ts1, time_current)
                    GAccessor.value(slider_ts2, time_current)
                end
                draw(can)
            end
        elseif k == 44 # ,
            if GAccessor.value(slider_time) > obj1.time_pts[1] + 1
                time_current = GAccessor.value(slider_time) - 1
                Gtk.@sigatom begin
                    GAccessor.value(slider_time, time_current)
                    GAccessor.value(slider_ts1, time_current)
                    GAccessor.value(slider_ts2, time_current)
                end
            else
                Gtk.@sigatom begin
                    GAccessor.value(slider_time, obj1.time_pts[1])
                    GAccessor.value(slider_ts1, obj1.time_pts[1])
                    GAccessor.value(slider_ts2, obj1.time_pts[1])
                end
            end
        elseif k == 46 # .
            if GAccessor.value(slider_time) < obj1.time_pts[end] - 1
                time_current = GAccessor.value(slider_time) + 1
                Gtk.@sigatom begin
                    GAccessor.value(slider_time, time_current)
                    GAccessor.value(slider_ts1, time_current)
                    GAccessor.value(slider_ts2, time_current)
                end
            end
        elseif k == 60 # <
            time_current = GAccessor.value(slider_time)
            if time_current > obj1.time_pts[1]
                time_current = time_current - 10
                Gtk.@sigatom begin
                    GAccessor.value(slider_time, time_current)
                    GAccessor.value(slider_ts1, time_current)
                    GAccessor.value(slider_ts2, time_current)
                end
            end
        elseif k == 62 # >
            if obj1.time_pts[end] % 10 == 0
                if GAccessor.value(slider_time) < obj1.time_pts[end] - 10
                    time_current = GAccessor.value(slider_time) + 10
                    Gtk.@sigatom begin
                        GAccessor.value(slider_time, time_current)
                        GAccessor.value(slider_ts1, time_current)
                        GAccessor.value(slider_ts2, time_current)
                    end
                end
            else
                if GAccessor.value(slider_time) < obj1.time_pts[end] - (obj1.time_pts[end] % 10)
                    time_current = GAccessor.value(slider_time) + 10
                    Gtk.@sigatom begin
                        GAccessor.value(slider_time, time_current)
                        GAccessor.value(slider_ts1, time_current)
                        GAccessor.value(slider_ts2, time_current)
                    end
                else
                    time_current = obj1.time_pts[end] - (obj1.time_pts[end] % 10)
                    Gtk.@sigatom begin
                        GAccessor.value(slider_time, time_current)
                        GAccessor.value(slider_ts1, time_current)
                        GAccessor.value(slider_ts2, time_current)
                    end
                end
            end
        elseif k == 65535 # DEL
            time_current = GAccessor.value(slider_time)
            time1 = GAccessor.value(slider_ts1)
            time2 = GAccessor.value(slider_ts2)
            if time1 < time2
                if ask_dialog("Delete segment $time1:$time2 ?", "No", "Yes")
                    trim!(obj1, seg=(time1, time2), remove_epochs=false)
                    trim!(obj2, seg=(time1, time2), remove_epochs=false)
                    _info("Deleted segment: $time1:$time2")
                    Gtk.@sigatom begin
                        GAccessor.range(slider_time, obj1.time_pts[1], obj1.time_pts[end])
                        if obj1.time_pts[end] % zoom_level == 0
                            time_current >= (obj1.time_pts[end] - zoom_level) && (time_current = obj1.time_pts[end] - zoom_level)
                        else
                            time_current >= obj1.time_pts[end] - (obj1.time_pts[end] % zoom_level) && (time_current = obj1.time_pts[end] - (obj1.time_pts[end] % zoom_level))
                        end
                        time_current < obj1.time_pts[1] && (time_current = obj1.time_pts[1])
                        GAccessor.value(slider_time, time_current)
                        GAccessor.value(slider_ts1, time_current)
                        GAccessor.value(slider_ts2, time_current)
                    end
                    draw(can)
                end
            else
                info_dialog("Time segment start must be < time segment end.")
            end
        end
    end

    g[1:15, 1] = can
    g[1, 2] = label_time
    g[2:15, 2] = slider_time
    g[1, 3] = label_ts1
    g[2:15, 3] = slider_ts1
    g[1, 4] = label_ts2
    g[2:15, 4] = slider_ts2
    g[1, 5] = separator1
    g[2, 5] = bt_start
    g[3, 5] = bt_prev10
    g[4, 5] = bt_prev
    g[5, 5] = bt_next
    g[6, 5] = bt_next10
    g[7, 5] = bt_end
    g[8, 5] = separator2
    g[9, 5] = bt_zoomin
    g[10, 5] = bt_zoomout
    g[11, 5] = separator3
    g[12, 5] = bt_delete
    g[13, 5] = separator4
    g[14, 5] = bt_help
    g[15, 5] = bt_close
    set_gtk_property!(g, :column_homogeneous, false)
    set_gtk_property!(g, :column_spacing, 10)  # introduce a 10-pixel gap between columns
    set_gtk_property!(win, :border_width, 20)
    set_gtk_property!(win, :resizable, false)
    set_gtk_property!(slider_time, :round_digits, 1)
    set_gtk_property!(slider_ts1, :round_digits, 1)
    set_gtk_property!(slider_ts2, :round_digits, 1)
    set_gtk_property!(slider_time, :digits, 1)
    set_gtk_property!(slider_ts1, :digits, 1)
    set_gtk_property!(slider_ts2, :digits, 1)
    push!(win, g)

    showall(win)

    return nothing

end
