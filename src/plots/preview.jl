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
    bt_help = GtkButton("🛈")
    bt_delete = GtkButton("DEL")
    bt_close = GtkButton("✖")

    @guarded draw(can) do widget
        time1 = GAccessor.value(slider_time)
        time2 = time1 + 10
        time2 > obj.time_pts[end] && (time2 = obj.time_pts[end])
        Gtk.@sigatom begin
            GAccessor.range(slider_ts1, 0, 10)
            GAccessor.range(slider_ts2, 0, 10)
            GAccessor.range(slider_ts1, time1, time2)
            GAccessor.range(slider_ts2, time1, time2)
        end
        ctx = getgc(can)
        show(io, MIME("image/png"), NeuroAnalyzer.plot(obj, ch=ch, seg=(time1, time2), mono=mono))
        img = read_from_png(io)
        set_source_surface(ctx, img, 0, 0)
        paint(ctx)
    end

    signal_connect((w) -> draw(can), slider_time, "value-changed")

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
    end

    signal_connect(bt_start, "clicked") do widget
        Gtk.@sigatom begin
            GAccessor.value(slider_time, obj.time_pts[1])
            GAccessor.value(slider_ts1, obj.time_pts[1])
            GAccessor.value(slider_ts2, obj.time_pts[1])
        end
    end

    signal_connect(bt_end, "clicked") do widget
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
                    if obj.time_pts[end] % 10 == 0
                        time_current >= (obj.time_pts[end] - 10) && (time_current = obj.time_pts[end] - 10)
                    else
                        time_current >= obj.time_pts[end] - (obj.time_pts[end] % 10) && (time_current = obj.time_pts[end] - (obj.time_pts[end] % 10))
                    end
                    time_current < obj.time_pts[1] && (time_current = obj.time_pts[1])
                    GAccessor.value(slider_time, time_current)
                    GAccessor.value(slider_ts1, time_current)
                    GAccessor.value(slider_ts2, time_current)
                    ctx = getgc(can)
                    if time_current + 10 < obj.time_pts[end]
                        show(io, MIME("image/png"), NeuroAnalyzer.plot(obj, ch=ch, seg=(time_current, time_current + 10), mono=mono))
                    else
                        show(io, MIME("image/png"), NeuroAnalyzer.plot(obj, ch=ch, seg=(time_current, time_current + (obj.time_pts[end] % 10)), mono=mono))
                    end
                    img = read_from_png(io)
                    set_source_surface(ctx, img, 0, 0)
                    paint(ctx)
                end
            end
        else
            info_dialog("Time segment start must be < time segment end.")
        end
    end

    signal_connect(bt_close, "clicked") do widget
        Gtk.destroy(win)
    end

    signal_connect(bt_help, "clicked") do widget
        info_dialog("Keyboard shortcuts:\n\nHOME\tgo to the signal beginning\nEND\t\tgo to the signal end\n,\t\tgo back by 1 second\n.\t\tgo forward by 1 second\n<\t\tgo back by 10 seconds\n>\t\tgo forward by 10 seconds\n\nDEL\t\tdelete current segment\n\nh\t\tthis info\nq\t\texit\n")
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
            info_dialog("Keyboard shortcuts:\n\nHOME\tgo to the signal beginning\nEND\t\tgo to the signal end\n,\t\tgo back by 1 second\n.\t\tgo forward by 1 second\n<\t\tgo back by 10 seconds\n>\t\tgo forward by 10 seconds\n\nDEL\t\tdelete current segment\n\nh\t\tthis info\nq\t\texit\n")
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
                        time_current > obj.time_pts[end] && (time_current = obj.time_pts[end] - 10)
                        time_current < obj.time_pts[1] && (time_current = obj.time_pts[1])
                        GAccessor.value(slider_time, time_current)
                        GAccessor.value(slider_ts1, time_current)
                        GAccessor.value(slider_ts2, time_current)
                        ctx = getgc(can)
                        show(io, MIME("image/png"), NeuroAnalyzer.plot(obj, ch=ch, seg=(time_current, time_current+10), mono=mono))
                        img = read_from_png(io)
                        set_source_surface(ctx, img, 0, 0)
                        paint(ctx)
                    end
                end
            else
                info_dialog("Time segment start must be < time segment end.")
            end
        end
    end

    g[1:9, 1] = can         # spans all columns
    g[1, 2] = label_time
    g[2:9, 2] = slider_time # spans all columns
    g[1, 3] = label_ts1
    g[2:9, 3] = slider_ts1  # spans all columns
    g[1, 4] = label_ts2
    g[2:9, 4] = slider_ts2  # spans all columns
    g[1, 5] = bt_start      # Cartesian coordinates, g[x,y]
    g[2, 5] = bt_prev10
    g[3, 5] = bt_prev
    g[4, 5] = bt_next
    g[5, 5] = bt_next10
    g[6, 5] = bt_end
    g[7, 5] = bt_delete
    g[8, 5] = bt_help
    g[9, 5] = bt_close
    set_gtk_property!(g, :column_homogeneous, true)
    set_gtk_property!(g, :column_spacing, 10)  # introduce a 5-pixel gap between columns
    set_gtk_property!(win, :border_width, 20)
    push!(win, g)

    showall(win)
    # show(can)

end
