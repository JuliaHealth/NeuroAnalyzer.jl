export iedit
export iedit_ep
export iedit_cont

"""
    iedit(obj; <keyword arguments>)

Interactive edit of continuous or epoched signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: channel(s) to plot, default is all channels
- `mono::Bool=true`: use color or grey palette
- `zoom::Int64=5`: how many seconds are displayed in one segment

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function iedit(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=NeuroAnalyzer._c(channel_n(obj)), mono::Bool=true, zoom::Int64=5)

    if epoch_n(obj) == 1
        iedit_cont(obj, ch=ch, mono=mono, zoom=zoom)
    else
        iedit_ep(obj, ch=ch, mono=mono)
    end

end

"""
    iedit(obj1, obj2; <keyword arguments>)

Interactive edit of two continuous or epoched signals.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object (before) - drawn in black
- `obj2::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object (after) - drawn in red
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj1))`: channel(s) to plot, default is all channels
- `zoom::Int64=5`: how many seconds are displayed in one segment

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function iedit(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=NeuroAnalyzer._c(channel_n(obj1)), zoom::Int64=5)

    if epoch_n(obj1) == 1
        iedit_cont(obj1, obj2, ch=ch, zoom=zoom)
    else
        iedit_ep(obj1, obj2, ch=ch)
    end

end

"""
    iedit_ep(obj; <keyword arguments>)

Interactive edit of epoched signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: channel(s) to plot, default is all channels
- `mono::Bool=true`: use color or grey palette

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function iedit_ep(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=NeuroAnalyzer._c(channel_n(obj)), mono::Bool=true)

    @assert epoch_n(obj) > 1 "iedit_cont() should be used for continuous object."
    _check_channels(obj, ch)

    p = NeuroAnalyzer.plot(obj, ch=ch, ep=1, mono=mono, title="")
    win = GtkWindow("NeuroAnalyzer: iedit_ep()", 1200, (p.attr[:size][2] + 40))
    set_gtk_property!(win, :border_width, 20)
    set_gtk_property!(win, :resizable, false)
    can = GtkCanvas(Int32(p.attr[:size][1]), Int32(p.attr[:size][2]))
    g = GtkGrid()
    set_gtk_property!(g, :column_homogeneous, false)
    set_gtk_property!(g, :column_spacing, 10)
    set_gtk_property!(g, :row_spacing, 10)
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
    g[1:10, 1] = can
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
        ep = parse(Int64, get_gtk_property(entry_epoch, :label, String))
        ctx = getgc(can)
        show(io, MIME("image/png"), NeuroAnalyzer.plot(obj, ch=ch, ep=ep, mono=mono, title=""))
        img = read_from_png(io)
        set_source_surface(ctx, img, 0, 0)
        paint(ctx)
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

    signal_connect(bt_delete, "clicked") do widget
        ep = parse(Int64, get_gtk_property(entry_epoch, :label, String))
        if ask_dialog("Delete epoch $ep ?", "No", "Yes")
            delete_epoch!(obj, ep=ep)
            _info("Deleted epoch: $ep")
            ep = ep > 1 ? ep -= 1 : ep = 1
            set_gtk_property!(entry_epoch, :label, string(ep))
            draw(can)
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
        elseif k == 65535 # DEL
            ep = parse(Int64, get_gtk_property(entry_epoch, :label, String))
            if ask_dialog("Delete epoch $ep ?", "No", "Yes")
                delete_epoch!(obj, ep=ep)
                _info("Deleted epoch: $ep")
                ep = ep > 1 ? ep -= 1 : ep = 1
                set_gtk_property!(entry_epoch, :label, string(ep))
                draw(can)
            end
        end
    end

    return nothing

end

"""
    iedit_ep(obj1, obj2; <keyword arguments>)

Interactive edit of epoched signal.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object (before) - drawn in black
- `obj2::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object (after) - drawn in red
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj1))`: channel(s) to plot, default is all channels

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function iedit_ep(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=NeuroAnalyzer._c(channel_n(obj1)))

    @assert epoch_n(obj1) > 1 "iedit_cont() should be used for continuous object."
    _check_channels(obj1, ch)

    p = NeuroAnalyzer.plot(obj1, obj2, ch=ch)
    win = GtkWindow("NeuroAnalyzer: iedit_ep()", 1200, (p.attr[:size][2] + 40))
    set_gtk_property!(win, :border_width, 20)
    set_gtk_property!(win, :resizable, false)
    can = GtkCanvas(Int32(p.attr[:size][1]), Int32(p.attr[:size][2]))
    g = GtkGrid()
    set_gtk_property!(g, :column_homogeneous, false)
    set_gtk_property!(g, :column_spacing, 10)
    set_gtk_property!(g, :row_spacing, 10)
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
    bt_close = GtkButton("✖")
    set_gtk_property!(bt_close, :tooltip_text, "Close this window")
    g[1:8, 1] = can
    g[1, 2] = bt_start
    g[2, 2] = bt_prev
    g[3, 2] = entry_epoch
    g[4, 2] = bt_next
    g[5, 2] = bt_end
    g[6, 2] = GtkLabel("")
    g[7, 2] = bt_help
    g[8, 2] = bt_close
    push!(win, g)

    showall(win)

    @guarded draw(can) do widget
        ep = parse(Int64, get_gtk_property(entry_epoch, :label, String))
        ctx = getgc(can)
        show(io, MIME("image/png"), NeuroAnalyzer.plot(obj1, obj2, ch=ch, ep=ep, title=""))
        img = read_from_png(io)
        set_source_surface(ctx, img, 0, 0)
        paint(ctx)
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
        if ep < epoch_n(obj1)
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
            set_gtk_property!(entry_epoch, :label, string(epoch_n(obj1)))
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
                elseif v > epoch_n(obj1)
                    warn_dialog("Value must be ≤ $(epoch_n(obj1)).")
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
                    elseif v > epoch_n(obj1)
                        warn_dialog("Value must be ≤ $(epoch_n(obj1)).")
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
        info_dialog("Keyboard shortcuts:\n\nHOME\tgo to first epoch\nEND\t\tgo to last epoch\n,\t\tprevious epoch\n.\t\tnext epoch\n\nDEL\t\tdelete current epoch\n\nh\t\tthis info\nq\t\texit\n")
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        if k == 113 # q
            Gtk.destroy(win)
        elseif k == 104 # h
            info_dialog("Keyboard shortcuts:\n\nHOME\tgo to first epoch\nEND\t\tgo to last epoch\n,\t\tprevious epoch\n.\t\tnext epoch\n\nDEL\t\tdelete current epoch\n\nh\t\tthis info\nq\t\texit\n")
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
                    elseif v > epoch_n(obj1)
                        warn_dialog("Value must be ≤ $(epoch_n(obj1)).")
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
                        elseif v > epoch_n(obj1)
                            warn_dialog("Value must be ≤ $(epoch_n(obj1)).")
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
                set_gtk_property!(entry_epoch, :label, string(epoch_n(obj1)))
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
            if ep < epoch_n(obj1)
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

"""
    iedit_cont(obj; <keyword arguments>)

Interactive edit of continuous signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: channel(s) to plot, default is all channels
- `mono::Bool=true`: use color or grey palette
- `zoom::Int64=5`: how many seconds are displayed in one segment

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function iedit_cont(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=NeuroAnalyzer._c(channel_n(obj)), mono::Bool=true, zoom::Int64=5)

    @assert epoch_n(obj) == 1 "iedit_ep() should be used for epoched object."
    _check_channels(obj, ch)

    @assert zoom >= 1 "zoom must be ≥ 1."
    @assert zoom <= signal_len(obj) / sr(obj) "zoom must be ≤ $(signal_len(obj) / sr(obj))."

    p = NeuroAnalyzer.plot(obj, ch=ch, mono=mono, title="")
    win = GtkWindow("NeuroAnalyzer: iedit_cont()", 1200, (p.attr[:size][2] + 40))
    set_gtk_property!(win, :border_width, 20)
    set_gtk_property!(win, :resizable, false)
    can = GtkCanvas(Int32(p.attr[:size][1]), Int32(p.attr[:size][2]))
    g = GtkGrid()
    set_gtk_property!(g, :column_homogeneous, false)
    set_gtk_property!(g, :column_spacing, 10)
    set_gtk_property!(g, :row_spacing, 10)
    entry_time = GtkButton(string(obj.time_pts[1]))
    set_gtk_property!(entry_time, :tooltip_text, "Time position [s]")
    entry_ts1 = GtkButton(string(obj.time_pts[1]))
    set_gtk_property!(entry_ts1, :tooltip_text, "Segment start [s]")
    entry_ts2 = GtkButton(string(obj.time_pts[1]))
    set_gtk_property!(entry_ts2, :tooltip_text, "Segment end [s]")
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
    g[1:16, 1] = can
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
    g[13, 2] = bt_delete
    g[14, 2] = GtkLabel("")
    g[15, 2] = bt_help
    g[16, 2] = bt_close
    push!(win, g)

    showall(win)

    @guarded draw(can) do widget
        time1 = parse(Float64, get_gtk_property(entry_time, :label, String))
        time2 = time1 + zoom
        time2 > obj.time_pts[end] && (time2 = obj.time_pts[end])
        ts1 = parse(Float64, get_gtk_property(entry_ts1, :label, String))
        ts2 = parse(Float64, get_gtk_property(entry_ts2, :label, String))
        tseg = (ts1, ts2)
        ctx = getgc(can)
        show(io, MIME("image/png"), NeuroAnalyzer.plot(obj, ch=ch, seg=(time1, time2), s_pos=tseg, mono=mono, title=""))
        img = read_from_png(io)
        set_source_surface(ctx, img, 0, 0)
        paint(ctx)
    end

    can.mouse.button1press = @guarded (widget, event) -> begin
        time_current = parse(Float64, get_gtk_property(entry_time, :label, String))
        x_pos = event.x
        x_pos < 52 && (x_pos = 52)
        x_pos > 1182 && (x_pos = 1182)
        if time_current + zoom < obj.time_pts[end]
            ts1 = time_current + round((x_pos - 52) / (1130 / zoom), digits=3)
        else
            ts1 = time_current + round((x_pos - 52) / (1130 / (obj.time_pts[end] - time_current)), digits=3)
        end
        Gtk.@sigatom begin
            set_gtk_property!(entry_ts1, :label, string(round(ts1, digits=3)))
        end
        draw(can)
    end

    can.mouse.button3press = @guarded (widget, event) -> begin
        time_current = parse(Float64, get_gtk_property(entry_time, :label, String))
        x_pos = event.x
        x_pos < 52 && (x_pos = 52)
        x_pos > 1182 && (x_pos = 1182)
        if time_current + zoom < obj.time_pts[end]
            ts2 = time_current + round((x_pos - 52) / (1130 / zoom), digits=3)
        else
            ts2 = time_current + round((x_pos - 52) / (1130 / (obj.time_pts[end] - time_current)), digits=3)
        end
        Gtk.@sigatom begin
            set_gtk_property!(entry_ts2, :label, string(round(ts2, digits=3)))
        end
        draw(can)
    end

    signal_connect(bt_prev, "clicked") do widget
        time_current = parse(Float64, get_gtk_property(entry_time, :label, String))
        if time_current >= obj.time_pts[1] + 1
            time_current -= 1
            Gtk.@sigatom begin
                set_gtk_property!(entry_time, :label, string(time_current))
                set_gtk_property!(entry_ts1, :label, string(time_current))
                set_gtk_property!(entry_ts2, :label, string(time_current))
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
                set_gtk_property!(entry_ts1, :label, string(time_current))
                set_gtk_property!(entry_ts2, :label, string(time_current))
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
                set_gtk_property!(entry_ts1, :label, string(time_current))
                set_gtk_property!(entry_ts2, :label, string(time_current))
            end
        else
            time_current = obj.time_pts[end] - zoom
            Gtk.@sigatom begin
                set_gtk_property!(entry_time, :label, string(time_current))
                set_gtk_property!(entry_ts1, :label, string(time_current))
                set_gtk_property!(entry_ts2, :label, string(time_current))
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
                set_gtk_property!(entry_ts1, :label, string(time_current))
                set_gtk_property!(entry_ts2, :label, string(time_current))
            end
        else
            time_current = obj.time_pts[end] - zoom
            Gtk.@sigatom begin
                set_gtk_property!(entry_time, :label, string(time_current))
                set_gtk_property!(entry_ts1, :label, string(time_current))
                set_gtk_property!(entry_ts2, :label, string(time_current))
            end
        end
        draw(can)
    end

    signal_connect(bt_start, "clicked") do widget
        Gtk.@sigatom begin
            set_gtk_property!(entry_time, :label, string(obj.time_pts[1]))
            set_gtk_property!(entry_ts1, :label, string(obj.time_pts[1]))
            set_gtk_property!(entry_ts1, :label, string(obj.time_pts[1]))
        end
        draw(can)
    end

    signal_connect(bt_end, "clicked") do widget
        time_current = obj.time_pts[end] - zoom
        Gtk.@sigatom begin
            set_gtk_property!(entry_time, :label, string(time_current))
            set_gtk_property!(entry_ts1, :label, string(time_current))
            set_gtk_property!(entry_ts2, :label, string(time_current))
        end
        draw(can)
    end

    signal_connect(bt_delete, "clicked") do widget
        time_current = parse(Float64, get_gtk_property(entry_time, :label, String))
        time1 = parse(Float64, get_gtk_property(entry_ts1, :label, String))
        time2 = parse(Float64, get_gtk_property(entry_ts2, :label, String))
        if time1 < time2
            if ask_dialog("Delete segment $time1:$time2 ?", "No", "Yes")
                trim!(obj, seg=(time1, time2), remove_epochs=false)
                _info("Deleted segment: $time1:$time2")
                if time1 == time_current && time2 > obj.time_pts[end]
                    Gtk.@sigatom begin
                        time_current = obj.time_pts[end] - zoom
                        time_current < obj.time_pts[1] && (time_current = obj.time_pts[1])
                    end
                else
                    Gtk.@sigatom begin
                        if obj.time_pts[end] % zoom == 0
                            time_current >= (obj.time_pts[end] - zoom) && (time_current = obj.time_pts[end] - zoom)
                        else
                            time_current >= obj.time_pts[end] - (obj.time_pts[end] % zoom) && (time_current = obj.time_pts[end] - (obj.time_pts[end] % zoom))
                        end
                        time_current < obj.time_pts[1] && (time_current = obj.time_pts[1])
                    end
                end
                set_gtk_property!(entry_time, :label, string(time_current))
                set_gtk_property!(entry_ts1, :label, string(time_current))
                set_gtk_property!(entry_ts2, :label, string(time_current))
                draw(can)
            end
        elseif time1 > time2
            warn_dialog("Cannot delete!\nSegment start is larger than segment end.")
        else
            warn_dialog("Cannot delete!\nSegment start must be different from segment end.")
        end
    end

    signal_connect(entry_time, "clicked") do widget
        value = parse(Float64, get_gtk_property(entry_time, :label, String))
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
            if _check_sfloat(value_s)
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
                if _check_sfloat(value_s)
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

    signal_connect(entry_ts1, "clicked") do widget
        value = parse(Float64, get_gtk_property(entry_ts1, :label, String))
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
            if _check_sfloat(value_s)
                v = parse(Float64, value_s)
                current_time = parse(Float64, get_gtk_property(entry_time, :label, String))
                if v < current_time
                    warn_dialog("Value must be ≥ $(current_time).")
                elseif v > current_time + zoom
                    warn_dialog("Value must be ≤ $(current_time + zoom).")
                else
                    value_s = string(obj.time_pts[vsearch(parse(Float64, value_s), obj.time_pts)])
                    set_gtk_property!(entry_ts1, :label, value_s)
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
                if _check_sfloat(value_s)
                    v = parse(Float64, value_s)
                    current_time = parse(Float64, get_gtk_property(entry_time, :label, String))
                    if v < current_time
                        warn_dialog("Value must be ≥ $(current_time).")
                    elseif v > current_time + zoom
                        warn_dialog("Value must be ≤ $(current_time + zoom).")
                    else
                        value_s = string(obj.time_pts[vsearch(parse(Float64, value_s), obj.time_pts)])
                        set_gtk_property!(entry_ts1, :label, value_s)
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

    signal_connect(entry_ts2, "clicked") do widget
        value = parse(Float64, get_gtk_property(entry_ts2, :label, String))
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
            if _check_sfloat(value_s)
                v = parse(Float64, value_s)
                current_time = parse(Float64, get_gtk_property(entry_time, :label, String))
                if v < current_time
                    warn_dialog("Value must be ≥ $(current_time).")
                elseif v > current_time + zoom
                    warn_dialog("Value must be ≤ $(current_time + zoom).")
                else
                    value_s = string(obj.time_pts[vsearch(parse(Float64, value_s), obj.time_pts)])
                    set_gtk_property!(entry_ts2, :label, value_s)
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
                if _check_sfloat(value_s)
                    v = parse(Float64, value_s)
                    current_time = parse(Float64, get_gtk_property(entry_time, :label, String))
                    if v < current_time
                        warn_dialog("Value must be ≥ $(current_time).")
                    elseif v > current_time + zoom
                        warn_dialog("Value must be ≤ $(current_time + zoom).")
                    else
                        value_s = string(obj.time_pts[vsearch(parse(Float64, value_s), obj.time_pts)])
                        set_gtk_property!(entry_ts2, :label, value_s)
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
        info_dialog("Keyboard shortcuts:\n\ng\t\tgo to time point\nHOME\tgo to the signal beginning\nEND\t\tgo to the signal end\n,\t\tgo back by 1 second\n.\t\tgo forward by 1 second\n<\t\tgo back by $zoom seconds\n>\t\tgo forward by $zoom seconds\n\n[\t\tedit segment start point\n]\t\tedit segment end point\nDEL\t\tdelete current segment\n\nh\t\tthis info\nq\t\texit\n")
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        if k == 113 # q
            Gtk.destroy(win)
        elseif k == 104 # h
            info_dialog("Keyboard shortcuts:\n\ng\t\tgo to time point\nHOME\tgo to the signal beginning\nEND\t\tgo to the signal end\n,\t\tgo back by 1 second\n.\t\tgo forward by 1 second\n<\t\tgo back by $zoom seconds\n>\t\tgo forward by $zoom seconds\n\n[\t\tedit segment start point\n]\t\tedit segment end point\nDEL\t\tdelete current segment\n\nh\t\tthis info\nq\t\texit\n")
        elseif k == 103 # g
            value = parse(Float64, get_gtk_property(entry_time, :label, String))
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
                if _check_sfloat(value_s)
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
                    if _check_sfloat(value_s)
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
        elseif k == 91 # [
            value = parse(Float64, get_gtk_property(entry_ts1, :label, String))
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
                if _check_sfloat(value_s)
                    v = parse(Float64, value_s)
                    current_time = parse(Float64, get_gtk_property(entry_time, :label, String))
                    if v < current_time
                        warn_dialog("Value must be ≥ $(current_time).")
                    elseif v > current_time + zoom
                        warn_dialog("Value must be ≤ $(current_time + zoom).")
                    else
                        value_s = string(obj.time_pts[vsearch(parse(Float64, value_s), obj.time_pts)])
                        set_gtk_property!(entry_ts1, :label, value_s)
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
                    if _check_sfloat(value_s)
                        v = parse(Float64, value_s)
                        current_time = parse(Float64, get_gtk_property(entry_time, :label, String))
                        if v < current_time
                            warn_dialog("Value must be ≥ $(current_time).")
                        elseif v > current_time + zoom
                            warn_dialog("Value must be ≤ $(current_time + zoom).")
                        else
                            value_s = string(obj.time_pts[vsearch(parse(Float64, value_s), obj.time_pts)])
                            set_gtk_property!(entry_ts1, :label, value_s)
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
        elseif k == 93 # ]
            value = parse(Float64, get_gtk_property(entry_ts2, :label, String))
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
                if _check_sfloat(value_s)
                    v = parse(Float64, value_s)
                    current_time = parse(Float64, get_gtk_property(entry_time, :label, String))
                    if v < current_time
                        warn_dialog("Value must be ≥ $(current_time).")
                    elseif v > current_time + zoom
                        warn_dialog("Value must be ≤ $(current_time + zoom).")
                    else
                        value_s = string(obj.time_pts[vsearch(parse(Float64, value_s), obj.time_pts)])
                        set_gtk_property!(entry_ts2, :label, value_s)
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
                    if _check_sfloat(value_s)
                        v = parse(Float64, value_s)
                        current_time = parse(Float64, get_gtk_property(entry_time, :label, String))
                        if v < current_time
                            warn_dialog("Value must be ≥ $(current_time).")
                        elseif v > current_time + zoom
                            warn_dialog("Value must be ≤ $(current_time + zoom).")
                        else
                            value_s = string(obj.time_pts[vsearch(parse(Float64, value_s), obj.time_pts)])
                            set_gtk_property!(entry_ts2, :label, value_s)
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
                set_gtk_property!(entry_ts1, :label, string(obj.time_pts[1]))
                set_gtk_property!(entry_ts1, :label, string(obj.time_pts[1]))
            end
            draw(can)
        elseif k == 65367 # END
            time_current = obj.time_pts[end] - zoom
            Gtk.@sigatom begin
                set_gtk_property!(entry_time, :label, string(time_current))
                set_gtk_property!(entry_ts1, :label, string(time_current))
                set_gtk_property!(entry_ts2, :label, string(time_current))
            end
            draw(can)
        elseif k == 44 # ,
            time_current = parse(Float64, get_gtk_property(entry_time, :label, String))
            if time_current >= obj.time_pts[1] + 1
                time_current -= 1
                Gtk.@sigatom begin
                    set_gtk_property!(entry_time, :label, string(time_current))
                    set_gtk_property!(entry_ts1, :label, string(time_current))
                    set_gtk_property!(entry_ts2, :label, string(time_current))
                end
            end
            draw(can)
        elseif k == 60 # <
            time_current = parse(Float64, get_gtk_property(entry_time, :label, String))
            if time_current >= obj.time_pts[1] + zoom
                time_current = time_current - zoom
                Gtk.@sigatom begin
                    set_gtk_property!(entry_time, :label, string(time_current))
                    set_gtk_property!(entry_ts1, :label, string(time_current))
                    set_gtk_property!(entry_ts2, :label, string(time_current))
                end
            end
            draw(can)
        elseif k == 46 # .
            time_current = parse(Float64, get_gtk_property(entry_time, :label, String))
            if time_current < obj.time_pts[end] - zoom
                time_current += 1
                Gtk.@sigatom begin
                    set_gtk_property!(entry_time, :label, string(time_current))
                    set_gtk_property!(entry_ts1, :label, string(time_current))
                    set_gtk_property!(entry_ts2, :label, string(time_current))
                end
            else
                time_current = obj.time_pts[end] - zoom
                Gtk.@sigatom begin
                    set_gtk_property!(entry_time, :label, string(time_current))
                    set_gtk_property!(entry_ts1, :label, string(time_current))
                    set_gtk_property!(entry_ts2, :label, string(time_current))
                end
            end
            draw(can)
        elseif k == 62 # >
            time_current = parse(Float64, get_gtk_property(entry_time, :label, String))
            if time_current < obj.time_pts[end] - zoom
                time_current += zoom
                Gtk.@sigatom begin
                    set_gtk_property!(entry_time, :label, string(time_current))
                    set_gtk_property!(entry_ts1, :label, string(time_current))
                    set_gtk_property!(entry_ts2, :label, string(time_current))
                end
            else
                time_current = obj.time_pts[end] - zoom
                Gtk.@sigatom begin
                    set_gtk_property!(entry_time, :label, string(time_current))
                    set_gtk_property!(entry_ts1, :label, string(time_current))
                    set_gtk_property!(entry_ts2, :label, string(time_current))
                end
            end
            draw(can)
        elseif k == 65535 # DEL
            time_current = parse(Float64, get_gtk_property(entry_time, :label, String))
            time1 = parse(Float64, get_gtk_property(entry_ts1, :label, String))
            time2 = parse(Float64, get_gtk_property(entry_ts2, :label, String))
            if time1 < time2
                if ask_dialog("Delete segment $time1:$time2 ?", "No", "Yes")
                    trim!(obj, seg=(time1, time2), remove_epochs=false)
                    _info("Deleted segment: $time1:$time2")
                    if time1 == time_current && time2 > obj.time_pts[end]
                        Gtk.@sigatom begin
                            time_current = obj.time_pts[end] - zoom
                            time_current < obj.time_pts[1] && (time_current = obj.time_pts[1])
                        end
                    else
                        Gtk.@sigatom begin
                            if obj.time_pts[end] % zoom == 0
                                time_current >= (obj.time_pts[end] - zoom) && (time_current = obj.time_pts[end] - zoom)
                            else
                                time_current >= obj.time_pts[end] - (obj.time_pts[end] % zoom) && (time_current = obj.time_pts[end] - (obj.time_pts[end] % zoom))
                            end
                            time_current < obj.time_pts[1] && (time_current = obj.time_pts[1])
                        end
                    end
                    set_gtk_property!(entry_time, :label, string(time_current))
                    set_gtk_property!(entry_ts1, :label, string(time_current))
                    set_gtk_property!(entry_ts2, :label, string(time_current))
                    draw(can)
                end
            elseif time1 > time2
                warn_dialog("Cannot delete!\nSegment start is larger than segment end.")
            else
                warn_dialog("Cannot delete!\nSegment start must be different from segment end.")
            end
        end
    end

    return nothing

end

"""
    iedit_cont(obj1, obj2; <keyword arguments>)

Interactive preview of two continuous signal.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object (before) - drawn in black
- `obj2::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object (after) - drawn in red
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj1))`: channel(s) to plot, default is all channels
- `zoom::Int64=5`: how many seconds are displayed in one segment

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function iedit_cont(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=NeuroAnalyzer._c(channel_n(obj1)), zoom::Int64=5)

    @assert epoch_n(obj1) == 1 "iedit_ep() should be used for epoched object."
    @assert epoch_n(obj2) == 1 "iedit_ep() should be used for epoched object."
    _check_channels(obj1, ch)
    _check_channels(obj2, ch)

    @assert zoom >= 1 "zoom must be ≥ 1."
    @assert zoom <= signal_len(obj1) / sr(obj1) "zoom must be ≤ $(signal_len(obj1) / sr(obj1))."
    @assert zoom <= signal_len(obj2) / sr(obj2) "zoom must be ≤ $(signal_len(obj2) / sr(obj2))."

    p = NeuroAnalyzer.plot(obj1, obj2, ch=ch, title="")
    win = GtkWindow("NeuroAnalyzer: iedit_cont()", 1200, (p.attr[:size][2] + 40))
    set_gtk_property!(win, :border_width, 20)
    set_gtk_property!(win, :resizable, false)
    can = GtkCanvas(Int32(p.attr[:size][1]), Int32(p.attr[:size][2]))
    g = GtkGrid()
    set_gtk_property!(g, :column_homogeneous, false)
    set_gtk_property!(g, :column_spacing, 10)
    set_gtk_property!(g, :row_spacing, 10)
    entry_time = GtkButton(string(obj1.time_pts[1]))
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
    bt_close = GtkButton("✖")
    set_gtk_property!(bt_close, :tooltip_text, "Close this window")
    g[1:10, 1] = can
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
        time1 = parse(Float64, get_gtk_property(entry_time, :label, String))
        time2 = time1 + zoom
        time2 > obj1.time_pts[end] && (time2 = obj1.time_pts[end])
        ctx = getgc(can)
        show(io, MIME("image/png"), NeuroAnalyzer.plot(obj1, obj2, ch=ch, seg=(time1, time2), title=""))
        img = read_from_png(io)
        set_source_surface(ctx, img, 0, 0)
        paint(ctx)
    end

    signal_connect(bt_prev, "clicked") do widget
        time_current = parse(Float64, get_gtk_property(entry_time, :label, String))
        if time_current >= obj1.time_pts[1] + 1
            time_current -= 1
            Gtk.@sigatom begin
                set_gtk_property!(entry_time, :label, string(time_current))
            end
        end
        draw(can)
    end

    signal_connect(bt_prev5, "clicked") do widget
        time_current = parse(Float64, get_gtk_property(entry_time, :label, String))
        if time_current >= obj1.time_pts[1] + zoom
            time_current = time_current - zoom
            Gtk.@sigatom begin
                set_gtk_property!(entry_time, :label, string(time_current))
            end
        end
        draw(can)
    end

    signal_connect(bt_next, "clicked") do widget
        time_current = parse(Float64, get_gtk_property(entry_time, :label, String))
        if time_current < obj1.time_pts[end] - zoom
            time_current += 1
            Gtk.@sigatom begin
                set_gtk_property!(entry_time, :label, string(time_current))
            end
        else
            time_current = obj1.time_pts[end] - zoom
            Gtk.@sigatom begin
                set_gtk_property!(entry_time, :label, string(time_current))
            end
        end
        draw(can)
    end

    signal_connect(bt_next5, "clicked") do widget
        time_current = parse(Float64, get_gtk_property(entry_time, :label, String))
        if time_current < obj1.time_pts[end] - zoom
            time_current += zoom
            Gtk.@sigatom begin
                set_gtk_property!(entry_time, :label, string(time_current))
            end
        else
            time_current = obj1.time_pts[end] - zoom
            Gtk.@sigatom begin
                set_gtk_property!(entry_time, :label, string(time_current))
            end
        end
        draw(can)
    end

    signal_connect(bt_start, "clicked") do widget
        Gtk.@sigatom begin
            set_gtk_property!(entry_time, :label, string(obj1.time_pts[1]))
        end
        draw(can)
    end

    signal_connect(bt_end, "clicked") do widget
        time_current = obj1.time_pts[end] - zoom
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
        set_gtk_property!(d_g, :column_spacing, 10)
        set_gtk_property!(d_g, :row_spacing, 10)
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
            if _check_sfloat(value_s)
                v = parse(Float64, value_s)
                if v < obj1.time_pts[1]
                    warn_dialog("Value must be ≥ $(obj1.time_pts[1]).")
                elseif v > obj1.time_pts[end] - zoom
                    warn_dialog("Value must be ≤ $(obj1.time_pts[end] - zoom).")
                else
                    value_s = string(obj1.time_pts[vsearch(parse(Float64, value_s), obj1.time_pts)])
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
                if _check_sfloat(value_s)
                    v = parse(Float64, value_s)
                    if v < obj1.time_pts[1]
                        warn_dialog("Value must be ≥ $(obj.time_pts[1]).")
                    elseif v > obj1.time_pts[end] - zoom
                        warn_dialog("Value must be ≤ $(obj1.time_pts[end] - zoom).")
                    else
                        value_s = string(obj1.time_pts[vsearch(parse(Float64, value_s), obj1.time_pts)])
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
        info_dialog("Keyboard shortcuts:\n\ng\t\tgo to time point\nHOME\tgo to the signal beginning\nEND\t\tgo to the signal end\n,\t\tgo back by 1 second\n.\t\tgo forward by 1 second\n<\t\tgo back by $zoom seconds\n>\t\tgo forward by $zoom seconds\n\n[\t\tedit segment start point\n]\t\tedit segment end point\nDEL\t\tdelete current segment\n\nh\t\tthis info\nq\t\texit\n")
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        if k == 113 # q
            Gtk.destroy(win)
        elseif k == 104 # h
            info_dialog("Keyboard shortcuts:\n\ng\t\tgo to time point\nHOME\tgo to the signal beginning\nEND\t\tgo to the signal end\n,\t\tgo back by 1 second\n.\t\tgo forward by 1 second\n<\t\tgo back by $zoom seconds\n>\t\tgo forward by $zoom seconds\n\n[\t\tedit segment start point\n]\t\tedit segment end point\nDEL\t\tdelete current segment\n\nh\t\tthis info\nq\t\texit\n")
        elseif k == 103 # g
            value = parse(Float64, get_gtk_property(entry_time, :label, String))
            d_w = GtkWindow("Enter value", 200, 100)
            set_gtk_property!(d_w, :border_width, 20)
            set_gtk_property!(d_w, :resizable, true)
            d_g = GtkGrid()
            set_gtk_property!(d_g, :column_homogeneous, true)
            set_gtk_property!(d_g, :column_spacing, 10)
            set_gtk_property!(d_g, :row_spacing, 10)
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
                if _check_sfloat(value_s)
                    v = parse(Float64, value_s)
                    if v < obj1.time_pts[1]
                        warn_dialog("Value must be ≥ $(obj1.time_pts[1]).")
                    elseif v > obj1.time_pts[end] - zoom
                        warn_dialog("Value must be ≤ $(obj1.time_pts[end] - zoom).")
                    else
                        value_s = string(obj1.time_pts[vsearch(parse(Float64, value_s), obj1.time_pts)])
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
                    if _check_sfloat(value_s)
                        v = parse(Float64, value_s)
                        if v < obj1.time_pts[1]
                            warn_dialog("Value must be ≥ $(obj1.time_pts[1]).")
                        elseif v > obj1.time_pts[end] - zoom
                            warn_dialog("Value must be ≤ $(obj1.time_pts[end] - zoom).")
                        else
                            value_s = string(obj1.time_pts[vsearch(parse(Float64, value_s), obj1.time_pts)])
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
                set_gtk_property!(entry_time, :label, string(obj1.time_pts[1]))
            end
            draw(can)
        elseif k == 65367 # END
            time_current = obj1.time_pts[end] - zoom
            Gtk.@sigatom begin
                set_gtk_property!(entry_time, :label, string(time_current))
            end
            draw(can)
        elseif k == 44 # ,
            time_current = parse(Float64, get_gtk_property(entry_time, :label, String))
            if time_current >= obj1.time_pts[1] + 1
                time_current -= 1
                Gtk.@sigatom begin
                    set_gtk_property!(entry_time, :label, string(time_current))
                end
            end
            draw(can)
        elseif k == 60 # <
            time_current = parse(Float64, get_gtk_property(entry_time, :label, String))
            if time_current >= obj1.time_pts[1] + zoom
                time_current = time_current - zoom
                Gtk.@sigatom begin
                    set_gtk_property!(entry_time, :label, string(time_current))
                end
            end
            draw(can)
        elseif k == 46 # .
            time_current = parse(Float64, get_gtk_property(entry_time, :label, String))
            if time_current < obj1.time_pts[end] - zoom
                time_current += 1
                Gtk.@sigatom begin
                    set_gtk_property!(entry_time, :label, string(time_current))
                end
            else
                time_current = obj1.time_pts[end] - zoom
                Gtk.@sigatom begin
                    set_gtk_property!(entry_time, :label, string(time_current))
                end
            end
            draw(can)
        elseif k == 62 # >
            time_current = parse(Float64, get_gtk_property(entry_time, :label, String))
            if time_current < obj1.time_pts[end] - zoom
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
