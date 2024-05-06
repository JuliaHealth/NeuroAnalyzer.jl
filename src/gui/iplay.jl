export iplay

"""
    iplay(obj; <keyword arguments>)

Interactive play channel signal as audio

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Int64=1`: channel number
- `ep::Int64=1`: epoch number
- `mono::Bool=true`: play mono or stereo
"""
function iplay(obj::NeuroAnalyzer.NEURO; ch::Int64=1, ep::Int64=1, mono::Bool=true, maxvol::Bool=false)

    _check_channels(obj, ch)
    _check_epochs(obj, ep)

    p = NeuroAnalyzer.plot(obj, ch=ch, ep=ep, mono=true, title="", scale=false)

    win = GtkWindow("NeuroAnalyzer: iplay()", Int32(p.attr[:size][1]), Int32(p.attr[:size][2]) + 40)

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
    set_gtk_property!(bt_chup, :tooltip_text, "Previous channel")
    ch == 1 && set_gtk_property!(bt_chup, :sensitive, false)
    bt_chdown = GtkButton("▽")
    set_gtk_property!(bt_chdown, :tooltip_text, "Next channel")
    ch == nchannels(obj) && set_gtk_property!(bt_chdown, :sensitive, false)
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
    bt_play = GtkButton("▷")
    set_gtk_property!(bt_play, :tooltip_text, "Play the signal")
    bt_close = GtkButton("✖")
    set_gtk_property!(bt_close, :tooltip_text, "Close this window")
    lab_ch = GtkLabel("$(string(ch))")
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
    g[11, 2] = bt_play
    g[12, 2] = GtkLabel("")
    g[13, 2] = bt_help
    g[14, 2] = bt_close
    push!(win, g)

    showall(win)

    @guarded draw(can) do widget
        ep = get_gtk_property(entry_epoch, :value, Int64)
        ctx = getgc(can)
        show(NeuroAnalyzer.io, MIME("image/png"), NeuroAnalyzer.plot(obj,
                                                                        ch=ch,
                                                                        ep=ep,
                                                                        mono=true,
                                                                        scale=false,
                                                                        title=""))
        img = read_from_png(NeuroAnalyzer.io)
        set_source_surface(ctx, img, 0, 0)
        paint(ctx)
    end

    can.mouse.scroll = @guarded (widget, event) -> begin
        if event.direction == 1 # down
            if ch < nchannels(obj)
                ch += 1
                Gtk.@sigatom begin
                    set_gtk_property!(lab_ch, :label, "$(string(ch))")
                    ch > 1 && set_gtk_property!(bt_chup, :sensitive, true)
                    ch == nchannels(obj) && set_gtk_property!(bt_chdown, :sensitive, false)
                end
                draw(can)
            end
        elseif event.direction == 0 # up
            if ch > 1
                ch -= 1
                Gtk.@sigatom begin
                    set_gtk_property!(lab_ch, :label, "$(string(ch))")
                    ch == 1 && set_gtk_property!(bt_chup, :sensitive, false)
                    ch < nchannels(obj) && set_gtk_property!(bt_chdown, :sensitive, true)
                end
                draw(can)
            end
        end
    end

    signal_connect(bt_chdown, "clicked") do widget
        if ch < nchannels(obj)
            ch += 1
            Gtk.@sigatom begin
                set_gtk_property!(lab_ch, :label, "$(string(ch))")
                ch > 1 && set_gtk_property!(bt_chup, :sensitive, true)
                ch == nchannels(obj) && set_gtk_property!(bt_chdown, :sensitive, false)
            end
            draw(can)
        end
    end

    signal_connect(bt_chup, "clicked") do widget
        if ch > 1
            ch -= 1
            Gtk.@sigatom begin
                set_gtk_property!(lab_ch, :label, "$(string(ch))")
                ch == 1 && set_gtk_property!(bt_chup, :sensitive, false)
                ch < nchannels(obj) && set_gtk_property!(bt_chdown, :sensitive, true)
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

    signal_connect(bt_play, "clicked") do widget
        ep = get_gtk_property(entry_epoch, :value, Int64)
        s = @views obj.data[ch, :, ep]
        fs = sr(obj)
        mono == false && (s = [s s])
        wavplay(s, fs)
    end

    signal_connect(bt_close, "clicked") do widget
        Gtk.destroy(win)
    end

    signal_connect(bt_help, "clicked") do widgete
        info_dialog("Keyboard shortcuts:\n\nctrl-b\tslide channels up\nctrl-n\tslide channels down\n\nctrl-a\tgo to first epoch\nctrl-s\tgo to last epoch\nctrl-z\tprevious epoch\nctrl-x\tnext epoch\n\nctrl-SPC\tplay current epoch\n\nctrl-h\tthis info\nctrl-q\texit\n")
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        s = event.state
        if s == 4
            if k == 113 # q
                Gtk.destroy(win)
            elseif k == 104 # h
                info_dialog("Keyboard shortcuts:\n\nctrl-b\tslide channels up\nctrl-n\tslide channels down\n\nctrl-a\tgo to first epoch\nctrl-s\tgo to last epoch\nctrl-z\tprevious epoch\nctrl-x\tnext epoch\n\nctrl-SPC\tplay current epoch\n\nctrl-h\tthis info\nctrl-q\texit\n")
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
            elseif k == 115 # s
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
            elseif k == 32 # space
                ep = get_gtk_property(entry_epoch, :value, Int64)
                s = @views obj.data[ch, :, ep]
                fs = sr(obj)
                if maxvol == true
                    s = normalize_minmax(s)
                end
                if mono == false
                    s = [s s]
                end
                wavplay(s, fs)
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