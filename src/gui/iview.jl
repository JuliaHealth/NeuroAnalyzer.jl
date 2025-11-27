export iview
export iview_ep

"""
    iview(obj; <keyword arguments>)

Interactive view of continuous signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `mch::Bool=true`: draw multichannel signal (up to 20 channels in one plot)
- `zoom::Real=10`: how many seconds are displayed in one segment
- `bad::Bool=true`: list of bad channels; if not false - plot bad channels using this list
- `snap::Bool=true`: snap region markers to grid at 0.0, 0.25, 0.5 and 0.75 time points

# Returns

- `seg::Union{Nothing, Tuple{Float64, Float64}}`
"""
function iview(obj::NeuroAnalyzer.NEURO; mch::Bool=true, zoom::Real=10, bad::Bool=true, snap::Bool=true)::Union{Nothing, Tuple{Float64, Float64}}

    @assert nepochs(obj) == 1 "For epoched object iview_ep() must be used."

    obj.time_pts[end] < zoom && (zoom = round(obj.time_pts[end]) / 2)

    @assert zoom > 0 "zoom must be > 0."
    @assert zoom <= signal_len(obj) / sr(obj) "zoom must be ≤ $(signal_len(obj) / sr(obj))."

    ch_order = obj.header.recording[:channel_order]
    cl = labels(obj)[ch_order]
    ch = 1:nchannels(obj)

    mono = false
    quit = true
    scale = true
    ts1 = 0
    ts2 = 0
    k = nothing

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

    function _activate(app)

        win = GtkApplicationWindow(app, "NeuroAnalyzer: iview()")
        Gtk4.default_size(win, p.attr[:size][1] + 40, p.attr[:size][2] + 40)

        can = GtkCanvas()
        can.content_width = p.attr[:size][1]
        can.content_height = p.attr[:size][2]
        g = GtkGrid()
        g.column_homogeneous = false
        g.column_spacing = 5
        g.row_spacing = 5
        g.margin_start = 5
        g.margin_end = 5
        g.margin_top = 5
        g.margin_bottom = 5

        entry_time = GtkSpinButton(obj.time_pts[1], obj.time_pts[end] - zoom, 1)
        entry_time.digits = 2
        entry_time.value = obj.time_pts[1]
        entry_time.tooltip_text = "Time position [s]"
        bt_start = GtkButton("⇤")
        bt_start.tooltip_text = "Go to the signal start"
        bt_prev = GtkButton("←")
        bt_prev.tooltip_text = "Go back by $(round(zoom)) seconds"
        bt_next = GtkButton("→")
        bt_next.tooltip_text = "Go forward by $(round(zoom)) seconds"
        bt_end = GtkButton("⇥")
        bt_end.tooltip_text = "Go to the signal end"
        entry_ts1 = GtkSpinButton(obj.time_pts[1], obj.time_pts[end], 0.5)
        bt_help = GtkButton("Help")
        bt_help.tooltip_text = "Show help"
        bt_close = GtkButton("Close")
        bt_close.tooltip_text = "Close this window"
        entry_ts1.tooltip_text = "Segment start [s]"
        entry_ts1.digits = 3
        entry_ts2 = GtkSpinButton(obj.time_pts[1], obj.time_pts[end], 0.5)
        entry_ts2.digits = 3
        entry_ts2.tooltip_text = "Segment end [s]"
        bt_ts = GtkButton("Return TS")
        bt_ts.tooltip_text = "Return selected time segment"
        bt_delete = GtkButton("Delete TS")
        bt_delete.tooltip_text = "Delete selected time segment"

        combo_ch = GtkComboBoxText()
        ch_types = uppercase.(unique(obj.header.recording[:channel_type]))
        for idx in ch_types
            length(get_channel(obj, type=lowercase(idx))) > 1 && push!(combo_ch, idx)
        end
        combo_ch.active = 0
        combo_ch.sensitive = false
        combo_ch.tooltip_text = "Channels"

        if !mch
            ch_slider = GtkScale(:v, 1:nchannels(obj))
            ch_slider.draw_value =  false
        else
            if length(ch) > 20
                ch_slider = GtkScale(:v, ch[ch_first]:(ch[end] - 19))
                ch_slider.draw_value =  false
            else
                ch_slider = GtkScale(:v, ch[1]:ch[end])
                ch_slider.draw_value =  false
                ch_slider.sensitive =  false
            end
        end
        ch_slider.tooltip_text = "Scroll channels"

        signal_slider = GtkScale(:h, obj.time_pts[1]:(obj.time_pts[end] - zoom))
        signal_slider.draw_value = false
        signal_slider.tooltip_text = "Time position"

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

        Gtk4.show(win)

        @guarded draw(can) do widget
            time1 = entry_time.value
            time2 = time1 + zoom
            time2 > obj.time_pts[end] && (time2 = obj.time_pts[end])
            ts1 = entry_ts1.value
            ts2 = entry_ts2.value
            ctx = getgc(can)
            channel_type = lowercase.(ch_types)[combo_ch.active + 1]
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
            can.content_width = p.attr[:size][1]
            can.content_height = p.attr[:size][2]
            Gtk4.default_size(win, p.attr[:size][1] + 40, p.attr[:size][2] + 40)
            img = read_from_png(io)
            set_source_surface(ctx, img, 0, 0)
            paint(ctx)
        end

        function _lmb_click(_, _, x, y)
            if x < 82
                if mch
                    ch_y = collect(50:39:793)[1:length(ch_first:ch_last)]
                    ch_idx = nothing
                    for idx in eachindex(ch_y)
                        if y > ch_y[idx] && y < ch_y[idx] + 15 && x > 0 && x < 82
                            ch_idx = idx + ch_first - 1
                        end
                    end
                    !isnothing(ch_idx) && channel_info(obj, ch=obj.header.recording[:channel_order][ch_idx])
                end
            else
                time_current = entry_time.value
                x > 1172 && (x = 1172)
                if time_current + zoom < obj.time_pts[end]
                    ts1 = time_current + round((x - 82) / (1090 / zoom), digits=3)
                else
                    ts1 = time_current + round((x - 82) / (1090 / (obj.time_pts[end] - time_current)), digits=3)
                end
                snap && (ts1 = round(ts1 * 4) / 4)
                @idle_add entry_ts1.value = round(ts1, digits=3)
            end
        end
        ggc_l = GtkGestureClick()
        ggc_l.button = 1
        push!(can, ggc_l)
        signal_connect(_lmb_click, ggc_l, "pressed")

        function _rmb_click(_, _, x, y)
            if x < 82
                if mch
                    ch_y = collect(50:39:793)[1:length(ch_first:ch_last)]
                    ch_idx = nothing
                    for idx in eachindex(ch_y)
                        if y > ch_y[idx] && y < ch_y[idx] + 15 && x > 0 && x < 82
                            ch_idx = idx + ch_first - 1
                        end
                    end
                    !isnothing(ch_idx) && (obj.header.recording[:bad_channel][obj.header.recording[:channel_order][ch_idx, 1]] = !obj.header.recording[:bad_channel][obj.header.recording[:channel_order][ch_idx, 1]])
                    draw(can)
                end
            else
                time_current = entry_time.value
                x > 1172 && (x = 1172)
                if time_current + zoom < obj.time_pts[end]
                    ts2 = time_current + ((x - 82) / (1090 / zoom))
                else
                    ts2 = time_current + ((x - 82) / (1090 / (obj.time_pts[end] - time_current)))
                end
                snap && (ts2 = round(ts2 * 4) / 4)
                @idle_add entry_ts2.value = round(ts2, digits=3)
            end
        end
        ggc_r = GtkGestureClick()
        ggc_r.button = 3
        push!(can, ggc_r)
        signal_connect(_rmb_click, ggc_r, "pressed")

        function _mwheel_scroll(_, dx, dy)
            if dy > 0.5
                if k == 0x0000ffe1 # shift
                    time_current = entry_time.value
                    if time_current < obj.time_pts[end] - zoom
                        time_current += 1
                    else
                        time_current = obj.time_pts[end] - zoom
                    end
                    @idle_add entry_time.value = time_current
                elseif k == 0x0000ffe9 # alt
                    time_current = entry_time.value
                    if time_current < obj.time_pts[end] - zoom
                        time_current += zoom
                    else
                        time_current = obj.time_pts[end] - zoom
                    end
                    @idle_add entry_time.value = time_current
                elseif k == 0x0000ffe3 # ctrl
                    if mch
                        if ch_last < length(ch)
                            ch_first += 1
                            ch_last += 1
                            @idle_add Gtk4.value(ch_slider, ch_first)
                        end
                    else
                        if ch_first < length(ch)
                            ch_first += 1
                            @idle_add Gtk4.value(ch_slider, ch_first)
                        end
                    end
                end
            elseif dy < -0.5
                if k == 0x0000ffe1 # shift
                    time_current = entry_time.value
                    if time_current >= obj.time_pts[1] + 1
                        time_current -= 1
                        @idle_add entry_time.value = time_current
                    end
                elseif k == 0x0000ffe9 # alt
                    time_current = entry_time.value
                    if time_current >= obj.time_pts[1] + zoom
                        time_current = time_current - zoom
                        @idle_add entry_time.value = time_current
                    end
                elseif k == 0x0000ffe3 # ctrl
                    if ch_first > 1
                        ch_first -= 1
                        ch_last -= 1
                        @idle_add Gtk4.value(ch_slider, ch_first)
                    end
                end
            end
        end
        ecsf = Gtk4.GtkEventControllerScroll(Gtk4.EventControllerScrollFlags_VERTICAL)
        push!(can, ecsf)
        signal_connect(_mwheel_scroll, ecsf, "scroll")

        signal_connect(signal_slider, "value-changed") do widget
            @idle_add entry_time.value = round(Gtk4.value(signal_slider))
        end

        signal_connect(ch_slider, "value-changed") do widget
            ch_first = round(Int64, Gtk4.value(ch_slider))
            length(ch) > 20 && (ch_last = ch_first + 19)
            draw(can)
        end

        signal_connect(entry_time, "value-changed") do widget
            @idle_add Gtk4.value(signal_slider, entry_time.value)
            draw(can)
        end

        signal_connect(bt_start, "clicked") do widget
            @idle_add Gtk4.value(entry_time, obj.time_pts[1])
        end

        signal_connect(bt_end, "clicked") do widget
            time_current = obj.time_pts[end] - zoom
            @idle_add Gtk4.value(entry_time, time_current)
        end

        signal_connect(entry_ts1, "value-changed") do widget
            entry_ts1.value = obj.time_pts[vsearch(entry_ts1.value, obj.time_pts)]
            ts1 = entry_ts1.value
            draw(can)
        end

        signal_connect(entry_ts2, "value-changed") do widget
            entry_ts2.value = obj.time_pts[vsearch(entry_ts2.value, obj.time_pts)]
            ts2 = entry_ts2.value
            draw(can)
        end

        signal_connect(bt_prev, "clicked") do widget
            time_current = entry_time.value
            if time_current >= obj.time_pts[1] + zoom
                time_current = time_current - zoom
                @idle_add Gtk4.value(entry_time, time_current)
            end
        end

        signal_connect(bt_next, "clicked") do widget
            time_current = entry_time.value
            if time_current < obj.time_pts[end] - zoom
                time_current += zoom
            else
                time_current = obj.time_pts[end] - zoom
            end
            @idle_add Gtk4.value(entry_time, time_current)
        end

        signal_connect(bt_delete, "clicked") do widget
            time_current = entry_time.value
            time1 = obj.time_pts[vsearch(entry_ts1.value, obj.time_pts)]
            time2 = obj.time_pts[vsearch(entry_ts2.value, obj.time_pts)]
            if time1 > time2
                warn_dialog(_nill, "Cannot delete!\nSegment start is larger than segment end.", win)
            elseif time1 == time2
                warn_dialog(_nill, "Cannot delete!\nSegment start must be different from segment end.", win)
            elseif time1 < time2
                ask_dialog("Delete segment $time1:$time2 ?", win) do ans
                    if ans
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

                        @idle_add entry_time.value = time_current
                        @idle_add entry_ts1.value = time_current
                        @idle_add entry_ts2.value = time_current

                        signal_adj = GtkAdjustment(signal_slider)
                        signal_adj.lower = obj.time_pts[1]
                        signal_adj.upper = obj.time_pts[end] - zoom
                        @idle_add Gtk4.adjustment(signal_slider, signal_adj)

                        time_adj = GtkAdjustment(entry_time)
                        time_adj.lower = obj.time_pts[1]
                        time_adj.upper = obj.time_pts[end] - zoom
                        @idle_add Gtk4.adjustment(entry_time, time_adj)

                        ts1_adj = GtkAdjustment(entry_ts1)
                        ts1_adj.lower = obj.time_pts[1]
                        ts1_adj.upper = obj.time_pts[end]
                        @idle_add Gtk4.adjustment(entry_ts1, ts1_adj)

                        ts2_adj = GtkAdjustment(entry_ts2)
                        ts2_adj.lower = obj.time_pts[1]
                        ts2_adj.upper = obj.time_pts[end]
                        @idle_add Gtk4.adjustment(entry_ts2, ts2_adj)
                    end
                end
            end
        end

        signal_connect(bt_close, "clicked") do widget
            close(win)
        end

        signal_connect(bt_ts, "clicked") do widget
            quit = false
            close(win)
        end

        signal_connect(combo_ch, "changed") do widget
            draw(can)
        end

        if mch
            help = "Keyboard shortcuts:\n\nCtrl + b\t\t\tToggle butterfly plot\nCtrl + m\t\t\tToggle mean plot\n\nCtrl + z\t\t\tScroll channels up\nCtrl + x\t\t\tScroll channels down\n\nCtrl + ,\t\t\tGo back by 1 second\nCtrl + .\t\t\tGo forward by 1 second\nAlt + ,\t\t\tGo back by $(round(zoom)) seconds\nAlt + .\t\t\tGo forward by $(round(zoom)) seconds\n\n[\t\t\t\tZoom in\n]\t\t\t\tZoom out\n\nCtrl + Enter\t\tReturn selected time segment\nCtrl + d\t\t\tDelete selected time segment\n\nCtrl + s\t\t\tToggle snapping\nAlt + s\t\t\tToggle scales\nAlt + m\t\t\tToggle monochromatic mode\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tClose\n"
        else
            help = "Keyboard shortcuts:\n\nCtrl + z\t\t\tScroll channels up\nCtrl + x\t\t\tScroll channels down\n\nCtrl + ,\t\t\tGo back by 1 second\nCtrl + .\t\t\tGo forward by 1 second\nAlt + ,\t\t\tGo back by $(round(zoom)) seconds\nAlt + .\t\t\tGo forward by $(round(zoom)) seconds\n\n[\t\t\t\tZoom in\n]\t\t\t\tZoom out\n\nCtrl + Enter\t\tReturn selected time segment\nCtrl + d\t\t\tDelete selected time segment\n\nCtrl + s\t\t\tToggle snapping\nAlt + s\t\t\tToggle scales\nAlt + m\t\t\tToggle monochromatic mode\n\nCtrl + p\t\t\tPlay current segment as audio\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tClose\n"
        end

        signal_connect(bt_help, "clicked") do widget
            info_dialog(_nill, help, win)
        end

        win_key = Gtk4.GtkEventControllerKey(win)

        signal_connect(win_key, "key-released") do widget, keyval, keycode, state
            k = nothing
        end

        signal_connect(win_key, "key-pressed") do widget, keyval, keycode, state
            k = keyval
            if keyval == UInt('[')
                if zoom > 1
                    zoom -= 1
                    bt_next.tooltip_text = "Go forward by $(round(zoom)) seconds"
                    bt_prev.tooltip_text = "Go back by $(round(zoom)) seconds"

                    signal_adj = GtkAdjustment(signal_slider)
                    signal_adj.lower = obj.time_pts[1]
                    signal_adj.upper = obj.time_pts[end] - zoom
                    @idle_add Gtk4.adjustment(signal_slider, signal_adj)

                    time_adj = GtkAdjustment(entry_time)
                    time_adj.lower = obj.time_pts[1]
                    time_adj.upper = obj.time_pts[end] - zoom
                    @idle_add Gtk4.adjustment(entry_time, time_adj)
                    draw(can)
                end
                if mch
                    help = "Keyboard shortcuts:\n\nCtrl + b\t\t\tToggle butterfly plot\nCtrl + m\t\t\tToggle mean plot\n\nCtrl + z\t\t\tScroll channels up\nCtrl + x\t\t\tScroll channels down\n\nCtrl + ,\t\t\tGo back by 1 second\nCtrl + .\t\t\tGo forward by 1 second\nAlt + ,\t\t\tGo back by $(round(zoom)) seconds\nAlt + .\t\t\tGo forward by $(round(zoom)) seconds\n\n[\t\t\t\tZoom in\n]\t\t\t\tZoom out\n\nCtrl + Enter\t\tReturn selected time segment\nCtrl + d\t\t\tDelete selected time segment\n\nCtrl + s\t\t\tToggle snapping\nAlt + s\t\t\tToggle scales\nAlt + m\t\t\tToggle monochromatic mode\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tClose\n"
                else
                    help = "Keyboard shortcuts:\n\nCtrl + z\t\t\tScroll channels up\nCtrl + x\t\t\tScroll channels down\n\nCtrl + ,\t\t\tGo back by 1 second\nCtrl + .\t\t\tGo forward by 1 second\nAlt + ,\t\t\tGo back by $(round(zoom)) seconds\nAlt + .\t\t\tGo forward by $(round(zoom)) seconds\n\n[\t\t\t\tZoom in\n]\t\t\t\tZoom out\n\nCtrl + Enter\t\tReturn selected time segment\nCtrl + d\t\t\tDelete selected time segment\n\nCtrl + s\t\t\tToggle snapping\nAlt + s\t\t\tToggle scales\nAlt + m\t\t\tToggle monochromatic mode\n\nCtrl + p\t\t\tPlay current segment as audio\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tClose\n"
                end
            elseif keyval == UInt(']')
                if zoom < 30 && zoom < obj.time_pts[end] - 1
                    zoom += 1
                    bt_next.tooltip_text = "Go forward by $(round(zoom)) seconds"
                    bt_prev.tooltip_text = "Go back by $(round(zoom)) seconds"

                    signal_adj = GtkAdjustment(signal_slider)
                    signal_adj.lower = obj.time_pts[1]
                    signal_adj.upper = obj.time_pts[end] - zoom
                    @idle_add Gtk4.adjustment(signal_slider, signal_adj)

                    time_adj = GtkAdjustment(entry_time)
                    time_adj.lower = obj.time_pts[1]
                    time_adj.upper = obj.time_pts[end] - zoom
                    @idle_add Gtk4.adjustment(entry_time, time_adj)
                    draw(can)
                else
                    zoom = obj.time_pts[end]
                    bt_next.tooltip_text = "Go forward by $(round(zoom)) seconds"
                    bt_prev.tooltip_text = "Go back by $(round(zoom)) seconds"

                    signal_adj = GtkAdjustment(signal_slider)
                    signal_adj.lower = obj.time_pts[1]
                    signal_adj.upper = obj.time_pts[end] - zoom
                    @idle_add Gtk4.adjustment(signal_slider, signal_adj)

                    time_adj = GtkAdjustment(entry_time)
                    time_adj.lower = obj.time_pts[1]
                    time_adj.upper = obj.time_pts[end] - zoom
                    @idle_add Gtk4.adjustment(entry_time, time_adj)
                    draw(can)
                end
                if mch
                    help = "Keyboard shortcuts:\n\nCtrl + b\t\t\tToggle butterfly plot\nCtrl + m\t\t\tToggle mean plot\n\nCtrl + z\t\t\tScroll channels up\nCtrl + x\t\t\tScroll channels down\n\nCtrl + ,\t\t\tGo back by 1 second\nCtrl + .\t\t\tGo forward by 1 second\nAlt + ,\t\t\tGo back by $(round(zoom)) seconds\nAlt + .\t\t\tGo forward by $(round(zoom)) seconds\n\n[\t\t\t\tZoom in\n]\t\t\t\tZoom out\n\nCtrl + Enter\t\tReturn selected time segment\nCtrl + d\t\t\tDelete selected time segment\n\nCtrl + s\t\t\tToggle snapping\nAlt + s\t\t\tToggle scales\nAlt + m\t\t\tToggle monochromatic mode\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tClose\n"
                else
                    help = "Keyboard shortcuts:\n\nCtrl + z\t\t\tScroll channels up\nCtrl + x\t\t\tScroll channels down\n\nCtrl + ,\t\t\tGo back by 1 second\nCtrl + .\t\t\tGo forward by 1 second\nAlt + ,\t\t\tGo back by $(round(zoom)) seconds\nAlt + .\t\t\tGo forward by $(round(zoom)) seconds\n\n[\t\t\t\tZoom in\n]\t\t\t\tZoom out\n\nCtrl + Enter\t\tReturn selected time segment\nCtrl + d\t\t\tDelete selected time segment\n\nCtrl + s\t\t\tToggle snapping\nAlt + s\t\t\tToggle scales\nAlt + m\t\t\tToggle monochromatic mode\n\nCtrl + p\t\t\tPlay current segment as audio\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tClose\n"
                end
            end

            # ALT
            if ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_alt == mask_alt) && keyval == UInt(','))
                time_current = entry_time.value
                if time_current >= obj.time_pts[1] + zoom
                    time_current = time_current - zoom
                    @idle_add entry_time.value = time_current
                end
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_alt == mask_alt) && keyval == UInt('.'))
                time_current = entry_time.value
                if time_current < obj.time_pts[end] - zoom
                    time_current += zoom
                    @idle_add entry_time.value = time_current
                else
                    time_current = obj.time_pts[end] - zoom
                    @idle_add entry_time.value = time_current
                end
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_alt == mask_alt) && keyval == UInt('m'))
                mono = !mono
                draw(can)
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_alt == mask_alt) && keyval == UInt('s'))
                scale = !scale
                draw(can)
            end

            # CONTROL
            if ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('q'))
                close(win)
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('h'))
                info_dialog(_nill, help, win)
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('b'))
                if mch
                    if plot_type != 1
                        ch_slider.sensitive = false
                        combo_ch.sensitive = true
                        plot_type = 1
                    else
                        ch_slider.sensitive = true
                        combo_ch.sensitive = false
                        plot_type = 0
                    end
                    draw(can)
                end
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('m'))
                if mch
                    if plot_type != 2
                        ch_slider.sensitive = false
                        combo_ch.sensitive = true
                        plot_type = 2
                    else
                        ch_slider.sensitive = true
                        combo_ch.sensitive = false
                        plot_type = 0
                    end
                    draw(can)
                end
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('p'))
                if !mch
                    _info("Playing current segment as audio")
                    time_current = entry_time.value
                    play(obj, ch=cl[ch_first], seg=(time_current, time_current+zoom), ep=1)
                end
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == 0x0000ff0d) # Enter
                quit = false
                close(win)
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('s'))
                snap = !snap
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt(','))
                time_current = entry_time.value
                if time_current >= obj.time_pts[1] + 1
                    time_current = time_current - 1
                    @idle_add entry_time.value = time_current
                end
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('.'))
                time_current = entry_time.value
                if time_current < obj.time_pts[end] - 1
                    time_current += 1
                    @idle_add entry_time.value = time_current
                else
                    time_current = obj.time_pts[end] - 1
                    @idle_add entry_time.value = time_current
                end
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('z'))
                if ch_first > 1
                    ch_first -= 1
                    ch_last -= 1
                    @idle_add Gtk4.value(ch_slider, ch_first)
                end
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('x'))
                if mch
                    if ch_last < length(ch)
                        ch_first += 1
                        ch_last += 1
                        @idle_add Gtk4.value(ch_slider, ch_first)
                    end
                else
                    if ch_first < length(ch)
                        ch_first += 1
                        @idle_add Gtk4.value(ch_slider, ch_first)
                    end
                end
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('d'))
                time_current = entry_time.value
                time1 = obj.time_pts[vsearch(entry_ts1.value, obj.time_pts)]
                time2 = obj.time_pts[vsearch(entry_ts2.value, obj.time_pts)]
                if time1 > time2
                    warn_dialog(_nill, "Cannot delete!\nSegment start is larger than segment end.", win)
                elseif time1 == time2
                    warn_dialog(_nill, "Cannot delete!\nSegment start must be different from segment end.", win)
                elseif time1 < time2
                    ask_dialog("Delete segment $time1:$time2 ?", win) do ans
                        if ans
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

                            @idle_add entry_time.value = time_current
                            @idle_add entry_ts1.value = time_current
                            @idle_add entry_ts2.value = time_current

                            signal_adj = GtkAdjustment(signal_slider)
                            signal_adj.lower = obj.time_pts[1]
                            signal_adj.upper = obj.time_pts[end] - zoom
                            @idle_add Gtk4.adjustment(signal_slider, signal_adj)

                            time_adj = GtkAdjustment(entry_time)
                            time_adj.lower = obj.time_pts[1]
                            time_adj.upper = obj.time_pts[end] - zoom
                            @idle_add Gtk4.adjustment(entry_time, time_adj)

                            ts1_adj = GtkAdjustment(entry_ts1)
                            ts1_adj.lower = obj.time_pts[1]
                            ts1_adj.upper = obj.time_pts[end]
                            @idle_add Gtk4.adjustment(entry_ts1, ts1_adj)

                            ts2_adj = GtkAdjustment(entry_ts2)
                            ts2_adj.lower = obj.time_pts[1]
                            ts2_adj.upper = obj.time_pts[end]
                            @idle_add Gtk4.adjustment(entry_ts2, ts2_adj)
                        end
                    end
                end
            end
        end
    end

    app = GtkApplication("org.neuroanalyzer.iview")
    Gtk4.signal_connect(_activate, app, :activate)
    Gtk4.GLib.stop_main_loop()
    Gtk4.run(app)

    if !quit
        time1 = obj.time_pts[vsearch(ts1, obj.epoch_time)]
        time2 = obj.time_pts[vsearch(ts2, obj.epoch_time)]
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
- `mch::Bool=true`: draw multichannel signal (up to 20 channels in one plot)
- `ep::Int64=1`: initial epoch to display
- `bad::Bool=true`: list of bad channels; if not false - plot bad channels using this list
- `snap::Bool=true`: snap region markers to grid at 0.0, 0.25, 0.5 and 0.75 time points

# Returns

- `seg::Union{Nothing, Tuple{Float64, Float64}}`
"""
function iview_ep(obj::NeuroAnalyzer.NEURO; mch::Bool=true, ep::Int64=1, bad::Bool=true, snap::Bool=true)::Union{Nothing, Tuple{Float64, Float64}}

    @assert nepochs(obj) > 1 "For continuous object iview() must be used."
    _check_epochs(obj, ep)

    ch_order = obj.header.recording[:channel_order]
    cl = labels(obj)[ch_order]
    ch = 1:nchannels(obj)

    mono = false
    quit = true
    scale = true
    k = nothing

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

    function _activate(app)

        win = GtkApplicationWindow(app, "NeuroAnalyzer: iview_ep()")
        Gtk4.default_size(win, p.attr[:size][1] + 40, p.attr[:size][2] + 40)

        can = GtkCanvas()
        can.content_height = p.attr[:size][1]
        can.content_width = p.attr[:size][2]
        g = GtkGrid()
        g.column_homogeneous = false
        g.column_spacing = 5
        g.row_spacing = 5
        g.margin_start = 5
        g.margin_end = 5
        g.margin_top = 5
        g.margin_bottom = 5

        entry_epoch = GtkSpinButton(1, nepochs(obj), 1)
        entry_epoch.digits = 0
        entry_epoch.value = ep
        entry_epoch.tooltip_text = "Epoch"
        bt_start = GtkButton("⇤")
        bt_start.tooltip_text = "Go to the first epoch"
        bt_end = GtkButton("⇥")
        bt_end.tooltip_text = "Go to the last epoch"
        bt_help = GtkButton("Help")
        bt_help.tooltip_text = "Show help"
        bt_close = GtkButton("Close")
        bt_close.tooltip_text = "Close this window"
        entry_ts1 = GtkSpinButton(obj.time_pts[1], obj.time_pts[end], 0.5)
        entry_ts1.tooltip_text = "Segment start [s]"
        entry_ts1.digits = 3
        entry_ts2 = GtkSpinButton(obj.time_pts[1], obj.time_pts[end], 0.5)
        entry_ts2.digits = 3
        entry_ts2.tooltip_text = "Segment end [s]"
        bt_ts = GtkButton("Return TS")
        bt_ts.tooltip_text = "Return selected time segment"
        bt_delete = GtkButton("Delete epoch")
        bt_delete.tooltip_text = "Delete current epoch"
        if !mch
            ch_slider = GtkScale(:v, 1:nchannels(obj))
            ch_slider.draw_value = false
        else
            if length(ch) > 20
                ch_slider = GtkScale(:v, ch[ch_first]:(ch[end] - 19))
                ch_slider.draw_value = false
            else
                ch_slider = GtkScale(:v, ch[1]:ch[end])
                ch_slider.draw_value = false
                ch_slider.sensitive = false
            end
        end

        combo_ch = GtkComboBoxText()
        ch_types = uppercase.(unique(obj.header.recording[:channel_type]))
        for idx in ch_types
            length(get_channel(obj, type=lowercase(idx))) > 1 && push!(combo_ch, idx)
        end
        combo_ch.active = 0
        combo_ch.sensitive = false
        combo_ch.tooltip_text = "Channels"

        ch_slider.tooltip_text = "Scroll channels"
        ch_slider.vexpand = true

        signal_slider = GtkScale(:h, 1:nepochs(obj))
        signal_slider.draw_value = false
        signal_slider.tooltip_text = "Current epoch"

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

        Gtk4.show(win)

        @guarded draw(can) do widget
            ep = Int64(entry_epoch.value)
            ts1 = entry_ts1.value
            ts2 = entry_ts2.value
            ctx = getgc(can)
            channel_type = lowercase.(ch_types)[combo_ch.active + 1]
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
            can.content_width = p.attr[:size][1]
            can.content_height = p.attr[:size][2]
            Gtk4.default_size(win, p.attr[:size][1] + 40, p.attr[:size][2] + 40)
            img = read_from_png(io)
            set_source_surface(ctx, img, 0, 0)
            paint(ctx)
        end

        function _lmb_click(_, _, x, y)
            if x < 82
                if mch
                    ch_y = collect(50:39:793)[1:length(ch_first:ch_last)]
                    ch_idx = nothing
                    for idx in eachindex(ch_y)
                        if y > ch_y[idx] && y < ch_y[idx] + 15 && x > 0 && x < 82
                            ch_idx = idx + ch_first - 1
                        end
                    end
                    !isnothing(ch_idx) && channel_info(obj, ch=obj.header.recording[:channel_order][ch_idx])
                end
            else
                x > 1172 && (x = 1172)
                ep = Int64(entry_epoch.value)
                ts1 = ((epoch_len(obj) / sr(obj)) * (ep - 1)) + obj.epoch_time[1] + (x - 82) / (1090 / (obj.epoch_time[end] - obj.epoch_time[1]))
                snap && (ts1 = round(ts1 * 4) / 4)
                @idle_add entry_ts1.value = round(ts1, digits=3)
            end
        end
        ggc_l = GtkGestureClick()
        ggc_l.button = 1
        push!(can, ggc_l)
        signal_connect(_lmb_click, ggc_l, "pressed")

        function _rmb_click(_, _, x, y)
            if x < 82
                if mch
                    ep = Int64(entry_epoch.value)
                    ch_y = collect(50:39:793)[1:length(ch_first:ch_last)]
                    ch_idx = nothing
                    for idx in eachindex(ch_y)
                        if y > ch_y[idx] && y < ch_y[idx] + 15 && x > 0 && x < 82
                            ch_idx = idx + ch_first - 1
                        end
                    end
                    !isnothing(ch_idx) && (obj.header.recording[:bad_channel][obj.header.recording[:channel_order][ch_idx], ep] = !obj.header.recording[:bad_channel][obj.header.recording[:channel_order][ch_idx], ep])
                    draw(can)
                end
            else
                x > 1172 && (x = 1172)
                ep = Int64(entry_epoch.value)
                ts2 = ((epoch_len(obj) / sr(obj)) * (ep - 1)) + obj.epoch_time[1] + ((x - 82) / (1090 / (obj.epoch_time[end] - obj.epoch_time[1])))
                snap && (ts2 = round(ts2 * 4) / 4)
                @idle_add entry_ts2.value = round(ts2, digits=3)
            end
        end
        ggc_r = GtkGestureClick()
        ggc_r.button = 3
        push!(can, ggc_r)
        signal_connect(_rmb_click, ggc_r, "pressed")

        function _mwheel_scroll(_, dx, dy)
            if dy > 0.5
                if k == 0x0000ffe1 # shift
                    ep = Int64(entry_epoch.value)
                    if ep < nepochs(obj)
                        ep += 1
                        @idle_add entry_epoch.value = ep
                    end
                elseif k == 0x0000ffe3 # ctrl
                    if mch
                        if ch_last < length(ch)
                            ch_first += 1
                            ch_last += 1
                            @idle_add Gtk4.value(ch_slider, ch_first)
                        end
                    else
                        if ch_first < length(ch)
                            ch_first += 1
                            @idle_add Gtk4.value(ch_slider, ch_first)
                        end
                    end
                end
            elseif dy < 0.5
                if k == 0x0000ffe1 # shift
                    ep = Int64(entry_epoch.value)
                    if ep > 1
                        ep -= 1
                        @idle_add entry_epoch.value = ep
                    end
                elseif k == 0x0000ffe3 # ctrl
                    if ch_first > 1
                        ch_first -= 1
                        ch_last -= 1
                        @idle_add Gtk4.value(ch_slider, ch_first)
                    end
                end
            end
        end
        ecsf = Gtk4.GtkEventControllerScroll(Gtk4.EventControllerScrollFlags_VERTICAL)
        push!(can, ecsf)
        signal_connect(_mwheel_scroll, ecsf, "scroll")

        signal_connect(signal_slider, "value-changed") do widget
            @idle_add entry_epoch.value = round(Int64, Gtk4.value(signal_slider))
            draw(can)
        end

        signal_connect(ch_slider, "value-changed") do widget
            ch_first = round(Int64, Gtk4.value(ch_slider))
            length(ch) > 20 && (ch_last = ch_first + 19)
            draw(can)
        end

        signal_connect(entry_epoch, "value-changed") do widget
            @idle_add Gtk4.value(signal_slider, Int64(entry_epoch.value))
            draw(can)
        end

        signal_connect(bt_start, "clicked") do widget
            @idle_add Gtk4.value(entry_epoch, 1)
        end

        signal_connect(bt_end, "clicked") do widget
            @idle_add Gtk4.value(entry_epoch, nepochs(obj))
        end

        signal_connect(entry_ts1, "value-changed") do widget
            entry_ts1.value = obj.time_pts[vsearch(entry_ts1.value, obj.time_pts)]
            ts1 = entry_ts1.value
            draw(can)
        end

        signal_connect(entry_ts2, "value-changed") do widget
            entry_ts2.value = obj.time_pts[vsearch(entry_ts2.value, obj.time_pts)]
            ts2 = entry_ts2.value
            draw(can)
        end

        signal_connect(bt_delete, "clicked") do widget
            if nepochs(obj) > 1
                ep = Int64(entry_epoch.value)
                ask_dialog("Delete epoch: $ep", win) do ans
                    if ans
                        delete_epoch!(obj, ep=ep)
                        _info("Deleted epoch: $ep")
                        @idle_add entry_epoch.value = ep

                        epoch_adj = GtkAdjustment(entry_epoch)
                        epoch_adj.lower = 1
                        epoch_adj.upper = nepochs(obj)
                        @idle_add Gtk4.adjustment(entry_epoch, epoch_adj)

                        signal_adj = GtkAdjustment(signal_slider)
                        signal_adj.lower = 1
                        signal_adj.upper = nepochs(obj)
                        @idle_add Gtk4.adjustment(signal_slider, signal_adj)

                        draw(can)
                    end
                end
            else
                warn_dialog(_nill, "You cannot delete the last epoch.", win)
            end
        end

        signal_connect(bt_close, "clicked") do widget
            close(win)
        end

        signal_connect(bt_ts, "clicked") do widget
            quit = false
            close(win)
        end

        signal_connect(combo_ch, "changed") do widget
            draw(can)
        end

        if mch
            help = "Keyboard shortcuts:\n\nCtrl + b\t\t\tToggle butterfly plot\nCtrl + m\t\t\tToggle mean plot\n\nCtrl + z\t\t\tScroll channels up\nCtrl + x\t\t\tScroll channels down\n\nCtrl + ,\t\t\tPrevious epoch\nCtrl + .\t\t\tNext epoch\n\nCtrl + Enter\t\tReturn selected time segment\nCtrl + d\t\t\tDelete current epoch\n\nCtrl + s\t\t\tToggle snapping\nAlt + s\t\t\tToggle scales\nAlt + m\t\t\tToggle monochromatic mode\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tClose\n"
        else
            help = "Keyboard shortcuts:\n\nCtrl + b\t\t\tToggle butterfly plot\nCtrl + m\t\t\tToggle mean plot\n\nCtrl + z\t\t\tScroll channels up\nCtrl + x\t\t\tScroll channels down\n\nCtrl + ,\t\t\tPrevious epoch\nCtrl + .\t\t\tNext epoch\n\nCtrl + Enter\t\tReturn selected time segment\nCtrl + d\t\t\tDelete current epoch\n\nCtrl + s\t\t\tToggle snapping\nAlt + s\t\t\tToggle scales\nAlt + m\t\t\tToggle monochromatic mode\n\nplay current epoch as audio\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tClose\n"
        end

        signal_connect(bt_help, "clicked") do widget
            info_dialog(_nill, help, win)
        end

        win_key = Gtk4.GtkEventControllerKey(win)

        signal_connect(win_key, "key-released") do widget, keyval, keycode, state
            k = nothing
        end

        signal_connect(win_key, "key-pressed") do widget, keyval, keycode, state
            k = keyval

            # ALT
            if ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_alt == mask_alt) && keyval == UInt('m'))
                mono = !mono
                draw(can)
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_alt == mask_alt) && keyval == UInt('s'))
                scale = !scale
                draw(can)
            end

            # CONTROL
            if ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('q'))
                close(win)
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('h'))
                info_dialog(_nill, help, win)
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('b'))
                if mch
                    if plot_type != 1
                        ch_slider.sensitive = false
                        combo_ch.sensitive = true
                        plot_type = 1
                    else
                        ch_slider.sensitive = true
                        combo_ch.sensitive = false
                        plot_type = 0
                    end
                    draw(can)
                end
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('m'))
                if mch
                    if plot_type != 2
                        ch_slider.sensitive = false
                        combo_ch.sensitive = true
                        plot_type = 2
                    else
                        ch_slider.sensitive = true
                        combo_ch.sensitive = false
                        plot_type = 0
                    end
                    draw(can)
                end
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('p'))
                if !mch
                    _info("Playing current segment as audio")
                    time_current = entry_time.value
                    play(obj, ch=cl[ch_first], seg=(time_current, time_current + ep_len(obj)), ep=1)
                end
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == 0x0000ff0d) # Enter
                quit = false
                close(win)
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('s'))
                snap = !snap
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt(','))
                ep = Int64(entry_epoch.value)
                if ep > 1
                    ep -= 1
                    @idle_add entry_epoch.value = ep
                end
                draw(can)
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('.'))
                ep = Int64(entry_epoch.value)
                if ep > 1
                    ep -= 1
                    @idle_add entry_epoch.value = ep
                end
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('z'))
                if ch_first > 1
                    ch_first -= 1
                    ch_last -= 1
                    @idle_add Gtk4.value(ch_slider, ch_first)
                end
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('x'))
                if mch
                    if ch_last < length(ch)
                        ch_first += 1
                        ch_last += 1
                        @idle_add Gtk4.value(ch_slider, ch_first)
                    end
                else
                    if ch_first < length(ch)
                        ch_first += 1
                        @idle_add Gtk4.value(ch_slider, ch_first)
                    end
                end
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('d'))
                ep = Int64(entry_epoch.value)
                ask_dialog("Delete epoch $ep ?", win) do ans
                    if ans
                        delete_epoch!(obj, ep=ep)
                        _info("Deleted epoch: $ep")
                        ep = ep > 1 ? ep -= 1 : ep = 1
                        @idle_add entry_epoch.value = ep

                        epoch_adj = GtkAdjustment(entry_epoch)
                        epoch_adj.lower = 1
                        epoch_adj.upper = nepochs(obj)
                        @idle_add Gtk4.adjustment(entry_epoch, epoch_adj)

                        signal_adj = GtkAdjustment(signal_slider)
                        signal_adj.lower = 1
                        signal_adj.upper = nepochs(obj)
                        @idle_add Gtk4.adjustment(signal_slider, signal_adj)

                        draw(can)
                    end
                end
            end
        end
    end

    app = GtkApplication("org.neuroanalyzer.iview_ep")
    Gtk4.signal_connect(_activate, app, :activate)
    Gtk4.GLib.stop_main_loop()
    Gtk4.run(app)

    if !quit
        time1 = obj.time_pts[vsearch(ts1, obj.epoch_time)]
        time2 = obj.time_pts[vsearch(ts2, obj.epoch_time)]
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

    @assert nepochs(obj) == 1 "For epoched object iview_ep() must be used."

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
    k = nothing

    if length(ch) > 20
        ch_first = 1
        ch_last = ch_first + 19
    else
        ch_first = 1
        ch_last = length(ch)
    end

    p = NeuroAnalyzer.plot(obj1, obj2, ch=cl[ch_first:ch_last], title="")

    function _activate(app)

        win = GtkApplicationWindow(app, "NeuroAnalyzer: iview()")
        Gtk4.default_size(win, p.attr[:size][1] + 40, p.attr[:size][2] + 40)

        can = GtkCanvas()
        can.content_width = p.attr[:size][1]
        can.content_height = p.attr[:size][2]
        g = GtkGrid()
        g.column_homogeneous = false
        g.column_spacing = 5
        g.row_spacing = 5
        g.margin_start = 5
        g.margin_end = 5
        g.margin_top = 5
        g.margin_bottom = 5

        entry_time = GtkSpinButton(obj1.time_pts[1], obj1.time_pts[end] - zoom, 1)
        entry_time.digits = 2
        entry_time.value = obj1.time_pts[1]
        entry_time.tooltip_text = "Time position [s]"
        bt_start = GtkButton("⇤")
        bt_start.tooltip_text = "Go to the signal start"
        bt_prev = GtkButton("←")
        bt_prev.tooltip_text = "Go back by $(round(zoom)) seconds"
        bt_next = GtkButton("→")
        bt_next.tooltip_text = "Go forward by $(round(zoom)) seconds"
        bt_end = GtkButton("⇥")
        bt_end.tooltip_text = "Go to the signal end"
        bt_help = GtkButton("Help")
        bt_help.tooltip_text = "Show help"
        bt_close = GtkButton("Close")
        bt_close.tooltip_text = "Close this window"
        if length(ch) > 20
            ch_slider = GtkScale(:v, ch[ch_first]:(ch[end] - 19))
            ch_slider.draw_value = false
        else
            ch_slider = GtkScale(:v, ch[ch_first]:ch[end])
            ch_slider.draw_value = false
            ch_slider.sensitive = false
        end
        ch_slider.tooltip_text = "Scroll channels"

        signal_slider = GtkScale(:h, 1:obj1.time_pts[end] - zoom)
        signal_slider.draw_value = false
        signal_slider.tooltip_text = "Time position"

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

        Gtk4.show(win)

        @guarded draw(can) do widget
            time1 = entry_time.value
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

        function _mwheel_scroll(_, dx, dy)
            if dy > 0.5
                if k == 0x0000ffe1 # shift
                    time_current = entry_time.value
                    if time_current < obj1.time_pts[end] - zoom
                        time_current += 1
                    else
                        time_current = obj1.time_pts[end] - zoom
                    end
                    @idle_add entry_time.value = time_current
                elseif k == 0x0000ffe9 # alt
                    time_current = entry_time.value
                    if time_current < obj1.time_pts[end] - zoom
                        time_current += zoom
                    else
                        time_current = obj1.time_pts[end] - zoom
                    end
                    @idle_add entry_time.value = time_current
                elseif k == 0x0000ffe3 # ctrl
                    if ch_last < length(ch)
                        ch_first += 1
                        ch_last += 1
                        @idle_add Gtk4.value(ch_slider, ch_first)
                    end
                end
            elseif dy < -0.5
                if k == 0x0000ffe1 # shift
                    time_current = entry_time.value
                    if time_current >= obj1.time_pts[1] + 1
                        time_current -= 1
                        @idle_add entry_time.value = time_current
                    end
                elseif k == 0x0000ffe9 # alt
                    time_current = entry_time.value
                    if time_current >= obj1.time_pts[1] + zoom
                        time_current = time_current - zoom
                        @idle_add entry_time.value = time_current
                    end
                elseif k == 0x0000ffe3 # ctrl
                    if ch_first > 1
                        ch_first -= 1
                        ch_last -= 1
                        @idle_add Gtk4.value(ch_slider, ch_first)
                    end
                end
            end
        end
        ecsf = Gtk4.GtkEventControllerScroll(Gtk4.EventControllerScrollFlags_VERTICAL)
        push!(can, ecsf)
        signal_connect(_mwheel_scroll, ecsf, "scroll")

        signal_connect(signal_slider, "value-changed") do widget
            @idle_add entry_time.value = round(Gtk4.value(signal_slider))
        end

        signal_connect(ch_slider, "value-changed") do widget
            ch_first = round(Int64, Gtk4.value(ch_slider))
            length(ch) > 20 && (ch_last = ch_first + 19)
            draw(can)
        end

        signal_connect(entry_time, "value-changed") do widget
            @idle_add Gtk4.value(signal_slider, entry_time.value)
            draw(can)
        end

        signal_connect(bt_prev, "clicked") do widget
            time_current = entry_time.value
            if time_current >= obj1.time_pts[1] + zoom
                time_current = time_current - zoom
                @idle_add Gtk4.value(entry_time, time_current)
            end
        end

        signal_connect(bt_next, "clicked") do widget
            time_current = entry_time.value
            if time_current < obj1.time_pts[end] - zoom
                time_current += zoom
            else
                time_current = obj1.time_pts[end] - zoom
            end
            @idle_add Gtk4.value(entry_time, time_current)
        end

        signal_connect(bt_start, "clicked") do widget
            @idle_add Gtk4.value(entry_time, obj1.time_pts[1])
        end

        signal_connect(bt_end, "clicked") do widget
            time_current = obj1.time_pts[end] - zoom
            @idle_add Gtk4.value(entry_time, time_current)
        end

        help = "Keyboard shortcuts:\n\nCtrl + z\t\t\tScroll channels up\nCtrl + x\t\t\tScroll channels down\n\nCtrl + ,\t\t\tGo back by 1 second\nCtrl + .\t\t\tGo forward by 1 second\nAlt + ,\t\t\tGo back by $(round(zoom)) seconds\nAlt + .\t\t\tGo forward by $(round(zoom)) seconds\n\n[\t\t\t\tZoom in\n]\t\t\t\tZoom out\n\nAlt + s\t\t\tToggle scales\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tClose\n"

        signal_connect(bt_help, "clicked") do widget
            info_dialog(_nill, help, win)
        end

        signal_connect(bt_close, "clicked") do widget
            close(win)
        end

        win_key = Gtk4.GtkEventControllerKey(win)

        signal_connect(win_key, "key-released") do widget, keyval, keycode, state
            k = nothing
        end

        signal_connect(win_key, "key-pressed") do widget, keyval, keycode, state
            k = keyval
            if keyval == UInt('[')
                if zoom > 1
                    zoom -= 1
                    bt_next.tooltip_text = "Go forward by $(round(zoom)) seconds"
                    bt_prev.tooltip_text = "Go back by $(round(zoom)) seconds"

                    signal_adj = GtkAdjustment(signal_slider)
                    signal_adj.lower = obj1.time_pts[1]
                    signal_adj.upper = obj1.time_pts[end] - zoom
                    @idle_add Gtk4.adjustment(signal_slider, signal_adj)

                    time_adj = GtkAdjustment(entry_time)
                    time_adj.lower = obj1.time_pts[1]
                    time_adj.upper = obj1.time_pts[end] - zoom
                    @idle_add Gtk4.adjustment(entry_time, time_adj)
                    draw(can)
                end
                help = "Keyboard shortcuts:\n\nCtrl + z\t\t\tScroll channels up\nCtrl + x\t\t\tScroll channels down\n\nCtrl + ,\t\t\tGo back by 1 second\nCtrl + .\t\t\tGo forward by 1 second\nAlt + ,\t\t\tGo back by $(round(zoom)) seconds\nAlt + .\t\t\tGo forward by $(round(zoom)) seconds\n\n[\t\t\t\tZoom in\n]\t\t\t\tZoom out\n\nAlt + s\t\t\tToggle scales\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tClose\n"
            elseif keyval == UInt(']')
                if zoom < 30 && zoom < obj1.time_pts[end] - 1
                    zoom += 1
                    bt_next.tooltip_text = "Go forward by $(round(zoom)) seconds"
                    bt_prev.tooltip_text = "Go back by $(round(zoom)) seconds"

                    signal_adj = GtkAdjustment(signal_slider)
                    signal_adj.lower = obj1.time_pts[1]
                    signal_adj.upper = obj1.time_pts[end] - zoom
                    @idle_add Gtk4.adjustment(signal_slider, signal_adj)

                    time_adj = GtkAdjustment(entry_time)
                    time_adj.lower = obj1.time_pts[1]
                    time_adj.upper = obj1.time_pts[end] - zoom
                    @idle_add Gtk4.adjustment(entry_time, time_adj)
                    draw(can)
                else
                    zoom = obj1.time_pts[end]
                    bt_next.tooltip_text = "Go forward by $(round(zoom)) seconds"
                    bt_prev.tooltip_text = "Go back by $(round(zoom)) seconds"

                    signal_adj = GtkAdjustment(signal_slider)
                    signal_adj.lower = obj1.time_pts[1]
                    signal_adj.upper = obj1.time_pts[end] - zoom
                    @idle_add Gtk4.adjustment(signal_slider, signal_adj)

                    time_adj = GtkAdjustment(entry_time)
                    time_adj.lower = obj1.time_pts[1]
                    time_adj.upper = obj1.time_pts[end] - zoom
                    @idle_add Gtk4.adjustment(entry_time, time_adj)
                    draw(can)
                end
                help = "Keyboard shortcuts:\n\nCtrl + z\t\t\tScroll channels up\nCtrl + x\t\t\tScroll channels down\n\nCtrl + ,\t\t\tGo back by 1 second\nCtrl + .\t\t\tGo forward by 1 second\nAlt + ,\t\t\tGo back by $(round(zoom)) seconds\nAlt + .\t\t\tGo forward by $(round(zoom)) seconds\n\n[\t\t\t\tZoom in\n]\t\t\t\tZoom out\n\nAlt + s\t\t\tToggle scales\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tClose\n"
            end

            # ALT
            if ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_alt == mask_alt) && keyval == UInt(','))
                time_current = entry_time.value
                if time_current >= obj1.time_pts[1] + zoom
                    time_current = time_current - zoom
                    @idle_add entry_time.value = time_current
                end
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_alt == mask_alt) && keyval == UInt('.'))
                time_current = entry_time.value
                if time_current < obj1.time_pts[end] - zoom
                    time_current += zoom
                    @idle_add entry_time.value = time_current
                else
                    time_current = obj1.time_pts[end] - zoom
                    @idle_add entry_time.value = time_current
                end
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_alt == mask_alt) && keyval == UInt('s'))
                scale = !scale
                draw(can)
            end

            # CONTROL
            if ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('q'))
                close(win)
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('h'))
                info_dialog(_nill, help, win)
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt(','))
                time_current = entry_time.value
                if time_current >= obj1.time_pts[1] + 1
                    time_current = time_current - 1
                    @idle_add entry_time.value = time_current
                end
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('.'))
                time_current = entry_time.value
                if time_current < obj1.time_pts[end] - 1
                    time_current += 1
                    @idle_add entry_time.value = time_current
                else
                    time_current = obj1.time_pts[end] - 1
                    @idle_add entry_time.value = time_current
                end
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('z'))
                if ch_first > 1
                    ch_first -= 1
                    ch_last -= 1
                    @idle_add Gtk4.value(ch_slider, ch_first)
                end
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('x'))
                if ch_last < length(ch)
                    ch_first += 1
                    ch_last += 1
                    @idle_add Gtk4.value(ch_slider, ch_first)
                end
            end
        end
    end

    app = GtkApplication("org.neuroanalyzer.iview")
    Gtk4.signal_connect(_activate, app, :activate)
    Gtk4.GLib.stop_main_loop()
    Gtk4.run(app)

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
    @assert nepochs(obj1) == nepochs(obj2) "Both signals must have the same number of epochs."
    @assert nepochs(obj1) > 1 "For continuous object iview() must be used."

    ch_order = obj1.header.recording[:channel_order]
    cl = labels(obj1)[ch_order]
    ch = 1:nchannels(obj1)
    _check_epochs(obj1, ep)

    scale = true
    k = nothing

    if length(ch) > 20
        ch_first = 1
        ch_last = ch_first + 19
    else
        ch_first = 1
        ch_last = length(ch)
    end

    p = NeuroAnalyzer.plot(obj1, obj2, ch=cl[ch_first:ch_last], ep=ep, title="")

    function _activate(app)

        win = GtkApplicationWindow(app, "NeuroAnalyzer: iview_ep()")
        Gtk4.default_size(win, p.attr[:size][1] + 40, p.attr[:size][2] + 40)

        can = GtkCanvas()
        can.content_width = p.attr[:size][1]
        can.content_height = p.attr[:size][2]
        g = GtkGrid()
        g.column_homogeneous = false
        g.column_spacing = 5
        g.row_spacing = 5
        g.margin_start = 5
        g.margin_end = 5
        g.margin_top = 5
        g.margin_bottom = 5

        entry_epoch = GtkSpinButton(1, nepochs(obj1), 1)
        entry_epoch.value = ep
        entry_epoch.tooltip_text = "Epoch"
        bt_start = GtkButton("⇤")
        bt_start.tooltip_text = "Go to the first epoch"
        bt_end = GtkButton("⇥")
        bt_end.tooltip_text = "Go to the last epoch"
        bt_help = GtkButton("Help")
        bt_help.tooltip_text = "Show help"
        bt_close = GtkButton("Close")
        bt_close.tooltip_text = "Close this window"

        if length(ch) > 20
            ch_slider = GtkScale(:v, ch[ch_first]:(ch[end] - 19))
            ch_slider.draw_value = false
        else
            ch_slider = GtkScale(:v, ch[ch_first]:ch[end])
            ch_slider.draw_value = false
            ch_slider.sensitive = false
        end
        ch_slider.tooltip_text = "Scroll channels"

        signal_slider = GtkScale(:h, 1:nepochs(obj1))
        signal_slider.draw_value = false
        signal_slider.tooltip_text = "Current epoch"

        g[1:5, 1] = can
        g[6, 1] = ch_slider
        g[1:5, 2] = signal_slider
        g[1, 3] = bt_start
        g[2, 3] = entry_epoch
        g[3, 3] = bt_end
        g[4, 3] = bt_help
        g[5, 3] = bt_close
        push!(win, g)

        Gtk4.show(win)

        @guarded draw(can) do widget
            ep = Int64(entry_epoch.value)
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

        function _mwheel_scroll(_, dx, dy)
            if dy > 0.5
                if k == 0x0000ffe1 # shift
                    ep = Int64(entry_epoch.value)
                    if ep < nepochs(obj1)
                        ep += 1
                        @idle_add entry_epoch.value = ep
                    end
                elseif k == 0x0000ffe3 # ctrl
                    if mch
                        if ch_last < length(ch)
                            ch_first += 1
                            ch_last += 1
                            @idle_add Gtk4.value(ch_slider, ch_first)
                        end
                    else
                        if ch_first < length(ch)
                            ch_first += 1
                            @idle_add Gtk4.value(ch_slider, ch_first)
                        end
                    end
                end
            elseif dy < 0.5
                if k == 0x0000ffe1 # shift
                    ep = Int64(entry_epoch.value)
                    if ep > 1
                        ep -= 1
                        @idle_add entry_epoch.value = ep
                    end
                elseif k == 0x0000ffe3 # ctrl
                    if ch_first > 1
                        ch_first -= 1
                        ch_last -= 1
                        @idle_add Gtk4.value(ch_slider, ch_first)
                    end
                end
            end
        end
        ecsf = Gtk4.GtkEventControllerScroll(Gtk4.EventControllerScrollFlags_VERTICAL)
        push!(can, ecsf)
        signal_connect(_mwheel_scroll, ecsf, "scroll")

        signal_connect(signal_slider, "value-changed") do widget
            @idle_add entry_epoch.value = round(Int64, Gtk4.value(signal_slider))
            draw(can)
        end

        signal_connect(ch_slider, "value-changed") do widget
            ch_first = round(Int64, Gtk4.value(ch_slider))
            length(ch) > 20 && (ch_last = ch_first + 19)
            draw(can)
        end

        signal_connect(entry_epoch, "value-changed") do widget
            @idle_add Gtk4.value(signal_slider, Int64(entry_epoch.value))
            draw(can)
        end

        signal_connect(bt_start, "clicked") do widget
            @idle_add Gtk4.value(entry_epoch, 1)
        end

        signal_connect(bt_end, "clicked") do widget
            @idle_add Gtk4.value(entry_epoch, nepochs(obj1))
        end

        signal_connect(bt_close, "clicked") do widget
            close(win)
        end

        help = "Keyboard shortcuts:\n\nCtrl + z\t\t\tScroll channels up\nCtrl + x\t\t\tScroll channels down\n\nCtrl + ,\t\t\tPrevious epoch\nCtrl + .\t\t\tNext epoch\n\nAlt + s\t\t\tToggle scales\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tClose\n"

        signal_connect(bt_help, "clicked") do widget
            info_dialog(_nill, help, win)
        end

        win_key = Gtk4.GtkEventControllerKey(win)

        signal_connect(win_key, "key-released") do widget, keyval, keycode, state
            k = nothing
        end

        signal_connect(win_key, "key-pressed") do widget, keyval, keycode, state
            k = keyval

            # ALT
            if ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_alt == mask_alt) && keyval == UInt('m'))
                mono = !mono
                draw(can)
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_alt == mask_alt) && keyval == UInt('s'))
                scale = !scale
                draw(can)
            end

            # CONTROL
            if ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('q'))
                close(win)
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('h'))
                info_dialog(_nill, help, win)
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt(','))
                ep = Int64(entry_epoch.value)
                if ep > 1
                    ep -= 1
                    @idle_add entry_epoch.value = ep
                end
                draw(can)
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('.'))
                ep = Int64(entry_epoch.value)
                if ep > 1
                    ep -= 1
                    @idle_add entry_epoch.value = ep
                end
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('z'))
                if ch_first > 1
                    ch_first -= 1
                    ch_last -= 1
                    @idle_add Gtk4.value(ch_slider, ch_first)
                end
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('x'))
                if mch
                    if ch_last < length(ch)
                        ch_first += 1
                        ch_last += 1
                        @idle_add Gtk4.value(ch_slider, ch_first)
                    end
                else
                    if ch_first < length(ch)
                        ch_first += 1
                        @idle_add Gtk4.value(ch_slider, ch_first)
                    end
                end
            end
        end
    end

    app = GtkApplication("org.neuroanalyzer.iview_ep")
    Gtk4.signal_connect(_activate, app, :activate)
    Gtk4.GLib.stop_main_loop()
    Gtk4.run(app)

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

    @assert nepochs(obj) == 1 "For epoched object iview_ep() must be used."

    @assert zoom > 0 "zoom must be > 0."
    @assert zoom <= signal_len(obj) / sr(obj) "zoom must be ≤ $(signal_len(obj) / sr(obj))."

    c isa Symbol && (c = _get_component(obj, c))
    ch = axes(c, 1)

    mono = false
    quit = false
    scale = true
    k = nothing

    if length(ch) > 20
        ch_first = 1
        ch_last = 20
    else
        ch_first = 1
        ch_last = length(ch)
    end

    length(ch) == 1 && (ch = [ch])

    p = NeuroAnalyzer.plot(obj, c, c_idx=ch[ch_first:ch_last], mono=mono, scale=scale, title="")

    function _activate(app)

        win = GtkApplicationWindow(app, "NeuroAnalyzer: iview()")
        Gtk4.default_size(win, p.attr[:size][1] + 40, p.attr[:size][2] + 40)

        can = GtkCanvas()
        can.content_width = p.attr[:size][1]
        can.content_height = p.attr[:size][2]
        g = GtkGrid()
        g.column_homogeneous = false
        g.column_spacing = 5
        g.row_spacing = 5

        entry_time = GtkSpinButton(obj.time_pts[1], obj.time_pts[end] - zoom, 1)
        entry_time.digits = 2
        entry_time.value = obj.time_pts[1]
        entry_time.tooltip_text = "Time position [s]"
        bt_start = GtkButton("⇤")
        bt_start.tooltip_text = "Go to the signal start"
        bt_prev = GtkButton("←")
        bt_prev.tooltip_text = "Go back by $(round(zoom)) seconds"
        bt_next = GtkButton("→")
        bt_next.tooltip_text = "Go forward by $(round(zoom)) seconds"
        bt_end = GtkButton("⇥")
        bt_end.tooltip_text = "Go to the signal end"
        bt_help = GtkButton("Help")
        bt_help.tooltip_text = "Show help"
        bt_close = GtkButton("Close")
        bt_close.tooltip_text = "Close this window"
        if length(ch) > 20
            ch_slider = GtkScale(:v, ch[ch_first]:(ch[end] - 19))
            ch_slider.draw_value = false
        else
            ch_slider = GtkScale(:v, ch[ch_first]:ch[end])
            ch_slider.draw_value = false
            ch_slider.sensitive = false
        end
        ch_slider.tooltip_text = "Scroll components"

        signal_slider = GtkScale(:h, obj.time_pts[1]:obj.time_pts[end] - zoom)
        signal_slider.draw_value = false
        signal_slider.tooltip_text = "Time position"

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

        Gtk4.show(win)

        @guarded draw(can) do widget
            time1 = entry_time.value
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

        function _mwheel_scroll(_, dx, dy)
            if dy > 0.5
                if k == 0x0000ffe1 # shift
                    time_current = entry_time.value
                    if time_current < obj.time_pts[end] - zoom
                        time_current += 1
                    else
                        time_current = obj.time_pts[end] - zoom
                    end
                    @idle_add entry_time.value = time_current
                elseif k == 0x0000ffe9 # alt
                    time_current = entry_time.value
                    if time_current < obj.time_pts[end] - zoom
                        time_current += zoom
                    else
                        time_current = obj.time_pts[end] - zoom
                    end
                    @idle_add entry_time.value = time_current
                elseif k == 0x0000ffe3 # ctrl
                    if mch
                        if ch_last < length(ch)
                            ch_first += 1
                            ch_last += 1
                            @idle_add Gtk4.value(ch_slider, ch_first)
                        end
                    else
                        if ch_first < length(ch)
                            ch_first += 1
                            @idle_add Gtk4.value(ch_slider, ch_first)
                        end
                    end
                end
            elseif dy < -0.5
                if k == 0x0000ffe1 # shift
                    time_current = entry_time.value
                    if time_current >= obj.time_pts[1] + 1
                        time_current -= 1
                        @idle_add entry_time.value = time_current
                    end
                elseif k == 0x0000ffe9 # alt
                    time_current = entry_time.value
                    if time_current >= obj.time_pts[1] + zoom
                        time_current = time_current - zoom
                        @idle_add entry_time.value = time_current
                    end
                elseif k == 0x0000ffe3 # ctrl
                    if ch_first > 1
                        ch_first -= 1
                        ch_last -= 1
                        @idle_add Gtk4.value(ch_slider, ch_first)
                    end
                end
            end
        end
        ecsf = Gtk4.GtkEventControllerScroll(Gtk4.EventControllerScrollFlags_VERTICAL)
        push!(can, ecsf)
        signal_connect(_mwheel_scroll, ecsf, "scroll")

        signal_connect(signal_slider, "value-changed") do widget
            @idle_add entry_time.value = round(Gtk4.value(signal_slider))
        end

        signal_connect(ch_slider, "value-changed") do widget
            ch_first = round(Int64, Gtk4.value(ch_slider))
            length(ch) > 20 && (ch_last = ch_first + 19)
            draw(can)
        end

        signal_connect(entry_time, "value-changed") do widget
            @idle_add Gtk4.value(signal_slider, entry_time.value)
            draw(can)
        end

        signal_connect(bt_start, "clicked") do widget
            @idle_add Gtk4.value(entry_time, obj.time_pts[1])
        end

        signal_connect(bt_end, "clicked") do widget
            time_current = obj.time_pts[end] - zoom
            @idle_add Gtk4.value(entry_time, time_current)
        end

        signal_connect(bt_prev, "clicked") do widget
            time_current = entry_time.value
            if time_current >= obj.time_pts[1] + zoom
                time_current = time_current - zoom
                @idle_add Gtk4.value(entry_time, time_current)
            end
        end

        signal_connect(bt_next, "clicked") do widget
            time_current = entry_time.value
            if time_current < obj.time_pts[end] - zoom
                time_current += zoom
            else
                time_current = obj.time_pts[end] - zoom
            end
            @idle_add Gtk4.value(entry_time, time_current)
        end

        signal_connect(bt_close, "clicked") do widget
            close(win)
        end

        help = "Keyboard shortcuts:\n\nCtrl + z\t\t\tScroll channels up\nCtrl + x\t\t\tScroll channels down\n\nCtrl + ,\t\t\tGo back by 1 second\nCtrl + .\t\t\tGo forward by 1 second\nAlt + ,\t\t\tGo back by $(round(zoom)) seconds\nAlt + .\t\t\tGo forward by $(round(zoom)) seconds\n\n[\t\t\t\tZoom in\n]\t\t\t\tZoom out\n\nToggle scales\nAlt + m\t\t\tToggle monochromatic mode\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tClose\n"

        signal_connect(bt_help, "clicked") do widget
            info_dialog(_nill, help, win)
        end

        win_key = Gtk4.GtkEventControllerKey(win)

        signal_connect(win_key, "key-released") do widget, keyval, keycode, state
            k = nothing
        end

        signal_connect(win_key, "key-pressed") do widget, keyval, keycode, state
            k = keyval
            if keyval == UInt('[')
                if zoom > 1
                    zoom -= 1
                    bt_next.tooltip_text = "Go forward by $(round(zoom)) seconds"
                    bt_prev.tooltip_text = "Go back by $(round(zoom)) seconds"

                    signal_adj = GtkAdjustment(signal_slider)
                    signal_adj.lower = obj.time_pts[1]
                    signal_adj.upper = obj.time_pts[end] - zoom
                    @idle_add Gtk4.adjustment(signal_slider, signal_adj)

                    time_adj = GtkAdjustment(entry_time)
                    time_adj.lower = obj.time_pts[1]
                    time_adj.upper = obj.time_pts[end] - zoom
                    @idle_add Gtk4.adjustment(entry_time, time_adj)
                    draw(can)
                end
                if mch
                    help = "Keyboard shortcuts:\n\nCtrl + b\t\t\tToggle butterfly plot\nCtrl + m\t\t\tToggle mean plot\n\nCtrl + z\t\t\tScroll channels up\nCtrl + x\t\t\tScroll channels down\n\nCtrl + ,\t\t\tGo back by 1 second\nCtrl + .\t\t\tGo forward by 1 second\nAlt + ,\t\t\tGo back by $(round(zoom)) seconds\nAlt + .\t\t\tGo forward by $(round(zoom)) seconds\n\n[\t\t\t\tZoom in\n]\t\t\t\tZoom out\n\nCtrl + Enter\t\tReturn selected time segment\nCtrl + d\t\t\tDelete selected time segment\n\nCtrl + s\t\t\tToggle snapping\nAlt + s\t\t\tToggle scales\nAlt + m\t\t\tToggle monochromatic mode\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tClose\n"
                else
                    help = "Keyboard shortcuts:\n\nCtrl + z\t\t\tScroll channels up\nCtrl + x\t\t\tScroll channels down\n\nCtrl + ,\t\t\tGo back by 1 second\nCtrl + .\t\t\tGo forward by 1 second\nAlt + ,\t\t\tGo back by $(round(zoom)) seconds\nAlt + .\t\t\tGo forward by $(round(zoom)) seconds\n\n[\t\t\t\tZoom in\n]\t\t\t\tZoom out\n\nCtrl + Enter\t\tReturn selected time segment\nCtrl + d\t\t\tDelete selected time segment\n\nCtrl + s\t\t\tToggle snapping\nAlt + s\t\t\tToggle scales\nAlt + m\t\t\tToggle monochromatic mode\n\nCtrl + p\t\t\tPlay current segment as audio\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tClose\n"
                end
            elseif keyval == UInt(']')
                if zoom < 30 && zoom < obj.time_pts[end] - 1
                    zoom += 1
                    bt_next.tooltip_text = "Go forward by $(round(zoom)) seconds"
                    bt_prev.tooltip_text = "Go back by $(round(zoom)) seconds"

                    signal_adj = GtkAdjustment(signal_slider)
                    signal_adj.lower = obj.time_pts[1]
                    signal_adj.upper = obj.time_pts[end] - zoom
                    @idle_add Gtk4.adjustment(signal_slider, signal_adj)

                    time_adj = GtkAdjustment(entry_time)
                    time_adj.lower = obj.time_pts[1]
                    time_adj.upper = obj.time_pts[end] - zoom
                    @idle_add Gtk4.adjustment(entry_time, time_adj)
                    draw(can)
                else
                    zoom = obj.time_pts[end]
                    bt_next.tooltip_text = "Go forward by $(round(zoom)) seconds"
                    bt_prev.tooltip_text = "Go back by $(round(zoom)) seconds"

                    signal_adj = GtkAdjustment(signal_slider)
                    signal_adj.lower = obj.time_pts[1]
                    signal_adj.upper = obj.time_pts[end] - zoom
                    @idle_add Gtk4.adjustment(signal_slider, signal_adj)

                    time_adj = GtkAdjustment(entry_time)
                    time_adj.lower = obj.time_pts[1]
                    time_adj.upper = obj.time_pts[end] - zoom
                    @idle_add Gtk4.adjustment(entry_time, time_adj)
                    draw(can)
                end
                if mch
                    help = "Keyboard shortcuts:\n\nCtrl + b\t\t\tToggle butterfly plot\nCtrl + m\t\t\tToggle mean plot\n\nCtrl + z\t\t\tScroll channels up\nCtrl + x\t\t\tScroll channels down\n\nCtrl + ,\t\t\tGo back by 1 second\nCtrl + .\t\t\tGo forward by 1 second\nAlt + ,\t\t\tGo back by $(round(zoom)) seconds\nAlt + .\t\t\tGo forward by $(round(zoom)) seconds\n\n[\t\t\t\tZoom in\n]\t\t\t\tZoom out\n\nCtrl + Enter\t\tReturn selected time segment\nCtrl + d\t\t\tDelete selected time segment\n\nCtrl + s\t\t\tToggle snapping\nAlt + s\t\t\tToggle scales\nAlt + m\t\t\tToggle monochromatic mode\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tClose\n"
                else
                    help = "Keyboard shortcuts:\n\nCtrl + z\t\t\tScroll channels up\nCtrl + x\t\t\tScroll channels down\n\nCtrl + ,\t\t\tGo back by 1 second\nCtrl + .\t\t\tGo forward by 1 second\nAlt + ,\t\t\tGo back by $(round(zoom)) seconds\nAlt + .\t\t\tGo forward by $(round(zoom)) seconds\n\n[\t\t\t\tZoom in\n]\t\t\t\tZoom out\n\nCtrl + Enter\t\tReturn selected time segment\nCtrl + d\t\t\tDelete selected time segment\n\nCtrl + s\t\t\tToggle snapping\nAlt + s\t\t\tToggle scales\nAlt + m\t\t\tToggle monochromatic mode\n\nCtrl + p\t\t\tPlay current segment as audio\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tClose\n"
                end
            end

            # ALT
            if ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_alt == mask_alt) && keyval == UInt(','))
                time_current = entry_time.value
                if time_current >= obj.time_pts[1] + zoom
                    time_current = time_current - zoom
                    @idle_add entry_time.value = time_current
                end
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_alt == mask_alt) && keyval == UInt('.'))
                time_current = entry_time.value
                if time_current < obj.time_pts[end] - zoom
                    time_current += zoom
                    @idle_add entry_time.value = time_current
                else
                    time_current = obj.time_pts[end] - zoom
                    @idle_add entry_time.value = time_current
                end
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_alt == mask_alt) && keyval == UInt('m'))
                mono = !mono
                draw(can)
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_alt == mask_alt) && keyval == UInt('s'))
                scale = !scale
                draw(can)
            end

            # CONTROL
            if ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('q'))
                close(win)
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('h'))
                info_dialog(_nill, help, win)
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt(','))
                time_current = entry_time.value
                if time_current >= obj.time_pts[1] + 1
                    time_current = time_current - 1
                    @idle_add entry_time.value = time_current
                end
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('.'))
                time_current = entry_time.value
                if time_current < obj.time_pts[end] - 1
                    time_current += 1
                    @idle_add entry_time.value = time_current
                else
                    time_current = obj.time_pts[end] - 1
                    @idle_add entry_time.value = time_current
                end
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('z'))
                if ch_first > 1
                    ch_first -= 1
                    ch_last -= 1
                    @idle_add Gtk4.value(ch_slider, ch_first)
                end
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('x'))
                if mch
                    if ch_last < length(ch)
                        ch_first += 1
                        ch_last += 1
                        @idle_add Gtk4.value(ch_slider, ch_first)
                    end
                else
                    if ch_first < length(ch)
                        ch_first += 1
                        @idle_add Gtk4.value(ch_slider, ch_first)
                    end
                end
            end
        end
    end

    app = GtkApplication("org.neuroanalyzer.iview")
    Gtk4.signal_connect(_activate, app, :activate)
    Gtk4.GLib.stop_main_loop()
    Gtk4.run(app)

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

    @assert nepochs(obj) > 1 "For continuous object iview() must be used."

    c isa Symbol && (c = _get_component(obj, c))
    ch = axes(c, 1)

    _check_epochs(obj, ep)

    mono = false
    scale = true
    k = nothing

    if length(ch) > 20
        ch_first = 1
        ch_last = 20
    else
        ch_first = 1
        ch_last = length(ch)
    end

    function _activate(app)

        win = GtkApplicationWindow(app, "NeuroAnalyzer: iview_ep()")
        Gtk4.default_size(win, p.attr[:size][1] + 40, p.attr[:size][2] + 40)

        can = GtkCanvas()
        can.content_width = p.attr[:size][1]
        can.content_height = p.attr[:size][2]
        g = GtkGrid()
        g.column_homogeneous = false
        g.column_spacing = 5
        g.row_spacing = 5

        entry_epoch = GtkSpinButton(1, nepochs(obj), 1)
        entry_epoch.value = ep
        entry_epoch.tooltip_text = "Epoch"
        bt_start = GtkButton("⇤")
        bt_start.tooltip_text = "Go to the first epoch"
        bt_end = GtkButton("⇥")
        bt_end.tooltip_text = "Go to the last epoch"
        bt_help = GtkButton("Help")
        bt_help.tooltip_text = "Show help"
        bt_close = GtkButton("Close")
        bt_close.tooltip_text = "Close this window"

        if length(ch) > 20
            ch_slider = GtkScale(:v, ch[ch_first]:(ch[end] - 19))
            ch_slider.draw_value = false
        else
            ch_slider = GtkScale(:v, ch[ch_first]:ch[end])
            ch_slider.draw_value = false
            ch_slider.sensitive = false
        end
        ch_slider.tooltip_text = "Scroll components"

        signal_slider = GtkScale(:h, 1:nepochs(obj))
        signal_slider.draw_value = false
        signal_slider.tooltip_text = "Current epoch"

        g[1:5, 1] = can
        g[6, 1] = ch_slider
        g[1:5, 2] = signal_slider
        g[1, 3] = bt_start
        g[2, 3] = entry_epoch
        g[3, 3] = bt_end
        g[4, 3] = bt_help
        g[5, 3] = bt_close
        push!(win, g)

        Gtk4.show(win)

        @guarded draw(can) do widget
            ep = Int64(entry_epoch)
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

        function _mwheel_scroll(_, dx, dy)
            if dy > 0.5
                if k == 0x0000ffe1 # shift
                    ep = Int64(entry_epoch.value)
                    if ep < nepochs(obj)
                        ep += 1
                        @idle_add entry_epoch.value = ep
                    end
                elseif k == 0x0000ffe3 # ctrl
                    if mch
                        if ch_last < length(ch)
                            ch_first += 1
                            ch_last += 1
                            @idle_add Gtk4.value(ch_slider, ch_first)
                        end
                    else
                        if ch_first < length(ch)
                            ch_first += 1
                            @idle_add Gtk4.value(ch_slider, ch_first)
                        end
                    end
                end
            elseif dy < 0.5
                if k == 0x0000ffe1 # shift
                    ep = Int64(entry_epoch.value)
                    if ep > 1
                        ep -= 1
                        @idle_add entry_epoch.value = ep
                    end
                elseif k == 0x0000ffe3 # ctrl
                    if ch_first > 1
                        ch_first -= 1
                        ch_last -= 1
                        @idle_add Gtk4.value(ch_slider, ch_first)
                    end
                end
            end
        end
        ecsf = Gtk4.GtkEventControllerScroll(Gtk4.EventControllerScrollFlags_VERTICAL)
        push!(can, ecsf)
        signal_connect(_mwheel_scroll, ecsf, "scroll")

        signal_connect(signal_slider, "value-changed") do widget
            @idle_add entry_epoch.value = round(Int64, Gtk4.value(signal_slider))
            draw(can)
        end

        signal_connect(ch_slider, "value-changed") do widget
            ch_first = round(Int64, Gtk4.value(ch_slider))
            length(ch) > 20 && (ch_last = ch_first + 19)
            draw(can)
        end

        signal_connect(entry_epoch, "value-changed") do widget
            @idle_add Gtk4.value(signal_slider, Int64(entry_epoch.value))
            draw(can)
        end

        signal_connect(bt_start, "clicked") do widget
            @idle_add Gtk4.value(entry_epoch, 1)
        end

        signal_connect(bt_end, "clicked") do widget
            @idle_add Gtk4.value(entry_epoch, nepochs(obj))
        end

        signal_connect(bt_close, "clicked") do widget
            close(win)
        end

        help = "Keyboard shortcuts:\n\nCtrl + z\t\t\tScroll channels up\nCtrl + x\t\t\tScroll channels down\n\nCtrl + ,\t\t\tPrevious epoch\nCtrl + .\t\t\tNext epoch\n\nAlt + s\t\t\tToggle scales\nAlt + m\t\t\tToggle monochromatic mode\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tClose\n"

        signal_connect(bt_help, "clicked") do widget
            info_dialog(_nill, help, win)
        end

        win_key = Gtk4.GtkEventControllerKey(win)

        signal_connect(win_key, "key-released") do widget, keyval, keycode, state
            k = nothing
        end

        signal_connect(win_key, "key-pressed") do widget, keyval, keycode, state
            k = keyval

            # ALT
            if ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_alt == mask_alt) && keyval == UInt('m'))
                mono = !mono
                draw(can)
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_alt == mask_alt) && keyval == UInt('s'))
                scale = !scale
                draw(can)
            end

            # CONTROL
            if ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('q'))
                close(win)
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('h'))
                info_dialog(_nill, help, win)
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt(','))
                ep = Int64(entry_epoch.value)
                if ep > 1
                    ep -= 1
                    @idle_add entry_epoch.value = ep
                end
                draw(can)
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('.'))
                ep = Int64(entry_epoch.value)
                if ep > 1
                    ep -= 1
                    @idle_add entry_epoch.value = ep
                end
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('z'))
                if ch_first > 1
                    ch_first -= 1
                    ch_last -= 1
                    @idle_add Gtk4.value(ch_slider, ch_first)
                end
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('x'))
                if mch
                    if ch_last < length(ch)
                        ch_first += 1
                        ch_last += 1
                        @idle_add Gtk4.value(ch_slider, ch_first)
                    end
                else
                    if ch_first < length(ch)
                        ch_first += 1
                        @idle_add Gtk4.value(ch_slider, ch_first)
                    end
                end
            end
        end
    end

    app = GtkApplication("org.neuroanalyzer.iview_ep")
    Gtk4.signal_connect(_activate, app, :activate)
    Gtk4.GLib.stop_main_loop()
    Gtk4.run(app)

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

    function _activate(app)
        win = GtkApplicationWindow(app, "NeuroAnalyzer: iview()")
        Gtk4.default_size(win, p.attr[:size][1] + 2, p.attr[:size][2] + 2)

        can = GtkCanvas()
        can.content_width = p.attr[:size][1] + 2
        can.content_height = p.attr[:size][2] + 2
        push!(win, can)

        Gtk4.show(win)

        @guarded draw(can) do widget
            show(io, MIME("image/png"), p)
            img = read_from_png(io)
            w = img.width
            h = img.height
            ctx = getgc(can)
            Cairo.set_source_surface(ctx, img, 1, 2)
            Cairo.paint(ctx)
        end

        win_key = Gtk4.GtkEventControllerKey(win)

        signal_connect(win_key, "key-pressed") do widget, keyval, keycode, state
            # CONTROL
            if ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('s'))
                save_dialog("Pick an image file", win, ["*.png"]) do file_name
                    if file_name != ""
                        surface_buf = Gtk4.cairo_surface(can)
                        if Cairo.write_to_png(surface_buf, file_name) == Cairo.STATUS_SUCCESS
                            _info("Plot saved as: $file_name")
                        else
                            warn_dialog(_nill, "File cannot be saved!", win)
                        end
                    end
                end
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('q'))
                close(win)
            end
        end
    end

    app = GtkApplication("org.neuroanalyzer.iview")
    Gtk4.signal_connect(_activate, app, :activate)
    Gtk4.GLib.stop_main_loop()
    Gtk4.run(app)

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
        _error("Incorrect filename!")
        return nothing
    end

    img = read_from_png(file_name)

    function _activate(app)
        win = GtkApplicationWindow(app, "NeuroAnalyzer: iview()")
        Gtk4.default_size(win, Int64(img.width), Int64(img.height))

        can = GtkCanvas()
        can.content_width = Int64(img.width)
        can.content_height = Int64(img.height)
        push!(win, can)
        Gtk4.show(win)

        @guarded draw(can) do widget
            ctx = getgc(can)
            Cairo.set_source_surface(ctx, img, 0, 0)
            Cairo.paint(ctx)
        end

        win_key = Gtk4.GtkEventControllerKey(win)

        signal_connect(win_key, "key-pressed") do widget, keyval, keycode, state
            # CONTROL
            if ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('q'))
                close(win)
            end
        end
    end

    app = GtkApplication("org.neuroanalyzer.iview")
    Gtk4.signal_connect(_activate, app, :activate)
    Gtk4.GLib.stop_main_loop()
    Gtk4.run(app)

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

    function _activate(app)
        win = GtkApplicationWindow(app, "NeuroAnalyzer: iview()")
        Gtk4.default_size(win, Int64(c.width), Int64(c.height))

        can = GtkCanvas()
        can.content_width = Int64(c.width)
        can.content_height = Int64(c.height)
        push!(win, can)
        Gtk4.show(win)

        @guarded draw(can) do widget
            ctx = getgc(can)
            Cairo.set_source_rgb(ctx, 255, 255, 255)
            Cairo.set_source_surface(ctx, c, 1, 1)
            Cairo.paint(ctx)
        end

        win_key = Gtk4.GtkEventControllerKey(win)

        signal_connect(win_key, "key-pressed") do widget, keyval, keycode, state
            # CONTROL
            if ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('s'))
                save_dialog("Pick an image file", win, ["*.png"]) do file_name
                    if file_name != ""
                        surface_buf = Gtk4.cairo_surface(can)
                        if Cairo.write_to_png(surface_buf, file_name) == Cairo.STATUS_SUCCESS
                            _info("Plot saved as: $file_name")
                        else
                            warn_dialog(_nill, "File cannot be saved!", win)
                        end
                    end
                end
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('q'))
                close(win)
            end
        end
    end

    app = GtkApplication("org.neuroanalyzer.iview")
    Gtk4.signal_connect(_activate, app, :activate)
    Gtk4.GLib.stop_main_loop()
    Gtk4.run(app)

    return nothing

end
