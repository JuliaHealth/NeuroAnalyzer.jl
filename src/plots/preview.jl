export preview

"""
    preview(obj; <keyword arguments>)

Interactive preview by epochs.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: channel(s) to plot, default is all channels
- `mono::Bool=false`: use color or grey palette

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function preview(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=NeuroAnalyzer._c(channel_n(obj)), mono::Bool=false)

    win = GtkWindow("NeuroAnalyzer", 1200, 900)
    can = GtkCanvas(1200, 800)
    g = GtkGrid()
    slider_ep = GtkScale(false, 1:epoch_n(obj))
    bt_start = GtkButton("First epoch [<]")
    bt_prev = GtkButton("Prev epoch [,]")
    bt_next = GtkButton("Next epoch [.]")
    bt_end = GtkButton("Last epoch [>]")
    bt_delete = GtkButton("Delete epoch [DEL]")
    bt_close = GtkButton("Close [q]")

    @guarded draw(can) do widget
        ctx = getgc(can)
        ep = Int(GAccessor.value(slider_ep))
        # refresh(obj, ch, ep)
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

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        # println(k)
        if k == 113 # Q
            Gtk.destroy(win)
        elseif k == 60 # <
            Gtk.@sigatom begin
                GAccessor.value(slider_ep, 1)
            end
        elseif k == 62 # >
            Gtk.@sigatom begin
                GAccessor.value(slider_ep, epoch_n(obj))
            end
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
                if ask_dialog("Delete epoch $ep?", "No", "Yes")
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

    g[1:6, 1] = can         # spans all columns
    g[1:6, 2] = slider_ep   # spans all columns
    g[1, 3] = bt_start      # Cartesian coordinates, g[x,y]
    g[2, 3] = bt_prev
    g[3, 3] = bt_next
    g[4, 3] = bt_end
    g[5, 3] = bt_delete
    g[6, 3] = bt_close
    set_gtk_property!(g, :column_homogeneous, true)
    set_gtk_property!(g, :column_spacing, 10)  # introduce a 5-pixel gap between columns
    push!(win, g)

    showall(win)
    show(can)

end
