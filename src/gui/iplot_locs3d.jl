export iplot_locs3d

"""
    iplot_locs3d(locs; <keyword arguments>)

3D interactive preview of channel locations.

# Arguments

- `locs::DataFrame`: columns: channel, labels, loc_radius, loc_theta, loc_x, loc_y, loc_z, loc_radius_sph, loc_theta_sph, loc_phi_sph
- `ch::Union{Int64, Vector{Int64}}=1:DataFrames.nrow(locs)`: channel(s) to plot, default is all channels
- `selected::Union{Int64, Vector{Int64}, AbstractRange}=0`: selected channel(s) to plot
- `ch_labels::Bool=true`: plot channel labels
- `head_labels::Bool=true`: plot head labels
- `mono::Bool=false`: use color or gray palette
- `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use spherical coordinates
- `camera::Tuple{Real, Real}=(20, 45)`: camera position -- (XY plane angle, XZ plane angle)

# Returns

Nothing
"""
function iplot_locs3d(locs::DataFrame; ch::Union{Int64, Vector{Int64}}=1:DataFrames.nrow(locs), selected::Union{Int64, Vector{Int64}, AbstractRange}=0, ch_labels::Bool=true, head_labels::Bool=true, mono::Bool=false, cart::Bool=false, camera::Tuple{Real, Real}=(20, 45))::Nothing

    p = NeuroAnalyzer.plot_locs3d(locs, ch=ch, selected=selected, ch_labels=ch_labels, head_labels=head_labels, cart=cart, camera=camera)
    win = GtkWindow("NeuroAnalyzer: iplot_locs3d()", Int32(p.attr[:size][1]), Int32(p.attr[:size][2]))
    can = GtkCanvas(Int32(p.attr[:size][1]), Int32(p.attr[:size][2]))
    set_gtk_property!(win, :border_width, 0)
    set_gtk_property!(win, :resizable, false)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    set_gtk_property!(win, :startup_id, "org.neuroanalyzer")
    push!(win, can)
    showall(win)

    camera_pos = camera
    x_pos_last = 0
    y_pos_last = 0

    @guarded draw(can) do widget
        p = NeuroAnalyzer.plot_locs3d(locs, camera=camera_pos, ch=ch, selected=selected, ch_labels=ch_labels, head_labels=head_labels, cart=cart, mono=mono);
        img = read_from_png(io)
        ctx = getgc(can)
        show(io, MIME("image/png"), p)
        img = read_from_png(io)
        set_source_surface(ctx, img, 0, 0)
        paint(ctx)
    end

    can.mouse.button1motion = @guarded (widget, event) -> begin
        x_pos = round(Int64, event.x)
        y_pos = round(Int64, event.y)
        x_pos > x_pos_last && (camera_pos = (camera_pos[1] + 5, camera_pos[2]))
        x_pos < x_pos_last && (camera_pos = (camera_pos[1] - 5, camera_pos[2]))
        y_pos > y_pos_last && (camera_pos = (camera_pos[1], camera_pos[2] + 5))
        y_pos < y_pos_last && (camera_pos = (camera_pos[1], camera_pos[2] - 5))
        x_pos_last = x_pos
        y_pos_last = y_pos
        draw(can)
    end

    can.mouse.button3press = @guarded (widget, event) -> begin
        camera_pos = camera
        draw(can)
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        s = event.state

        if s == 0x00000004 || s == 0x00000014 # ctrl
            if k == 0x00000071 # q
                Gtk.destroy(win)
            end
        end

        if  k == 116 # t
            camera_pos = (0, 90)
            draw(can)
        elseif k == 115 # s
            camera_pos = (90, 0)
            draw(can)
        elseif k == 102 # f
            camera_pos = (180, 0)
            draw(can)
        end
    end

    return nothing

end

"""
    iplot_locs3d(obj; <keyword arguments>)

3D interactive preview of channel locations.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel(s) to plot
- `selected::Union{String, Vector{String}}=""`: selected channel(s) to plot
- `ch_labels::Bool=true`: plot channel labels
- `head_labels::Bool=true`: plot head labels
- `mono::Bool=false`: use color or gray palette
- `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use spherical coordinates
- `camera::Tuple{Real, Real}=(20, 45)`: camera position -- (XY plane angle, XZ plane angle)

# Returns

Nothing
"""
function iplot_locs3d(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, selected::Union{String, Vector{String}}="", ch_labels::Bool=true, head_labels::Bool=true, mono::Bool=false, cart::Bool=false, camera::Tuple{Real, Real}=(20, 45))::Nothing

    # select channels, default is all channels
    ch = get_channel(obj, ch=ch)

    # get selected channels
    if selected == ""
        selected = 0
    else
        selected = get_channel(obj, ch=selected)
        selected = intersect(locs[!, :label], labels(obj)[selected])
        selected = _find_bylabel(locs, selected)
    end

    p = NeuroAnalyzer.plot_locs3d(obj.locs, ch=ch, selected=selected, ch_labels=ch_labels, head_labels=head_labels, cart=cart, camera=camera)
    win = GtkWindow("NeuroAnalyzer: iplot_locs3d()", Int32(p.attr[:size][1]), Int32(p.attr[:size][2]))
    can = GtkCanvas(Int32(p.attr[:size][1]), Int32(p.attr[:size][2]))
    set_gtk_property!(win, :border_width, 0)
    set_gtk_property!(win, :resizable, false)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    set_gtk_property!(win, :startup_id, "org.neuroanalyzer")
    push!(win, can)
    showall(win)

    camera_pos = camera
    x_pos_last = 0
    y_pos_last = 0

    @guarded draw(can) do widget
        p = plot_locs3d(obj.locs, camera=camera_pos, ch=ch, selected=selected, ch_labels=ch_labels, head_labels=head_labels, cart=cart, mono=mono);
        img = read_from_png(io)
        ctx = getgc(can)
        show(io, MIME("image/png"), p)
        img = read_from_png(io)
        set_source_surface(ctx, img, 0, 0)
        paint(ctx)
    end

    can.mouse.button1motion = @guarded (widget, event) -> begin
        x_pos = round(Int64, event.x)
        y_pos = round(Int64, event.y)
        x_pos > x_pos_last && (camera_pos = (camera_pos[1] + 5, camera_pos[2]))
        x_pos < x_pos_last && (camera_pos = (camera_pos[1] - 5, camera_pos[2]))
        y_pos > y_pos_last && (camera_pos = (camera_pos[1], camera_pos[2] + 5))
        y_pos < y_pos_last && (camera_pos = (camera_pos[1], camera_pos[2] - 5))
        x_pos_last = x_pos
        y_pos_last = y_pos
        draw(can)
    end

    can.mouse.button3press = @guarded (widget, event) -> begin
        camera_pos = camera
        draw(can)
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        s = event.state

        if s == 0x00000004 || s == 0x00000014 # ctrl
            if k == 0x00000071 # q
                Gtk.destroy(win)
            end
        end

        if  k == 116 # t
            camera_pos = (0, 90)
            draw(can)
        elseif k == 115 # s
            camera_pos = (90, 0)
            draw(can)
        elseif k == 102 # f
            camera_pos = (180, 0)
            draw(can)
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
