export iplot_locs3d

"""
    iplot_locs3d(locs; <keyword arguments>)

3D interactive preview of channel locations.

# Arguments

- `locs::DataFrame`: columns: channel, labels, loc_radius, loc_theta, loc_x, loc_y, loc_z, loc_radius_sph, loc_theta_sph, loc_phi_sph
- `ch::Union{Int64, Vector{Int64}, AbstractRange}=1:DataFrames.nrow(locs)`: channel(s) to plot, default is all channels
- `selected::Union{Int64, Vector{Int64}, AbstractRange}=0`: selected channel(s) to plot
- `ch_labels::Bool=true`: plot channel labels
- `head_labels::Bool=true`: plot head labels
- `mono::Bool=false`: use color or gray palette
- `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use spherical coordinates
- `camera::Tuple{Real, Real}=(20, 45)`: camera position - (XY plane angle, XZ plane angle)

# Returns

Nothing
"""
function iplot_locs3d(locs::DataFrame; ch::Union{Int64, Vector{Int64}, AbstractRange}=1:DataFrames.nrow(locs), selected::Union{Int64, Vector{Int64}, AbstractRange}=0, ch_labels::Bool=true, head_labels::Bool=true, mono::Bool=false, cart::Bool=false, camera::Tuple{Real, Real}=(20, 45))::Nothing

    camera_pos = camera
    x_pos_last = 0
    y_pos_last = 0

    p = NeuroAnalyzer.plot_locs3d(locs, ch=ch, selected=selected, ch_labels=ch_labels, head_labels=head_labels, cart=cart, camera=camera)

    function _activate(app)

        win = GtkApplicationWindow(app, "NeuroAnalyzer: iplot_locs3d()")
        win.width_request = p.attr[:size][1]
        win.height_request = p.attr[:size][2]

        can = GtkCanvas(p.attr[:size][1], p.attr[:size][2])
        push!(win, can)

        Gtk4.show(win)

        @guarded draw(can) do widget
            p = NeuroAnalyzer.plot_locs3d(locs, camera=camera_pos, ch=ch, selected=selected, ch_labels=ch_labels, head_labels=head_labels, cart=cart, mono=mono)
            img = read_from_png(io)
            ctx = getgc(can)
            show(io, MIME("image/png"), p)
            img = read_from_png(io)
            set_source_surface(ctx, img, 0, 0)
            paint(ctx)
        end

        function _lmb_click(_, x, y)
            x_pos = round(Int64, x)
            y_pos = round(Int64, y)
            x_pos > x_pos_last && (camera_pos = (camera_pos[1] - 5, camera_pos[2]))
            x_pos < x_pos_last && (camera_pos = (camera_pos[1] + 5, camera_pos[2]))
            y_pos > y_pos_last && (camera_pos = (camera_pos[1], camera_pos[2] + 5))
            y_pos < y_pos_last && (camera_pos = (camera_pos[1], camera_pos[2] - 5))
            x_pos_last = x_pos
            y_pos_last = y_pos
            draw(can)
        end
        ggc_l = GtkGestureDrag()
        ggc_l.button = 1
        push!(can, ggc_l)
        signal_connect(_lmb_click, ggc_l, "drag-begin")
        signal_connect(_lmb_click, ggc_l, "drag-update")

        function _rmb_click(_, _, x, y)
            camera_pos = camera
            draw(can)
        end
        ggc_r = GtkGestureClick()
        ggc_r.button = 3
        push!(can, ggc_r)
        signal_connect(_rmb_click, ggc_r, "pressed")

        help = "Keyboard shortcuts:\n\nLeft click\t\tRotate view\nRight click\t\tReset view\n\nCtrl + t\t\t\tTop view\nCtrl + s\t\t\tSide view\nCtrl + f\t\t\tFront view\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tClose\n"

        win_key = Gtk4.GtkEventControllerKey(win)
        signal_connect(win_key, "key-pressed") do widget, keyval, keycode, state
            if ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('q'))
                close(win)
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('h'))
                info_dialog(help, win) do
                    nothing
                end
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('t'))
                camera_pos = (0, 90)
                draw(can)
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('s'))
                camera_pos = (90, 0)
                draw(can)
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('f'))
                camera_pos = (180, 0)
                draw(can)
            end
        end
    end

    app = GtkApplication("org.neuroanalyzer.iplot_locs3d")
    Gtk4.signal_connect(_activate, app, :activate)
    Gtk4.GLib.stop_main_loop()
    Gtk4.run(app)

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
- `camera::Tuple{Real, Real}=(20, 45)`: camera position - (XY plane angle, XZ plane angle)

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

    iplot_locs3d(obj.locs, ch=ch, selected=selected, ch_labels=ch_labels, head_labels=head_labels, mono=mono, cart=cart, camera=camera)

    return nothing

end
