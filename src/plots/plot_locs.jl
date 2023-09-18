export plot_locs
export plot_locs3d
export iplot_locs3d

"""
    plot_locs(locs; <keyword arguments>)

Preview channel locations.

# Arguments

- `locs::DataFrame`: columns: channel, labels, loc_theta, loc_radius, loc_x, loc_y, loc_z, loc_radius_sph, loc_theta_sph, loc_phi_sph
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1:nrow(locs)`: channel(s) to plot, default is all channels
- `selected::Union{Int64, Vector{Int64}, <:AbstractRange}=0`: selected channel(s) to plot
- `ch_labels::Bool=true`: plot channel labels
- `head::Bool=true`: draw head
- `head_labels::Bool=true`: plot head labels
- `mono::Bool=false`: Use color or gray palette
- `head_details::Bool=true`: draw nose and ears
- `grid::Bool=false`: draw grid, useful for locating positions
- `plot_size::Int64=400`: plot dimensions in pixels (size × size)
- `cart::Bool=false`: if true, use Cartesian x and y coordinates, otherwise use polar radius and theta coordinates

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_locs(locs::DataFrame; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1:nrow(locs), selected::Union{Int64, Vector{Int64}, <:AbstractRange}=0, ch_labels::Bool=true, head::Bool=true, head_labels::Bool=true, mono::Bool=false, head_details::Bool=true, grid::Bool=false, plot_size::Int64=400, cart::Bool=false)

    pal = mono == true ? :grays : :darktest

    if cart == false
        loc_x = zeros(size(locs, 1))
        loc_y = zeros(size(locs, 1))
        for idx in 1:size(locs, 1)
            loc_x[idx], loc_y[idx] = pol2cart(locs[!, :loc_radius][idx], locs[!, :loc_theta][idx])
        end
    else
        loc_x = locs[!, :loc_x]
        loc_y = locs[!, :loc_y]
    end
    loc_x = _s2v(loc_x)
    loc_y = _s2v(loc_y)

    if plot_size > 300
        marker_size = plot_size ÷ 75
        font_size = plot_size ÷ 75
    else
        marker_size = plot_size ÷ 50
        font_size = plot_size ÷ 50
        ch_labels = false
    end

    length(ch) > 64 && (font_size = plot_size ÷ 100)

    if grid == false
        p = Plots.plot(grid=false,
                       framestyle=:none,
                       palette=pal,
                       size=(plot_size, plot_size),
                       border=:none,
                       aspect_ratio=1,
                       right_margin=-30*Plots.px,
                       bottom_margin=-20*Plots.px,
                       top_margin=-30*Plots.px,
                       left_margin=-50*Plots.px,
                       xlim=(-1.22, 1.23),
                       ylim=(-1.1, 1.2))
    else
        p = Plots.plot(grid=true,
                       palette=pal,
                       size=(plot_size, plot_size),
                       aspect_ratio=1,
                       right_margin=-30*Plots.px,
                       bottom_margin=-50*Plots.px,
                       top_margin=-50*Plots.px,
                       left_margin=-5*Plots.px,
                       xticks=-1:0.1:1,
                       yticks=-1:0.1:1,
                       xtickfontsize=4,
                       ytickfontsize=4;
                       xlim=(-1.22, 1.23),
                       ylim=(-1.1, 1.2))
    end

    if head == true
        hd = _draw_head(p, head_labels=head_labels, head_details=head_details)
        p = Plots.plot!(hd)
    end

    for idx in eachindex(locs[!, :labels])
        if idx in ch
            if selected != 0
                p = Plots.scatter!((loc_x[idx], loc_y[idx]),
                                color=:lightgrey,
                                markerstrokecolor = Colors.RGBA(255/255, 255/255, 255/255, 0/255),
                                grid=true,
                                label="",
                                markershape=:circle,
                                markersize=marker_size,
                                markerstrokewidth=0,
                                markerstrokealpha=0)
            else
                p = Plots.scatter!((loc_x[idx], loc_y[idx]),
                                color=:lightgrey,
                                markerstrokecolor = Colors.RGBA(255/255, 255/255, 255/255, 0/255),
                                grid=true,
                                label="",
                                markershape=:circle,
                                markersize=marker_size,
                                markerstrokewidth=0,
                                markerstrokealpha=0)
            end
        end
        if idx in selected
            if mono != true
                p = Plots.scatter!((loc_x[idx], loc_y[idx]),
                                color=idx,
                                markerstrokecolor = Colors.RGBA(255/255, 255/255, 255/255, 0/255),
                                grid=true,
                                label="",
                                markershape=:circle,
                                markersize=marker_size,
                                markerstrokewidth=0,
                                markerstrokealpha=0)
            else
                #p = Plots.plot!((loc_x[idx], loc_y[idx]),
                p = Plots.scatter!((loc_x[idx], loc_y[idx]),
                                color=:lightgrey,
                                markerstrokecolor = Colors.RGBA(255/255, 255/255, 255/255, 0/255),
                                grid=true,
                                label="",
                                markershape=:circle,
                                markersize=marker_size,
                                markerstrokewidth=0,
                                markerstrokealpha=0)
            end
        end
    end
    if ch_labels
        for idx in eachindex(locs[!, :labels])
            if idx in ch
                Plots.plot!(annotations=(loc_x[idx], loc_y[idx] + 0.075, Plots.text(locs[!, :labels][idx], pointsize=font_size)))
                Plots.plot!(annotations=(loc_x[idx], loc_y[idx] + 0.075, Plots.text(locs[!, :labels][idx], pointsize=font_size)))
            end
            if idx in selected
                Plots.plot!(annotations=(loc_x[idx], loc_y[idx] + 0.075, Plots.text(locs[!, :labels][idx], pointsize=font_size)))
            end
        end
    end

    Plots.plot(p)

    return p

end

"""
    plot_locs3d(locs; <keyword arguments>)

3D preview of channel locations.

# Arguments

- `locs::DataFrame`: columns: channel, labels, loc_theta, loc_radius, loc_x, loc_y, loc_z, loc_radius_sph, loc_theta_sph, loc_phi_sph
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1:nrow(locs)`: channel(s) to plot, default is all channels
- `selected::Union{Int64, Vector{Int64}, <:AbstractRange}=0`: selected channel(s) to plot
- `ch_labels::Bool=true`: plot channel labels
- `head_labels::Bool=true`: plot head labels
- `mono::Bool=false`: Use color or gray palette
- `plot_size::Int64=800`: plot dimensions in pixels (plot_size×plot_size)
- `cart::Bool=false`: if true, use Cartesian x, y and z coordinates, otherwise use spherical radius, theta and phi coordinates
- `camera::Tuple{Real, Real}=(20, 45)`: camera position -- (X-Y plane angle, X-Z plane angle)

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_locs3d(locs::DataFrame; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1:nrow(locs), selected::Union{Int64, Vector{Int64}, <:AbstractRange}=0, ch_labels::Bool=true, head_labels::Bool=true, mono::Bool=false, plot_size::Int64=800, cart::Bool=false, camera::Tuple{Real, Real}=(20, 45))

    pal = mono == true ? :grays : :darktest

    if cart == false
        loc_x = zeros(nrow(locs))
        loc_y = zeros(nrow(locs))
        loc_z = zeros(nrow(locs))
        for idx in 1:nrow(locs)
            loc_x[idx], loc_y[idx], loc_z[idx] = sph2cart(locs[idx, :loc_radius_sph], locs[idx, :loc_theta_sph], locs[idx, :loc_phi_sph])
        end
    else
        loc_x = locs[!, :loc_x]
        loc_y = locs[!, :loc_y]
        loc_z = locs[!, :loc_z]
    end

    x_lim = (-1.1, 1.1)
    y_lim = (-1.1, 1.1)
    z_lim = (-1.1, 1.1)

    marker_size = plot_size ÷ 200
    font_size = plot_size ÷ 75

    p = Plots.scatter3d(grid=true,
                        palette=pal,
                        size=(plot_size, plot_size),
                        aspect_ratios=:equal,
                        right_margin=-20*Plots.px,
                        bottom_margin=-20*Plots.px,
                        top_margin=-20*Plots.px,
                        left_margin=-20*Plots.px,
                        legend=false,
                        camera=camera,
                        xticks=([-1, 0, 1]),
                        yticks=([-1, 0, 1]),
                        zticks=([-1, 0, 1]),
                        xlabel="X",
                        ylabel="Y",
                        zlabel="Z",
                        xlim=x_lim,
                        ylim=y_lim,
                        zlim=z_lim)

    p = Plots.scatter3d!(loc_x[ch], loc_y[ch], loc_z[ch], ms=marker_size, mc=:gray)

    if selected != 0
        if length(selected) > 1
            if mono == true
                p = Plots.scatter3d!(loc_x[selected], loc_y[selected], loc_z[selected], ms=marker_size, mc=:gray)
            else
                p = Plots.scatter3d!(loc_x[selected], loc_y[selected], loc_z[selected], ms=marker_size, mc=:red)
            end
        else
            if mono == true
                p = Plots.scatter3d!((loc_x[selected], loc_y[selected], loc_z[selected]), ms=marker_size, mc=:gray)
            else
                p = Plots.scatter3d!((loc_x[selected], loc_y[selected], loc_z[selected]), ms=marker_size, mc=:red)
            end
        end
    end

    if ch_labels == true
        for idx in eachindex(locs[!, :labels])
            if idx in ch
                Plots.annotate!(loc_x[idx] * 1.1, loc_y[idx] * 1.1, loc_z[idx] * 1.1, Plots.text(locs[!, :labels][idx], font_size))
            end
            if idx in selected
                Plots.annotate!(loc_x[idx] * 1.1, loc_y[idx] * 1.1, loc_z[idx] * 1.1, Plots.text(locs[!, :labels][idx], font_size))
            end
        end
    end

    if head_labels == true
        Plots.annotate!(0, 1.5, 0, Plots.text("Nz", font_size))
        Plots.annotate!(0, -1.5, 0, Plots.text("In", font_size))
        Plots.annotate!(-1.5, 0, 0, Plots.text("LPA", font_size))
        Plots.annotate!(1.5, 0, 0, Plots.text("RPA", font_size))
        Plots.annotate!(0, 0, 1.5, Plots.text("top", font_size))
    end
    
    Plots.plot!(p)

    return p

end

"""
    plot_locs(obj; <keyword arguments>)

Preview of channel locations.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `selected::Union{Int64, Vector{Int64}, <:AbstractRange}=0`: which channel should be highlighted
- `ch_labels::Bool=true`: plot channel labels
- `src_labels::Bool=false`: plot source labels
- `det_labels::Bool=false`: plot detector labels
- `opt_labels::Bool=false`: plot optode type (S for source, D for detector) and number
- `head::Bool=true`: draw head
- `head_labels::Bool=false`: plot head labels
- `threed::Bool=false`: 3-dimensional plot
- `plot_size::Int64=threed ? 800 : 400`: plot dimensions in pixels (plot_size×plot_size)
- `head_details::Bool=true`: draw nose and ears
- `mono::Bool=false`: Use color or gray palette
- `grid::Bool=false`: draw grid, useful for locating positions
- `cart::Bool=false`: if true, use polar coordinates, otherwise use Cartesian spherical x and y coordinates
- `interactive::Bool=true`: if true, use interactive 3-dimensional plot
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_locs(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), selected::Union{Int64, Vector{Int64}, <:AbstractRange}=0, ch_labels::Bool=true, src_labels::Bool=false, det_labels::Bool=false, opt_labels::Bool=false, head::Bool=true, head_labels::Bool=false, threed::Bool=false, plot_size::Int64=threed ? 800 : 400, head_details::Bool=true, mono::Bool=false, grid::Bool=false, cart::Bool=false, interactive::Bool=true, kwargs...)

    # select channels, default is all channels
    _check_channels(obj, ch, Symbol(obj.header.recording[:data_type]))
    selected != 0 && _check_channels(obj, selected)

    if obj.header.recording[:data_type] == "ecog"
        @error "Use plot_locs_ecog() for ECoG data."
    elseif obj.header.recording[:data_type] == "nirs"
        @error "Use plot_locs_nirs() for NIRS data."
    elseif threed == false
        if obj.header.recording[:data_type] == "nirs"
            ch_pairs = obj.header.recording[:channel_pairs]
            src_n = length(unique(ch_pairs[:, 1]))
            det_n = length(unique(ch_pairs[:, 2]))
            p = plot_locs_nirs(obj.locs, ch_pairs, src_n, det_n; src_labels=src_labels, det_labels=det_labels, opt_labels=opt_labels, head=head, head_labels=head_labels, head_details=head_details, plot_size=plot_size, grid=grid, mono=mono)
        else
            p = plot_locs(obj.locs, ch=ch, selected=selected, ch_labels=ch_labels, head=head, head_labels=head_labels, head_details=head_details, plot_size=plot_size, grid=grid, mono=mono, cart=cart)
        end
    else
        if interactive
            iplot_locs3d(obj.locs, ch=ch, selected=selected, ch_labels=ch_labels, head_labels=head_labels, mono=mono, plot_size=plot_size)
            return
        else
            p = plot_locs3d(obj.locs, ch=ch, selected=selected, ch_labels=ch_labels, head_labels=head_labels, mono=mono, plot_size=plot_size)
        end
    end

    return p
    
end

"""
    iplot_locs3d(locs; <keyword arguments>)

3D interactive preview of channel locations.

# Arguments

- `locs::DataFrame`: columns: channel, labels, loc_theta, loc_radius, loc_x, loc_y, loc_z, loc_radius_sph, loc_theta_sph, loc_phi_sph
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1:nrow(locs)`: channel(s) to plot, default is all channels
- `selected::Union{Int64, Vector{Int64}, <:AbstractRange}=0`: selected channel(s) to plot
- `ch_labels::Bool=true`: plot channel labels
- `head_labels::Bool=true`: plot head labels
- `mono::Bool=false`: Use color or gray palette
- `plot_size::Int64=800`: plot dimensions in pixels (plot_size×plot_size)
- `cart::Bool=false`: if true, use Cartesian x, y and z coordinates, otherwise use spherical radius, theta and phi coordinates
- `camera::Tuple{Real, Real}=(20, 45)`: camera position -- (X-Y plane angle, X-Z plane angle)
"""
function iplot_locs3d(locs::DataFrame; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1:nrow(locs), selected::Union{Int64, Vector{Int64}, <:AbstractRange}=0, ch_labels::Bool=true, head_labels::Bool=true, mono::Bool=false, plot_size::Int64=800, cart::Bool=false, camera::Tuple{Real, Real}=(20, 45))

    p = NeuroAnalyzer.plot_locs3d(locs, ch=ch, selected=selected, ch_labels=ch_labels, head_labels=head_labels, plot_size=plot_size, cart=cart, camera=camera)
    win = GtkWindow("NeuroAnalyzer: iplot_locs3d()", Int32(p.attr[:size][1]), Int32(p.attr[:size][2]))
    can = GtkCanvas(Int32(p.attr[:size][1]), Int32(p.attr[:size][2]))
    set_gtk_property!(win, :border_width, 0)
    set_gtk_property!(win, :resizable, false)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    push!(win, can)
    showall(win)

    camera_pos = camera
    x_pos_last = 0
    y_pos_last = 0

    @guarded draw(can) do widget
        p = NeuroAnalyzer.plot_locs3d(locs, camera=camera_pos, ch=ch, selected=selected, ch_labels=ch_labels, head_labels=head_labels, plot_size=plot_size, cart=cart);
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

    return nothing

end


"""
    iplot_locs3d(obj; <keyword arguments>)

3D interactive preview of channel locations.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1:nrow(locs)`: channel(s) to plot, default is all channels
- `selected::Union{Int64, Vector{Int64}, <:AbstractRange}=0`: selected channel(s) to plot
- `ch_labels::Bool=true`: plot channel labels
- `head_labels::Bool=true`: plot head labels
- `mono::Bool=false`: Use color or gray palette
- `plot_size::Int64=800`: plot dimensions in pixels (plot_size×plot_size)
- `cart::Bool=false`: if true, use Cartesian x, y and z coordinates, otherwise use spherical radius, theta and phi coordinates
- `camera::Tuple{Real, Real}=(20, 45)`: camera position -- (X-Y plane angle, X-Z plane angle)
"""
function iplot_locs3d(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1:nrow(obj.locs), selected::Union{Int64, Vector{Int64}, <:AbstractRange}=0, ch_labels::Bool=true, head_labels::Bool=true, mono::Bool=false, plot_size::Int64=800, cart::Bool=false, camera::Tuple{Real, Real}=(20, 45))

    # select channels, default is all channels
    _check_channels(obj, ch, Symbol(obj.header.recording[:data_type]))
    selected != 0 && _check_channels(obj, selected)

    p = NeuroAnalyzer.plot_locs3d(obj.locs, ch=ch, selected=selected, ch_labels=ch_labels, head_labels=head_labels, plot_size=plot_size, cart=cart, camera=camera)
    win = GtkWindow("NeuroAnalyzer: iplot_locs3d()", Int32(p.attr[:size][1]), Int32(p.attr[:size][2]))
    can = GtkCanvas(Int32(p.attr[:size][1]), Int32(p.attr[:size][2]))
    set_gtk_property!(win, :border_width, 0)
    set_gtk_property!(win, :resizable, false)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    push!(win, can)
    showall(win)

    camera_pos = camera
    x_pos_last = 0
    y_pos_last = 0

    @guarded draw(can) do widget
        p = NeuroAnalyzer.plot_locs3d(obj.locs, camera=camera_pos, ch=ch, selected=selected, ch_labels=ch_labels, head_labels=head_labels, plot_size=plot_size, cart=cart);
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

    return nothing
end
