export plot_locs3d

"""
    plot_locs3d_mesh(locs; <keyword arguments>)

3D interactive preview of channel locations with optional brain/head mesh visualization.

# Arguments

- `locs::DataFrame`: channel location data
- `ch::Union{Int64, Vector{Int64}}=1:DataFrames.nrow(locs)`: channels to plot, default is all channels
- `sch::Union{Int64, Vector{Int64}, AbstractRange}=0`: significant channels to highlight
- `ch_labels::Bool=true`: if `true`, plot channel labels
- `head_labels::Bool=true`: if `true`, plot head labels
- `mono::Bool=false`: if `true`, use a monochrome palette
- `cart::Bool=false`: if `true`, use Cartesian coordinates, otherwise use spherical coordinates
- `cam::Tuple{Real, Real}=(20, 45)`: camera position as (XY plane angle, XZ plane angle) in degrees
- `mesh_type::Symbol=:disabled`: type of mesh to plot:
    - `:disabled`: no mesh
    - `:brain`: brain mesh
    - `:head`: head mesh
- `mesh_alpha::Float64=0.95`: mesh opacity (0.0-1.0, higher = more opaque)
- `gui::Bool=true`: if `true`, keep window open and interactive

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_locs3d(
    locs::DataFrame;
    ch::Union{Int64, Vector{Int64}, AbstractRange} = 1:DataFrames.nrow(locs),
    sch::Union{Int64, Vector{Int64}, AbstractRange} = 0,
    ch_labels::Bool = true,
    head_labels::Bool = true,
    mono::Bool = false,
    cart::Bool = false,
    cam::Tuple{Real, Real} = (20, 45),
    mesh_type::Symbol = :disabled,
    mesh_alpha::Float64 = 0.95,
    gui::Bool = true,
)::GLMakie.Figure
    # validate
    _check_var(mesh_type, [:disabled, :brain, :head], "mesh_type")
    _in(mesh_alpha, (0.0, 1.0), "mesh_alpha")

    msh = nothing

    # load meshes
    if mesh_type !== :disabled
        if mesh_type === :brain
            msh = FileIO.load(joinpath(res_path, "mesh/brain_hires.stl"))
        else
            msh = FileIO.load(joinpath(res_path, "mesh/head.stl"))
        end
        if mesh_type === :brain
            msh.position ./= _mesh_normalize_xyz(msh)
            msh.position .*= 0.95
        else
            msh.position ./= _mesh_normalize_xy(msh)
            msh.position .*= 1.1
        end
    end

    # set color palette
    pal = mono ? :grays : :darktest
    effective_alpha = mesh_type === :head ? 1.0 : mesh_alpha

    if !cart
        loc_x = zeros(DataFrames.nrow(locs))
        loc_y = zeros(DataFrames.nrow(locs))
        loc_z = zeros(DataFrames.nrow(locs))
        for idx in axes(locs, 1)
            loc_x[idx], loc_y[idx], loc_z[idx] = sph2cart(
                locs.loc_radius_sph,
                locs.loc_theta_sph,
                locs.loc_phi_sph,
            )
        end
    else
        loc_x = locs.loc_x
        loc_y = locs.loc_y
        loc_z = locs.loc_z
    end

    if maximum(abs.(locs.loc_x)) <= 1.2 &&
       maximum(abs.(locs.loc_y)) <= 1.2 &&
       maximum(abs.(locs.loc_z)) <= 1.5
        x_lim = (-1.5, 1.5)
        y_lim = (-1.5, 1.5)
        z_lim = (-1.5, 1.5)
    else
        x_lim = (-2.0, 2.0)
        y_lim = (-2.0, 2.0)
        z_lim = (-2.0, 2.0)
    end

    plot_size   = (850, 850)
    marker_size = length(ch) > 64 ? 8 : 16
    font_size   = 14

    # prepare plot
    GLMakie.activate!(; title = "plot_locs3d()")
    fig = GLMakie.Figure(; size = plot_size)

    # create axis with customizable properties
    ax = GLMakie.Axis3(
        fig[1, 1];
        xlabel = "X",
        ylabel = "Y",
        zlabel = "Z",
        limits = (x_lim, y_lim, z_lim),
        title = "",
        aspect = (1, 1, 1),
        xticks = [-1, 0, 1],
        yticks = [-1, 0, 1],
        zticks = [-1, 0, 1],
        elevation = deg2rad(cam[1]),
        azimuth = deg2rad(cam[2]),
    )

    # draw mesh
    mesh_type !== :disabled && GLMakie.mesh!(msh; alpha = effective_alpha, color = :gray)

    ch_n = length(ch)
    cmap = GLMakie.resample_cmap(pal, ch_n)
    ch = setdiff(ch, sch)

    sch_set = Set(sch)
    for (i, idx) in enumerate(ch)
        if idx in sch_set
            GLMakie.scatter!(
                loc_x[idx], loc_y[idx], loc_z[idx];
                markersize  = marker_size,
                color       = mono ? :gray : cmap[i],
                colormap    = pal,
                colorrange  = 1:ch_n,
                strokewidth = 1,
                strokecolor = :black,
            )
        else
            GLMakie.scatter!(
                loc_x[idx], loc_y[idx], loc_z[idx];
                markersize  = marker_size,
                color       = :gray,
                strokewidth = 1,
                strokecolor = :black,
            )
        end
    end

    if ch_labels
        ch_vec  = collect(ch)
        sch_vec = collect(sch)
        GLMakie.text!(
            loc_x[ch_vec] .* 1.15,
            loc_y[ch_vec] .* 1.15,
            loc_z[ch_vec] .* 1.15;
            text     = locs[ch_vec, :label],
            fontsize = font_size,
            align    = (:center, :center),
        )
        if !isempty(sch_vec)
            GLMakie.text!(
                loc_x[sch_vec] .* 1.15,
                loc_y[sch_vec] .* 1.15,
                loc_z[sch_vec] .* 1.15;
                text     = locs[sch_vec, :label],
                fontsize = font_size,
                align    = (:center, :center),
            )
        end
    end

    head_labels && _draw_head_labels!(ax; font_size = font_size)

    # mouse events
    on(events(fig).keyboardbutton) do event
        if event.action == Keyboard.press
            new_pov = nothing
            event.key == Keyboard.home && (new_pov = (20, 45))
            event.key == Keyboard.r && (new_pov = (10, 10))
            event.key == Keyboard.l && (new_pov = (10, 190))
            event.key == Keyboard.f && (new_pov = (10, 100))
            event.key == Keyboard.b && (new_pov = (10, 280))
            event.key == Keyboard.t && (new_pov = (90, 270))
            if event.key == Keyboard.s
                save_dialog("Pick an image file", nothing, ["*.png"]) do file_name
                    if file_name != ""
                        GLMakie.save(file_name, fig.scene)
                        @info "Image saved as $file_name"
                    end
                end
            end
            if !isnothing(new_pov)
                ax.elevation[] = deg2rad(new_pov[1])
                ax.azimuth[]   = deg2rad(new_pov[2])
            end
        end
    end

    gui && wait(display(fig))

    return fig
end

"""
    plot_locs3d(obj; <keyword arguments>)

3D interactive preview of channel locations from a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `sch::Union{String, Vector{String}, Regex}`: significant channels to highlight
- `ch_labels::Bool=true`: if `true`, plot channel labels
- `head_labels::Bool=true`: if `true`, plot head labels
- `cart::Bool=false`: if `true`, use Cartesian coordinates, otherwise use polar coordinates for XY plane and spherical coordinates for XZ and YZ planes
- `cam::Tuple{Real, Real}=(20, 45)`: camera position as (XY plane angle, XZ plane angle) in degrees
- `mesh_type::Symbol=:disabled`: type of mesh to plot:
    - `:disabled`: no mesh
    - `:brain`: brain mesh
    - `:head`: head mesh
- `mesh_alpha::Float64=0.95`: mesh opacity (0.0-1.0, higher = more opaque)
- `gui::Bool=true`: if `true`, keep window open and interactive

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_locs3d(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    sch::Union{String, Vector{String}, Regex} = "",
    ch_labels::Bool = true,
    head_labels::Bool = false,
    cart::Bool = false,
    mono::Bool = false,
    cam::Tuple{Real, Real} = (20, 45),
    mesh_type::Symbol = :disabled,
    mesh_alpha::Float64 = 0.95,
    gui::Bool = true,
)::GLMakie.Figure
    # validate
    datatype(obj) in ["eeg"] ||
        throw(ArgumentError("Currently plot_locs3d() works for EEG objects only."))

    # resolve channel names to integer indices
    ch   = get_channel(obj; ch = ch)
    chs  = intersect(obj.locs[!, :label], labels(obj)[ch])
    locs = Base.filter(:label => in(chs), obj.locs)
    ch   = collect(1:DataFrames.nrow(locs))

    sch_resolved = if sch == ""
        Int64[]
    else
        sch_idx = get_channel(obj; ch = sch)
        sch_chs = intersect(locs[!, :label], labels(obj)[sch_idx])
        _find_bylabel(locs, sch_chs)
    end

    return plot_locs3d(
        locs;
        ch          = ch,
        sch         = sch_resolved,
        ch_labels   = ch_labels,
        head_labels = head_labels,
        mono        = mono,
        cart        = cart,
        cam         = cam,
        mesh_type   = mesh_type,
        mesh_alpha  = mesh_alpha,
        gui         = gui,
    )
end
