export iplot_locs

"""
    mplot_locs(obj; <keyword arguments>)

Preview of channel locations.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `selected::Union{String, Vector{String}, Regex}`: which channels should be highlighted
- `ch_labels::Bool=true`: plot channel labels
- `head_labels::Bool=false`: plot head labels
- `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use polar coordinates for XY plane and spherical coordinates for XZ and YZ planes
- `cam::Tuple{Real, Real}=(20, 45)`: camera position - (XY plane angle, XZ plane angle)
- `mesh_type::Symbol=:disabled`: type of mesh to plot (`:disabled`, `:brain` or `:head`)
- `mesh_alpha::Float64=0.95`: mesh opacity, from 1 (no opacity) to 0 (complete opacity)

# Returns

- `Nothing`
"""
function iplot_locs(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, selected::Union{String, Vector{String}, Regex}="", ch_labels::Bool=true, head_labels::Bool=false, cart::Bool=false, mono::Bool=false, mesh_type::Symbol=:disabled, mesh_alpha::Float64=0.95, cam::Tuple{Real, Real}=(20, 45))::Nothing

    @assert datatype(obj) in ["eeg"] "Currently iplot_locs() works for EEG objects only."
    _check_var(mesh_type, [:disabled, :brain, :head], "mesh_type")
    _in(mesh_alpha, (0.0, 1.0), "mesh_alpha")

    ch = get_channel(obj, ch=ch)
    ch_n = length(ch)
    chs = intersect(obj.locs[!, :label], labels(obj)[ch])
    locs = Base.filter(:label => in(chs), obj.locs)
    ch = collect(1:DataFrames.nrow(locs))

    if selected == ""
        selected = 0
    else
        selected = get_channel(obj, ch=selected)
        selected = intersect(locs[!, :label], labels(obj)[selected])
        selected = _find_bylabel(locs, selected)
    end

    msh = nothing

    if mesh_type !== :disabled
        if mesh_type === :brain
            msh = FileIO.load(joinpath(res_path, "mesh/brain_hires.stl"))
        else
            msh = FileIO.load(joinpath(res_path, "mesh/head.stl"))
            mesh_alpha = 1.0
        end
        # scale mesh
        if mesh_type === :brain
            msh.position ./= _mesh_normalize_xyz(msh)
            msh.position .*= 0.95
        else
            msh.position ./= _mesh_normalize_xy(msh)
            msh.position .*= 1.1
        end
    end

    pal = mono ? :grays : :darktest

    if !cart
        loc_x = zeros(DataFrames.nrow(locs))
        loc_y = zeros(DataFrames.nrow(locs))
        loc_z = zeros(DataFrames.nrow(locs))
        for idx in 1:DataFrames.nrow(locs)
            loc_x[idx], loc_y[idx], loc_z[idx] = sph2cart(locs[idx, :loc_radius_sph], locs[idx, :loc_theta_sph], locs[idx, :loc_phi_sph])
        end
    else
        loc_x = locs[!, :loc_x]
        loc_y = locs[!, :loc_y]
        loc_z = locs[!, :loc_z]
    end

    if maximum(locs[:, :loc_x]) <= 1.2 && maximum(locs[:, :loc_y]) <= 1.2 && maximum(locs[:, :loc_z]) <= 1.5
        x_lim = (-1.5, 1.5)
        y_lim = (-1.5, 1.5)
        z_lim = (-1.5, 1.5)
    else
        x_lim = (-2.0, 2.0)
        y_lim = (-2.0, 2.0)
        z_lim = (-0.5, 1.0)
    end

    marker_size = length(ch) > 64 ? 8 : 16
    font_size = 15

    p = Figure(size=(850, 900))

    but_r = Button(p[2, 1], label="Right [r]", tellwidth=false)
    but_l = Button(p[3, 1], label="Left [l]", tellwidth=false)
    but_f = Button(p[2, 2], label="Front [f]", tellwidth=false)
    but_b = Button(p[3, 2], label="Back [b]", tellwidth=false)
    but_t = Button(p[2, 3], label="Top [t]", tellwidth=false)
    but_reset = Button(p[3, 3], label="Reset [Home]", tellwidth=false)
    but_png = Button(p[2:3, 4], label="Save as PNG", tellwidth=false)

    ax = Axis3(p[1, 1:4],
               xlabel="X",
               ylabel="Y",
               zlabel="Z",
               limits=(x_lim, y_lim, z_lim),
               title="",
               aspect=(1, 1, 1),
               xticks=[-1, 0, 1],
               yticks=[-1, 0, 1],
               zticks=[-1, 0, 1],
               elevation=deg2rad(cam[1]),
               azimuth=deg2rad(cam[2]))

#    on(sl_xy.value) do val
#        ax.elevation[] = val
#    end
#    on(sl_yz.value) do val
#        ax.azimuth[] = val
#    end

    if mesh_type !== :disabled
        GLMakie.mesh!(msh,
                      alpha=mesh_alpha,
                      color=:gray)
    end

    ch_n = length(ch)
    cmap = GLMakie.resample_cmap(pal, ch_n)
    ch = setdiff(ch, selected)

    for idx in 1:ch_n
        if idx in selected
            if mono
                GLMakie.scatter!(loc_x[idx],
                                 loc_y[idx],
                                 loc_z[idx],
                                 markersize=marker_size,
                                 color=:gray,
                                 strokewidth=1,
                                 strokecolor=:black)

            else
                GLMakie.scatter!(loc_x[idx],
                                 loc_y[idx],
                                 loc_z[idx],
                                 markersize=marker_size,
                                 color=cmap[idx],
                                 colormap=pal,
                                 colorrange=1:ch_n,
                                 strokewidth=1,
                                 strokecolor=:black)
            end
        else
            GLMakie.scatter!(loc_x[idx],
                             loc_y[idx],
                             loc_z[idx],
                             markersize=marker_size,
                             color=:gray,
                             strokewidth=1,
                             strokecolor=:black)
        end
    end

    if ch_labels
        GLMakie.text!(loc_x[ch] * 1.15,
                      loc_y[ch] * 1.15,
                      loc_z[ch] * 1.15,
                      text=locs[ch, :label],
                      fontsize=font_size,
                      align=(:center, :center))
        if selected != 0
            GLMakie.text!(loc_x[selected] * 1.15,
                          loc_y[selected] * 1.15,
                          loc_z[selected] * 1.15,
                          text=locs[selected, :label],
                          fontsize=font_size,
                          align=(:center, :center))
        end
    end

    if head_labels
        fid_names = ["NAS", "IN", "LPA", "RPA"]
        for idx in 1:length(NeuroAnalyzer.fiducial_points)
            GLMakie.text!(NeuroAnalyzer.fiducial_points[idx][1],
                          NeuroAnalyzer.fiducial_points[idx][2],
                          NeuroAnalyzer.fiducial_points[idx][3],
                          text=fid_names[idx],
                          fontsize=font_size,
                          align=(:center, :center))
        end
    end

    screen = display(GLMakie.Screen(), p)

    on(events(p).keyboardbutton) do event
        new_pov = (0, 0)
        event.key == Keyboard.home && (new_pov = (20, 45))
        event.key == Keyboard.r && (new_pov = (10, 10))
        event.key == Keyboard.l && (new_pov = (10, 190))
        event.key == Keyboard.f && (new_pov = (10, 100))
        event.key == Keyboard.b && (new_pov = (10, 280))
        event.key == Keyboard.t && (new_pov = (90, 270))
        if new_pov != (0, 0)
            ax.elevation[] = deg2rad(new_pov[1])
            ax.azimuth[] = deg2rad(new_pov[2])
        end
    end

    on(but_png.clicks) do _
        save_dialog("Pick an image file", nothing, ["*.png"]) do file_name
            if file_name != ""
                GLMakie.save(file_name, p.scene)
                @info("Image saved as $file_name")
            end
        end
    end

    on(but_reset.clicks) do _
        ax.elevation[] = deg2rad(cam[1])
        ax.azimuth[] = deg2rad(cam[2])
    end

    on(but_r.clicks) do _
        ax.elevation[] = deg2rad(10)
        ax.azimuth[] = deg2rad(10)
    end

    on(but_l.clicks) do _
        ax.elevation[] = deg2rad(10)
        ax.azimuth[] = deg2rad(190)
    end

    on(but_t.clicks) do _
        ax.elevation[] = deg2rad(90)
        ax.azimuth[] = deg2rad(270)
    end

    on(but_f.clicks) do _
        ax.elevation[] = deg2rad(10)
        ax.azimuth[] = deg2rad(100)
    end

    on(but_b.clicks) do _
        ax.elevation[] = deg2rad(10)
        ax.azimuth[] = deg2rad(280)
    end

#    on(sl_xy.value) do val
#        ax.elevation[] = val
#    end
#    on(sl_yz.value) do val
#        ax.azimuth[] = val
#    end

    wait(display(p))

    return nothing

end
