export plot_dipole3d

"""
    plot_dipole3d(d; <keyword arguments>)

Plot a 3D dipole inside a schematic brain model.

# Arguments

- `d::NeuroAnalyzer.DIPOLE`: dipole object with `pos` (position) and `q` (moment) fields
- `t::Int64=1`: time point to plot (column of `q`; row of `pos` for a moving dipole)
- `project::Bool=true`: if `true`, plot lines projecting the dipole onto the x, y, and z axes

# Returns

- `GLMakie.Figure`: the plotted figure

# Notes

- Brain volume is within -1.0 to +1.0 for all axes (x, y, z).
- The dipole position is marked with a red dot.
- If `project=true`, dashed lines show the dipole's projection onto the coordinate planes.
- Position and moment must be normalized to the unit head; for a fitted dipole (in `locs` units), use `(pos - sphere_center) / sphere_radius`.
"""
function plot_dipole3d(
    d::NeuroAnalyzer.DIPOLE;
    t::Int64 = 1,
    project::Bool = true,
)::GLMakie.Figure
    _wip()

    # validate
    pos, mag = _dipole_at(d, t)
    all(-1.0 .≤ pos .≤ 1.0) ||
        throw(ArgumentError("Position must be within [-1.0, 1.0]."))
    all(-1.0 .≤ mag .≤ 1.0) ||
        throw(ArgumentError("Magnitude must be within [-1.0, 1.0]."))

    # define texture file paths (adjust as needed)
    brain_top_texture_path = joinpath(NeuroAnalyzer.res_path, "brain_top.png")
    brain_side_texture_path = joinpath(NeuroAnalyzer.res_path, "brain_side.png")
    brain_front_texture_path = joinpath(NeuroAnalyzer.res_path, "brain_front.png")

    # load textures (fallback to simple colors if files are missing)
    brain_top_texture = nothing
    brain_side_texture = nothing
    brain_front_texture = nothing
    try
        brain_top_texture = FileIO.load(brain_top_texture_path)
        brain_side_texture = FileIO.load(brain_side_texture_path)
        brain_front_texture = FileIO.load(brain_front_texture_path)
        brain_side_texture = brain_side_texture[:, end:-1:1]
        brain_front_texture = rotr90(brain_front_texture)
    catch e
        @warn "Failed to load brain textures: $e. Using fallback colors."
        brain_top_texture = rand(RGB, 100, 100)
        brain_side_texture = rand(RGB, 100, 100)
        brain_front_texture = rand(RGB, 100, 100)
    end

    # prepare meshes for brain surfaces (top, side, front)
    brain_top_vertices =
        Point3f[(-1.2, -1.2, -0.1), (1.2, -1.2, -0.1), (1.2, 1.2, -0.1), (-1.2, 1.2, -0.1)]
    brain_top_uvs = Vec2f[(0, 0), (1, 0), (1, 1), (0, 1)]
    brain_top_faces = [TriangleFace(1, 2, 3), TriangleFace(1, 3, 4)]
    brain_top_mesh =
        GeometryBasics.Mesh(brain_top_vertices, brain_top_faces; uv = brain_top_uvs)

    brain_side_vertices =
        Point3f[(-1.2, -1.2, -0.1), (-1.2, 1.2, -0.1), (-1.2, 1.2, 1.2), (-1.2, -1.2, 1.2)]
    brain_side_uvs = Vec2f[(0, 0), (1, 0), (1, 1), (0, 1)]
    brain_side_faces = [TriangleFace(1, 2, 3), TriangleFace(1, 3, 4)]
    brain_side_mesh =
        GeometryBasics.Mesh(brain_side_vertices, brain_side_faces; uv = brain_side_uvs)

    brain_front_vertices =
        Point3f[(-1.2, 1.2, -0.1), (-1.2, 1.2, 1.2), (1.2, 1.2, 1.2), (1.2, 1.2, -0.1)]
    brain_front_uvs = Vec2f[(0, 0), (1, 0), (1, 1), (0, 1)]
    brain_front_faces = [TriangleFace(1, 2, 3), TriangleFace(1, 3, 4)]
    brain_front_mesh =
        GeometryBasics.Mesh(brain_front_vertices, brain_front_faces; uv = brain_front_uvs)

    # extract position and magnitude
    x, y, z = pos
    mx, my, mz = mag

    # prepare plot
    GLMakie.activate!(; title = "plot_dipole_3d()")
    plot_size = (800, 800)
    fig = Figure(;
        backgroundcolor = :black,
        size = plot_size,
    )
    ax = Axis3(fig[1, 1])
    hidedecorations!(ax)

    # draw brain surfaces with textures
    GLMakie.mesh!(ax, brain_top_mesh; color = brain_top_texture, shading = NoShading)
    GLMakie.mesh!(ax, brain_side_mesh; color = brain_side_texture, shading = NoShading)
    GLMakie.mesh!(ax, brain_front_mesh; color = brain_front_texture, shading = NoShading)

    # draw dipole position
    GLMakie.scatter!(ax, x, y, z; markersize = 20, color = :red)

    # draw projection lines if requested
    if project
        # project onto top plane (z=-0.1)
        GLMakie.lines!(ax, [x, x], [y, y], [z, -0.1]; linestyle = :dash, color = :red)
        # project onto side plane (x=-1.2)
        GLMakie.lines!(ax, [x, -1.2], [y, y], [z, z]; linestyle = :dash, color = :red)
        # project onto front plane (y=1.2)
        GLMakie.lines!(ax, [x, x], [y, 1.2], [z, z]; linestyle = :dash, color = :red)
    end

    return fig
end
