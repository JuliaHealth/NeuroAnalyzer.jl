export plot_dipole3d

"""
    plot_dipole3d(d; <keyword arguments>)

Plot dipole in 3D.

# Arguments

  - `d::NeuroAnalyzer.DIPOLE`
  - `project::Bool=true`: plot lines projected onto X, Y and Z axes

# Returns

  - `p::GLMakie.Figure`

# Notes

Brain volume is within -1.0 to +1.0 (X-, Y- and Z-axis)
"""
function plot_dipole3d(d::NeuroAnalyzer.DIPOLE; project::Bool = true)

    _wip()

    # load textures
    brain_top_texture = FileIO.load(joinpath(res_path, "brain_top.png"))
    brain_side_texture = FileIO.load(joinpath(res_path, "brain_side.png"))
    brain_front_texture = FileIO.load(joinpath(res_path, "brain_front.png"))
    brain_side_texture = brain_side_texture[:, end:-1:1]
    brain_front_texture = rotr90(brain_front_texture)

    # prepare meshes
    brain_top_vertices = Point3f[(-1.2, -1.2, -0.1), (1.2, -1.2, -0.1), (1.2, 1.2, -0.1), (-1.2, 1.2, -0.1)]
    brain_top_uvs = Vec2f[(0, 0), (1, 0), (1, 1), (0, 1)]
    brain_top_faces = [TriangleFace(1, 2, 3), TriangleFace(1, 3, 4)]
    brain_top_mesh = GeometryBasics.Mesh(brain_top_vertices, brain_top_faces, uv = brain_top_uvs)

    brain_side_vertices = Point3f[(-1.2, -1.2, -0.1), (-1.2, 1.2, -0.1), (-1.2, 1.2, 1.2), (-1.2, -1.2, 1.2)]
    brain_side_uvs = Vec2f[(0, 0), (1, 0), (1, 1), (0, 1)]
    brain_side_faces = [TriangleFace(1, 2, 3), TriangleFace(1, 3, 4)]
    brain_side_mesh = GeometryBasics.Mesh(brain_side_vertices, brain_side_faces, uv = brain_side_uvs)

    brain_front_vertices = Point3f[(-1.2, 1.2, -0.1), (-1.2, 1.2, 1.2), (1.2, 1.2, 1.2), (1.2, 1.2, -0.1)]
    brain_front_uvs = Vec2f[(0, 0), (1, 0), (1, 1), (0, 1)]
    brain_front_faces = [TriangleFace(1, 2, 3), TriangleFace(1, 3, 4)]
    brain_front_mesh = GeometryBasics.Mesh(brain_front_vertices, brain_front_faces, uv = brain_front_uvs)

    # get dipole position
    x = d.pos[1]
    y = d.pos[2]
    z = d.pos[3]
    pos = (x, y, z)

    # get dipole magnitude
    mx = d.mag[1]
    my = d.mag[2]
    mz = d.mag[3]

    # prepare plot
    GLMakie.activate!(title = "plot_dipole_3d()")
    plot_size = (800, 800)
    p = Figure(
        backgroundcolor = :black,
        size = plot_size,
    )
    ax = Axis3(p[1, 1])
    hidedecorations!(ax)

    GLMakie.mesh!(ax, brain_top_mesh, color = brain_top_texture, shading = NoShading)
    GLMakie.mesh!(ax, brain_side_mesh, color = brain_side_texture, shading = NoShading)
    GLMakie.mesh!(ax, brain_front_mesh, color = brain_front_texture, shading = NoShading)

    # draw dipole
    GLMakie.scatter!(ax, x, y, z, markersize = 20, color = :red)

    if project == true
        # project at top-plane
        GLMakie.lines!(ax, [x, x], [y, y], [z, 0], linestyle = :dash, color = :red)
        # project at side-axis
        GLMakie.lines!(ax, [x, -1.2], [y, y], [z, z], linestyle = :dash, color = :red)
        # project at front-axis
        GLMakie.lines!(ax, [x, x], [y, 1.2], [z, z], linestyle = :dash, color = :red)
    end

    return p

end
