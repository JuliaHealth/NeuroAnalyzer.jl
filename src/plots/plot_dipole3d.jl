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
function plot_dipole3d(d::NeuroAnalyzer.DIPOLE; project::Bool=true)

    # prepare meshes
    brain_top = Point3f[[-1.2, -1.2, 0.0], [1.2, -1.2, 0.0], [1.2, 1.2, 0.0], [-1.2, 1.2, 0.0]]
    brain_top_uvs = Vec2f[(0, 0), (1, 0), (1, 1), (0, 1)]
    brain_top_fs = GLTriangleFace[(1, 2, 3), (1, 3, 4)]
    brain_top_mesh = GeometryBasics.Mesh(GeometryBasics.meta(brain_top, uv = brain_top_uvs, normals = normals(brain_top, brain_top_fs)), brain_top_fs)

    brain_side = Point3f[[-1.2, -1.2, 0],[-1.2, 1.2, 0], [-1.2, 1.2, 1.2], [-1.2, -1.2, 1.2]]
    brain_side_uvs = Vec2f[(0, 0), (1, 0), (1, 1), (0, 1)]
    brain_side_fs = GLTriangleFace[(1, 2, 3), (1, 3, 4)]
    brain_side_mesh = GeometryBasics.Mesh(GeometryBasics.meta(brain_side, uv = brain_side_uvs, normals = normals(brain_side, brain_side_fs)), brain_side_fs)

    brain_front = Point3f[[-1.2, 1.2, 0],[-1.2, 1.2, 1.2], [1.2, 1.2, 1.2], [1.2, 1.2, 0]]
    brain_front_uvs = Vec2f[(0, 0), (1, 0), (1, 1), (0, 1)]
    brain_front_fs = GLTriangleFace[(1, 2, 3), (1, 3, 4)]
    brain_front_mesh = GeometryBasics.Mesh(GeometryBasics.meta(brain_front, uv = brain_front_uvs, normals = normals(brain_front, brain_front_fs)), brain_front_fs)

    # load textures
    brain_top_texture = FileIO.load("resources/brain_top.png")
    brain_side_texture = FileIO.load("resources/brain_side.png")
    brain_front_texture = FileIO.load("resources/brain_front.png")

    # get dipole position
    x = d.pos[1]
    y = d.pos[2]
    z = d.pos[3]

    # get dipole magnitude
    mx = d.mag[1]
    my = d.mag[2]
    mz = d.mag[3]

    # prepare figure
    p = Figure(backgroundcolor=:black, resolution=(800, 800))
    ax = Axis3(p[1, 1])
    hidedecorations!(ax)
    mesh!(ax, brain_side_mesh, color=brain_side_texture)
    mesh!(ax, brain_top_mesh, color=brain_top_texture)
    mesh!(ax, brain_front_mesh, color=brain_front_texture)

    # draw dipole
    GLMakie.scatter!(ax, x, y, z, markersize=20, color=:red)

    if project == true
        # project at top-plane
        GLMakie.lines!(ax, [x, x], [y, y], [z, 0], linestyle=:dash, color=:red)
        # project at side-axis
        GLMakie.lines!(ax, [x, -1.2], [y, y], [z, z], linestyle=:dash, color=:red)
        # project at front-axis
        GLMakie.lines!(ax, [x, x], [y, 1.2], [z, z], linestyle=:dash, color=:red)
    end

    GLMakie.show(p)
    
    return p

end
