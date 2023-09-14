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

    _wip()

#=
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
    brain_top_texture = FileIO.load(joinpath(res_path, "brain_t.png"))
    brain_side_texture = FileIO.load(joinpath(res_path, "brain_s.png"))
    brain_front_texture = FileIO.load(joinpath(res_path, "brain_f.png"))

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
=#
  
    # load textures
    brain_top_texture = FileIO.load(joinpath(NeuroAnalyzer.res_path, "brain_t512.png"))
    brain_side_texture = FileIO.load(joinpath(NeuroAnalyzer.res_path, "brain_s512.png"))
    brain_front_texture = FileIO.load(joinpath(NeuroAnalyzer.res_path, "brain_f512.png"))

    # get dipole position
    x = d.pos[1]
    y = d.pos[2]
    z = d.pos[3]
    p = (x, y, z) .* 43

    # get dipole magnitude
    mx = d.mag[1]
    my = d.mag[2]
    mz = d.mag[3]
    m = (mx, my, mz) .* 10

    # origins
    # XY: x: 14..100
    # XY: y: 14..100
    # YZ: y: 9..93
    # YZ: z: 8..58
    # XZ: x: 14..98
    # XZ: z: 8..58
    oxy = reverse(size(brain_top_texture) .÷ 2 .+ 1)
    oxy = (oxy[1], oxy[2] + 3)
    oyz = reverse(size(brain_side_texture) .÷ 2 .+ 1)
    oyz = (oyz[1] - 2, oyz[2] + 5)
    oxz = reverse(size(brain_front_texture) .÷ 2 .+ 1)
    oxz = (oxz[1], oxz[2] + 5)

    # convert to image size; y- and z-coordinate has to be flipped
    # origin location
    lxy = (oxy[1] + p[1], oxy[2] - p[2])
    lyz = (oyz[1] + p[2], oyz[2] - p[3])
    lxz = (oxz[1] + p[1], oxz[2] - p[3])

    # magnitude
    mxy = (lxy[1] + m[1], lxy[2] - m[2])
    myz = (lyz[1] + m[2], lyz[2] - m[3])
    mxz = (lxz[1] + m[1], lxz[2] - m[3])

    # mx_xy = lxy[1] + (mx * oxy[1] * 0.25)
    # my_xy = lxy[2] + -(my * oxy[2] * 0.25)

    # mx_xz = lxz[1] + (mx * oxz[1] * 0.25)
    # mz_xz = lxz[2] + -(mz * oxz[2] * 0.25)

    # my_yz = lyz[1] + (my * oyz[1] * 0.25)
    # mz_yz = lyz[2] + -(mz * oyz[2] * 0.25)

    # prepare figure
    pxy = Plots.plot(brain_top_texture, border=:none, framestyle=:none, aspect_ratio=1, margins=0Plots.px)
    pyz = Plots.plot(brain_side_texture, border=:none, framestyle=:none, aspect_ratio=1, margins=0Plots.px)
    pxz = Plots.plot(rotl90(brain_front_texture), border=:none, framestyle=:none, aspect_ratio=1, margins=0Plots.px)

    # draw dipole
    pxy = Plots.scatter!(pxy, lxy, label="", c=:red, msc=:red, msa=0, ms=5)
    pyz = Plots.scatter!(pyz, lyz, label="", c=:red, msc=:red, msa=0, ms=5)
    pxz = Plots.scatter!(pxz, lxz, label="", c=:red, msc=:red, msa=0, ms=5)

    pxy = Plots.plot!(pxy, lw=3, [lxy[1], mxy[1]], [lxy[2], mxy[2]], label="", c=:red)
    pyz = Plots.plot!(pyz, lw=3, [lyz[1], myz[1]], [lyz[2], myz[2]], label="", c=:red)
    pxz = Plots.plot!(pxz, lw=3, [lxz[1], mxz[1]], [lxz[2], mxz[2]], label="", c=:red)

    p3d = Plots.scatter3d((x, y, z), label="", c=:red, msc=:red, msa=0, ms=5, xlims=(-1.5, 1.5), ylims=(-1.5, 1.5), zlims=(-1.5, 1.5), xticks=[-1, 0, 1], yticks=[-1, 0, 1], zticks=[-1, 0, 1], legend=false, foreground_color=:white, background_color=:black, foreground_color_grid=:white)
        xaxis = Plots.get_axis(Plots.get_subplot(p3d,1),:x)
        yaxis = Plots.get_axis(Plots.get_subplot(p3d,1),:y)
        zaxis = Plots.get_axis(Plots.get_subplot(p3d,1),:z)
        xaxis[:gridalpha] = 0.6
        yaxis[:gridalpha] = 0.6
        zaxis[:gridalpha] = 0.6
        xaxis[:foreground_color_grid] = colorant"white"
        yaxis[:foreground_color_grid] = colorant"white"
        zaxis[:foreground_color_grid] = colorant"white"
    p3d = Plots.plot3d!(p3d, [x, x + mx], [y, y + my], [z, z + mz], color=:red, lw=3)
    p = Plots.plot(pxy, pxz, pyz, layout=(1, 3), margins=0Plots.px)
    p = Plots.plot(p, p3d, layout=(2, 1), size=(800, 600), margins=0Plots.px)

    return p

end
