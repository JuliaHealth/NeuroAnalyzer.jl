export plot_dipole3d

"""
    plot_dipole3d(d; <keyword arguments>)

Plot dipole in 3D.

# Arguments

- `d::NeuroAnalyzer.DIPOLE`
- `project::Bool=true`: plot lines projected onto X, Y and Z axes

# Returns

- `p::Plots.Plot{Plots.GRBackend}`

# Notes

Brain volume is within -1.0 to +1.0 (X-, Y- and Z-axis)
"""
function plot_dipole3d(d::NeuroAnalyzer.DIPOLE; project::Bool=true)

    _wip()

    # load textures
    brain_top_texture = FileIO.load(joinpath(NeuroAnalyzer.res_path, "brain_top.png"))
    brain_side_texture = FileIO.load(joinpath(NeuroAnalyzer.res_path, "brain_side.png"))
    brain_front_texture = FileIO.load(joinpath(NeuroAnalyzer.res_path, "brain_front.png"))

    # get dipole position
    x = d.pos[1]
    y = d.pos[2]
    z = d.pos[3]

    # get dipole magnitude
    mx = d.mag[1]
    my = d.mag[2]
    mz = d.mag[3]
    m = round(Int64, sqrt(mx^2 + my^2 + mz^2) * 5)

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
    lxy = (oxy[1] + x * 10, oxy[2] - y * 10)
    lyz = (oyz[1] + y * 10, oyz[2] - z * 10)
    lxz = (oxz[1] + x * 10, oxz[2] - z * 10)

    # prepare figure
    pxy = Plots.plot(brain_top_texture, border=:none, framestyle=:none, aspect_ratio=1, margins=0Plots.px)
    pyz = Plots.plot(brain_side_texture, border=:none, framestyle=:none, aspect_ratio=1, margins=0Plots.px)
    pxz = Plots.plot(brain_front_texture, border=:none, framestyle=:none, aspect_ratio=1, margins=0Plots.px)

    # draw dipole
    pxy = Plots.scatter!(pxy, lxy, label="", c=:red, msc=:red, msa=0, ms=m)
    pyz = Plots.scatter!(pyz, lyz, label="", c=:red, msc=:red, msa=0, ms=m)
    pxz = Plots.scatter!(pxz, lxz, label="", c=:red, msc=:red, msa=0, ms=m)

    p3d = Plots.scatter3d((x, y, z),
                          label="",
                          c=:red,
                          msc=:red,
                          msa=0,
                          ms=m,
                          xlims=(-1.5, 1.5),
                          ylims=(-1.5, 1.5),
                          zlims=(-1.5, 1.5),
                          xticks=[-1, 0, 1],
                          yticks=[-1, 0, 1],
                          zticks=[-1, 0, 1],
                          xlabel="horizontal",
                          ylabel="lateral",
                          zlabel="coronal",
                          legend=false,
                          foreground_color=:white,
                          background_color=:black,
                          foreground_color_text=:white)
    xaxis = Plots.get_axis(Plots.get_subplot(p3d,1),:x)
    yaxis = Plots.get_axis(Plots.get_subplot(p3d,1),:y)
    zaxis = Plots.get_axis(Plots.get_subplot(p3d,1),:z)
    xaxis[:gridalpha] = 0.6
    yaxis[:gridalpha] = 0.6
    zaxis[:gridalpha] = 0.6
    xaxis[:foreground_color_grid] = colorant"white"
    xaxis[:foreground_color_guide] = colorant"white"
    yaxis[:foreground_color_grid] = colorant"white"
    yaxis[:foreground_color_guide] = colorant"white"
    zaxis[:foreground_color_grid] = colorant"white"
    zaxis[:foreground_color_guide] = colorant"white"
    if project
        p3d = Plots.plot3d!(p3d, [x, -1.1], [y, y], [z, z], color=:red, lw=1, ls=:dash)
        p3d = Plots.plot3d!(p3d, [x, x], [y, 1.1], [z, z], color=:red, lw=1, ls=:dash)
        p3d = Plots.plot3d!(p3d, [x, x], [y, y], [z, -1.1], color=:red, lw=1, ls=:dash)
    end
    p = Plots.plot(pxy, pxz, pyz, p3d, layout=(2, 2), size=(800, 600), margins=0Plots.px, background_color=:black, foreground_color_grid=:white)

    return p

end
