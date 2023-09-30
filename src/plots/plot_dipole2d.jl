export plot_dipole2d

"""
    plot_dipole2d(d; <keyword arguments>)

Plot dipole in 2D.

# Arguments

- `d::NeuroAnalyzer.DIPOLE`

# Returns

- `p::Plots.Plot{Plots.GRBackend}`

# Notes

Brain volume is within -1.0 to +1.0 (X-, Y- and Z-axis)
"""
function plot_dipole2d(d::NeuroAnalyzer.DIPOLE)

    _wip()

    # load textures
    head_top_texture = FileIO.load(joinpath(res_path, "head_t.png"))
    head_side_texture = FileIO.load(joinpath(res_path, "head_s.png"))
    head_front_texture = FileIO.load(joinpath(res_path, "head_f.png"))

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
    # 
    oxy = reverse(size(head_top_texture) .÷ 2 .+ 1)
    oxy = (oxy[1], oxy[2] + 3)
    oyz = reverse(size(head_side_texture) .÷ 2 .+ 1)
    oyz = (oyz[1] - 2, oyz[2] + 5)
    oxz = reverse(size(head_front_texture) .÷ 2 .+ 1)
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
    pxy = Plots.plot(head_top_texture, border=:none, framestyle=:none, aspect_ratio=1)
    pyz = Plots.plot(head_side_texture, border=:none, framestyle=:none, aspect_ratio=1);
    pxz = Plots.plot(head_front_texture, border=:none, framestyle=:none, aspect_ratio=1);

    # draw dipole
    pxy = Plots.plot!(pxy, lt=:scatter, lxy, label="", c=:red, msc=:red, msa=0, ms=3)
    pyz = Plots.plot!(pyz, lt=:scatter, lyz, label="", c=:red, msc=:red, msa=0, ms=3)
    pxz = Plots.plot!(pxz, lt=:scatter, lxz, label="", c=:red, msc=:red, msa=0, ms=3)

    pxy = Plots.plot!(pxy, lw=3, [lxy[1], mxy[1]], [lxy[2], mxy[2]], label="", c=:red)
    pyz = Plots.plot!(pyz, lw=3, [lyz[1], myz[1]], [lyz[2], myz[2]], label="", c=:red)
    pxz = Plots.plot!(pxz, lw=3, [lxz[1], mxz[1]], [lxz[2], mxz[2]], label="", c=:red)

    # pxy = Plots.plot!(pxy, lw=4, [lxy[1], mx_xy], [lxy[2], my_xy], label="", c=:red)
    # pyz = Plots.plot!(pyz, lw=4, [lyz[1], my_yz], [lyz[2], mz_yz], label="", c=:red)
    # pxz = Plots.plot!(pxz, lw=4, [lxz[1], mx_xz], [lxz[2], mz_xz], label="", c=:red)

    p = Plots.plot(pxy, pxz, pyz, layout=(1, 3), size=(600, 200), left_margin=-30Plots.px, top_margin=-20Plots.px, bottom_margin=-15Plots.px)

    return p

end
