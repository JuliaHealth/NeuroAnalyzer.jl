export plot_dipole2d

"""
    plot_dipole2d(d; <keyword arguments>)

Plot dipole in 2D.

# Arguments

- `d::NeuroAnalyzer.DIPOLE`

# Returns

- `p::GLMakie.Figure`

# Notes

Brain volume is within -1.0 to +1.0 (X-, Y- and Z-axis)
"""
function plot_dipole2d(d::NeuroAnalyzer.DIPOLE)

    # load textures
    head_top_texture = FileIO.load("resources/head_t.png")
    head_side_texture = FileIO.load("resources/head_s.png")
    head_front_texture = FileIO.load("resources/head_f.png")

    # get dipole position
    x = d.loc[1]
    y = d.loc[2]
    z = d.loc[3]

    # get dipole magnitude
    mx = d.mag[1]
    my = d.mag[2]
    mz = d.mag[3]

    # origins
    oxy = reverse(size(head_top_texture) .÷ 2 .+ 1)
    oyz = reverse(size(head_side_texture) .÷ 2 .+ 1)
    oxz = reverse(size(head_front_texture) .÷ 2 .+ 1)

    # oyz = (oyz[1], oyz[2] - (oyz[2] * 0.25))
    # oxz = (oxz[1], oxz[2] - (oxz[2] * 0.25))

    x = 0
    y = 0
    z = -0.2

    mx = 0
    my = 0
    mz = -1

    # convert to image size; y- and z-coordinate has to be flipped
    # origin location
    lxy = (oxy[1] + (x * oxy[1] * 0.25), oxy[2] + -(y * oxy[2] * 0.25))
    lyz = (oyz[1] + (y * oyz[1] * 0.25), oyz[2] + -(z * oyz[2] * 0.85))
    lxz = (oxz[1] + (x * oxz[1] * 0.25), oxz[2] + -(z * oxz[2] * 0.85))

    # magnitude
    mx_xy = lxy[1] + (mx * oxy[1] * 0.25)
    my_xy = lxy[2] + -(my * oxy[2] * 0.25)

    mx_xz = lxz[1] + (mx * oxz[1] * 0.25)
    mz_xz = lxz[2] + -(mz * oxz[2] * 0.25)

    my_yz = lyz[1] + (my * oyz[1] * 0.25)
    mz_yz = lyz[2] + -(mz * oyz[2] * 0.25)

    # prepare figure
    pxy = Plots.plot(head_top_texture, border=:none, framestyle=:none, aspect_ratio=1)
    pyz = Plots.plot(head_side_texture, border=:none, framestyle=:none, aspect_ratio=1);
    pxz = Plots.plot(head_front_texture, border=:none, framestyle=:none, aspect_ratio=1);

    # draw dipole
    pxy = Plots.plot!(pxy, lt=:scatter, lxy, label="", c=:red, msc=:red, msa=0)
    pyz = Plots.plot!(pyz, lt=:scatter, lyz, label="", c=:red, msc=:red, msa=0)
    pxz = Plots.plot!(pxz, lt=:scatter, lxz, label="", c=:red, msc=:red, msa=0)

    pxy = Plots.plot!(pxy, lw=4, [lxy[1], mx_xy], [lxy[2], my_xy], label="", c=:red)
    pyz = Plots.plot!(pyz, lw=4, [lyz[1], my_yz], [lyz[2], mz_yz], label="", c=:red)
    pxz = Plots.plot!(pxz, lw=4, [lxz[1], mx_xz], [lxz[2], mz_xz], label="", c=:red)

    p = Plots.plot(pxy, pxz, pyz, layout=(1, 3), size=(600, 200), left_margin=-30Plots.px)

    return p

end
