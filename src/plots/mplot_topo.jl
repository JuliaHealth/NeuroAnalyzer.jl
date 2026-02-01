export mplot_topo

"""
    mplot_topo(c; <keyword arguments>)

Plot topographical view.

# Arguments

- `s::Vector{<:Real}`: values to plot (one value per channel)
- `locs::DataFrame`: columns: channel, labels, loc_radius, loc_theta, loc_x, loc_y, loc_z, loc_radius_sph, loc_theta_sph, loc_phi_sph
- `ch::Union{Int64, Vector{Int64}}=1:DataFrames.nrow(locs)`: list of channels, default is all channels
- `cb::Bool=true`: plot color bar
- `cb_title::String="[A.U.]"`: color bar title
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `imethod::Symbol=:sh`: interpolation method:
    - `:sh`: Shepard
    - `:mq`: Multiquadratic
    - `:imq`: InverseMultiquadratic
    - `:tp`: ThinPlate
    - `:nn`: NearestNeighbour
    - `:ga`: Gaussian
- `nmethod::Symbol=:minmax`: method for normalization, see `normalize()`
- `plot_contours::Bools=true`: plot contours over topo plot
- `plot_electrodes::Bools=true`: plot electrodes over topo plot
- `ps::Symbol=:l`: plot size (`:l`: large (800×800 px), `:m`: medium (300×300 px), `:s`: small (100×100 px))
- `head::Bool=true`: draw head
- `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use polar coordinates for XY plane and spherical coordinates for XZ and YZ planes
- `threshold::Union{Nothing, Real}=nothing`: if set, use threshold to mark a region
- `threshold_type::Symbol=:neq`: rule for thresholding:
    - `:eq`: draw region is values are equal to threshold
    - `:neq`: draw region is values are not equal to threshold
    - `:geq`: draw region is values are ≥ to threshold
    - `:leq`: draw region is values are ≤ to threshold
    - `:g`: draw region is values are > to threshold
    - `:l`: draw region is values are < to threshold

# Returns

- `p::GLMakie.Figure`
"""
function mplot_topo(s::Vector{<:Real}; locs::DataFrame, ch::Union{Int64, Vector{Int64}}=1:DataFrames.nrow(locs), cb::Bool=true, cb_title::String="[A.U.]", title::String="default", mono::Bool=false, imethod::Symbol=:sh, nmethod::Symbol=:minmax, plot_contours::Bool=true, plot_electrodes::Bool=true, ps::Symbol=:l, head::Bool=true, cart::Bool=false, threshold::Union{Nothing, Real}=nothing, threshold_type::Symbol=:neq)::GLMakie.Figure

    pal = mono ? :grays : :bluesreds
    _check_var(imethod, [:sh, :mq, :imq, :tp, :nn, :ga], "imethod")
    _check_var(ps, [:l, :m, :s], "ps")

    ch_n = length(ch)
    locs = locs[ch, :]

    if cart == false
        loc_x = zeros(length(ch))
        loc_y = zeros(length(ch))
        for idx in eachindex(ch)
            loc_x[idx], loc_y[idx] = pol2cart(locs[!, :loc_radius][idx], locs[!, :loc_theta][idx])
        end
    else
        loc_x = locs[ch, :loc_x]
        loc_y = locs[ch, :loc_y]
    end

    loc_x = _s2v(loc_x)
    loc_y = _s2v(loc_y)

    if ps === :l
        plot_size = (800, 800)
        marker_size = length(ch) > 64 ? 8 : 16
        length(ch) > 64 && (ch_labels = false)
        iter = 512
        font_size = 20
        !occursin("\n", title) && title !== "" && (title *= "\n")
    elseif ps === :m
        plot_size = (300, 300)
        marker_size = length(ch) > 64 ? 2 : 4
        font_size = 8
        iter = 256
        title=""
        cb_title=""
        font_size = 10
    elseif ps === :s
        plot_size = (100, 100)
        marker_size = length(ch) > 64 ? 1 : 2
        iter = 128
        title=""
        cb = false
        plot_contours = false
        cb_title=""
        font_size = 5
    end

    s_interpolated, interpolated_x, interpolated_y = _interpolate2d(s, loc_x, loc_y, iter, imethod, nmethod)
    s_interpolated = s_interpolated'[:, end:-1:1]

    # we need to calculate thresholded region before removing peripherals
    if !isnothing(threshold)
        _, bm = seg_extract(s_interpolated, threshold=threshold, threshold_type=threshold_type)
        # _, bm = seg_extract(s, threshold=1.0, threshold_type=:geq)
        # reg = ones(size(s_interpolated)) .* minimum(s_interpolated)
        reg = zeros(size(s_interpolated))
        # reg = zeros(size(s))
        # reg[bm] .= maximum(s_interpolated)
        reg[bm] .= 1.0
        # reg = reg[:, end:-1:1]
    end

    head12 = false
    maximum(abs.(locs[:, :loc_x])) <= 1.2 && maximum(abs.(locs[:, :loc_y])) <= 1.2 && maximum(abs.(locs[:, :loc_z])) <= 1.5 && (head12 = true)

    if head12
        xl = (-1.2, 1.2)
        yl = (-1.2, 1.2)
        r = 1.2
    else
        xl = (-1.6, 1.6)
        yl = (-1.6, 1.6)
        r = 1.6
    end
    xt = (round.(linspace(xl[1], xl[2], 25), digits=1))
    yt = (round.(linspace(yl[1], yl[2], 25), digits=1))
    interpolated_x = round.(linspace(xl[1], xl[2], length(interpolated_x)), digits=2)
    interpolated_y = round.(linspace(yl[1], yl[2], length(interpolated_y)), digits=2)

    d = zeros(length(interpolated_x), length(interpolated_y))
    for idx1 in eachindex(interpolated_x)
        for idx2 in eachindex(interpolated_y)
            d[idx1, idx2] = distance((0, 0), (interpolated_x[idx1], interpolated_y[idx2]))
        end
    end
    s_interpolated[d .>= xl[2]] .= NaN
    !isnothing(threshold) && (reg[d .>= xl[2]] .= NaN)

    # prepare plot
    p = GLMakie.Figure(size=plot_size,
                       figure_padding=10)
    ax = GLMakie.Axis(p[1, 1],
                      aspect=1,
                      xlabel="",
                      ylabel="",
                      title=title,
                      xautolimitmargin=(0, 0),
                      yautolimitmargin=(0, 0),
                      backgroundcolor=:transparent,
                      titlesize=font_size)
    hidedecorations!(ax)
    hidespines!(ax)
    GLMakie.xlims!(ax, xl)
    GLMakie.ylims!(ax, yl)

    hm = GLMakie.heatmap!(ax,
                          interpolated_x,
                          interpolated_y,
                          s_interpolated,
                          colormap=pal)

    # draw contours
    if plot_contours && isnothing(threshold)
        GLMakie.contour!(ax,
                         interpolated_x,
                         interpolated_y,
                         s_interpolated,
                         linestyle=:dash,
                         levels=5,
                         linewidth=0.5,
                         color=:black)
    end

    # draw head
    if head
        ps === :l && (lw = 3)
        ps === :m && (lw = 2)
        ps === :s && (lw = 1)
        # nose
        GLMakie.lines!(ax, [-0.1, 0], [0.995, 1.1], linewidth=lw, color=:black)
        GLMakie.lines!(ax, [0, 0.1], [1.1, 0.995], linewidth=lw, color=:black)

        # ears
        # left
        GLMakie.lines!(ax, [-0.995, -1.03], [0.1, 0.15], linewidth=lw, color=:black)
        GLMakie.lines!(ax, [-1.03, -1.06], [0.15, 0.16], linewidth=lw, color=:black)
        GLMakie.lines!(ax, [-1.06, -1.1], [0.16, 0.14], linewidth=lw, color=:black)
        GLMakie.lines!(ax, [-1.1, -1.12], [0.14, 0.05], linewidth=lw, color=:black)
        GLMakie.lines!(ax, [-1.12, -1.10], [0.05, -0.1], linewidth=lw, color=:black)
        GLMakie.lines!(ax, [-1.10, -1.13], [-0.1, -0.3], linewidth=lw, color=:black)
        GLMakie.lines!(ax, [-1.13, -1.09], [-0.3, -0.37], linewidth=lw, color=:black)
        GLMakie.lines!(ax, [-1.09, -1.02], [-0.37, -0.39], linewidth=lw, color=:black)
        GLMakie.lines!(ax, [-1.02, -0.98], [-0.39, -0.33], linewidth=lw, color=:black)
        GLMakie.lines!(ax, [-0.98, -0.975], [-0.33, -0.22], linewidth=lw, color=:black)
        # right
        GLMakie.lines!(ax, [0.995, 1.03], [0.1, 0.15], linewidth=lw, color=:black)
        GLMakie.lines!(ax, [1.03, 1.06], [0.15, 0.16], linewidth=lw, color=:black)
        GLMakie.lines!(ax, [1.06, 1.1], [0.16, 0.14], linewidth=lw, color=:black)
        GLMakie.lines!(ax, [1.1, 1.12], [0.14, 0.05], linewidth=lw, color=:black)
        GLMakie.lines!(ax, [1.12, 1.10], [0.05, -0.1], linewidth=lw, color=:black)
        GLMakie.lines!(ax, [1.10, 1.13], [-0.1, -0.3], linewidth=lw, color=:black)
        GLMakie.lines!(ax, [1.13, 1.09], [-0.3, -0.37], linewidth=lw, color=:black)
        GLMakie.lines!(ax, [1.09, 1.02], [-0.37, -0.39], linewidth=lw, color=:black)
        GLMakie.lines!(ax, [1.02, 0.98], [-0.39, -0.33], linewidth=lw, color=:black)
        GLMakie.lines!(ax, [0.98, 0.975], [-0.33, -0.22], linewidth=lw, color=:black)

        # head
        GLMakie.arc!(ax,(0, 0), 1, 0, 2pi, linewidth=lw, color=:black)
    end

    # draw electrodes
    if plot_electrodes
        ps === :l && (sw = 2)
        ps === :m && (sw = 1)
        ps === :s && (sw = 0.0)

        for idx in 1:ch_n
            GLMakie.scatter!(loc_x[idx],
                             loc_y[idx],
                             markersize=marker_size,
                             color=:black)
        end
    end

    # draw thresholded region
    if !isnothing(threshold)
        GLMakie.contour!(interpolated_x,
                         interpolated_y,
                         reg,
                         levels=1,
                         linestyle=:dash,
                         color=:black,
                         linewidth=1)
    end

    # draw mask
    GLMakie.arc!(ax,
                 Point2f(0),
                 r,
                 -pi, pi,
                 linewidth=5,
                 color=:white)

    # draw colorbar
    if cb
        GLMakie.Colorbar(p[1, 2],
                         hm,
                         label=cb_title,
                         labelsize=font_size,
                         ticklabelsize=font_size,
                         width=ps === :l ? 25 : 10,
                         tellheight=true)
        rowsize!(p.layout, 1, ax.scene.viewport[].widths[2])
        colgap!(p.layout, 10)
    end

    return p

end

"""
    mplot_topo(obj; <keyword arguments>)

Topographical plot.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ep::Union{Int64, AbstractRange}=0`: epoch to display
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
- `title::String="default"`: plot title, default is Amplitude topographical plot [channels: 1:19, epoch: 1, time window: 0 ms:20 s]
- `mono::Bool=false`: use color or gray palette
- `cb::Bool=true`: plot color bar
- `cb_title::String="[A.U.]"`: color bar title
- `amethod::Symbol=:mean`: averaging method:
    - `:mean`
    - `:median`
- `imethod::Symbol=:sh`: interpolation method:
    - `:sh`: Shepard
    - `:mq`: Multiquadratic
    - `:imq`: InverseMultiquadratic
    - `:tp`: ThinPlate
    - `:nn`: NearestNeighbour
    - `:ga`: Gaussian
- `nmethod::Symbol=:minmax`: method for normalization, see `normalize()`
- `plot_contours::Bools=true`: plot contours over topo plot
- `plot_electrodes::Bools=true`: plot electrodes over topo plot
- `ps::Symbol=:l`: plot size (`:l`: large (800×800 px), `:m`: medium (300×300 px), `:s`: small (100×100 px))
- `head::Bool=true`: draw head
- `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use polar coordinates for XY plane and spherical coordinates for XZ and YZ planes
- `threshold::Union{Nothing, Real}=nothing`: if set, use threshold to mark a region
- `threshold_type::Symbol=:neq`: rule for thresholding:
    - `:eq`: draw region is values are equal to threshold
    - `:neq`: draw region is values are not equal to threshold
    - `:geq`: draw region is values are ≥ to threshold
    - `:leq`: draw region is values are ≤ to threshold
    - `:g`: draw region is values are > to threshold
    - `:l`: draw region is values are < to threshold

# Returns

- `p::GLMakie.Figure`
"""
function mplot_topo(obj::NeuroAnalyzer.NEURO; ep::Union{Int64, AbstractRange}=0, ch::Union{String, Vector{String}, Regex}, seg::Tuple{Real, Real}=(0, 10), title::String="default", mono::Bool=false, cb::Bool=true, cb_title::String="default", amethod::Symbol=:mean, imethod::Symbol=:sh, nmethod::Symbol=:minmax, plot_contours::Bool=true, plot_electrodes::Bool=true, ps::Symbol=:l, head::Bool=true, cart::Bool=false, threshold::Union{Nothing, Real}=nothing, threshold_type::Symbol=:neq)::GLMakie.Figure

    if obj.time_pts[end] < 10 && seg == (0, 10)
        seg = (0, obj.time_pts[end])
    else
        _check_segment(obj, seg)
    end
    seg = (vsearch(seg[1], obj.time_pts), vsearch(seg[2], obj.time_pts))

    _check_var(imethod, [:sh, :mq, :imq, :tp, :nn, :ga], "imethod")
    _check_var(amethod, [:mean, :median], "amethod")

    if ep != 0
        _check_epochs(obj, ep)
        if nepochs(obj) == 1
            ep = 0
        else
            seg = (((ep[1] - 1) * epoch_len(obj) + 1), seg[2])
            if ep isa Int64
                seg = (seg[1], (seg[1] + epoch_len(obj) - 1))
            else
                seg = (seg[1], (ep[end] * epoch_len(obj)))
            end
            ep = 0
        end
    end

    ch = get_channel(obj, ch=ch)
    @assert length(ch) >= 2 "plot_topo() requires ≥ 2 channels."
    chs = intersect(obj.locs[!, :label], labels(obj)[ch])
    locs = Base.filter(:label => in(chs), obj.locs)
    _check_ch_locs(ch, labels(obj), obj.locs[!, :label])

    # get time vector
    if seg[2] <= epoch_len(obj)
        s = obj.data[ch, seg[1]:seg[2], 1]
    else
        s = epoch(obj, ep_n=1).data[ch, seg[1]:seg[2], 1]
    end
    # t = _get_t(seg[1], seg[2], sr(obj))
    t = obj.time_pts[seg[1]:seg[2]]
    _, t_s1, _, t_s2 = _convert_t(t[1], t[end])
    ep = _s2epoch(obj, seg[1], seg[2])

    # average signal and convert to vector
    if size(s, 2) > 1
        if amethod === :mean
            s = vec(mean(s, dims=2))
        elseif amethod === :median
            s = vec(median(s, dims=2))
        end
    else
        s = vec(s)
    end

    if seg[2] != seg[1]
        title == "default" && (title = "Amplitude\n[$(string(amethod)) over time window: $t_s1:$t_s2]")
    else
        title == "default" && (title = "Amplitude\n[time point: $t_s1]")
    end
    cb_title == "default" && (cb_title = "[A.U.]")

    p = mplot_topo(s,
                  ch=collect(1:DataFrames.nrow(locs)),
                  locs=locs,
                  cb=cb,
                  cb_title=cb_title,
                  title=title,
                  mono=mono,
                  imethod=imethod,
                  nmethod=nmethod,
                  plot_contours=plot_contours,
                  plot_electrodes=plot_electrodes,
                  ps=ps,
                  head=head,
                  cart=cart,
                  threshold=threshold,
                  threshold_type=threshold_type)

    return p

end

"""
    mplot_topo(obj; <keyword arguments>)

Topographical plot of embedded or external component.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `c::Union{Symbol, AbstractArray}`: component to plot
- `ep::Union{Int64, AbstractRange}=0`: epoch to display
- `c_idx::Union{Int64, Vector{Int64}, AbstractRange}=0`: component channel to display, default is all component channels
- `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
- `title::String="default"`: plot title, default is Amplitude topographical plot [channels: 1:19, epoch: 1, time window: 0 ms:20 s]
- `mono::Bool=false`: use color or gray palette
- `cb::Bool=true`: plot color bar
- `cb_title::String="[A.U.]"`: color bar title
- `amethod::Symbol=:mean`: averaging method:
    - `:mean`
    - `:median`
- `imethod::Symbol=:sh`: interpolation method:
    - `:sh`: Shepard
    - `:mq`: Multiquadratic
    - `:imq`: InverseMultiquadratic
    - `:tp`: ThinPlate
    - `:nn`: NearestNeighbour
    - `:ga`: Gaussian
- `nmethod::Symbol=:minmax`: method for normalization, see `normalize()`
- `plot_contours::Bools=true`: plot contours over topo plot
- `plot_electrodes::Bools=true`: plot electrodes over topo plot
- `ps::Symbol=:l`: plot size (`:l`: large (800×800 px), `:m`: medium (300×300 px), `:s`: small (100×100 px))
- `head::Bool=true`: draw head
- `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use polar coordinates for XY plane and spherical coordinates for XZ and YZ planes
- `threshold::Union{Nothing, Real}=nothing`: if set, use threshold to mark a region
- `threshold_type::Symbol=:neq`: rule for thresholding:
    - `:eq`: draw region is values are equal to threshold
    - `:neq`: draw region is values are not equal to threshold
    - `:geq`: draw region is values are ≥ to threshold
    - `:leq`: draw region is values are ≤ to threshold
    - `:g`: draw region is values are > to threshold
    - `:l`: draw region is values are < to threshold

# Returns

- `p::GLMakie.Figure`
"""
function mplot_topo(obj::NeuroAnalyzer.NEURO, c::Union{Symbol, AbstractArray}; ep::Union{Int64, AbstractRange}=0, c_idx::Union{Int64, Vector{Int64}, AbstractRange}=0, seg::Tuple{Real, Real}=(0, 10), title::String="default", mono::Bool=false, cb::Bool=true, cb_title::String="default", amethod::Symbol=:mean, imethod::Symbol=:sh, nmethod::Symbol=:minmax, plot_contours::Bool=true, plot_electrodes::Bool=true, ps::Symbol=:l, head::Bool=true, cart::Bool=false, threshold::Union{Nothing, Real}=nothing, threshold_type::Symbol=:neq)::GLMakie.Figure

    if obj.time_pts[end] < 10 && seg == (0, 10)
        seg = (0, obj.time_pts[end])
    else
        _check_segment(obj, seg)
    end

    seg = (vsearch(seg[1], obj.time_pts), vsearch(seg[2], obj.time_pts))

    _has_locs(obj)
    _check_var(imethod, [:sh, :mq, :imq, :tp, :nn, :ga], "imethod")
    _check_var(amethod, [:mean, :median], "amethod")

    time_segment = true

    if c isa Matrix{Float64}
        c = reshape(c, size(c, 1), size(c, 2), 1)
        time_segment = false
    elseif c isa Array{Float64, 3} && size(c, 2) != epoch_len(obj)
        time_segment = false
    elseif c isa Vector{Float64}
        c = reshape(c, length(c), 1, 1)
        time_segment = false
    end

    if time_segment
        if ep != 0
            _check_epochs(obj, ep)
            if nepochs(obj) == 1
                ep = 0
            else
                seg = (((ep[1] - 1) * epoch_len(obj) + 1), seg[2])
                if ep isa Int64
                    seg = (seg[1], (seg[1] + epoch_len(obj) - 1))
                else
                    seg = (seg[1], (ep[end] * epoch_len(obj)))
                end
                ep = 0
            end
        end
    end

    # select component channels, default is all channels
    c isa Symbol && (c = _get_component(obj, c).c)
    c_idx == 0 && (c_idx = _select_cidx(c, c_idx))
    _check_cidx(c, c_idx)
    clabels = _gen_clabels(c)[c_idx]
    c_idx isa Int64 && (clabels = [clabels])

    @assert length(c_idx) >= 2 "plot_topo() requires ≥ 2 channels."
    chs = intersect(obj.locs[!, :label], labels(obj)[c_idx])
    locs = Base.filter(:label => in(chs), obj.locs)
    @assert length(chs) == DataFrames.nrow(locs) "Some channels do not have locations."

    # get time vector
    if time_segment
        if seg[2] <= epoch_len(obj)
            s = c[c_idx, seg[1]:seg[2], 1]
        else
            s = _make_epochs(c, ep_n=1)[c_idx, seg[1]:seg[2], 1]
        end
        if seg[1] != seg[2]
            t = _get_t(seg[1], seg[2], sr(obj))
        else
            t = _get_t(seg[1], seg[2] + 1, sr(obj))
        end
        _, t_s1, _, t_s2 = _convert_t(t[1], t[end])
        ep = _s2epoch(obj, seg[1], seg[2])
    else
        s = c
    end

    # average signal and convert to vector
    if size(s, 2) > 1
        if amethod === :mean
            s = vec(mean(s, dims=2))
        elseif amethod === :median
            s = vec(median(s, dims=2))
        end
    else
        s = vec(s)
    end

    if time_segment
        if seg[2] != seg[1]
          title == "default" && (title = "Amplitude\n[channel$(_pl(length(ch))): $(_channel2channel_name(ch)), epoch$(_pl(length(ep))): $ep, $(string(amethod)) over time window: $t_s1:$t_s2]")
        else
            title == "default" && (title = "Amplitude\n[channel$(_pl(length(ch))): $(_channel2channel_name(ch)), epoch$(_pl(length(ep))): $ep, time point: $t_s1]")
        end
    else
        title == "default" && (title = "Amplitude\n[component$(_pl(length(c_idx))): $(_channel2channel_name(c_idx))]")
    end

    cb_title == "default" && (cb_title = "[A.U.]")

    p = mplot_topo(s,
                  ch=c_idx,
                  locs=locs,
                  cb=cb,
                  cb_title=cb_title,
                  title=title,
                  mono=mono,
                  imethod=imethod,
                  nmethod=nmethod,
                  plot_contours=plot_contours,
                  plot_electrodes=plot_electrodes,
                  ps=ps,
                  head=head,
                  cart=cart,
                  threshold=threshold,
                  threshold_type=threshold_type)

    return p

end
