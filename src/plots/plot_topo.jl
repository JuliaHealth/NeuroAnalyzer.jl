export plot_topo

"""
    plot_topo(c; <keyword arguments>)

Plot topographical view.

# Arguments

- `s::Vector{<:Real}`: values to plot (one value per channel)
- `locs::DataFrame`: columns: channel, labels, loc_radius, loc_theta, loc_x, loc_y, loc_z, loc_radius_sph, loc_theta_sph, loc_phi_sph
- `ch::Union{Int64, Vector{Int64}}=1:nrow(locs)`: list of channels, default is all channels
- `cb::Bool=true`: plot color bar
- `cb_label::String="[A.U.]"`: color bar label
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
- `large::Bool=true`: draw large (size of electrodes area 600×600 px, more details) or small (size of electrodes area 240×240 px, less details) plot
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
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_topo(s::Vector{<:Real}; locs::DataFrame, ch::Union{Int64, Vector{Int64}}=1:nrow(locs), cb::Bool=true, cb_label::String="[A.U.]", title::String="default", mono::Bool=false, imethod::Symbol=:sh, nmethod::Symbol=:minmax, plot_contours::Bool=true, plot_electrodes::Bool=true, large::Bool=true, head::Bool=true, cart::Bool=false, threshold::Union{Nothing, Real}=nothing, threshold_type::Symbol=:neq, kwargs...)::Plots.Plot{Plots.GRBackend}

    pal = mono ? :grays : :bluesreds
    _check_var(imethod, [:sh, :mq, :imq, :tp, :nn, :ga], "imethod")

    locs = locs[ch, :]

    if large
        head_shape = FileIO.load(joinpath(res_path, "head_t_outline_large.png"))
        head_mask = FileIO.load(joinpath(res_path, "mask_large.png"))
    else
        head_shape = FileIO.load(joinpath(res_path, "head_t_small.png"))
        head_mask = FileIO.load(joinpath(res_path, "mask_small.png"))
    end

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

    s_interpolated, interpolated_x, interpolated_y = _interpolate2d(s, loc_x, loc_y, 100, imethod, nmethod)

    head12 = false
    maximum(abs.(locs[:, :loc_x])) <= 1.2 && maximum(abs.(locs[:, :loc_y])) <= 1.2 && maximum(abs.(locs[:, :loc_z])) <= 1.5 && (head12 = true)

    if head12
        xt = (round.(linspace(0, size(head_shape, 1), 25)), string.(round.(linspace(-1.2, 1.2, 25), digits=1)))
        yt = (round.(linspace(0, size(head_shape, 2), 25)), string.(round.(linspace(1.2, -1.2, 25), digits=1)))
        interpolated_x = round.(linspace(0, size(head_shape, 1), length(interpolated_x)), digits=2)
        interpolated_y = round.(linspace(0, size(head_shape, 2), length(interpolated_y)), digits=2)
        xl = (0, size(head_shape, 1))
        yl = (0, size(head_shape, 2))
        !head && (loc_y .= -loc_y)
        origin = size(head_shape) ./ 2 .+ 1
        if large
            marker_size = length(ch) > 64 ? 4 : 8
            loc_x = @. round(origin[1] + (loc_x * 250), digits=2)
            loc_y = @. round(origin[2] - (loc_y * 250), digits=2)
        else
            marker_size = length(ch) > 64 ? 2 : 4
            loc_x = @. round(origin[1] + (loc_x * 100), digits=2)
            loc_y = @. round(origin[2] - (loc_y * 100), digits=2)
        end
    else
        if large
            m = zeros(RGBA{FixedPointNumbers.N0f8}, size(head_shape) .+ 200)
            m[101:100+size(head_shape, 1), 101:100+size(head_shape, 2)] .= head_shape
            head_shape = m
            xt = (round.(linspace(0, size(head_shape, 1), 25)), string.(round.(linspace(-1.6, 1.6, 25), digits=1)))
            yt = (round.(linspace(0, size(head_shape, 2), 25)), string.(round.(linspace(1.6, -1.6, 25), digits=1)))
            interpolated_x = round.(linspace(0, size(head_shape, 1), length(interpolated_x)), digits=2)
            interpolated_y = round.(linspace(0, size(head_shape, 2), length(interpolated_y)), digits=2)
            xl = (0, size(head_shape, 1))
            yl = (0, size(head_shape, 2))
            !head && (loc_y .= -loc_y)
            origin = size(head_shape) ./ 2 .+ 1
            marker_size = length(ch) > 64 ? 4 : 8
            loc_x = @. round(origin[1] + (loc_x * 250), digits=2)
            loc_y = @. round(origin[2] - (loc_y * 250), digits=2)
            length(ch) > 64 && (ch_labels = false)
        else
            m = zeros(RGBA{FixedPointNumbers.N0f8}, size(head_shape) .+ 100)
            m[51:50+size(head_shape, 1), 51:50+size(head_shape, 2)] .= head_shape
            head_shape = m
            xt = (round.(linspace(0, size(head_shape, 1), 25)), string.(round.(linspace(-1.6, 1.6, 25), digits=1)))
            yt = (round.(linspace(0, size(head_shape, 2), 25)), string.(round.(linspace(1.6, -1.6, 25), digits=1)))
            interpolated_x = round.(linspace(0, size(head_shape, 1), length(interpolated_x)), digits=2)
            interpolated_y = round.(linspace(0, size(head_shape, 2), length(interpolated_y)), digits=2)
            xl = (0, size(head_shape, 1))
            yl = (0, size(head_shape, 2))
            !head && (loc_y .= -loc_y)
            origin = size(head_shape) ./ 2 .+ 1
            marker_size = length(ch) > 64 ? 2 : 4
            font_size = 4
            ch_labels = false
            sch_labels = false
            grid = false
            loc_x = @. round(origin[1] + (loc_x * 100), digits=2)
            loc_y = @. round(origin[2] - (loc_y * 100), digits=2)
        end
    end

    if large
        font_size = 10
        !occursin("\n", title) && title !== "" && (title *= "\n")
    else
        title=""
        cb_label=""
        marker_size = 2
        font_size = 2
    end

    if large
        if cb
            p = Plots.plot(grid=false,
                           framestyle=:none,
                           border=:none,
                           palette=pal,
                           aspect_ratio=1,
                           size=head12 ? size(head_shape) .+ 150 : size(head_shape) .+ 175,
                           right_margin=10*Plots.px,
                           bottom_margin=0*Plots.px,
                           top_margin=-50*Plots.px,
                           left_margin=0*Plots.px,
                           titlefontsize=font_size,
                           colorbar=cb,
                           colorbar_title=cb_label,
                           colorbar_tickfontsize=1,
                           colorbar_titlefontsize=8,
                           xlims=xl,
                           ylims=yl,
                           title=title;
                           kwargs...)
        else
            p = Plots.plot(grid=false,
                           framestyle=:none,
                           border=:none,
                           colorbar=cb,
                           palette=pal,
                           aspect_ratio=1,
                           size=head12 ? size(head_shape) .+ 105 : size(head_shape) .+ 110,
                           right_margin=-10*Plots.px,
                           bottom_margin=10*Plots.px,
                           top_margin=20*Plots.px,
                           left_margin=-30*Plots.px,
                           ticks_fontsize=font_size,
                           xticks=xt,
                           yticks=yt,
                           xlims=xl,
                           ylims=yl,
                           titlefontsize=font_size,
                           title=title;
                           kwargs...)
        end
    else
        if cb
            p = Plots.plot(grid=false,
                           framestyle=:none,
                           border=:none,
                           palette=pal,
                           aspect_ratio=1,
                           size=head12 ? size(head_shape) .+ 45 : (size(head_shape)[1] + 58, size(head_shape)[2] + 55),
                           right_margin=head12 ? -30*Plots.px : -20*Plots.px,
                           bottom_margin=head12 ? -30*Plots.px : -25*Plots.px,
                           top_margin=head12 ? -30*Plots.px : -20*Plots.px,
                           left_margin=head12 ? -20*Plots.px : -20*Plots.px,
                           titlefontsize=font_size,
                           colorbar=cb,
                           colorbar_titlefontsize=font_size - 2,
                           colorbar_title=cb_label,
                           colorbar_ticks=false,
                           xlims=xl,
                           ylims=yl,
                           title=title;
                           kwargs...)
        else
            p = Plots.plot(grid=false,
                           framestyle=:none,
                           border=:none,
                           colorbar=cb,
                           palette=pal,
                           aspect_ratio=1,
                           size=head12 ? size(head_shape) .+ 2 : size(head_shape) .+ 1,
                           right_margin=-10*Plots.px,
                           bottom_margin=-30*Plots.px,
                           top_margin=-20*Plots.px,
                           left_margin=-40*Plots.px,
                           xticks=xt,
                           yticks=yt,
                           xlims=xl,
                           ylims=yl,
                           foreground_color=:black)
        end
    end

    p = Plots.plot!(interpolated_x,
                    interpolated_y,
                    s_interpolated,
                    fill=:bluesreds,
                    seriestype=:heatmap,
                    seriescolor=pal,
                    levels=10,
                    linewidth=0)

    if plot_contours
        p = Plots.plot!(interpolated_x,
                        interpolated_y,
                        s_interpolated,
                        fill=:bluesreds,
                        seriestype=:contour,
                        seriescolor=pal,
                        colorbar_title=cb_label,
                        levels=5,
                        linecolor=:black,
                        linewidth=0.2)
    end

    # draw electrodes
    if plot_electrodes
        p = Plots.scatter!((loc_x, loc_y),
                            color=:black,
                            markerstrokecolor=Colors.RGBA(255/255, 255/255, 255/255, 0/255),
                            label="",
                            markershape=:circle,
                            markersize=marker_size,
                            markerstrokewidth=0,
                            markerstrokealpha=0)
    end

    # draw region
    if !isnothing(threshold)
        _, bm = seg_extract(s_interpolated, threshold=threshold, threshold_type=threshold_type)
        reg = ones(size(s_interpolated)) .* minimum(s_interpolated)
        reg[bm] .= maximum(s_interpolated)
        p = Plots.plot!(interpolated_x,
                        interpolated_y,
                        reg,
                        seriestype=:contour,
                        levels=1,
                        linecolor=:black,
                        colorbar_entry=false,
                        linewidth=2;
                        kwargs=kwargs)
    end

    # draw head
    if head
        if large
            cb ? head_mask = head_mask[158:end-140, 147:end-140] : head_mask = head_mask[158:end-140, 148:end-135]
            p = Plots.plot!(head_shape)
            if head12
                p = Plots.plot!(head_mask)
                p = Plots.plot!(Shape([0, size(head_shape, 1), size(head_shape, 1), 0], [0, 0, size(head_shape, 2), size(head_shape, 2)]), lc=:white, lw=1, fill=nothing, legend=false)
            end
        else
            head_mask = head_mask[80:end-80, 80:end-80]
            p = Plots.plot!(head_shape)
            if head12
                p = Plots.plot!(head_mask)
                p = Plots.plot!(Shape([0, size(head_shape, 1), size(head_shape, 1), 0], [0, 0, size(head_shape, 2), size(head_shape, 2)]), lc=:white, lw=2, fill=nothing, legend=false)
            end
        end
    end

    return p

end

"""
    plot_topo(obj; <keyword arguments>)

Topographical plot.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ep::Union{Int64, AbstractRange}=0`: epoch to display
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
- `title::String="default"`: plot title, default is Amplitude topographical plot [channels: 1:19, epoch: 1, time window: 0 ms:20 s]
- `mono::Bool=false`: use color or gray palette
- `cb::Bool=true`: plot color bar
- `cb_label::String="[A.U.]"`: color bar label
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
- `large::Bool=true`: draw large (size of electrodes area 600×600 px, more details) or small (size of electrodes area 240×240 px, less details) plot
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
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_topo(obj::NeuroAnalyzer.NEURO; ep::Union{Int64, AbstractRange}=0, ch::Union{String, Vector{String}, Regex}, seg::Tuple{Real, Real}=(0, 10), title::String="default", mono::Bool=false, cb::Bool=true, cb_label::String="default", amethod::Symbol=:mean, imethod::Symbol=:sh, nmethod::Symbol=:minmax, plot_contours::Bool=true, plot_electrodes::Bool=true, large::Bool=true, head::Bool=true, cart::Bool=false, threshold::Union{Nothing, Real}=nothing, threshold_type::Symbol=:neq, kwargs...)::Plots.Plot{Plots.GRBackend}

    if obj.time_pts[end] < 10 && seg == (0, 10)
        seg = (0, obj.time_pts[end])
    else
        _check_segment_topo(obj, seg)
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
    cb_label == "default" && (cb_label = "[A.U.]")

    p = plot_topo(s,
                  ch=collect(1:nrow(locs)),
                  locs=locs,
                  cb=cb,
                  cb_label=cb_label,
                  title=title,
                  mono=mono,
                  imethod=imethod,
                  nmethod=nmethod,
                  plot_contours=plot_contours,
                  plot_electrodes=plot_electrodes,
                  large=large,
                  head=head,
                  cart=cart,
                  threshold=threshold,
                  threshold_type=threshold_type,
                  kwargs=kwargs)

    return p

end

"""
    plot_topo(obj; <keyword arguments>)

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
- `cb_label::String="[A.U.]"`: color bar label
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
- `large::Bool=true`: draw large (size of electrodes area 600×600 px, more details) or small (size of electrodes area 240×240 px, less details) plot
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
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_topo(obj::NeuroAnalyzer.NEURO, c::Union{Symbol, AbstractArray}; ep::Union{Int64, AbstractRange}=0, c_idx::Union{Int64, Vector{Int64}, AbstractRange}=0, seg::Tuple{Real, Real}=(0, 10), title::String="default", mono::Bool=false, cb::Bool=true, cb_label::String="default", amethod::Symbol=:mean, imethod::Symbol=:sh, nmethod::Symbol=:minmax, plot_contours::Bool=true, plot_electrodes::Bool=true, large::Bool=true, head::Bool=true, cart::Bool=false, threshold::Union{Nothing, Real}=nothing, threshold_type::Symbol=:neq, kwargs...)::Plots.Plot{Plots.GRBackend}

    if obj.time_pts[end] < 10 && seg == (0, 10)
        seg = (0, obj.time_pts[end])
    else
        _check_segment_topo(obj, seg)
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
    @assert length(chs) == nrow(locs) "Some channels do not have locations."

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

    cb_label == "default" && (cb_label = "[A.U.]")

    p = plot_topo(s,
                  ch=c_idx,
                  locs=locs,
                  cb=cb,
                  cb_label=cb_label,
                  title=title,
                  mono=mono,
                  imethod=imethod,
                  nmethod=nmethod,
                  plot_contours=plot_contours,
                  plot_electrodes=plot_electrodes,
                  large=large,
                  head=head,
                  cart=cart,
                  threshold=threshold,
                  threshold_type=threshold_type,
                  kwargs=kwargs)

    return p

end
