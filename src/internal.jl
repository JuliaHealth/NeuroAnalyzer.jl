##################################
#                                #
#  Low-level internal functions  #
#                                #
##################################

_reflect(signal::AbstractArray) = vcat(signal[end:-1:1], signal, signal[end:-1:1])
_reflect(s1::AbstractArray, s2::AbstractArray, s3::AbstractArray) = vcat(s2[end:-1:1], s1, s3[end:-1:1])

_chop(signal::AbstractArray) = signal[(length(signal) ÷ 3 + 1):(length(signal) ÷ 3) * 2]
_chop(s1::AbstractArray, n::Int64) = s1[(n + 1):end - n]

_xlims(t::Union{Vector{<:Real}, AbstractRange}) = floor(t[1], digits=2), ceil(t[end], digits=2)

_ticks(t::Union{Vector{<:Real}, AbstractRange}) = floor(t[1], digits=2):((ceil(t[end]) - floor(t[1])) / 10):ceil(t[end], digits=2)
_ticks(t::Tuple{Real, Real}) = floor(t[1], digits=2):((ceil(t[2]) - floor(t[1])) / 10):ceil(t[2], digits=2)

_erpticks(t::Union{Vector{<:Real}, AbstractRange}) = vcat(collect(range(floor(t[1], digits=2), 0, 3)), collect(range(0, ceil(t[end], digits=2), 9))[2:end])
_erpticks(t::Tuple{Real, Real}) = vcat(collect(range(floor(t[1], digits=2), 0, 3)), collect(range(0, ceil(t[2], digits=2), 9))[2:end])

_pl(x::Union{AbstractRange, AbstractVector}) = length(collect(x)) > 1 ? "s" : ""
_pl(x::Real) = x > 1 ? "s" : ""

_get_range(signal::Union{AbstractVector, AbstractArray}) = round(abs(minimum(signal)) + abs(maximum(signal)), digits=0)

function _check_channels(eeg::NeuroAnalyzer.EEG, channel::Union{Int64, Vector{Int64}, AbstractRange})
    for idx in channel
        (idx < 1 || idx > eeg_channel_n(eeg)) && throw(ArgumentError("channel must be ≥ 1 and ≤ $(eeg_channel_n(eeg))."))
    end
end

function _check_channels(eeg::NeuroAnalyzer.EEG, channel::Union{Int64, Vector{Int64}, AbstractRange}, type::Symbol)
    channels = eeg_get_channel_bytype(eeg, type=type)
    for idx in channel
        idx in channels || throw(ArgumentError("Channel $idx does not match type: $(uppercase(string(type))) signal channels."))
        (idx < 1 || idx > eeg_channel_n(eeg)) && throw(ArgumentError("channel must be ≥ 1 and ≤ $(eeg_channel_n(eeg))."))
    end
end

function _check_channels(channels::Union{Int64, Vector{Int64}, AbstractRange}, channel::Union{Int64, Vector{Int64}, AbstractRange})
    for idx in channel
        idx in channels || throw(ArgumentError("Channel $idx does not match signal channels."))
        (idx < 1 || idx > length(channels)) && throw(ArgumentError("channel must be ≥ 1 and ≤ $(eeg_channel_n(eeg))."))
    end
end

function _check_epochs(eeg::NeuroAnalyzer.EEG, epoch::Union{Int64, Vector{Int64}, AbstractRange})
    for idx in epoch
        (idx < 1 || idx > eeg_epoch_n(eeg)) && throw(ArgumentError("epoch must be ≥ 1 and ≤ $(eeg_epoch_n(eeg))."))
    end
end

function _check_cidx(eeg::NeuroAnalyzer.EEG, c::Symbol, c_idx::Union{Int64, Vector{Int64}, AbstractRange})
    c, _ = _get_component(eeg, c)
    for idx in c_idx
        (idx < 1 || idx > size(c, 1)) && throw(ArgumentError("c_idx must be ≥ 1 and ≤ $(size(c, 1))."))
    end
end

function _check_cidx(c::Array{Float64, 3}, c_idx::Union{Int64, Vector{Int64}, AbstractRange})
    for idx in c_idx
        (idx < 1 || idx > size(c, 1)) && throw(ArgumentError("c_idx must be ≥ 1 and ≤ $(size(c, 1))."))
    end
end

function _gen_clabels(eeg::NeuroAnalyzer.EEG, c::Symbol)
    c, _ = _get_component(eeg, c)
    clabels = Vector{String}()
    for idx in 1:size(c, 1)
        push!(clabels, lpad(string(idx), length(string(size(c, 1))), "0"))
    end
    return clabels
end

function _gen_clabels(c::Array{Float64, 3})
    clabels = Vector{String}()
    for idx in 1:size(c, 1)
        push!(clabels, lpad(string(idx), length(string(size(c, 1))), "0"))
    end
    return clabels
end

function _select_channels(eeg::NeuroAnalyzer.EEG, channel::Union{Int64, Vector{Int64}, AbstractRange}, def_chn::Int64=0)
    # select channels, default is all or def_chn
    def_chn > eeg_channel_n(eeg) && (def_chn = eeg_channel_n(eeg))
    def_chn == 0 && (def_chn = eeg_channel_n(eeg))
    channel == 0 && (channel = 1:def_chn)
    typeof(channel) <: AbstractRange && (channel = collect(channel))
    length(channel) > 1 && sort!(channel)
    return channel
end

function _select_epochs(eeg::NeuroAnalyzer.EEG, epoch::Union{Int64, Vector{Int64}, AbstractRange}, def_ep::Int64=0)
    # select epochs, default is all or def_ep
    def_ep > eeg_epoch_n(eeg) && (def_ep = eeg_epoch_n(eeg))
    def_ep == 0 && (def_ep = eeg_epoch_n(eeg))
    epoch == 0 && (epoch = 1:def_ep)
    typeof(epoch) <: AbstractRange && (epoch = collect(epoch))
    length(epoch) > 1 && sort!(epoch)
    return epoch
end

function _select_cidx(eeg::NeuroAnalyzer.EEG, c::Symbol, c_idx::Union{Int64, Vector{Int64}, AbstractRange}, def_cidx::Int64=0)
    # select component channels, default is all or def_cidx
    c, _ = _get_component(eeg, c)
    def_cidx > size(c, 1) && (def_cidx = size(c, 1))
    def_cidx == 0 && (def_cidx = size(c, 1))
    c_idx == 0 && (c_idx = 1:def_cidx)
    typeof(c_idx) <: AbstractRange && (c_idx = collect(c_idx))
    length(c_idx) > 1 && sort!(c_idx)
    for idx in c_idx
        (idx < 1 || idx > size(c, 1)) && throw(ArgumentError("c_idx must be ≥ 1 and ≤ $(size(c, 1))."))
    end
    return c_idx
end

function _select_cidx(c::AbstractArray, c_idx::Union{Int64, Vector{Int64}, AbstractRange}, def_cidx::Int64=0)
    # select component channels, default is all or def_cidx
    def_cidx > size(c, 1) && (def_cidx = size(c, 1))
    def_cidx == 0 && (def_cidx = size(c, 1))
    c_idx == 0 && (c_idx = 1:def_cidx)
    typeof(c_idx) <: AbstractRange && (c_idx = collect(c_idx))
    length(c_idx) > 1 && sort!(c_idx)
    for idx in c_idx
        (idx < 1 || idx > size(c, 1)) && throw(ArgumentError("c_idx must be ≥ 1 and ≤ $(size(c, 1))."))
    end
    return c_idx
end

function _get_component(eeg::NeuroAnalyzer.EEG, c::Symbol)
    c in eeg.eeg_header[:components] || throw(ArgumentError("Component $c not found."))
    c_idx = findfirst(isequal(c), eeg.eeg_header[:components])
    c = eeg.eeg_components[c_idx]
    return (c=c, c_idx=c_idx)
end

function _set_defaults(xl::String, yl::String, tt::String, x::String, y::String, t::String)
    xl == "default" && (xl = x)
    yl == "default" && (yl = y)
    tt == "default" && (tt = t)
    return xl, yl, tt
end

function _len(eeg::NeuroAnalyzer.EEG, len::Int64, def_l::Int64)
    # return default length: one epoch (if epoch_len_seconds < def_l) or def_l seconds
    if len == 0
        if eeg_epoch_len(eeg) > def_l * eeg_sr(eeg)
            len = def_l * eeg_sr(eeg)
        else
            len = eeg_epoch_len(eeg)
        end
    end
    return len
end

function _draw_head(p::Plots.Plot{Plots.GRBackend}; head_labels::Bool=true, head_details::Bool=true, topo::Bool=false, kwargs...)
    # Draw head over a topographical plot `p`.
    # - `p::Plots.Plot{Plots.GRBackend}`: electrodes plot
    # - `loc_x::Vector{<:Real}`: vector of x electrode position
    # - `loc_y::Vector{<:Real}`: vector of y electrode position
    # - `head_labels::Bool=true`: add text labels to the plot
    # - `topo::Bool=false`: if true, perform peripheral erasing for topo plots
    # - `kwargs`: optional arguments for plot() function
    # loc_x, loc_y = loc_y, loc_x
    pts = Plots.partialcircle(0, 2π, 100, 1.1)
    x, y = Plots.unzip(pts)
    maxx = maximum(x)
    maxy = maximum(y)
    minx = minimum(x)
    miny = minimum(y)
    head = Plots.Shape(x, y)
    if head_details == true
        nose = Plots.Shape([(-0.2, maxy - 0.015),
                            (0, maxy + 0.08),
                            (0.2, maxy - 0.015),
                            (-0.01, maxy)
                           ])
        ear_r = Plots.Shape([(maxx, -0.05),
                             (maxx - 0.005, 0.09),
                             (maxx + 0.02, 0.125),
                             (maxx + 0.04, 0.13),
                             (maxx + 0.06, 0.115),
                             (maxx + 0.075, 0.085),
                             (maxx + 0.07, 0),
                             (maxx + 0.08, -0.155),
                             (maxx + 0.05, -0.215),
                             (maxx + 0.015, -0.225),
                             (maxx - 0.016, -0.19)
                            ])
        ear_l = Plots.Shape([(minx, -0.05),
                             (minx + 0.005, 0.09),
                             (minx - 0.02, 0.125),
                             (minx - 0.04, 0.13),
                             (minx - 0.06, 0.115),
                             (minx - 0.075, 0.085),
                             (minx - 0.07, 0),
                             (minx - 0.08, -0.155),
                             (minx - 0.05, -0.215),
                             (minx - 0.015, -0.225),
                             (minx + 0.016, -0.19)
                            ])
        p = Plots.plot!(nose, fill=nothing, label="")
        p = Plots.plot!(ear_l, fill=nothing, label="")
        p = Plots.plot!(ear_r, fill=nothing, label="")
    end

    p = Plots.plot!(p, head, fill=nothing, label="")

    if head_labels == true
        p = Plots.plot!(annotation=(0, 1.05, Plots.text("IN", pointsize=4, halign=:center, valign=:center)))
        p = Plots.plot!(annotation=(0, -1.05, Plots.text("NAS", pointsize=4, halign=:center, valign=:center)))
        p = Plots.plot!(annotation=(-1.05, 0, Plots.text("LPA", pointsize=4, halign=:center, valign=:center, rotation=90)))
        p = Plots.plot!(annotation=(1.05, 0, Plots.text("RPA", pointsize=4, halign=:center, valign=:center, rotation=-90)))
    end

    if topo == true
        pts = Plots.partialcircle(0, 2π, 100, 1.4)
        x, y = Plots.unzip(pts)
        for idx in 1:0.001:1.7
            peripheral = Shape(x .* idx, y .* idx)
            p = Plots.plot!(p, peripheral, label="", fill=nothing, lc=:white)
        end
        p = Plots.plot!(xlims=(-1.4, 1.4), ylims=(-1.4, 1.4); kwargs...)
    end

    return p
end

function _get_epoch_markers(eeg::NeuroAnalyzer.EEG)
    return round.(s2t.(collect(1:eeg_epoch_len(eeg):eeg_epoch_len(eeg) * eeg_epoch_n(eeg)), eeg_sr(eeg)), digits=2)
end

function _get_t(from::Int64, to::Int64, fs::Int64)
    t = collect((from / fs):(1 / fs):(to / fs))
    t = t[1:(end - 1)]
    t[1] = floor(t[1], digits=2)
    t[2:(end - 1)] = round.(t[2:(end - 1)], digits=3)
    t[end] = ceil(t[end], digits=2)
    return t
end

function _check_segment(eeg::NeuroAnalyzer.EEG, from::Int64, to::Int64)
    from < 1 && throw(ArgumentError("from must be > 0."))
    to < 1 && throw(ArgumentError("to must be > 0."))
    to < from && throw(ArgumentError("to must be ≥ $from."))
    (from > eeg_signal_len(eeg)) && throw(ArgumentError("from must be ≤ $(eeg_signal_len(eeg))."))
    (to > eeg_signal_len(eeg)) && throw(ArgumentError("to must be ≤ $(eeg_signal_len(eeg))."))
end

function _check_segment(signal::AbstractVector, from::Int64, to::Int64)
    from < 0 && throw(ArgumentError("from must be > 0."))
    to < 0 && throw(ArgumentError("to must be > 0."))
    to < from && throw(ArgumentError("to must be ≥ $from."))
    from > length(signal) && throw(ArgumentError("from must be ≤ $(length(signal))."))
    to > length(signal) && throw(ArgumentError("to must be ≤ $(length(signal))."))
end

function _convert_t(t1::Float64, t2::Float64)
    abs(t1) < 1.0 && (ts1 = string(floor(t1 * 1000, digits=2)) * " ms")
    abs(t1) >= 1.0 && (ts1 = string(floor(t1, digits=2)) * " s")
    abs(t2) < 1.0 && (ts2 = string(ceil(t2 * 1000, digits=2)) * " ms")
    abs(t2) >= 1.0 && (ts2 = string(ceil(t2, digits=2)) * " s")
    return t1, ts1, t2, ts2
end

function _tuple_max(t::Tuple{Real, Real})
    abs(t[1]) > abs(t[2]) && (t = (-abs(t[1]), abs(t[1])))
    abs(t[1]) < abs(t[2]) && (t = (-abs(t[2]), abs(t[2])))
    return t
end

function _channel2channel_name(channel::Union{Int64, Vector{Int64}, AbstractRange})
    if typeof(channel) == Int64
        return channel
    else
        if collect(channel[1]:channel[end]) == channel
            channel_name = string(channel[1]) * ":" * string(channel[end])
        else
            channel_name = "" 
            for idx in 1:(length(channel) - 1)
                channel_name *= string(channel[idx])
                channel_name *= ", "
            end
            channel_name *= string(channel[end])
        end
    end
    return channel_name
end

function _s2epoch(eeg::NeuroAnalyzer.EEG, from::Int64, to::Int64)
    epoch = floor(Int64, from / eeg_epoch_len(eeg)):ceil(Int64, to / eeg_epoch_len(eeg))
    from / eeg_epoch_len(eeg) > from ÷ eeg_epoch_len(eeg) && (epoch = epoch[1] + 1:epoch[end])
    epoch[1] == 0 && (epoch = 1:epoch[end])
    epoch[1] == epoch[end] && (epoch = epoch[1])
    return epoch
end

function _epoch2s(eeg::NeuroAnalyzer.EEG, epoch::Int64)
    t1 = (epoch - 1) * eeg_epoch_len(eeg) + 1
    t2 = epoch * eeg_epoch_len(eeg)
    return t1, t2
end

function _fir_response(f::Vector{<:Real}, w=range(0, stop=π, length=1024))
    # code based on Matti Pastell "FIR filter design with Julia"
    n = length(w)
    h = Array{ComplexF32}(undef, n)
    sw = 0
    for i = 1:n
        for j = 1:length(f)
            sw += f[j] * exp(-im * w[i])^-j
        end
        h[i] = sw
        sw = 0
    end
    return h
end

function _make_epochs(signal::Matrix{<:Real}; epoch_n::Union{Int64, Nothing}=nothing, epoch_len::Union{Int64, Nothing}=nothing)

    (epoch_len === nothing && epoch_n === nothing) && throw(ArgumentError("Either epoch_n or epoch_len must be specified."))
    (epoch_len !== nothing && epoch_n !== nothing) && throw(ArgumentError("Both epoch_n and epoch_len cannot be specified."))
    (epoch_len !== nothing && epoch_len < 1) && throw(ArgumentError("epoch_len must be ≥ 1."))
    (epoch_n !== nothing && epoch_n < 1) && throw(ArgumentError("epoch_n must be ≥ 1."))

    channel_n = size(signal, 1)
    if epoch_n === nothing
        epoch_n = size(signal, 2) ÷ epoch_len
    else
        epoch_len = size(signal, 2) ÷ epoch_n
    end

    epoch_len > size(signal, 2) && throw(ArgumentError("epoch_len must be ≤ $(size(signal, 2))."))
    epoch_len < 1 && throw(ArgumentError("epoch_len must be ≥ 1."))
    epoch_n < 1 && throw(ArgumentError("epoch_n must be ≥ 1."))

    epochs = reshape(signal[:, 1:(epoch_len * epoch_n)], channel_n, epoch_len, epoch_n)

    return epochs
end

function _make_epochs(signal::Array{T, 3}; epoch_n::Union{Int64, Nothing}=nothing, epoch_len::Union{Int64, Nothing}=nothing) where {T <: Real}

    (epoch_len === nothing && epoch_n === nothing) && throw(ArgumentError("Either epoch_n or epoch_len must be specified."))
    (epoch_len !== nothing && epoch_n !== nothing) && throw(ArgumentError("Both epoch_n and epoch_len cannot be specified."))
    (epoch_len !== nothing && epoch_len < 1) && throw(ArgumentError("epoch_len must be ≥ 1."))
    (epoch_n !== nothing && epoch_n < 1) && throw(ArgumentError("epoch_n must be ≥ 1."))

    channel_n = size(signal, 1)
    if epoch_n === nothing
        epoch_n = size(signal, 2) * size(signal, 3) ÷ epoch_len
    else
        epoch_len = size(signal, 2) * size(signal, 3) ÷ epoch_n
    end
    signal = reshape(signal, channel_n, (size(signal, 2) * size(signal, 3)), 1)
    epochs = reshape(signal[:, 1:(epoch_len * epoch_n), 1], channel_n, epoch_len, epoch_n)

    return epochs
end

function _make_epochs_bymarkers(signal::Array{<:Real, 3}; markers::DataFrame, marker_start::Vector{Int64}, epoch_offset::Int64, epoch_len::Int64)

    if size(signal, 3) > 1
        _info("EEG has already been epoched, parts of the signal might have been removed.")
        signal = reshape(signal, channel_n, (size(signal, 2) * size(signal, 3)), 1)
    end

    marker_n = length(marker_start)
    epoch_start = marker_start .- epoch_offset
    epoch_end = epoch_start .+ epoch_len .- 1

    # delete epochs outside signal limits
    for marker_idx in marker_n:-1:1
        if epoch_start[marker_idx] < 1
            deleteat!(epoch_start, marker_idx)
            deleteat!(epoch_end, marker_idx)
            marker_n -= 1
        elseif epoch_end[marker_idx] > size(signal, 2)
            deleteat!(epoch_start, marker_idx)
            deleteat!(epoch_end, marker_idx)
            marker_n -= 1
        end
    end

    (epoch_offset < 1) && throw(ArgumentError("epoch_offset must be ≥ 1."))
    (epoch_len !== nothing && epoch_len < 1) && throw(ArgumentError("epoch_len must be ≥ 1."))

    (epoch_offset + epoch_len > size(signal, 2)) && throw(ArgumentError("epoch_offset + epoch_len must be ≤ signal length $(size(signal, 2))."))

    epochs = zeros(size(signal, 1), epoch_len, marker_n)
    for marker_idx in 1:marker_n
        epochs[:, :, marker_idx] = reshape(signal[:, epoch_start[marker_idx]:epoch_end[marker_idx], :],
                                           size(signal, 1), 
                                           epoch_len)
    end

    # remove markers outside epoch limits
    for marker_idx in nrow(markers):-1:1
        within_epoch = false
        for epoch_idx in 1:marker_n
            markers[!, :start][marker_idx] in epoch_start[epoch_idx]:epoch_end[epoch_idx] && (within_epoch = true)
        end
        within_epoch == false && deleteat!(markers, marker_idx)
    end

    # calculate marker offsets
    for epoch_idx in 1:marker_n
        for marker_idx in 1:nrow(markers)
            if markers[!, :start][marker_idx] in epoch_start[epoch_idx]:epoch_end[epoch_idx]
                markers[!, :start][marker_idx] = (epoch_idx - 1) * epoch_len + markers[!, :start][marker_idx] .- epoch_start[epoch_idx]
            end
        end
    end

    return epochs, markers
end

function _get_channel_idx(labels::Vector{String}, channel::Union{String, Int64})
    if typeof(channel) == String
        channel_found = nothing
        for idx in eachindex(labels)
            if channel == labels[idx]
                channel_found = idx
            end
        end
        if channel_found === nothing
            throw(ArgumentError("channel name does not match signal labels."))
        end
    else
        channel < 1 || channel > length(labels) && throw(ArgumentError("channel index does not match signal channels."))
        channel_found = channel
    end

    return channel_found
end

function _angle_quadrant(a::Real)
    if a >= 0
        a = mod(a, 360)
        a <= 90 && (q = 1)
        (a > 90 && a <= 180) && (q = 2)
        (a > 180 && a <= 270) && (q = 3)
        (a > 270 && a < 360) && (q = 4)
    else
        a = mod(a, -360)
        a >= -90 && (q = 4)
        (a < -90 && a >= -180) && (q = 3)
        (a < -180 && a >= -270) && (q = 2)
        (a < -270 && a > -360) && (q = 1)
    end
    return q    
end 

function _locnorm(x::Union{AbstractVector, Real}, y::Union{AbstractVector, Real})
    xy = s_normalize_minmax(hcat(x, y))
    x = xy[:, 1]
    y = xy[:, 2]
    return x, y
end

function _locnorm(x::Union{AbstractVector, Real}, y::Union{AbstractVector, Real}, z::Union{AbstractVector, Real})
    xyz = s_normalize_minmax(hcat(x, y, z))
    x = xyz[:, 1]
    y = xyz[:, 2]
    z = xyz[:, 3]
    return x, y, z
end

function _round_locs(locs::DataFrame)
    locs[!, :loc_x] = round.(locs[!, :loc_x], digits=3)
    locs[!, :loc_y] = round.(locs[!, :loc_y], digits=3)
    locs[!, :loc_z] = round.(locs[!, :loc_z], digits=3)
    locs[!, :loc_radius] = round.(locs[!, :loc_radius], digits=3)
    locs[!, :loc_theta] = round.(locs[!, :loc_theta], digits=3)
    locs[!, :loc_radius_sph] = round.(locs[!, :loc_radius_sph], digits=3)
    locs[!, :loc_theta_sph] = round.(locs[!, :loc_theta_sph], digits=3)
    locs[!, :loc_phi_sph] = round.(locs[!, :loc_phi_sph], digits=3)
    return locs
end

function _free_gpumem(threshold::Real=0.95)
    m = CUDA.MemoryInfo()
    usedmem = m.total_bytes - m.free_bytes
    totalmem = m.total_bytes
    if usedmem / totalmem > threshold
        # CUDA.reclaim()
        GC.gc(true)
    end
end 

function _m2df(markers::Vector{String})
    # convert EDF/BDF markers to DataFrame
    markers = replace.(markers, "\x14\x14\0" => "|")
    markers = replace.(markers, "\x14\x14" => "|")
    markers = replace.(markers, "\x14" => "|")
    markers = replace.(markers, "\0" => "")
    a_start = Vector{Float64}()
    a_event = Vector{String}()
    # what about markers containing event duration?
    for idx in eachindex(markers)
        s = split(markers[idx], "|")
        if length(s) > 2
            push!(a_start, parse(Float64, strip(s[2])))
            push!(a_event, strip(s[3]))
        end
    end
    return DataFrame(:id => repeat([""], length(a_event)), :start => a_start, :length => zeros(Int64, length(a_event)), :description => a_event, :channel => zeros(Int64, length(a_event)))
end

function _clean_labels(labels::Vector{String})
    labels = replace.(labels, "EEG " => "")
    labels = replace.(labels, "EDF " => "")
    labels = replace.(labels, "BDF " => "")
    return labels
end

function _has_markers(channel_types::Vector{String})
    markers = false
    markers_channel = 0
    if "mrk" in channel_types
        markers = true
        markers_channel = nothing
        for channel_idx in eachindex(channel_types)
            channel_types[channel_idx] == "mrk" && (markers_channel = channel_idx)
        end
    end
    return markers, markers_channel
end

function _set_channel_types(labels::Vector{String})
    eeg_channels = ["af3", "af4", "af7", "af8", "afz", "c1", "c2", "c3", "c4", "c5", "c6", "cp1", "cp2", "cp3", "cp4", "cp5", "cp6", "cpz", "cz", "f1", "f10", "f2", "f3", "f4", "f5", "f6", "f7", "f8", "f9", "fc1", "fc2", "fc3", "fc4", "fc5", "fc6", "fcz", "fp1", "fp2", "fpz", "ft10", "ft7", "ft8", "ft9", "fz", "nz", "o1", "o2", "oz", "p1", "p10", "p2", "p3", "p4", "p5", "p6", "p7", "p8", "p9", "po3", "po4", "po7", "po8", "poz", "pz", "t10", "t3", "t4", "t5", "t6", "t7", "t8", "t9", "tp10", "tp7", "tp8", "tp9"]
    ref_channels = ["a1", "a2", "m1", "m2", "pg1", "pg2"]
    eog_channel = ["e", "e1", "e2"]
    channel_type = repeat(["???"], length(labels))
    for idx in eachindex(labels)
        occursin("meg", lowercase(labels[idx])) && (channel_type[idx] = "meg")
        occursin("ecg", lowercase(labels[idx])) && (channel_type[idx] = "ecg")
        occursin("ekg", lowercase(labels[idx])) && (channel_type[idx] = "ecg")
        occursin("eog", lowercase(labels[idx])) && (channel_type[idx] = "eog")
        occursin("rr", lowercase(labels[idx])) && (channel_type[idx] = "misc")
        occursin("mic", lowercase(labels[idx])) && (channel_type[idx] = "misc")
        occursin("flw", lowercase(labels[idx])) && (channel_type[idx] = "misc")
        occursin("tho", lowercase(labels[idx])) && (channel_type[idx] = "misc")
        occursin("abd", lowercase(labels[idx])) && (channel_type[idx] = "misc")
        occursin("sao2", lowercase(labels[idx])) && (channel_type[idx] = "misc")
        occursin("sa02", lowercase(labels[idx])) && (channel_type[idx] = "misc")
        occursin("plr", lowercase(labels[idx])) && (channel_type[idx] = "misc")
        occursin("body", lowercase(labels[idx])) && (channel_type[idx] = "misc")
        occursin("ux", lowercase(labels[idx])) && (channel_type[idx] = "misc")
        for idx2 in eachindex(eog_channel)
            occursin(eeg_channels[idx2], lowercase(labels[idx])) && (channel_type[idx] = "eog")
        end
        occursin("emg", lowercase(labels[idx])) && (channel_type[idx] = "emg")
        in(lowercase(labels[idx]), ref_channels) && (channel_type[idx] = "ref")
        for idx2 in eachindex(ref_channels)
            occursin(ref_channels[idx2], lowercase(labels[idx])) && (channel_type[idx] = "ref")
        end
        occursin("mark", lowercase(labels[idx])) && (channel_type[idx] = "mrk")
        occursin("marker", lowercase(labels[idx])) && (channel_type[idx] = "mrk")
        occursin("markers", lowercase(labels[idx])) && (channel_type[idx] = "mrk")
        occursin("event", lowercase(labels[idx])) && (channel_type[idx] = "mrk")
        occursin("events", lowercase(labels[idx])) && (channel_type[idx] = "mrk")
        occursin("annotation", lowercase(labels[idx])) && (channel_type[idx] = "mrk")
        occursin("annotations", lowercase(labels[idx])) && (channel_type[idx] = "mrk")
        occursin("status", lowercase(labels[idx])) && (channel_type[idx] = "mrk")
        # eeg channels should have priority, e.g. C3A1 (C3 referenced to A1 should be of eeg type, not ref)
        in(lowercase(labels[idx]), eeg_channels) && (channel_type[idx] = "eeg")
        for idx2 in eachindex(eeg_channels)
            occursin(eeg_channels[idx2], lowercase(labels[idx])) && (channel_type[idx] = "eeg")
        end
        lowercase(labels[idx])[1] == 'c' && (channel_type[idx] = "eeg")
        lowercase(labels[idx])[1] == 'f' && (channel_type[idx] = "eeg")
        lowercase(labels[idx])[1] == 'n' && (channel_type[idx] = "eeg")
        lowercase(labels[idx])[1] == 'o' && (channel_type[idx] = "eeg")
        lowercase(labels[idx])[1] == 'p' && (channel_type[idx] = "eeg")
        lowercase(labels[idx])[1] == 't' && (channel_type[idx] = "eeg")
        lowercase(labels[idx])[1] == 'i' && (channel_type[idx] = "eeg")
        (length(labels[idx]) > 1 && lowercase(labels[idx])[1:2] == "af") && (channel_type[idx] = "eeg")
    end
    return channel_type
end

function _check_var(s1::Symbol, s2::Vector{Symbol}, var::String)
    m = var * " must be "
    for idx in 1:(length(s2) - 2)
        m *= ":" * string(s2[idx]) * ", "
    end
    m *= ":" * string(s2[end - 1]) * " or :" * string(s2[end]) * "."
    s1 in s2 || throw(ArgumentError(m))
end

function _check_var(s1::String, s2::Vector{String}, var::String)
    m = var * " must be "
    for idx in 1:(length(s2) - 2)
        m *= ":" * s2[idx] * ", "
    end
    m *= ":" * s2[end - 1] * " or " * s2[end] * "."
    s1 in s2 || throw(ArgumentError(m))
end

function _interpolate(signal::AbstractVector, loc_x::Vector{Float64}, loc_y::Vector{Float64}, interpolation_factor::Int64=100, imethod::Symbol=:sh, nmethod::Symbol=:minmax)
    # `imethod::Symbol=:sh`: interpolation method Shepard (`:sh`), Multiquadratic (`:mq`), InverseMultiquadratic (`:imq`), ThinPlate (`:tp`), NearestNeighbour (`:nn`), Gaussian (`:ga`)
    _check_var(imethod, [:sh, :mq, :imq, :tp, :nn, :ga], "imethod")
    x_lim_int = (-1.4, 1.4)
    y_lim_int = (-1.4, 1.4)
    interpolated_x = linspace(x_lim_int[1], x_lim_int[2], interpolation_factor)
    interpolated_y = linspace(y_lim_int[1], y_lim_int[2], interpolation_factor)
    interpolated_x = round.(interpolated_x, digits=2)
    interpolated_y = round.(interpolated_y, digits=2)
    interpolation_m = Matrix{Tuple{Float64, Float64}}(undef, interpolation_factor, interpolation_factor)
    @inbounds @simd for idx1 in 1:interpolation_factor
        for idx2 in 1:interpolation_factor
            interpolation_m[idx1, idx2] = (interpolated_x[idx1], interpolated_y[idx2])
        end
    end
    signal_interpolated = zeros(interpolation_factor, interpolation_factor)
    electrode_locations = [loc_x loc_y]'
    imethod === :sh && (itp = ScatteredInterpolation.interpolate(Shepard(), electrode_locations, signal))
    imethod === :mq && (itp = ScatteredInterpolation.interpolate(Multiquadratic(), electrode_locations, signal))
    imethod === :imq && (itp = ScatteredInterpolation.interpolate(InverseMultiquadratic(), electrode_locations, signal))
    imethod === :tp && (itp = ScatteredInterpolation.interpolate(ThinPlate(), electrode_locations, signal))
    imethod === :nn && (itp = ScatteredInterpolation.interpolate(NearestNeighbor(), electrode_locations, signal))
    imethod === :ga && (itp = ScatteredInterpolation.interpolate(Gaussian(), electrode_locations, signal))
    @inbounds @simd for idx1 in 1:interpolation_factor
        for idx2 in 1:interpolation_factor
            signal_interpolated[idx1, idx2] = ScatteredInterpolation.evaluate(itp, [interpolation_m[idx1, idx2][1]; interpolation_m[idx1, idx2][2]])[1]
        end
    end
    signal_interpolated = signal_interpolated'
    
    return s_normalize(signal_interpolated, method=nmethod), interpolated_x, interpolated_y
end

function _dict2labeled_matrix(d::Dict; rev::Bool=true)
    l = Vector{String}()
    v = Vector{Vector{Float64}}()
    for (kk, vv) in d
        push!(l, kk)
        push!(v, vv)
    end
    rev == true && return reverse!(l), reverse!(c)
    rev == false && return l, c
end

function _labeled_matrix2dict(l::Vector{String}, v::Vector{Vector{Float64}})
    length(l) == length(v) || throw(ArgumentError("Length of labels and values do not match."))
    return Dict(zip(l, v))
end

function _map_channels(channel::Union{Int64, Vector{Int64}, AbstractRange}, channels=Vector{Int64})
    channel_orig = channel
    if typeof(channel) == Int64
        channel = vsearch(channel, channels)
    else
        for idx in eachindex(channel)
            channel[idx] = vsearch(channel[idx], channels)
        end
    end
    return channel, channel_orig
end

function _sort_channels(ch_t::Vector{String})
    replace!(ch_t, "eeg" => "1")
    replace!(ch_t, "meg" => "2")
    replace!(ch_t, "ref" => "3")
    replace!(ch_t, "eog" => "4")
    replace!(ch_t, "ecg" => "5")
    replace!(ch_t, "emg" => "6")
    replace!(ch_t, "other" => "7")
    replace!(ch_t, "markers" => "7")
    return sortperm(ch_t)
end

_c(n) = collect(1:n)

function _check_markers(markers::Vector{String}, marker::String)
    marker in markers || throw(ArgumentError("Marker: $marker not found in markers."))
end

function _check_markers(eeg::NeuroAnalyzer.EEG, marker::String)
    marker in unique(eeg.eeg_markers[!, :description]) || throw(ArgumentError("Marker: $marker not found in markers."))
end

function _delete_markers(markers::DataFrame, segment::Tuple{Int64, Int64})
    for marker_idx in nrow(markers):-1:1
        markers[marker_idx, :start] in segment[1]:segment[2] && deleteat!(markers, marker_idx)
    end
    return markers
end

function _shift_markers(m::DataFrame, pos::Int64, offset::Int64)
    markers = deepcopy(m)
    for marker_idx in 1:nrow(markers)
        markers[marker_idx, :start] > pos && (markers[marker_idx, :start] -= offset)
    end
    return markers
end

function _s2v(s::Union{<:Number, Vector{<:Number}})
    if typeof(s) <: Number
        return [s]
    else
        return s
    end
end

function _split(df::DataFrame, ratio::Float64=0.8)
    n = nrow(df)
    idx = shuffle(1:n)
    train_idx = view(idx, 1:floor(Int, ratio * n))
    test_idx = view(idx, (floor(Int, ratio * n) + 1):n)
    return df[train_idx, :], df[test_idx, :]
end

function _info(s::String)
    verbose == true && @info s
end

function _read_fiff_tag(fid::IOStream)
    tag = reinterpret(Int32, read(fid, sizeof(Int32) * 4))
    tag .= ntoh.(tag)
    tag_kind = tag[1]
    tag_type = tag[2]
    tag_size = tag[3]
    tag_next = tag[4]
    current_position = position(fid)

    if tag_next == 0
        seek(fid, current_position + tag_size)
    elseif tag_next > 0
        seek(fid, tag_next)
    end

    return tag_kind, tag_type, tag_size, tag_next
end

function _get_fiff_block_type(fid::IOStream, tags::Vector{NTuple{5, Int64}}, id)
    tag = tags[id]
    seek(fid, tag[1] + 16)
    buf = zeros(UInt8, tag[4])
    readbytes!(fid, buf, tag[4])
    return reinterpret(Int32, reverse(buf))
end

function _read_fiff_tag(fid::IOStream, fiff_blocks::Matrix{Int64}, tag_id::Int64)
    id = _find_fiff_tag(fiff_blocks, tag_id)
    if length(id) == 0
        return nothing
    elseif length(id) == 1
        out = _read_fiff_data(fid, fiff_blocks, id)
        if typeof(out) == String
            return out
        else
            return out[]
        end
    else
        out = Vector{Any}()
        for idx in 1:length(id)
            push!(out, _read_fiff_data(fid, fiff_blocks, id[idx]))
        end
        return out
    end
end

function _read_fiff_data(fid::IOStream, fiff_blocks::Matrix{Int64}, id::Union{Int64, Nothing})
    id === nothing && throw(ArgumentError("id does not contain valid tag id."))
    tag = fiff_blocks[id, :]
    seek(fid, tag[1] + 16)
    buf = zeros(UInt8, tag[4])
    readbytes!(fid, buf, tag[4])
    
    if tag[3] == 0
        return nothing
    elseif tag[3] == 1
        return reinterpret(Int8, buf)
    elseif tag[3] == 2
        return ntoh.(reinterpret(Int16, buf))
    elseif tag[3] == 3
        return ntoh.(reinterpret(Int32, buf))
    elseif tag[3] == 4
        return ntoh.(reinterpret(Float32, buf))
    elseif tag[3] == 5
        return ntoh.(reinterpret(Float64, buf))
    elseif tag[3] == 6
        return ntoh.(reinterpret(Float32, buf))
        # julian
    elseif tag[3] == 7
        return ntoh.(reinterpret(UInt16, buf))
    elseif tag[3] == 8
        return ntoh.(reinterpret(UInt32, buf))
    elseif tag[3] == 9
        return ntoh.(reinterpret(UInt64, buf))
    elseif tag[3] == 10
        return String(Char.(buf))
    elseif tag[3] == 11
        return ntoh.(reinterpret(Int64, buf))
    elseif tag[3] == 16
        return ntoh.(reinterpret(Int64, buf))
    elseif tag[3] == 20
        return ntoh.(reinterpret(ComplexF32, buf))
    elseif tag[3] == 21
        return ntoh.(reinterpret(ComplexF64, buf))
    elseif tag[3] == 21
        # old_pack
    elseif tag[3] == 30
        # ch_info_struct
        # recording order
        scan_no = ntoh.(reinterpret(Int32, buf[1:4]))[]
        # logical channels order
        log_no = ntoh.(reinterpret(Int32, buf[5:8]))[]
        # channel type
        kind = ntoh.(reinterpret(Int32, buf[9:12]))[]
        r = reinterpret(Float32, reverse(buf[13:16]))[]
        cal = reinterpret(Float32, reverse(buf[17:20]))[]
        coil_type = reinterpret(Int32, reverse(buf[21:24]))[]
        # coil coordinate system origin
        r01 = reinterpret(Float32, reverse(buf[25:28]))[]
        r02 = reinterpret(Float32, reverse(buf[29:32]))[]
        r03 = reinterpret(Float32, reverse(buf[33:36]))[]
        # coil coordinate system x-axis unit vector
        ex1 = reinterpret(Float32, reverse(buf[37:40]))[]
        ex2 = reinterpret(Float32, reverse(buf[41:44]))[]
        ex3 = reinterpret(Float32, reverse(buf[45:48]))[]
        # coil coordinate system y-axis unit vector
        ey1 = reinterpret(Float32, reverse(buf[49:52]))[]
        ey2 = reinterpret(Float32, reverse(buf[53:56]))[]
        ey3 = reinterpret(Float32, reverse(buf[57:60]))[]
        # coil coordinate system z-axis unit vector
        ez1 = reinterpret(Float32, reverse(buf[61:64]))[]
        ez2 = reinterpret(Float32, reverse(buf[65:68]))[]
        ez3 = reinterpret(Float32, reverse(buf[69:72]))[]
        unit = reinterpret(Int32, reverse(buf[73:76]))[]
        unit_mul = reinterpret(Int32, reverse(buf[77:80]))[]
        return(scan_no, log_no, kind, r, cal, coil_type, r01, r02, r03, ex1, ex2, ex2, ey1, ey2, ey2, ez1, ez2, ez2, unit, unit_mul)
    elseif tag[3] == 31
        # id_struct
        reverse(buf)
    elseif tag[3] == 32
        # dir_entry_struct
        reverse(buf)
    elseif tag[3] == 33
        # dig_point_struct
        reverse(buf)
    elseif tag[3] == 34
        # ch_pos_struct
        reverse(buf)
    elseif tag[3] == 35
        # coord_trans_struct
        reverse(buf)
    elseif tag[3] == 36
        # dig_string_struct
        reverse(buf)
    elseif tag[3] == 37
        # stream_segment_struct
        reverse(buf)
    else
        throw(ArgumentError("Unknown tag type $(tag[3])."))
    end
end

function _find_fiff_tag(fiff_blocks::Matrix{Int64}, id::Int64)
    ids = findall(x -> x == id, fiff_blocks[:, 2])
    if length(ids) == 1
        return ids[]
    else
        return ids
    end
end

function _extract_struct(s, id::Int64)
    id < 1 && throw(ArgumentError("id must be ≥ 1."))
    id > length(s[1]) && throw(ArgumentError("id must be ≤ $(length(s[1]))."))
    out = Vector{typeof(s[1][id])}()
    for channels_idx in 1:length(s)
        push!(out, s[channels_idx][id])
    end
    return out
end

function _fiff_channel_type(channel_types::Vector{Int32})
    channel_type = Vector{String}()
    for idx in 1:length(channel_types)
        if channel_types[idx] == 1
            push!(channel_type, "meg")
        elseif channel_types[idx] == 2
            push!(channel_type, "eeg")
        elseif channel_types[idx] == 3
            push!(channel_type, "stim")
        elseif channel_types[idx] == 102
            push!(channel_type, "bio")
        elseif channel_types[idx] == 201
            push!(channel_type, "mcg")
        elseif channel_types[idx] == 202
            push!(channel_type, "eog")
        elseif channel_types[idx] == 301
            push!(channel_type, "meg_ref")
        elseif channel_types[idx] == 302
            push!(channel_type, "emg")
        elseif channel_types[idx] == 402
            push!(channel_type, "ecg")
        elseif channel_types[idx] == 900
            push!(channel_type, "system")
        elseif channel_types[idx] == 910
            push!(channel_type, "ias")
        else
            push!(channel_type, "misc")
        end
    end
    return channel_type
end

function _create_fiff_block(file_name::String)

    fid = ""
    try
        fid = open(file_name, "r")
    catch
        error("File $file_name cannot be loaded.")
    end

    # read file_id tag
    tag_kind = nothing
    tag_type = nothing
    tag_size = nothing
    tag_next = nothing
    try
        tag_kind, tag_type, tag_size, tag_next = _read_fiff_tag(fid)
    catch
        error("File $file_name first tag cannot be read.")
    end

    # check file_id tag
    fiff_file_id = 100
    fiff_id_struct = 31
    fiff_next_seq = 0
    tag_kind != fiff_file_id && throw(ArgumentError("File $file_name is not a FIFF file."))
    tag_type != fiff_id_struct && throw(ArgumentError("File $file_name is not a FIFF file."))
    tag_size != 20 && throw(ArgumentError("File $file_name is not a FIFF file."))

    # read dir_pointer tag
    try
        tag_kind, tag_type, tag_size, tag_next = _read_fiff_tag(fid)
    catch
        error("File $file_name dir_pointer tag cannot be read.")
    end

    # check dir_pointer tag
    fiff_dir_pointer = 101
    tag_kind != fiff_dir_pointer && throw(ArgumentError("File $file_name has no dir_pointer tag."))

    # read tags
    seek(fid, 0)
    # tags: position in file, tag_id, data_type, data_size, next
    tags = Vector{Tuple{Int64, Int64, Int64, Int64, Int64}}()
    while tag_next != -1
        current_position = position(fid)
        tag_kind, tag_type, tag_size, tag_next = _read_fiff_tag(fid)
        push!(tags, (current_position, tag_kind, tag_type, tag_size, tag_next))
    end
    seek(fid, 0)

    # create list of tag IDs
    tag_pos = Vector{Int64}()
    for tag_idx in 1:length(tags)
        push!(tag_pos, tags[tag_idx][1])
    end
    tag_ids = Vector{Int64}()
    for tag_idx in 1:length(tags)
        push!(tag_ids, tags[tag_idx][2])
    end
    tag_type = Vector{Int64}()
    for tag_idx in 1:length(tags)
        push!(tag_type, tags[tag_idx][3])
    end
    tag_size = Vector{Int64}()
    for tag_idx in 1:length(tags)
        push!(tag_size, tags[tag_idx][4])
    end

    # create block structure
    fiff_block_start = 104
    fiff_block_end = 105
    fiff_block_nop = 108
    fiff_block_root = 999
    block_level = ones(Int64, length(tags))
    block_number = ones(Int64, length(tags))
    block_type = Vector{Int64}()
    block_type_current = 999
    for tag_idx in 1:length(tags)
        if tags[tag_idx][2] == fiff_block_start
            block_level[tag_idx:end] .+= 1
            block_type_current = _get_fiff_block_type(fid, tags, tag_idx)[]
            push!(block_type, block_type_current)
        elseif tags[tag_idx][2] == fiff_block_end
            push!(block_type, block_type_current)
            block_level[tag_idx:end] .-= 1
        else
            push!(block_type, block_type_current)
        end
        tags[tag_idx][2] == fiff_block_start && (block_number[tag_idx:end] .+= 1)
    end

    return fid, hcat(tag_pos, tag_ids, tag_type, tag_size, block_number, block_level, block_type)
end

function _view_fiff_block(fiff_block::Matrix{Int64})
    for tag_idx in 1:size(fiff_block, 1)
        indent = repeat(" ", (fiff_block[tag_idx, 6] - 1))
        println(indent * "$(fiff_block[tag_idx, 6]) [$(fiff_block[tag_idx, 5])] $(fiff_block[tag_idx, 7])")
    end
end