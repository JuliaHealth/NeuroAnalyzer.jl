##################################
#                                #
#  Low-level internal functions  #
#                                #
##################################

_reflect(signal::AbstractArray) = vcat(signal[end:-1:1], signal, signal[end:-1:1])

_chop(signal::AbstractArray) = signal[(length(signal) ÷ 3 + 1):(length(signal) ÷ 3) * 2]

_xlims(t::Vector{<:Real}) = (floor(t[1], digits=2), ceil(t[end], digits=2))

_xticks(t::Vector{<:Real}) = floor(t[1], digits=2):((ceil(t[end]) - floor(t[1])) / 10):ceil(t[end], digits=2)

_pl(x) = ((length(collect(x)) > 1) && return "s") || return ""

function _check_channels(eeg::NeuroAnalyzer.EEG, channel::Union{Int64, Vector{Int64}, AbstractRange})
    for idx in 1:length(channel)
        (channel[idx] < 1 || channel[idx] > eeg_channel_n(eeg)) && throw(ArgumentError("channel must be ≥ 1 and ≤ $(eeg_channel_n(eeg))."))
    end
end

function _check_channels(channels::Union{Int64, Vector{Int64}, AbstractRange}, channel::Union{Int64, Vector{Int64}, AbstractRange})
    for idx in 1:length(channel)
        channel[idx] in channels || throw(ArgumentError("Channel $(channel[idx]) does not match signal channels."))
        (channel[idx] < 1 || channel[idx] > length(channels)) && throw(ArgumentError("channel must be ≥ 1 and ≤ $(eeg_channel_n(eeg))."))
    end
end

function _check_epochs(eeg::NeuroAnalyzer.EEG, epoch::Union{Int64, Vector{Int64}, AbstractRange})
    for idx in 1:length(epoch)
        (epoch[idx] < 1 || epoch[idx] > eeg_epoch_n(eeg)) && throw(ArgumentError("epoch must be ≥ 1 and ≤ $(eeg_epoch_n(eeg))."))
    end
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
    c, _ = _get_component(eeg, c)
    # select channels, default is all or def_cidx
    def_cidx > size(c, 1) && (def_cidx = size(c, 1))
    def_cidx == 0 && (def_cidx = size(c, 1))
    channel == 0 && (channel = 1:def_cidx)
    typeof(c_idx) <: AbstractRange && (channel = collect(channel))
    length(c_idx) > 1 && sort!(c_idx)
    for idx in 1:length(c_idx)
        (c_idx[idx] < 1 || c_idx[idx] > size(c, 1)) && throw(ArgumentError("c_idx must be ≥ 1 and ≤ $(size(c, 1))."))
    end
    return c_idx
end

function _get_component(eeg::NeuroAnalyzer.EEG, c::Symbol)
    c in eeg.eeg_header[:components] || throw(ArgumentError("Component $c not found."))
    c_idx = findfirst(isequal(c), eeg.eeg_header[:components])
    c = eeg.eeg_components[c_idx]
    return (c=c, c_idx=c_idx)
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

function _draw_head(p::Plots.Plot{Plots.GRBackend}; head_labels::Bool=true, topo::Bool=false, kwargs...)
    # Draw head over a topographical plot `p`.
    # - `p::Plots.Plot{Plots.GRBackend}`: electrodes plot
    # - `loc_x::Vector{Float64}`: vector of x electrode position
    # - `loc_y::Vector{Float64}`: vector of y electrode position
    # - `head_labels::Bool=true`: add text labels to the plot
    # - `topo::Bool=false`: if true, perform peripheral erasing for topo plots
    # - `kwargs`: optional arguments for plot() function
    # loc_x, loc_y = loc_y, loc_x
    pts = Plots.partialcircle(0, 2π, 100, 1.1)
    x, y = Plots.unzip(pts)
    head = Plots.Shape(x, y)
    nose = Plots.Shape([(-0.05, maximum(y)), (0, maximum(y) + 0.1 * maximum(y)), (0.05, maximum(y))])
    ear_l = Plots.Shape([(minimum(x), -0.2), (minimum(x) + 0.05 * minimum(x), -0.2), (minimum(x) + 0.05 * minimum(x), 0.2), (minimum(x), 0.2)])
    ear_r = Plots.Shape([(maximum(x), -0.2), (maximum(x) + 0.05 * maximum(x), -0.2), (maximum(x) + 0.05 * maximum(x), 0.2), (maximum(x), 0.2)])
    p = Plots.plot!(p, head, fill=nothing, label="")
    p = Plots.plot!(nose, fill=nothing, label="")
    p = Plots.plot!(ear_l, fill=nothing, label="")
    p = Plots.plot!(ear_r, fill=nothing, label="")

    if head_labels == true
        p = Plots.plot!(annotation=(0, 1.05, Plots.text("IN", pointsize=4, halign=:center, valign=:center)))
        p = Plots.plot!(annotation=(0, -1.05, Plots.text("NAS", pointsize=4, halign=:center, valign=:center)))
        p = Plots.plot!(annotation=(-1.05, 0, Plots.text("LPA", pointsize=4, halign=:center, valign=:center, rotation=90)))
        p = Plots.plot!(annotation=(1.05, 0, Plots.text("RPA", pointsize=4, halign=:center, valign=:center, rotation=-90)))
    end

    if topo == true
        pts = Plots.partialcircle(0, 2π, 100, 1.3)
        x, y = Plots.unzip(pts)
        for idx in 1:0.001:1.7
            peripheral = Shape(x .* idx, y .* idx)
            p = Plots.plot!(p, peripheral, label="", fill=nothing, lc=:white)
        end
        p = Plots.plot!(xlims=(-1.4, 1.4), ylims=(-1.4, 1.4); kwargs...)
    end

    return p
end

function _check_epochs(eeg::NeuroAnalyzer.EEG, epoch)
    epoch[1] < 1 || epoch[end] > eeg_epoch_n(eeg) && throw(ArgumentError("epoch must be ≥ 1 and ≤ $(eeg_epoch_n(eeg))."))
    for idx in 1:length(epoch)
        (epoch[idx] < 1 || epoch[idx] > eeg_epoch_n(eeg)) && throw(ArgumentError("epoch must be ≥ 1 and ≤ $(eeg_epoch_n(eeg))."))
    end
end

function _get_epoch_markers(eeg::NeuroAnalyzer.EEG, offset, len)
    # get epochs markers for len > epoch_len
    epoch_markers = Vector{Int64}[]
    if len + offset > eeg_epoch_len(eeg) && eeg_epoch_n(eeg) > 1
        eeg_tmp = eeg_epochs(eeg, epoch_n=1)
        epoch_len = eeg_epoch_len(eeg)
        epoch_n = eeg_epoch_n(eeg)
        epoch_markers = collect(1:epoch_len:epoch_len * epoch_n)[2:end] 
        epoch_markers = floor.(Int64, (epoch_markers ./ eeg_sr(eeg)))
        epoch_markers = epoch_markers[epoch_markers .> floor(Int64, offset / eeg_sr(eeg))]
        epoch_markers = epoch_markers[epoch_markers .<= ceil(Int64, (offset + len) / eeg_sr(eeg))]
    else
        eeg_tmp = eeg
    end
    return eeg_tmp, epoch_markers
end

function _get_t(eeg::NeuroAnalyzer.EEG, offset, len)
    t = collect(0:(1 / eeg_sr(eeg)):(len / eeg_sr(eeg)))
    t = t .+ (offset / eeg_sr(eeg))
    t = t[1:(end - 1)]
    t[1] = floor(t[1], digits=2)
    t[end] = ceil(t[end], digits=2)
    return t
end

function _check_offset_len(eeg::NeuroAnalyzer.EEG, offset, len)
    (offset < 0 || offset > eeg_epoch_len(eeg)) && throw(ArgumentError("offset must be > 0 and ≤ $(eeg_epoch_len(eeg))."))
    (offset + len > eeg_epoch_len(eeg)) && throw(ArgumentError("offset + len must be ≤ $(eeg_epoch_len(eeg))."))
end

function _convert_t(t)
    t_1 = floor(t[1], digits=2)
    t_2 = ceil(t[end], digits=2)
    abs(t_1) < 1.0 && (t_s1 = string(floor(t_1 * 1000, digits=2)) * " ms")
    abs(t_1) >= 1.0 && (t_s1 = string(floor(t_1, digits=2)) * " s")
    abs(t_2) < 1.0 && (t_s2 = string(ceil(t_2 * 1000, digits=2)) * " ms")
    abs(t_2) >= 1.0 && (t_s2 = string(ceil(t_2, digits=2)) * " s")
    return t_1, t_s1, t_2, t_s2
end

function _check_channels(eeg::NeuroAnalyzer.EEG, channel)
    for idx in 1:length(channel)
        (channel[idx] < 1 || channel[idx] > eeg_channel_n(eeg)) && throw(ArgumentError("channel must be ≥ 1 and ≤ $(eeg_channel_n(eeg))."))
    end
end

function _check_cidx(eeg::NeuroAnalyzer.EEG, c::Symbol, c_idx)
    c, _ = _get_component(eeg, c)
    for idx in 1:length(c_idx)
        (c_idx[idx] < 1 || c_idx[idx] > size(c, 1)) && throw(ArgumentError("c_idx must be ≥ 1 and ≤ $(size(c, 1))."))
    end
end

function _tuple_max(t::Tuple{Real, Real})
    abs(t[1]) > abs(t[2]) && (t = (-abs(t[1]), abs(t[1])))
    abs(t[1]) < abs(t[2]) && (t = (-abs(t[2]), abs(t[2])))
    return t
end

function _channel2channel_name(channel)
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
    return channel_name
end

function _t2epoch(eeg::NeuroAnalyzer.EEG, offset, len, epoch_tmp)
    if (1 + offset) > eeg_epoch_len(eeg)
        if (floor(Int64, (1 + offset) / eeg_epoch_len(eeg)) + 1) < (ceil(Int64, (1 + offset + len) / eeg_epoch_len(eeg)) - 1)
            (epoch_tmp = (floor(Int64, (1 + offset) / eeg_epoch_len(eeg)) + 1):(ceil(Int64, (1 + offset + len) / eeg_epoch_len(eeg)) - 1))
        else
            (epoch_tmp = (floor(Int64, (1 + offset) / eeg_epoch_len(eeg)) + 1):(floor(Int64, (1 + offset) / eeg_epoch_len(eeg)) + 1))
        end
    end
    return epoch_tmp
end

function _fir_response(f::Vector{Float64}, w=range(0, stop=π, length=1024))
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

function _make_epochs(signal::Matrix{Float64}; epoch_n::Union{Int64, Nothing}=nothing, epoch_len::Union{Int64, Nothing}=nothing, average::Bool=false)
    (epoch_len === nothing && epoch_n === nothing) && throw(ArgumentError("Either epoch_n or epoch_len must be set."))
    (epoch_len !== nothing && epoch_n !== nothing) && throw(ArgumentError("Both epoch_n and epoch_len cannot be set."))
    (epoch_len !== nothing && epoch_len < 1) && throw(ArgumentError("epoch_len must be ≥ 1."))
    (epoch_n !== nothing && epoch_n < 1) && throw(ArgumentError("epoch_n must be ≥ 1."))

    channel_n, _ = size(signal)

    if epoch_n === nothing
        epoch_n = size(signal, 2) ÷ epoch_len
    else
        epoch_len = size(signal, 2) ÷ epoch_n
    end

    epochs = zeros(channel_n, epoch_len, epoch_n)

    idx1 = 1
    for idx2 in 1:epoch_len:(epoch_n * epoch_len - 1)
        epochs[:, :, idx1] = signal[:, idx2:(idx2 + epoch_len - 1), 1]
        idx1 += 1
    end

    average == true && (epochs = mean(epochs, dims=3)[:, :])

    return epochs
end

function _get_channel_idx(labels::Vector{String}, channel::Union{String, Int64})
    if typeof(channel) == String
        channel_found = nothing
        for idx in 1:length(labels)
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

function _free_gpumem(threshold::Real=0.95)
    m = CUDA.MemoryInfo()
    usedmem = m.total_bytes - m.free_bytes
    totalmem = m.total_bytes
    if usedmem / totalmem > threshold
        # CUDA.reclaim()
        GC.gc(true)
    end
end 

function _a2df(annotations::Vector{String})
    # convert EDF/BDF annotations to DataFrame
    annotations = replace.(annotations, "\x14\x14\0" => "|")
    annotations = replace.(annotations, "\x14\x14" => "|")
    annotations = replace.(annotations, "\x14" => "|")
    annotations = replace.(annotations, "\0" => "")
    a_onset = Vector{Float64}()
    a_event = Vector{String}()
    # what about annotations containing event duration?
    for ann_idx in 1:length(annotations)
        s = split(annotations[ann_idx], "|")
        if length(s) > 2
            push!(a_onset, parse(Float64, strip(s[2])))
            push!(a_event, strip(s[3]))
        end
    end
    return DataFrame(onset=a_onset, event=a_event)
end

function _clean_labels(labels::Vector{String})
    labels = replace.(labels, "EEG " => "")
    labels = replace.(labels, "EDF " => "")
    labels = replace.(labels, "BDF " => "")
    return labels
end

function _has_annotations(channel_type::Vector{String})
    if "annotations" in channel_type
        annotations = true
        for channel_idx in 1:length(channel_type)
            channel_type[channel_idx] == "annotations" && (annotations_channel = channel_idx)
        end
    else
        annotations = false
        annotations_channel = 0
    end
    return annotations, annotations_channel
end

function _set_channel_types(labels::Vector{String})
    eeg_channels = ["af3", "af4", "af7", "af8", "afz", "c1", "c2", "c3", "c4", "c5", "c6", "cp1", "cp2", "cp3", "cp4", "cp5", "cp6", "cpz", "cz", "f1", "f10", "f2", "f3", "f4", "f5", "f6", "f7", "f8", "f9", "fc1", "fc2", "fc3", "fc4", "fc5", "fc6", "fcz", "fp1", "fp2", "fpz", "ft10", "ft7", "ft8", "ft9", "fz", "nz", "o1", "o2", "oz", "p1", "p10", "p2", "p3", "p4", "p5", "p6", "p7", "p8", "p9", "po3", "po4", "po7", "po8", "poz", "pz", "t10", "t3", "t4", "t5", "t6", "t7", "t8", "t9", "tp10", "tp7", "tp8", "tp9"]
    ref_channels = ["a1", "a2", "m1", "m2", "pg1", "pg2"]
    eog_channel = ["e", "e1", "e2"]
    channel_type = repeat(["unknown"], length(labels))
    for idx in 1:length(labels)
        in(lowercase(labels[idx]), eeg_channels) && (channel_type[idx] = "eeg")
        for idx2 in 1:length(eeg_channels)
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
        occursin("meg", lowercase(labels[idx])) && (channel_type[idx] = "meg")
        occursin("ecg", lowercase(labels[idx])) && (channel_type[idx] = "ecg")
        occursin("ekg", lowercase(labels[idx])) && (channel_type[idx] = "ecg")
        occursin("eog", lowercase(labels[idx])) && (channel_type[idx] = "eog")
        for idx2 in 1:length(eog_channel)
            occursin(eeg_channels[idx2], lowercase(labels[idx])) && (channel_type[idx] = "eog")
        end
        occursin("emg", lowercase(labels[idx])) && (channel_type[idx] = "emg")
        occursin("temp", lowercase(labels[idx])) && (channel_type[idx] = "temp")
        in(lowercase(labels[idx]), ref_channels) && (channel_type[idx] = "ref")
        for idx2 in 1:length(ref_channels)
            occursin(ref_channels[idx2], lowercase(labels[idx])) && (channel_type[idx] = "ref")
        end
        occursin("mark", lowercase(labels[idx])) && (channel_type[idx] = "markers")
        occursin("event", lowercase(labels[idx])) && (channel_type[idx] = "events")
        occursin("annotations", lowercase(labels[idx])) && (channel_type[idx] = "annotations")
        occursin("status", lowercase(labels[idx])) && (channel_type[idx] = "annotations")
    end
    return channel_type
end