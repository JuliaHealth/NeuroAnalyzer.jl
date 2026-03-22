export erp_peaks
export amp_at
export avgamp_at
export maxamp_at
export minamp_at
export erp_auc

"""
    erp_peaks(obj)

Detect the positive and negative peak of each channel's ERP/ERF/MEP.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object

# Returns

- `p::Matrix{Int64}`: : shape `(channels, 2)` — column 1 is the positive peak sample index, column 2 is the negative peak sample index
"""
function erp_peaks(obj::NeuroAnalyzer.NEURO)::Matrix{Int64}

    _check_datatype(obj, ["erp", "erf", "mep"])

    # number of channels
    ch_n = size(obj, 1)

    #  pre-allocate output
    p = zeros(ch_n, ep_n)

    @inbounds Threads.@threads :dynamic for ch_idx in 1:ch_n
        s = @view obj.data[ch_idx, :, 1]
        # positive peak: sample index of maximum
        p[ch_idx, 1] = argmax(s)
        # negative peak: sample index of minimum
        p[ch_idx, 2] = argmin(s)
    end

    return p

end

"""
    amp_at(obj; <keyword arguments>)

Calculate amplitude at a given time point.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `t::Real`: time in seconds

# Returns

- `p::Matrix{Float64}`: amplitude for each channel per epoch, shape `(channels, epochs)` or `(channels, 1)` for continuous objects
"""
function amp_at(obj::NeuroAnalyzer.NEURO; t::Real)::Matrix{Float64}

    if datatype(obj) in ["erp", "erf", "mep"]

        !(t >= obj.epoch_time[1]) && throw(ArgumentError("t must be ≥ $(obj.epoch_time[1])."))
        !(t <= obj.epoch_time[end]) && throw(ArgumentError("t must be ≤ $(obj.epoch_time[end])."))

        t_idx = vsearch(t, obj.epoch_time)

        # number of channels
        ch_n = size(obj, 1)
        # number of epochs
        ep_n = size(obj, 3)

        #  pre-allocate output
        p = zeros(ch_n, ep_n)

        @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
            ch_idx, ep_idx = idx[1], idx[2]
            p[ch_idx, ep_idx] = obj.data[ch_idx, t_idx, ep_idx]
        end

    else

        !(t >= obj.time_pts[1]) && throw(ArgumentError("t must be ≥ $(obj.time_pts[1])."))
        !(t <= obj.time_pts[end]) && throw(ArgumentError("t must be ≤ $(obj.time_pts[end])."))

        t_idx = vsearch(t, obj.time_pts)

        # number of channels
        ch_n = size(obj, 1)

        # pre-allocate output
        p = zeros(ch_n, 1)

        @inbounds Threads.@threads :dynamic for ch_idx in 1:ch_n
            p[ch_idx, 1] = obj.data[ch_idx, t_idx, 1]
        end

    end

    return p

end

"""
    avgamp_at(obj; <keyword arguments>)

Calculate mean amplitude over a time segment.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `t::Tuple{Real, Real}`: time segment in seconds

# Returns

- `p::Matrix{Float64}`: mean amplitude for each channel per epoch, shape `(channels, epochs)` or `(channels, 1)` for continuous objects
"""
function avgamp_at(obj::NeuroAnalyzer.NEURO; t::Tuple{Real, Real})::Matrix{Float64}

    if datatype(obj) in ["erp", "erf", "mep"]

        _check_tuple(t, (obj.epoch_time[1], obj.epoch_time[end]), "seg")

        t_idx1 = vsearch(t[1], obj.epoch_time)
        t_idx2 = vsearch(t[2], obj.epoch_time)

        # number of channels
        ch_n = size(obj, 1)
        # number of epochs
        ep_n = size(obj, 3)

        # pre-allocate output
        p = zeros(ch_n, ep_n)

        @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
            ch_idx, ep_idx = idx[1], idx[2]
            p[ch_idx, ep_idx] = mean(@view(obj.data[ch_idx, t_idx1:t_idx2, ep_idx]))
        end

    else

        _check_tuple(t, (obj.time_pts[1], obj.time_pts[end]), "seg")

        t_idx1 = vsearch(t[1], obj.time_pts)
        t_idx2 = vsearch(t[2], obj.time_pts)

        # number of channels
        ch_n = size(obj, 1)

        # pre-allocate output
        p = zeros(ch_n, 1)

        @inbounds Threads.@threads :dynamic for ch_idx in 1:ch_n
            p[ch_idx, 1] = mean(@view(obj.data[ch_idx, t_idx1:t_idx2, 1]))
        end

    end

    return p

end

"""
    maxamp_at(obj; <keyword arguments>)

Calculate maximum amplitude over a time segment.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `t::Tuple{Real, Real}`: time segment in seconds

# Returns

- `p::Matrix{Float64}`: maximum amplitude for each channel per epoch, shape `(channels, epochs)` or `(channels, 1)` for continuous objects
"""
function maxamp_at(obj::NeuroAnalyzer.NEURO; t::Tuple{Real, Real})::Matrix{Float64}

    if datatype(obj) in ["erp", "erf", "mep"]

        _check_tuple(t, (obj.epoch_time[1], obj.epoch_time[end]), "seg")

        t_idx1 = vsearch(t[1], obj.epoch_time)
        t_idx2 = vsearch(t[2], obj.epoch_time)

        # number of channels
        ch_n = size(obj, 1)
        # number of epochs
        ep_n = size(obj, 3)

        # pre-allocate output
        p = zeros(ch_n, ep_n)

        @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
            ch_idx, ep_idx = idx[1], idx[2]
            p[ch_idx, ep_idx] = maximum(@view(obj.data[ch_idx, t_idx1:t_idx2, ep_idx]))
        end

    else

        _check_tuple(t, (obj.time_pts[1], obj.time_pts[end]), "seg")

        t_idx1 = vsearch(t[1], obj.time_pts)
        t_idx2 = vsearch(t[2], obj.time_pts)

        # number of channels
        ch_n = size(obj, 1)

        # pre-allocate output
        p = zeros(ch_n, 1)

        @inbounds Threads.@threads :dynamic for ch_idx in 1:ch_n
            p[ch_idx, 1] = maximum(@view(obj.data[ch_idx, t_idx1:t_idx2, 1]))
        end

    end

    return p

end

"""
    minamp_at(obj; <keyword arguments>)

Calculate minimum amplitude over a time segment.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `t::Tuple{Real, Real}`: time segment in seconds

# Returns

- `p::Matrix{Float64}`: minimum amplitude for each channel per epoch, shape `(channels, epochs)` or `(channels, 1)` for continuous objects
"""
function minamp_at(obj::NeuroAnalyzer.NEURO; t::Tuple{Real, Real})::Matrix{Float64}

    if datatype(obj) in ["erp", "erf", "mep"]

        _check_tuple(t, (obj.epoch_time[1], obj.epoch_time[end]), "seg")

        t_idx1 = vsearch(t[1], obj.epoch_time)
        t_idx2 = vsearch(t[2], obj.epoch_time)

        # number of channels
        ch_n = size(obj, 1)
        # number of epochs
        ep_n = size(obj, 3)

        # pre-allocate output
        p = zeros(ch_n, ep_n)

        @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
            ch_idx, ep_idx = idx[1], idx[2]
            p[ch_idx, ep_idx] = minimum(@view(obj.data[ch_idx, t_idx1:t_idx2, ep_idx]))
        end

    else

        _check_tuple(t, (obj.time_pts[1], obj.time_pts[end]), "seg")

        t_idx1 = vsearch(t[1], obj.time_pts)
        t_idx2 = vsearch(t[2], obj.time_pts)

        # number of channels
        ch_n = size(obj, 1)

        # pre-allocate output
        p = zeros(ch_n, 1)

        @inbounds Threads.@threads :dynamic for ch_idx in 1:ch_n
            p[ch_idx, 1] = minimum(@view(obj.data[ch_idx, t_idx1:t_idx2, 1]))
        end

    end

    return p

end

"""
    erp_auc(obj; <keyword arguments>)

Compute area under curve of an ERP/ERF/MEP (epoch 1).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `seg::Tuple{Real, Real}=(obj.epoch_time[1], obj.epoch_time[end])`: time segment in seconds
- `type::Symbol=:all`: which part of the signal to integrate:
    - `:all` - full waveform
    - `:pos` - positive values only
    - `:neg` - negative values only

# Returns

- `auc::Vector{Float64}`
"""
function erp_auc(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    seg::Tuple{Real, Real} = (obj.epoch_time[1], obj.epoch_time[end]),
    type::Symbol = :all
)::Vector{Float64}

    _check_datatype(obj, ["erp", "erf", "mep"])
    _check_var(type, [:all, :pos, :neg], "type")
    _check_tuple(seg, (obj.epoch_time[1], obj.epoch_time[end]), "seg")
    t1 = vsearch(seg[1], obj.epoch_time)
    t2 = vsearch(seg[2], obj.epoch_time)
    ch = get_channel(obj, ch = ch)

    auc = zeros(length(ch))

    # time vector and resolution for the selected segment.
    t = obj.epoch_time[t1:t2]
    dx = t[2] - t[1]

    @inbounds for ch_idx in eachindex(ch)

        # extract the signal segment for this channel (epoch 1 = ERP/ERF average).
        s = @view obj.data[ch_idx, t1:t2, 1]

        if type === :all

            auc[ch_idx] = simpson(s, t, dx = dx)

        elseif type === :pos

            mask = s .> 0
            !(any(mask)) && throw(ArgumentError("No positive values in channel $(ch[ch_idx]) segment, cannot compute AUC."))
            auc[ch_idx] = simpson(s[mask], t[mask], dx = dx)

        elseif type === :neg

            mask = s .< 0
            !(any(mask)) && throw(ArgumentError("No negative values in channel $(ch[ch_idx]) segment, cannot compute AUC."))
            auc[ch_idx] = simpson(s[mask], t[mask], dx = dx)

        end

    end

    return auc

end
