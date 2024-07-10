export erp_peaks
export amp_at
export avgamp_at
export maxamp_at
export minamp_at

"""
    erp_peaks(obj)

Detect a pair of positive and negative peaks of ERP.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `p::Array{Int64, 2}`: peaks: channels × positive peak position, negative peak position
"""
function erp_peaks(obj::NeuroAnalyzer.NEURO)

    _check_datatype(obj, "erp")

    ch_n = size(obj)[1]
    p = zeros(Int64, ch_n, 2)
    @inbounds for ch_idx in 1:ch_n
        pp_pos = @views maximum(obj.data[ch_idx, :, 1])
        pp_neg = @views minimum(obj.data[ch_idx, :, 1])
        p[ch_idx, :] = @views [vsearch(pp_pos, obj.data[ch_idx, :, 1]), vsearch(pp_neg, obj.data[ch_idx, :, 1])]
    end

    return p

end

"""
    amp_at(obj; <keyword arguments>)

Calculate amplitude at given time.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `t::Real`: time in seconds

# Returns

- `p::Matrix{Float64, 2}`: amplitude for each channel per epoch
"""
function amp_at(obj::NeuroAnalyzer.NEURO; t::Real)

    _check_datatype(obj, ["erp", "mep"])

    if datatype(obj) == "erp"
        @assert t >= obj.epoch_time[1] "t must be ≥ $(obj.epoch_time[1])."
        @assert t <= obj.epoch_time[end] "t must be ≤ $(obj.epoch_time[end])."

        t_idx = vsearch(t, obj.epoch_time)

        ch_n = size(obj)[1]
        ep_n = size(obj)[3]
        p = zeros(ch_n, ep_n)

        @inbounds for ep_idx in 1:ep_n
            Threads.@threads for ch_idx in 1:ch_n
                p[ch_idx, ep_idx] = @views obj.data[ch_idx, t_idx, ep_idx]
            end
        end
    else
        @assert t >= obj.time_pts[1] "t must be ≥ $(obj.time_pts[1])."
        @assert t <= obj.time_pts[end] "t must be ≤ $(obj.time_pts[end])."

        t_idx = vsearch(t, obj.time_pts)

        ch_n = size(obj)[1]
        p = zeros(ch_n)

        Threads.@threads for ch_idx in 1:ch_n
            @inbounds p[ch_idx] = @views obj.data[ch_idx, t_idx, 1]
        end
    end

    return p

end

"""
    avgamp_at(obj; <keyword arguments>)

Calculate average amplitude at given time segment.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `t::Tuple{Real, Real}`: time segment in seconds

# Returns

- `p::Matrix{Float64, 2}`: mean amplitude for each channel per epoch
"""
function avgamp_at(obj::NeuroAnalyzer.NEURO; t::Tuple{Real, Real})

    _check_datatype(obj, ["erp", "mep"])

    if datatype(obj) == "erp"
        @assert t[1] >= obj.epoch_time[1] "t[1] must be ≥ $(obj.epoch_time[1])."
        @assert t[2] <= obj.epoch_time[end] "t[2] must be ≤ $(obj.epoch_time[end])."
        @assert t[1] <= t[2] "t[1] must be < t[2]."

        t_idx1 = vsearch(t[1], obj.epoch_time)
        t_idx2 = vsearch(t[2], obj.epoch_time)

        ch_n = size(obj)[1]
        ep_n = size(obj)[3]
        p = zeros(ch_n, ep_n)

        @inbounds for ep_idx in 1:ep_n
            Threads.@threads for ch_idx in 1:ch_n
                p[ch_idx, ep_idx] = mean(obj.data[ch_idx, t_idx1:t_idx2, ep_idx])
            end
        end
    else
        @assert t[1] >= obj.time_pts[1] "t[1] must be ≥ $(obj.time_pts[1])."
        @assert t[2] <= obj.time_pts[end] "t[2] must be ≤ $(obj.time_pts[end])."
        @assert t[1] <= t[2] "t[1] must be < t[2]."

        t_idx1 = vsearch(t[1], obj.time_pts)
        t_idx2 = vsearch(t[2], obj.time_pts)

        ch_n = size(obj)[1]
        p = zeros(ch_n)

        Threads.@threads for ch_idx in 1:ch_n
            @inbounds p[ch_idx] = mean(obj.data[ch_idx, t_idx1:t_idx2, 1])
        end
    end

    return p

end

"""
    maxamp_at(obj; <keyword arguments>)

Calculate maximum amplitude at given time segment.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `t::Tuple{Real, Real}`: time segment in seconds

# Returns

- `p::Matrix{Float64, 2}`: maximum amplitude for each channel per epoch
"""
function maxamp_at(obj::NeuroAnalyzer.NEURO; t::Tuple{Real, Real})

    _check_datatype(obj, ["erp", "mep"])

    if datatype(obj) == "erp"
        @assert t[1] >= obj.epoch_time[1] "t[1] must be ≥ $(obj.epoch_time[1])."
        @assert t[2] <= obj.epoch_time[end] "t[2] must be ≤ $(obj.epoch_time[end])."
        @assert t[1] <= t[2] "t[1] must be < t[2]."

        t_idx1 = vsearch(t[1], obj.epoch_time)
        t_idx2 = vsearch(t[2], obj.epoch_time)

        ch_n = size(obj)[1]
        ep_n = size(obj)[3]
        p = zeros(ch_n, ep_n)

        @inbounds for ep_idx in 1:ep_n
            Threads.@threads for ch_idx in 1:ch_n
                p[ch_idx, ep_idx] = maximum(obj.data[ch_idx, t_idx1:t_idx2, ep_idx])
            end
        end
    else
        @assert t[1] >= obj.time_pts[1] "t[1] must be ≥ $(obj.time_pts[1])."
        @assert t[2] <= obj.time_pts[end] "t[2] must be ≤ $(obj.time_pts[end])."
        @assert t[1] <= t[2] "t[1] must be < t[2]."

        t_idx1 = vsearch(t[1], obj.time_pts)
        t_idx2 = vsearch(t[2], obj.time_pts)

        ch_n = size(obj)[1]
        p = zeros(ch_n)

        Threads.@threads for ch_idx in 1:ch_n
            @inbounds p[ch_idx] = maximum(obj.data[ch_idx, t_idx1:t_idx2, 1])
        end
    end

    return p

end

"""
    minamp_at(obj; <keyword arguments>)

Calculate minimum amplitude at given time segment.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `t::Tuple{Real, Real}`: time segment in seconds

# Returns

- `p::Matrix{Float64, 2}`: minimum amplitude for each channel per epoch
"""
function minamp_at(obj::NeuroAnalyzer.NEURO; t::Tuple{Real, Real})

    _check_datatype(obj, ["erp", "mep"])

    if datatype(obj) == "erp"
        @assert t[1] >= obj.epoch_time[1] "t[1] must be ≥ $(obj.epoch_time[1])."
        @assert t[2] <= obj.epoch_time[end] "t[2] must be ≤ $(obj.epoch_time[end])."
        @assert t[1] <= t[2] "t[1] must be < t[2]."

        t_idx1 = vsearch(t[1], obj.epoch_time)
        t_idx2 = vsearch(t[2], obj.epoch_time)

        ch_n = size(obj)[1]
        ep_n = size(obj)[3]
        p = zeros(ch_n, ep_n)

        @inbounds for ep_idx in 1:ep_n
            Threads.@threads for ch_idx in 1:ch_n
                p[ch_idx, ep_idx] = minimum(obj.data[ch_idx, t_idx1:t_idx2, 1])
            end
        end
    else
        @assert t[1] >= obj.time_pts[1] "t[1] must be ≥ $(obj.time_pts[1])."
        @assert t[2] <= obj.time_pts[end] "t[2] must be ≤ $(obj.time_pts[end])."
        @assert t[1] <= t[2] "t[1] must be < t[2]."

        t_idx1 = vsearch(t[1], obj.time_pts)
        t_idx2 = vsearch(t[2], obj.time_pts)

        ch_n = size(obj)[1]
        p = zeros(ch_n)

        Threads.@threads for ch_idx in 1:ch_n
            @inbounds p[ch_idx] = minimum(obj.data[ch_idx, t_idx1:t_idx2, 1])
        end
    end

    return p

end
