export remove_pops
export remove_pops!

"""
    remove_pops(s; r, repair)

Detect and repair electrode pop (rapid amplitude change). Signal is recovered within the segments starting and ending at zero-crossing. Only one pop is detected, signal length should be ≥2 seconds.

# Arguments

- `s::AbstractVector`
- `r::Int64=20`: detection segment length; pops are checked within `pop_location - r:pop_location + r` samples
- `repair::Bool=true`: recover the segment if `true`

# Returns

Named tuple containing:
- `s::Vector{Float64}`
- `pop_location::Int64`: sample number in the signal
- `left_segment::Int64`: length of segment before the pop that starts when signal crosses 0
- `right_segment::Int64`: length of segment after the pop that ends when signal crosses 0
"""
function remove_pops(s::AbstractVector; r::Int64=20, repair::Bool=true)

    length(s) < 2 * r + 1 && throw(ArgumentError("s length must be ≥ $(2 * r + 1)."))

    sdiff = diff(s)
    pushfirst!(sdiff, 0)
    sdiff_abs = abs.(sdiff)

    pop_location = 0
    pop_location = vsearch(maximum(sdiff_abs), sdiff_abs)

    # pop found too close to the segment border, cannot be analyzed
    (pop_location < r + 1 || pop_location >= length(s) - r) && return nothing

    if abs(s[pop_location]) > abs(s[pop_location + 1])
        if abs(s[pop_location + 1]) < abs(s[pop_location + 2])
            pop_location += 1
        end
    end

    # check if amplitude of the pop segment is significantly larger than the rest of the signal
    s_pre = s[1:pop_location - r]
    s_post = s[pop_location + r:end]
    s_pop = s[(pop_location - r):(pop_location + r)]
    s_pop_ci = (mean(s_pop) - 2 * std(s_pop), mean(s_pop) + 2 * std(s_pop))
    s_ci = extrema([mean(s_pre) - 2 * std(s_pre), mean(s_post) - 2 * std(s_post), mean(s_pre) + 2 * std(s_pre), mean(s_post) + 2 * std(s_post)])
    # if not, this is ignore this pop
    s_pop_ci[1] > s_ci[1] || s_pop_ci[2] < s_ci[2] && return nothing

    # there is no pop if it does not cross Y at 0
    sign(s[pop_location - 5]) == sign(s[pop_location + 5]) && return nothing

    zero_1 = 1
    zero_2 = length(s)
    for idx in  (pop_location - 5):-1:2
        if sign(s[idx]) != sign(s[idx - 1])
            zero_1 = idx
            break
        end
    end
    for idx in (pop_location + 5):(length(s) - 1)
        if sign(s[idx]) != sign(s[idx + 1])
            zero_2 = idx
            break
        end
    end

    left_seg = zero_1:(pop_location - 1)
    right_seg = (pop_location + 1):zero_2

    if repair
        s_pop = s[zero_1:zero_2]
        s_pop_min = vsearch(minimum(s_pop), s_pop)
        s_pop_max = vsearch(maximum(s_pop), s_pop)
    
        if s_pop_max < s_pop_min
            # /|
            #  |/
            m1 = mean(s_pop[1:s_pop_max])
            l1 = length(s_pop[1:s_pop_max])
            ll1 = linspace(0, 1, l1) .* (s_pop[s_pop_max] * 0.9) .* sign(m1)
            m2 = mean(s_pop[s_pop_min:end])
            l2 = length(s_pop[s_pop_min:end])
            ll2 = linspace(1, 0, l2) .* (s_pop[s_pop_min] * 0.9) .* sign(m2)
            s_pop[1:s_pop_max] -= ll1
            s_pop[s_pop_min:end] += ll2
            s_pop[s_pop_max + 1:s_pop_min - 1] = normalize_minmax(s_pop[s_pop_max + 1:s_pop_min - 1])
            s_pop[s_pop_max + 1:s_pop_min - 1][s_pop[s_pop_max + 1:s_pop_min - 1] .> 0] .*= abs(maximum(s_pop[1:s_pop_max]))
            s_pop[s_pop_max + 1:s_pop_min - 1][s_pop[s_pop_max + 1:s_pop_min - 1] .< 0] .*= abs(minimum(s_pop[s_pop_min:end]))
            s[zero_1:zero_2] = s_pop
        else
            #  |\
            # \|
            m1 = mean(s_pop[1:s_pop_min])
            l1 = length(s_pop[1:s_pop_min])
            ll1 = linspace(0, 1, l1) .* s_pop[s_pop_min] .* sign(m1)
            m2 = mean(s_pop[s_pop_max:end])
            l2 = length(s_pop[s_pop_max:end])
            ll2 = linspace(1, 0, l2) .* s_pop[s_pop_max] .* sign(m2)
            s_pop[1:s_pop_min] += ll1
            s_pop[s_pop_max:end] -= ll2
            s_pop[s_pop_min:s_pop_max] = normalize_minmax(s_pop[s_pop_min:s_pop_max])
            s_pop[s_pop_min:s_pop_max][s_pop[s_pop_min:s_pop_max] .> 0] .*= maximum(s[1:s_pop_max])
            s_pop[s_pop_min:s_pop_max][s_pop[s_pop_min:s_pop_max] .< 0] .*= minimum(s[1:s_pop_min])
            s[zero_1:zero_2] = s_pop
        end

    end

    if repair
        return (s=s, pop_location=pop_location, left_seg=left_seg[1], right_seg=right_seg[end])
    else
        return (pop_location=pop_location, left_seg=left_seg[1], right_seg=right_seg[end])
    end

end

"""
    remove_pops(obj; <keyword arguments>)

Detect and repair electrode pop (rapid amplitude change). Signal is recovered within the segments starting and ending at zero-crossing. Only one pop is detected per segment, signal length should be ≈2 seconds.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `repair::Bool=true`: recover the segment if `true`
- `window::Real=10.0`: window length (in seconds) in which the signal is scanned and repaired (windows are non-overlapping)
- `r::Int64=sr(obj)÷2`: detection segment length; pops are checked within `(pop_location - r):(pop_location + r)` samples

# Returns

- `obj_new::NeuroAnalyzer.NEURO`: returned if `repair=true`
- `pop_loc::Vector{Vector{Int64}}`: location of pops: channel, epoch and sample number in the signal
- `l_seg::Vector{Int64}`: length of segment before the pop that starts when signal crosses 0
- `r_seg::Vector{Int64}`: length of segment after the pop that ends when signal crosses 0
"""
function remove_pops(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), repair::Bool=true, window::Real=10.0, r::Int64=sr(obj)÷2)

    _check_channels(obj, ch)
    epoch_n(obj) > 1 && throw(ArgumentError("pop() should be applied to continuous (non-epoched) signal."))

    obj_new = deepcopy(obj)

    if length(ch) == 1
        s = @views reshape(obj_new.data[ch, :, 1], 1, :, size(obj_new.data[ch, :, 1], 2))
    else
        s = @views obj_new.data[ch, :, :]
    end

    pop_loc = Vector{Vector{Int64}}()
    l_seg = Vector{Int64}()
    r_seg = Vector{Int64}()

    window *= sr(obj)
    window > signal_len(obj) && throw(ArgumentError("window must be ≤ $(signal_len(obj) / sr(obj))."))
    ch_n = size(s, 1)

    @inbounds @simd for ch_idx in 1:ch_n
        for window_idx in Int64.(1:window:(signal_len(obj) - signal_len(obj) % window))
            p = @views remove_pops(obj_new.data[ch[ch_idx], Int64.(window_idx:(window_idx + window - 1)), 1], repair=repair, r=r)
            if p !== nothing
                if repair
                    s[ch_idx, Int64.(window_idx:(window_idx + window - 1)), 1] = p.s
                end
                push!(pop_loc, [ch[ch_idx], window_idx + p.pop_location])
                push!(l_seg, p.pop_location - p.left_seg)
                push!(r_seg, p.right_seg - p.pop_location)
            end
        end
    end

    if repair
        obj_new.data[ch, :, :] = s
        reset_components!(obj_new)
        push!(obj_new.history, "pops(OBJ, ch=$ch, repair=true, window=$window)")
        return obj_new, pop_loc, l_seg, r_seg
    else
        return pop_loc, l_seg, r_seg
    end

end

"""
    remove_pops!(obj; <keyword arguments>)

Detect and repair electrode pop (rapid amplitude change). Signal is recovered within the segments starting and ending at zero-crossing. Only one pop is detected, signal length should be ≈2 seconds.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `repair::Bool=true`: recover the segment if `true`
- `window::Real=20.0`: window length (in seconds) in which the signal is scanned and repaired (windows are non-overlapping)
- `r::Int64=sr(obj)÷2`: detection segment length; pops are checked within `pop_location - r:pop_location + r` samples

# Returns

- `pop_loc::Vector{Vector{Int64}}`: location of pops: channel, epoch and sample number in the signal
- `l_seg::Vector{Int64}`: length of segment before the pop that starts when signal crosses 0
- `r_seg::Vector{Int64}`: length of segment after the pop that ends when signal crosses 0
"""
function remove_pops!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), repair::Bool=true, window::Real=10.0, r::Int64=sr(obj)÷2)

    if repair
        obj_new, pop_loc, l_seg, r_seg = remove_pops(obj, ch=ch, repair=true, window=window, r=r)
        obj.data = obj_new.data
        obj.components = obj_new.components
        obj.history = obj_new.history
        return obj_new, pop_loc, l_seg, r_seg
    else
        pop_loc, l_seg, r_seg = remove_pops(obj, ch=ch, repair=false, window=window, r=r)
        return pop_loc, l_seg, r_seg
    end

end