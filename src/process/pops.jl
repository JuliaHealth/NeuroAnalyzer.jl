export pops
export pops!

"""
    pops(s; repair)

Detect and repair electrode pop (rapid amplitude change). Signal is recovered within the segments starting and ending at zero-crossing. Only one pop is detected, signal length should be ≥2 seconds.

# Arguments

- `s::AbstractVector`
- `r::Int64=20`: detection segment length; pops are checked within `pop_location - r:pop_location + r` samples
- `repair::Bool=true`: recover the segment if `true`

# Returns

Named tuple containing:
- `s_new::Vector{Float64}`
- `pop_location::Int64`: sample number in the signal
- `left_segment::UnitRange{Int64}`: range before the pop that starts when signal crosses 0
- `right_segment::UnitRange{Int64}`: range after the pop that ends when signal crosses 0
"""
function pops(s::AbstractVector; r::Int64=20, repair::Bool=true)

    length(s) < 2 * r + 1 && throw(ArgumentError("s length must be ≥ $(2 * r + 1)."))

    sdiff = diff(s)
    pushfirst!(sdiff, 0)
    sdiff_abs = abs.(sdiff)

    pop_location = 0
    pop_location = vsearch(maximum(sdiff_abs), sdiff_abs)

    if pop_location < r + 1 || pop_location >= length(s) - r
        # NeuroAnalyzer._info("Pop found too close to the segment border, cannot be analyzed. Use different window size.")
        return nothing
    end

    if abs(s[pop_location]) > abs(s[pop_location + 1])
        if abs(s[pop_location + 1]) < abs(s[pop_location + 2])
            pop_location += 1
        end
    end

    s1 = sign.(s[pop_location - r:(pop_location - 1)])
    e1 = extrema(s[pop_location - r:(pop_location - 1)])
    m1 = mean(s[pop_location - r:(pop_location - 1)])
    s2 = sign.(s[pop_location + 1:pop_location + r])
    e2 = extrema(s[pop_location + 1:pop_location + r])
    m2 = mean(s[pop_location + 1:pop_location + r])

    if !allequal(s1)
        return nothing
    elseif !allequal(s2)
        return nothing
    elseif sign(e1[1]) == sign(e2[1])
        return nothing
    elseif sign(e1[2]) == sign(e2[2])
        return nothing
    elseif m1 > mean(s) - 2 * std(s) && m1 < mean(s) + 2 * std(s)
        return nothing
    elseif m2 > mean(s) - 2 * std(s) && m2 < mean(s) + 2 * std(s)
        return nothing
    end

    zero_1 = 1
    zero_2 = length(s)
    for idx in  (pop_location - 1):-1:2
        if sign(s[idx]) != sign(s[idx - 1])
            zero_1 = idx
            break
        end
    end
    for idx in (pop_location + 1):(length(s) - 1)
        if sign(s[idx]) != sign(s[idx + 1])
            zero_2 = idx
            break
        end
    end

    left_seg = zero_1:(pop_location - 1)
    right_seg = (pop_location + 1):zero_2

    if repair
        m1 = mean(s[left_seg])
        l1 = length(s[left_seg])
        ll1 = linspace(0, 1, l1) .* s[pop_location - 1] .* sign(m1)
        m2 = mean(s[right_seg])
        l2 = length(s[right_seg])
        ll2 = linspace(1, 0, l2) .* s[pop_location + 1] .* sign(m2)

        s_new = deepcopy(s)

        if sign(m1) == -1
            s_new[pop_location - l1:pop_location - 1] += ll1
        else
            s_new[pop_location - l1:pop_location - 1] -= ll1
        end
        if sign(m2) == -1
            s_new[pop_location + 1:pop_location + l2] += ll2
        else
            s_new[pop_location + 1:pop_location + l2] -= ll2
        end

        s_new[pop_location] = mean(s_new[pop_location - 1:pop_location + 1])
    
    end

    if repair
        return (s_new=s_new, pop_location=pop_location, left_seg=left_seg[1], right_seg=right_seg[end])
    else
        return (pop_location=pop_location, left_seg=left_seg[1], right_seg=right_seg[end])
    end

end

"""
    pops(s; repair)

Detect and repair electrode pop (rapid amplitude change). Signal is recovered within the segments starting and ending at zero-crossing. Only one pop is detected, signal length should be ≈2 seconds.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `repair::Bool=true`: recover the segment if `true`
- `window::Real=10.0`: window length (in seconds) in which the signal is scanned and repaired (windows are non-overlapping)
- `r::Int64=20`: detection segment length; pops are checked within `pop_location - r:pop_location + r` samples

# Returns

Named tuple containing:
- `obj_new::NeuroAnalyzer.NEURO`
- `pop_location::Vector{Vector{Int64}}`: location of pops: channel, epoch and sample number in the signal
- `left_segment::UnitRange{Int64}`: range before the pop that starts when signal crosses 0
- `right_segment::UnitRange{Int64}`: range after the pop that ends when signal crosses 0
"""
function pops(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), repair::Bool=true, window::Real=10.0, r::Int64=20)

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
            p = @views pops(obj_new.data[ch[ch_idx], Int64.(window_idx:(window_idx + window - 1)), 1], repair=repair, r=r)
            if p !== nothing
                if repair
                    s[ch_idx, Int64.(window_idx:(window_idx + window - 1)), 1] = p.s_new
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
    pops!(s; repair)

Detect and repair electrode pop (rapid amplitude change). Signal is recovered within the segments starting and ending at zero-crossing. Only one pop is detected, signal length should be ≈2 seconds.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `repair::Bool=true`: recover the segment if `true`
- `window::Real=20.0`: window length (in seconds) in which the signal is scanned and repaired (windows are non-overlapping)
- `r::Int64=20`: detection segment length; pops are checked within `pop_location - r:pop_location + r` samples

# Returns

Named tuple containing:
- `pop_location::Vector{Vector{Int64}}`: location of pops: channel, epoch and sample number in the signal
- `left_segment::UnitRange{Int64}`: range before the pop that starts when signal crosses 0
- `right_segment::UnitRange{Int64}`: range after the pop that ends when signal crosses 0
"""
function pops!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), repair::Bool=true, window::Real=10.0, r::Int64=20)

    if repair
        obj_new, pop_loc, l_seg, r_seg = pops(obj, ch=ch, repair=true, window=window, r=r)
        obj.data = obj_new.data
        obj.components = obj_new.components
        obj.history = obj_new.history
        return obj_new, pop_loc, l_seg, r_seg
    else
        pop_loc, l_seg, r_seg = pops(obj, ch=ch, repair=false, window=window, r=r)
        return pop_loc, l_seg, r_seg
    end

    return nothing

end