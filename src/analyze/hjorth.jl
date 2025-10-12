export hjorth

_h_mob(s::AbstractVector)::Float64 = sqrt(var(derivative(s)) / var(s))

"""
    hjorth(s)

Calculate Hjorths parameters.

# Arguments

- `s::AbstractVector`

# Returns

Named tuple containing:
- `h_act::Float64`: activity
- `h_mob::Float64`: mobility
- `h_comp::Float64`: complexity

# Notes:

- Activity: the total power of the signal
- Mobility: an estimate of the mean frequency
- Complexity: indicates the similarity of the shape of the signal to a pure sine wave
"""
function hjorth(s::AbstractVector)::@NamedTuple{h_act::Float64, h_mob::Float64, h_comp::Float64}

    h_act = var(s)
    h_mob = _h_mob(s)
    h_comp = _h_mob(derivative(s)) / _h_mob(s)

    return (h_act=h_act, h_mob=h_mob, h_comp=h_comp)

end

"""
    hjorth(s)

Calculate Hjorths parameters.

# Arguments

- `s::AbstractArray`

# Returns

Named tuple containing:
- `h_act::Matrix{Float64}`: activity
- `h_mob::Matrix{Float64}`: mobility
- `h_comp::Matrix{Float64}`: complexity

# Notes:

- Activity: the total power of the signal
- Mobility: an estimate of the mean frequency
- Complexity: indicates the similarity of the shape of the signal to a pure sine wave
"""
function hjorth(s::AbstractArray)::@NamedTuple{h_act::Matrix{Float64}, h_mob::Matrix{Float64}, h_comp::Matrix{Float64}}

    _chk3d(s)
    ch_n = size(s, 1)
    ep_n = size(s, 3)

    h_act = zeros(ch_n, ep_n)
    h_mob = zeros(ch_n, ep_n)
    h_comp = zeros(ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            h_act[ch_idx, ep_idx], h_mob[ch_idx, ep_idx], h_comp[ch_idx, ep_idx] = @views hjorth(s[ch_idx, :, ep_idx])
        end
    end

    return (h_act=h_act, h_mob=h_mob, h_comp=h_comp)

end

"""
    hjorth(obj)

Calculate Hjorths parameters.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names

# Returns

Named tuple containing:
- `h_act::Matrix{Float64}`: activity
- `h_mob::Matrix{Float64}`: mobility
- `h_comp::Matrix{Float64}`: complexity

# Notes:

- Activity: the total power of the signal
- Mobility: an estimate of the mean frequency
- Complexity: indicates the similarity of the shape of the signal to a pure sine wave
"""
function hjorth(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::@NamedTuple{h_act::Matrix{Float64}, h_mob::Matrix{Float64}, h_comp::Matrix{Float64}}


    ch = exclude_bads ? get_channel(obj, ch=ch, exclude="bad") : get_channel(obj, ch=ch, exclude="")
    h_act, h_mob, h_comp = @views hjorth(obj.data[ch, :, :])

    return (h_act=h_act, h_mob=h_mob, h_comp=h_comp)

end
