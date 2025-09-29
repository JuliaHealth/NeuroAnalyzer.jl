export ghexp

"""
    ghexp(s; tau_range, q_range)

Calculate the Generalised Hurst Exponents (GHEs).

# Arguments

- `s::AbstractVector`
- `tau_range::UnitRange{Int64}`: GHEs are estimated the moments of the absolute increments over time delay
- `q_range::Union{Nothing, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}}=nothing`: moments at which GHE are estimated

# Returns

- `ghe::Matrix{Float64}`
"""
function ghexp(s::AbstractVector; tau_range::UnitRange{Int64}, q_range::Union{Nothing, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}}=nothing)::Matrix{Float64}

    @assert tau_range[end] < length(s) "End of tau_range ($(tau_range[end])) must be < length of s ($(length(s)))."

    if isnothing(q_range)
        ghe = hurst_exponent(s, tau_range)
    else
        ghe = generalised_hurst_range(s, tau_range, q_range)
    end

    return ghe

end

"""
    ghexp(s)

Calculate the Generalised Hurst Exponents (GHEs).

# Arguments

- `s::AbstractArray`
- `tau_range::UnitRange{Int64}`: GHEs are estimated the moments of the absolute increments over time delay
- `q_range::Union{Nothing, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}}=nothing`: moments at which GHE are estimated

# Returns

- `ghe::Array{Float64, 4}`
"""
function ghexp(s::AbstractArray; tau_range::UnitRange{Int64}, q_range::Union{Nothing, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}}=nothing)::Array{Float64, 4}

    _chk3d(s)
    ch_n = size(s, 1)
    ep_n = size(s, 3)

    if isnothing(q_range)
        ghe = zeros(ch_n, 1, 2, ep_n)
    else
        ghe = zeros(ch_n, length(q_range), 2, ep_n)
    end

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads :greedy for ch_idx in 1:ch_n
            ghe[ch_idx, :, :, ep_idx] = ghexp(s[ch_idx, :, ep_idx], tau_range=tau_range, q_range=q_range)
        end
    end

    return ghe

end

"""
    ghexp(obj; <keyword arguments>)

Calculate the Generalised Hurst Exponents (GHEs).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `tau_range::UnitRange{Int64}`: GHEs are estimated the moments of the absolute increments over time delay
- `q_range::Union{Nothing, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}}=nothing`: moments at which GHE are estimated

# Returns

- `ghe::Array{Float64, 4}`
"""
function ghexp(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, tau_range::UnitRange{Int64}, q_range::Union{Nothing, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}}=nothing)::Array{Float64, 4}

    ch = exclude_bads ? get_channel(obj, ch=ch, exclude="bad") : get_channel(obj, ch=ch, exclude="")

    ghe = @views ghexp(obj.data[ch, :, :], tau_range=tau_range, q_range=q_range)

    return ghe

end
