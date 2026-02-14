export gfp
export gfp_norm

"""
    gfp(s)

Calculate GFP (Global Field Power).

# Arguments

- `s::AbstractMatrix`

# Returns

- `g::Vector{Float64}`: GFP values over time points
"""
function gfp(s::AbstractMatrix)::Vector{Float64}

    g = std(s, dims=1)[:]

    return g

end

"""
    gfp_norm(s)

Calculate signal normalized for GFP (Global Field Power).

# Arguments

- `s::AbstractMatrix`

# Returns

- `gn::Matrix{Float64}`: normalized signal
"""
function gfp_norm(s::AbstractMatrix)::Matrix{Float64}

    g = gfp(s)
    gn = similar(s)
    Threads.@threads for ch_idx in axes(s, 1)
        @views gn[ch_idx, :] = s[ch_idx, :] ./ g
    end

    return gn

end

"""
    gfp(obj)

Calculate Global Field Power (GFP).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `norm::Bool=false`: if true, calculate signal normalized for GFP

# Returns

- `gfp::Union{Vector{Float64}, Matrix{Float64}}`

# Notes

Global field power is a measure of agreement of the signals picked up by all sensors across the entire scalp: if all sensors have the same value at a given time point, the GFP will be zero at that time point; if the signals differ, the GFP will be non-zero at that time point. GFP peaks may reflect “interesting” brain activity, warranting further investigation. Mathematically, the GFP is the population standard deviation across all sensors, calculated separately for every time point.
"""
function gfp(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::Union{Vector{Float64}, Matrix{Float64}}

    _check_datatype(obj, ["erp", "erf"])

    ch = exclude_bads ? get_channel(obj, ch=ch, exclude="bad") : get_channel(obj, ch=ch, exclude="")
    @assert length(ch) > 1 "More than 1 channel must be selected."
    s = @views obj.data[ch, :, 1]
    if norm
        return gfp(s)
    else
        return gfp_norm(s)
    end

end
