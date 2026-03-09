export erp_gfp
export erp_gfp_norm
export erp_gfp

"""
    erp_gfp(s)

Calculate GFP (Global Field Power).

GFP is the population standard deviation across all selected channels at each time point. It is zero when all channels agree and maximal when they diverge. GFP peaks may reflect "interesting" brain activity warranting further investigation.

# Arguments

- `s::AbstractMatrix`: signal matrix (channels × samples)

# Returns

- `g::Vector{Float64}`: GFP at each time point (length = number of samples)

# Notes

GFP is the population standard deviation across all channels at each time point. It is zero when all channels have identical values and maximal when channel activity is maximally divergent.

GFP(t) = std_channels( s[:, t] )

This is used internally by diss() to make amplitudes comparable.
"""
function erp_gfp(s::AbstractMatrix)::Vector{Float64}

    # GFP = population std across channels at each time point
    # dropdims removes the trailing singleton dimension left by std(..., dims=1) without allocating a copy
    g = dropdims(std(s, dims = 1), dims = 1)

    return g

end

"""
    erp_gfp_norm(s)

Normalize a signal matrix by its GFP (Global Field Power).

Each column (time point) is divided by the GFP value at that time, so that the resulting matrix has unit standard deviation across channels at every sample. Returns the original matrix unchanged at any time point where GFP = 0.

# Arguments

- `s::AbstractMatrix`: signal matrix (channels × samples)

# Returns

- `gn::Matrix{Float64}`: GFP-normalised signal (channels × samples)
"""
function erp_gfp_norm(s::AbstractMatrix)::Matrix{Float64}

    g = erp_gfp(s)
    gn = similar(s)

    # divide each channel row by the GFP vector element-wise
    # g' broadcasts across rows (channels)
    @. gn = s / g'

    return gn

end

"""
    erp_gfp(obj)

Calculate GFP (Global Field Power).

GFP is the population standard deviation across all selected channels at each time point. It is zero when all channels agree and maximal when they diverge. GFP peaks may reflect "interesting" brain activity warranting further investigation.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `norm::Bool=false`: if true, return the GFP-normalised signal instead of GFP values

# Returns

- `erp_gfp::Union{Vector{Float64}, Matrix{Float64}}`: GFP values over time (norm=false) or the GFP-normalized signal matrix (norm=true)

# Notes

GFP is the population standard deviation across all channels at each time point. It is zero when all channels have identical values and maximal when channel activity is maximally divergent.

GFP(t) = std_channels( s[:, t] )
"""
function erp_gfp(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    norm::Bool = false
)::Union{Vector{Float64}, Matrix{Float64}}

    _check_datatype(obj, ["erp", "erf"])

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")
    @assert length(ch) > 1 "More than 1 channel must be selected."

    s = @view obj.data[ch, :, 1]

    return norm ? erp_gfp_norm(s) : erp_gfp(s)

end
