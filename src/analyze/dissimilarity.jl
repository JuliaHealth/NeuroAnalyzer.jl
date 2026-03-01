export topo_var
export diss

"""
    topo_var(obj; <keyword arguments>)

Calculate topographical variance (variance calculated at each time point across all channels).

# Arguments

  - `obj::NeuroAnalyzer.NEURO`
  - `ch::Union{String, Vector{String}}`: channels to analyze

# Returns

  - `tv::Vector{Float64}`: topographical variance
"""
function topo_var(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}})::Vector{Float64}

    @assert datatype(obj) in ["erp", "erf"] "topo_var() should be applied for ERP or ERF object only."
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")

    tv = @views var(obj.data[ch, :, 1], dims = 1)[:]

    return tv

end

"""
    diss(s1, s2)

Calculate DISS (global dissimilarity) and spatial correlation (channels vs channels).

# Arguments

  - `s1::AbstractMatrix`
  - `s2::AbstractMatrix`

# Returns

Named tuple containing:

  - `gd::Vector{Float64}`: global dissimilarity
  - `sc::Vector{Float64}`: spatial correlation
"""
function diss(s1::AbstractMatrix, s2::AbstractMatrix)::@NamedTuple{gd::Vector{Float64}, sc::Vector{Float64}}

    @assert size(s1) == size(s2) "s1 and s2 must have the same size."

    n_ch = size(s1, 1)
    g1 = erp_gfp_norm(s1)
    g2 = erp_gfp_norm(s2)

    gd = zeros(size(s1, 2))
    sc = zeros(size(s1, 2))
    Threads.@threads for idx in axes(s1, 2)
        gd[idx] = sqrt(sum((g1[:, idx] .- g2[:, idx]) .^ 2) / n_ch)
        sc[idx] = 0.5 * (2 - gd[idx]^2)
    end

    return (gd = gd, sc = sc)

end

"""
    diss(obj1, obj2; <keyword arguments>)

Calculate DISS (global dissimilarity) and spatial correlation (channels vs channels). This works for ERP/ERF objects only and calculates DISS for the first epoch (ERP/ERF) only.

# Arguments

  - `obj::NeuroAnalyzer.NEURO`
  - `ch1::Union{String, Vector{String}}: list of channels
  - `ch2::Union{String, Vector{String}}: list of channels
  - `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1))`: default use all epochs
  - `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2))`: default use all epochs

# Returns

Named tuple containing:

  - `gd::Vector{Float64}`: global dissimilarity
  - `sc::Vector{Float64}`: spatial correlation
"""
function diss(
        obj1::NeuroAnalyzer.NEURO,
        obj2::NeuroAnalyzer.NEURO;
        ch1::Union{String, Vector{String}},
        ch2::Union{String, Vector{String}},
    )::@NamedTuple{gd::Vector{Float64}, sc::Vector{Float64}}

    @assert datatype(obj1) in ["erp", "erf"] "diss() should be applied for ERP or ERF object only."
    @assert datatype(obj2) in ["erp", "erf"] "diss() should be applied for ERP or ERF object only."
    @assert sr(obj1) == sr(obj2) "OBJ1 and OBJ2 must have the same sampling rate."
    @assert length(ch1) == length(ch2) "Lengths of ch1 ($(length(ch1)) and ch2 ($(length(ch2)) must be equal."
    @assert epoch_len(obj1) == epoch_len(obj2) "OBJ1 and OBJ2 must have the same epoch lengths."

    ch1 = exclude_bads ? get_channel(obj1, ch = ch1, exclude = "bad") : get_channel(obj1, ch = ch1, exclude = "")
    ch2 = exclude_bads ? get_channel(obj2, ch = ch2, exclude = "bad") : get_channel(obj2, ch = ch2, exclude = "")

    gd, sc = @views diss(obj1.data[ch1, :, 1], obj2.data[ch2, :, 1])

    return (gd = gd, sc = sc)

end
