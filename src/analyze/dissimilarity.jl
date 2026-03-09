export topo_var
export diss

"""
    topo_var(obj; <keyword arguments>)

Calculate the variance across channels at each time point of an ERP/ERF object (epoch 1 = the trial-averaged waveform), returning a time series of spatial variability.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: must be an ERP or ERF object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)

# Returns

- `tv::Vector{Float64}`: topographical variance
"""
function topo_var(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::Vector{Float64}

    @assert datatype(obj) in ["erp", "erf"] "topo_var() should be applied for ERP or ERF object only."
    
    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")

    tv = dropdims(var(@view(obj.data[ch, :, 1]), dims = 1), dims = 1)

    return tv

end

"""
    diss(s1, s2)

Calculate Global Dissimilarity and Spatial Correlation. Measure how different two GFP-normalised scalp topographies are at each time point:

- gd[t] = √( Σ_ch (g1[ch,t] - g2[ch,t])² / n_ch ) ∈ [0, 2]
- sc[t] = 1 - gd[t]² / 2  =  0.5 · (2 - gd[t]²) ∈ [-1, 1]

DISS = 0 (gd = 0) means identical topographies; DISS = 2 means opposite.

Spatial correlation is the linear rescaling of DISS onto [-1, 1].

# Arguments

- `s1::AbstractMatrix`: signal matrix (channels × samples)
- `s2::AbstractMatrix`: signal matrix (channels × samples)

# Returns

Named tuple containing:

- `gd::Vector{Float64}`: global dissimilarity ∈ [0, 2], one value per time point
- `sc::Vector{Float64}`: spatial correlation ∈ [-1, 1], one value per time point
"""
function diss(s1::AbstractMatrix, s2::AbstractMatrix)::@NamedTuple{gd::Vector{Float64}, sc::Vector{Float64}}

    @assert size(s1) == size(s2) "s1 and s2 must have the same size."

    # number of channels
    ch_n = size(s1, 1)

    g1 = erp_gfp_norm(s1)
    g2 = erp_gfp_norm(s2)

    # pre-allocate outputs
    gd = zeros(size(s1, 2))
    sc = zeros(size(s1, 2))

    Threads.@threads :dynamic for idx in axes(s1, 2)
        # Global dissimilarity: RMS difference of GFP-normalized topographies
        # sum(abs2, a .- b) = Σ(a-b)² avoids allocating (a.-b).^2 separately
        gd[idx] = sqrt(sum(abs2, @view(g1[:, idx]) .- @view(g2[:, idx])) / ch_n)

        # Spatial correlation: linear rescaling of gd² onto [-1, 1].
        # sc = 1 - gd²/2 = 0.5*(2 - gd²). Computed from the already-stored gd[idx] so no re-calculation is needed
        sc[idx] = 0.5 * (2 - gd[idx]^2)
    end

    return (gd = gd, sc = sc)

end

"""
    diss(obj1, obj2; <keyword arguments>)

Calculate DISS (global dissimilarity) and spatial correlation. Operates on ERP/ERF objects only and uses epoch 1 (the trial-averaged waveform).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: must be an ERP or ERF object
- `ch1::Union{String, Vector{String}, Regex}`: channel name(s)
- `ch2::Union{String, Vector{String}, Regex}`: channel name(s)

# Returns

Named tuple containing:

- `gd::Vector{Float64}`: global dissimilarity ∈ [0, 2], one value per time point
- `sc::Vector{Float64}`: spatial correlation ∈ [-1, 1], one value per time point
"""
function diss(
    obj1::NeuroAnalyzer.NEURO,
    obj2::NeuroAnalyzer.NEURO;
    ch1::Union{String, Vector{String}, Regex},
    ch2::Union{String, Vector{String}, Regex},
)::@NamedTuple{gd::Vector{Float64}, sc::Vector{Float64}}

    @assert datatype(obj1) in ["erp", "erf"] "diss() must be applied to ERP or ERF object only."
    @assert datatype(obj2) in ["erp", "erf"] "diss() must be applied to ERP or ERF object only."

    @assert sr(obj1) == sr(obj2) "OBJ1 and OBJ2 must have the same sampling rate."
    @assert epoch_len(obj1) == epoch_len(obj2) "OBJ1 and OBJ2 must have the same epoch lengths."

    # resolve channel names to integer indices, optionally skipping bad channels
    ch1 = exclude_bads ? get_channel(obj1, ch = ch1, exclude = "bad") : get_channel(obj1, ch = ch1, exclude = "")
    ch2 = exclude_bads ? get_channel(obj2, ch = ch2, exclude = "bad") : get_channel(obj2, ch = ch2, exclude = "")
    @assert length(ch1) == length(ch2) "Lengths of ch1 ($(length(ch1))) and ch2 ($(length(ch2))) must be equal."

    result = diss(
        @view(obj1.data[ch1, :, 1]),
        @view(obj2.data[ch2, :, 1])
    )

    return result

end
