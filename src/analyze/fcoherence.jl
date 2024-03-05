export fcoherence

"""
    fcoherence(s; fs, frq_lim)

Calculate coherence and MSC (magnitude-squared coherence).

# Arguments

- `s::AbstractMatrix`
- `fs::Int64`: sampling rate
- `frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing`: return coherence only for the given frequency range

# Returns

Named tuple containing:
- `c::Array{Float64, 3}`: coherence
- `msc::Array{Float64, 3}`: magnitude-squared coherence
- `f::Vector{Float64}`: frequencies
"""
function fcoherence(s::AbstractMatrix; fs::Int64, frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing)

    @assert fs >= 1 "fs must be ≥ 1."

    c = mt_coherence(s, fs=fs)
    f = Vector(c.freq)
    c = c.coherence

    if frq_lim !== nothing
        _check_tuple(frq_lim, "frq_lim", (0, fs / 2))
        idx1 = vsearch(frq_lim[1], f)
        idx2 = vsearch(frq_lim[2], f)
        c = c[:, :, idx1:idx2]
        f = f[idx1:idx2]
    end

    return (c=c, msc=c.^2, f=f)

end

"""
    fcoherence(s1, s2; fs, frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing)

Calculate coherence and MSC (magnitude-squared coherence).

# Arguments

- `s1::AbstractMatrix`
- `s2::AbstractMatrix`
- `fs::Int64`
- `frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing`: return coherence only for the given frequency range

# Returns

- `c::Array{Float64, 3}`: coherence
- `msc::Array{Float64, 3}`: magnitude-squared coherence
- `f::Vector{Float64}`: frequencies
"""
function fcoherence(s1::AbstractMatrix, s2::AbstractMatrix; fs::Int64, frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing)

    @assert length(s1) == length(s2) "s1 and s2 must have the same length."
    @assert fs >= 1 "fs must be ≥ 1."

    s = vcat(s1, s2)

    c = mt_coherence(s, fs=fs)
    f = Vector(c.freq)
    c = c.coherence
    msc = @. abs(c)^2

    if frq_lim !== nothing
        _check_tuple(frq_lim, "frq_lim", (0, fs / 2))
        idx1 = vsearch(frq_lim[1], f)
        idx2 = vsearch(frq_lim[2], f)
        c = c[:, :, idx1:idx2]
        msc = msc[:, :, idx1:idx2]
        f = f[idx1:idx2]
    end
    
    return (c=c, msc=msc, f=f)

end

"""
    fcoherence(s1, s2; fs, frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing)

Calculate coherence and MSC (magnitude-squared coherence).

# Arguments

- `s1::AbstractArray`
- `s2::AbstractArray`
- `fs::Int64`
- `frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing`: return coherence only for the given frequency range

# Returns

- `c::Array{Float64, 3}`: coherence
- `msc::Array{Float64, 3}`: magnitude-squared coherence
- `f::Vector{Float64}`: frequencies
"""
function fcoherence(s1::AbstractArray, s2::AbstractArray; fs::Int64, frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing)

    @assert size(s1) == size(s2) "s1 and s2 must have the same length."

    ep_n = size(s1, 3)
    
    c, msc, f = fcoherence(s1[:, :, 1], s2[:, :, 1], fs=fs, frq_lim=frq_lim)

    c = zeros(size(c, 1), size(c, 2), size(c, 3), ep_n)
    msc = zeros(size(msc, 1), size(msc, 2), size(msc, 3), ep_n)

    @inbounds for ep_idx in 1:ep_n
        c[:, :, :, ep_idx], msc[:, :, :, ep_idx], _ = fcoherence(s1[:, :, ep_idx], s2[:, :, ep_idx], fs=fs, frq_lim=frq_lim)
    end

    return (c=c, msc=msc, f=f)

end

"""
    fcoherence(obj1, obj2; ch1, ch2, ep1, ep2, frq_lim)

Calculate coherence and MSC (magnitude-squared coherence).

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `ch1::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj1)`: index of channels, default is all signal channels
- `ch2::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj2)`: index of channels, default is all signal channels
- `ep1::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj1))`: default use all epochs
- `ep2::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj2))`: default use all epochs
- `frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing`: return coherence only for the given frequency range

# Returns

Named tuple containing:
- `c::Array{Float64, 3}`: coherence
- `msc::Array{Float64, 3}`: magnitude-squared coherence
- `f::Vector{Float64}`: frequencies
"""
function fcoherence(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch1::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj1), ch2::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj2), ep1::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj1)), ep2::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj2)), frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing)

    _check_channels(obj1, ch1)
    _check_channels(obj2, ch2)
    @assert length(ch1) == length(ch2) "ch1 and ch2 lengths must be equal."
    
    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    @assert length(ep1) == length(ep2) "ep1 and ep2 lengths must be equal."
    @assert epoch_len(obj1) == epoch_len(obj2) "OBJ1 and OBJ2 must have the same epoch lengths."

    @assert sr(obj1) == sr(obj2) "OBJ1 and OBJ2 must have the same sampling rate."

    size(ch1) == () && (ch1 = [ch1])
    size(ch2) == () && (ch2 = [ch2])
    size(ep1) == () && (ep1 = [ep1])
    size(ep2) == () && (ep2 = [ep2])

    c, msc, f = @views fcoherence(obj1.data[ch1, :, ep1], obj2.data[ch2, :, ep2], fs=sr(obj1), frq_lim=frq_lim)

    return (c=c, msc=msc, f=f)

end
