export tcoherence

"""
    tcoherence(s1, s2; pad)

Calculate coherence, IC (imaginary coherence) and MSC (magnitude-squared coherence). If signals have different lengths, the shorter signal will be padded with zeros to match their lengths.

# Arguments

- `s1::AbstractVector`
- `s2::AbstractVector`
- `pad::Int64=0`: number of zeros to add
- `frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing`: return coherence only for the given frequency range

# Returns

Named tuple containing:
- `c::Vector{Float64}`: coherence
- `msc::Vector{Float64}`: magnitude-squared coherence
- `ic::Vector{Float64}`: imaginary part of coherence
- `f::Vector{Float64}`: frequencies
"""
function tcoherence(s1::AbstractVector, s2::AbstractVector; pad::Int64=0, fs::Int64, frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing)

    s1, s2 = _veqlen(s1, s2)

    s1_fft = rfft0(s1, pad) ./ length(s1)
    s2_fft = rfft0(s2, pad) ./ length(s2)

    coh = @. (abs((s1_fft) * conj.(s2_fft))^2) / (s1_fft * s2_fft)

    c = real.(coh)
    msc = @. abs(coh)^2
    ic = imag.(coh)
    f, _ = freqs(s1, fs)

    if frq_lim !== nothing
        _check_tuple(frq_lim, "frq_lim", (0, fs / 2))
        idx1 = vsearch(frq_lim[1], f)
        idx2 = vsearch(frq_lim[2], f)
        c = c[idx1:idx2]
        msc = msc[idx1:idx2]
        ic = ic[idx1:idx2]
        f = f[idx1:idx2]
    end

    return (c=c, msc=msc, ic=ic, f=f)
end

"""
    tcoherence(s1, s2; ch1, ch2, ep1, ep2)

Calculate coherence, IC (imaginary coherence) and MSC (magnitude-squared coherence).

# Arguments

- `s1::AbstractArray`
- `s2::AbstractArray`
- `pad::Int64=0`: number of zeros to add signal for FFT
- `frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing`: return coherence only for the given frequency range

# Returns

Named tuple containing:
- `c::Array{Float64, 3}`: coherence
- `msc::Array{Float64, 3}`: MSC
- `ic::Array{Float64, 3}`: imaginary part of coherence
- `f::Vector{Float64}`: frequencies
"""
function tcoherence(s1::AbstractArray, s2::AbstractArray; pad::Int64=0, fs::Int64, frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing)

    @assert size(s1) == size(s2) "s1 and s2 must have the same size."
    
    ch_n = size(s1, 1)
    ep_n = size(s1, 3)

    f, _ = freqs(s1[1, :, 1], fs)
    if frq_lim !== nothing
        _check_tuple(frq_lim, "frq_lim", (0, fs / 2))
        idx1 = vsearch(frq_lim[1], f)
        idx2 = vsearch(frq_lim[2], f)
        f = f[idx1:idx2]
    end

    c = zeros(ch_n, length(f), ep_n)
    msc = zeros(ch_n, length(f), ep_n)
    ic = zeros(ch_n, length(f), ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            c[ch_idx, :, ep_idx], msc[ch_idx, :, ep_idx], ic[ch_idx, :, ep_idx], _ = @views tcoherence(s1[ch_idx, :, ep_idx], s2[ch_idx, :, ep_idx], pad=pad, fs=fs, frq_lim=frq_lim)
        end
    end

    return (c=c, msc=msc, ic=ic, f=f)
end

"""
    tcoherence(obj1, obj2; ch1, ch2, ep1, ep2)

Calculate coherence, IC (imaginary coherence) and MSC (magnitude-squared coherence).

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `ch1::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj1)`: index of channels, default is all signal channels
- `ch2::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj2)`: index of channels, default is all signal channels
- `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1))`: default use all epochs
- `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2))`: default use all epochs
- `pad::Int64=0`: number of zeros to add signal for FFT
- `frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing`: return coherence only for the given frequency range

# Returns

Named tuple containing:
- `c::Array{Float64, 3}`: coherence
- `msc::Array{Float64, 3}`: MSC
- `ic::Array{Float64, 3}`: imaginary part of coherence
- `f::Vector{Float64}`: frequencies
"""
function tcoherence(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch1::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj1), ch2::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj2), ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1)), ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2)), pad::Int64=0, frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing)

    _check_channels(obj1, ch1)
    _check_channels(obj2, ch2)
    @assert length(ch1) == length(ch2) "ch1 and ch2 must have the same length."

    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    @assert length(ep1) == length(ep2) "ep1 and ep2 must have the same length."
    @assert epoch_len(obj1) == epoch_len(obj2) "OBJ1 and OBJ2 must have the same epoch lengths."

    @assert sr(obj1) == sr(obj2) "OBJ1 and OBJ2 must have the same sampling rate."

    size(ch1) == () && (ch1 = [ch1])
    size(ch2) == () && (ch2 = [ch2])
    size(ep1) == () && (ep1 = [ep1])
    size(ep2) == () && (ep2 = [ep2])

    c, msc, ic, f = @views tcoherence(obj1.data[ch1, :, ep1], obj2.data[ch2, :, ep2], pad=pad, fs=sr(obj1), frq_lim=frq_lim)

    return (c=c, msc=msc, ic=ic, f=f)

end
