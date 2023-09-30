export tcoherence

"""
    tcoherence(s1, s2; pad)

Calculate coherence (mean over time), IC (imaginary coherence) and MSC (magnitude-squared coherence).

# Arguments

- `s1::AbstractVector`
- `s2::AbstractVector`
- `pad::Int64=0`: number of zeros to add

# Returns

Named tuple containing:
- `c::Vector{Float64}`: coherence
- `msc::Vector{Float64}`: magnitude-squares coherence
- `ic::Vector{Float64}`: imaginary part of coherence
"""
function tcoherence(s1::AbstractVector, s2::AbstractVector; pad::Int64=0)

    @assert length(s1) == length(s2) "Both signals must have the same length."

    s1_fft = fft0(s1, pad) ./ length(s1)
    s2_fft = fft0(s2, pad) ./ length(s2)

    coh = @. (abs((s1_fft) * conj.(s2_fft))^2) / (s1_fft * s2_fft)

    c=real.(coh)
    msc = @. abs(coh)^2
    ic = imag.(coh)

    return (c=c, msc=msc, ic=ic)
end

"""
    tcoherence(s1, s2; ch1, ch2, ep1, ep2)

Calculate coherence (mean over time), IC (imaginary coherence) and MSC (magnitude-squared coherence).

# Arguments

- `s1::AbstractArray`
- `s2::AbstractArray`
- `pad::Int64=0`: number of zeros to add signal for FFT

# Returns

Named tuple containing:
- `c::Array{Float64, 3}`: coherence
- `msc::Array{Float64, 3}`: MSC
- `ic::Array{Float64, 3}`: imaginary part of coherence
"""
function tcoherence(s1::AbstractArray, s2::AbstractArray; pad::Int64=0)

    @assert size(s1) == size(s2) "s1 and s2 must have the same size."
    
    ch_n = size(s1, 1)
    ep_n = size(s1, 3)

    c = similar(s1)
    msc = similar(s1)
    ic = similar(s1)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            c[ch_idx, :, ep_idx], msc[ch_idx, :, ep_idx], ic[ch_idx, :, ep_idx] = @views tcoherence(s1[ch_idx, :, ep_idx], s2[ch_idx, :, ep_idx], pad=pad)
        end
    end

    return (c=c, msc=msc, ic=ic)
end

"""
    tcoherence(obj1, obj2; ch1, ch2, ep1, ep2)

Calculate coherence (mean over time), IC (imaginary coherence) and MSC (magnitude-squared coherence).

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `ch1::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj1)`: index of channels, default is all signal channels
- `ch2::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj2)`: index of channels, default is all signal channels
- `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1))`: default use all epochs
- `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2))`: default use all epochs
- `pad::Int64=0`: number of zeros to add signal for FFT

# Returns

Named tuple containing:
- `c::Array{Float64, 3}`: coherence
- `msc::Array{Float64, 3}`: MSC
- `ic::Array{Float64, 3}`: imaginary part of coherence
"""
function tcoherence(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch1::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj1), ch2::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj2), ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1)), ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2)), pad::Int64=0)

    _check_channels(obj1, ch1)
    _check_channels(obj2, ch2)
    @assert length(ch1) == length(ch2) "ch1 and ch2 must have the same length."

    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    @assert length(ep1) == length(ep2) "ep1 and ep2 must have the same length."
    @assert epoch_len(obj1) == epoch_len(obj2) "OBJ1 and OBJ2 must have the same epoch lengths."

    c, msc, ic = @views tcoherence(reshape(obj1.data[ch1, :, ep1], length(ch1), :, length(ep1)), reshape(obj2.data[ch2, :, ep2], length(ch2), :, length(ep2)), pad=pad)

    return (c=c, msc=msc, ic=ic)
end
