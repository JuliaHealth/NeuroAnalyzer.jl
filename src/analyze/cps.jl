export cps

"""
    cps(s1, s2; fs)

Calculate cross power spectrum.

# Arguments

- `s1::AbstractVector`
- `s2::AbstractVector`
- `fs::Int64`: sampling rate

# Returns

Named tuple containing:
- `cps_pw::Vector{Float64}`: cross power spectrum power
- `cps_ph::Vector{Float64}`: cross power spectrum phase (in radians)
- `cps_fq::Vector{Float64}`: cross power spectrum frequencies
"""
function cps(s1::AbstractVector, s2::AbstractVector; fs::Int64)

    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))
    length(s1) == length(s2) || throw(ArgumentError("s1 and s2 must have the same length."))

    p = mt_cross_power_spectra(hcat(s1, s2)', fs=fs)

    cps_pw = real.(p.power)[1, 2, :]
    cps_ph = angle.(imag.(p.power))[1, 2, :]
    cps_fq = Vector(p.freq)

    return (cps_pw=cps_pw, cps_ph=cps_ph, cps_fq=cps_fq)

end

"""
    cps(s; fs)

Calculate cross power spectrum (channels vs channels).

# Arguments

- `s::AbstractArray`
- `fs::Int64`: sampling rate

# Returns

Named tuple containing:
- `cps_pw::Array{Float64, 4}`: cross power spectrum power
- `cps_ph::Array{Float64, 4}`: cross power spectrum phase (in radians)
- `cps_fq::Vector{Float64}`: cross power spectrum frequencies
"""
function cps(s::AbstractArray; fs::Int64)

    ch_n = size(s, 1)
    ep_n = size(s, 3)
    
    cps_pw, cps_ph, cps_fq = @views cps(s[1, :, 1], s[1, :, 1], fs=fs)
    cps_pw = zeros(ch_n, ch_n, length(cps_pw), ep_n)
    cps_ph = zeros(ch_n, ch_n, length(cps_ph), ep_n)

    # initialize progress bar
    progress_bar == true && (p = Progress(ep_n * ch_n, 1))

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx1 in 1:ch_n
           for ch_idx2 in 1:ch_idx1
                cps_pw[ch_idx1, ch_idx2, :, ep_idx], cps_ph[ch_idx1, ch_idx2, :, ep_idx], _ = @views cps(s[ch_idx1, :, ep_idx], s[ch_idx2, :, ep_idx], fs=fs)
            end

        # update progress bar
        progress_bar == true && next!(p)
        end
    end

    @inbounds @simd for cps_idx in 1:size(cps_pw, 3)
        Threads.@threads for ep_idx in 1:ep_n
            for ch_idx1 in 1:(ch_n - 1)
                for ch_idx2 in (ch_idx1 + 1):ch_n
                    cps_pw[ch_idx1, ch_idx2, cps_idx, ep_idx] = @views cps_pw[ch_idx2, ch_idx1, cps_idx, ep_idx]
                    cps_ph[ch_idx1, ch_idx2, cps_idx, ep_idx] = @views cps_ph[ch_idx2, ch_idx1, cps_idx, ep_idx]
                end
            end
        end
    end

    return (cps_pw=cps_pw, cps_ph=cps_ph, cps_fq=cps_fq)

end

"""
    cps(s1, s2; fs)

Calculate cross power spectrum (channels of `s1` vs channels of `s2`).

# Arguments

- `s1::AbstractArray`
- `s2::AbstractArray`
- `fs::Int64`: sampling rate

# Returns

Named tuple containing:
- `cps_pw::Array{Float64, 4}`: cross power spectrum power
- `cps_ph::Array{Float64, 4}`: cross power spectrum phase (in radians)
- `cps_fq::Vector{Float64}`: cross power spectrum frequencies
"""
function cps(s1::AbstractArray, s2::AbstractArray; fs::Int64)

    size(s1) == size(s2) || throw(ArgumentError("s1 and s2 must have the same size."))

    ch_n = size(s1, 1)
    ep_n = size(s1, 3)

    cps_pw, cps_ph, cps_fq = @views cps(s1[1, :, 1], s2[1, :, 1], fs=fs)

    cps_pw = zeros(ch_n, length(cps_pw), ep_n)
    cps_ph = zeros(ch_n, length(cps_ph), ep_n)
    cps_fq = zeros(ch_n, length(cps_fq), ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            cps_pw[ch_idx, :, ep_idx], cps_ph[ch_idx, :, ep_idx], cps_fq[ch_idx, :, ep_idx] = @views cps(s1[ch_idx, :, ep_idx], s2[ch_idx, :, ep_idx], fs=fs)
        end
    end

    return (cps_pw=cps_pw, cps_ph=cps_ph, cps_fq=cps_fq)

end

"""
    cps(obj; ch)

Calculate cross power spectrum (channels vs channels).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

# Returns

Named tuple containing:
- `cps_pw::Array{Float64, 4}`: cross power spectrum power
- `cps_ph::Array{Float64, 4}`: cross power spectrum phase (in radians)
- `cps_fq::Vector{Float64, 4}`: cross power spectrum frequencies
"""
function cps(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj))

    _check_channels(obj, ch)

    cps_pw, cps_ph, cps_fq = cps(obj.data[ch, :, :], fs=sr(obj))

    return (cps_pw=cps_pw, cps_ph=cps_ph, cps_fq=cps_fq)

end

"""
    cps(obj1, obj2; ch1, ch2, ep1, ep2, norm)

Calculate cross power spectrum (`ch1` of `obj1` vs `ch2` of `obj2`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch1::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj1)`: index of channels, default is all signal channels
- `ch2::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj2)`: index of channels, default is all signal channels
- `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1))`: default use all epochs
- `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2))`: default use all epochs

# Returns

Named tuple containing:
- `cps_pw::Array{Float64, 3}`: cross power spectrum power
- `cps_ph::Array{Float64, 3}`: cross power spectrum phase (in radians)
- `cps_fq::Vector{Float64, 3}`: cross power spectrum frequencies
"""
function cps(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch1::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj1), ch2::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj2), ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1)), ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2)))

    sr(obj1) == sr(obj2) || throw(ArgumentError("OBJ1 and OBJ2 must have the same sampling rate."))

    _check_channels(obj1, ch1)
    _check_channels(obj2, ch2)
    length(ch1) == length(ch2) || throw(ArgumentError("ch1 and ch2 must have the same length."))
    
    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    length(ep1) == length(ep2) || throw(ArgumentError("ep1 and ep2 must have the same length."))
    epoch_len(obj1) == epoch_len(obj2) || throw(ArgumentError("OBJ1 and OBJ2 must have the same epoch lengths."))

    cps_pw, cps_ph, cps_fq = @views cps(obj1.data[ch1, :, ep1], obj2.data[ch2, :, ep2], fs=sr(obj1))

    return (cps_pw=cps_pw, cps_ph=cps_ph, cps_fq=cps_fq)

end
