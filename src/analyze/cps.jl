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
- `pw::Vector{Float64}`: cross power spectrum power
- `ph::Vector{Float64}`: cross power spectrum phase (in radians)
- `f::Vector{Float64}`: cross power spectrum frequencies
"""
function cps(s1::AbstractVector, s2::AbstractVector; fs::Int64)

    @assert fs >= 1 "fs must be ≥ 1."
    @assert length(s1) == length(s2) "s1 and s2 must have the same length."

    p = mt_cross_power_spectra(hcat(s1, s2)', fs=fs)

    pw = real.(p.power)[1, 2, :]
    ph = angle.(imag.(p.power))[1, 2, :]
    f = Vector(p.freq)

    return (pw=pw, ph=ph, f=f)

end

"""
    cps(s; fs)

Calculate cross power spectrum (channels vs channels).

# Arguments

- `s::AbstractArray`
- `fs::Int64`: sampling rate

# Returns

Named tuple containing:
- `pw::Array{Float64, 4}`: cross power spectrum power
- `ph::Array{Float64, 4}`: cross power spectrum phase (in radians)
- `f::Vector{Float64}`: cross power spectrum frequencies
"""
function cps(s::AbstractArray; fs::Int64)

    ch_n = size(s, 1)
    ep_n = size(s, 3)
    
    pw, ph, f = @views cps(s[1, :, 1], s[1, :, 1], fs=fs)
    pw = zeros(ch_n, ch_n, length(pw), ep_n)
    ph = zeros(ch_n, ch_n, length(ph), ep_n)

    # initialize progress bar
    progress_bar == true && (progbar = Progress(ep_n * ch_n, dt=1, barlen=20, color=:white))

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx1 in 1:ch_n
           for ch_idx2 in 1:ch_idx1
                pw[ch_idx1, ch_idx2, :, ep_idx], ph[ch_idx1, ch_idx2, :, ep_idx], _ = @views cps(s[ch_idx1, :, ep_idx], s[ch_idx2, :, ep_idx], fs=fs)
            end

        # update progress bar
        progress_bar == true && next!(progbar)
        end
    end

    @inbounds for cps_idx in 1:size(pw, 3)
        Threads.@threads for ep_idx in 1:ep_n
            for ch_idx1 in 1:(ch_n - 1)
                for ch_idx2 in (ch_idx1 + 1):ch_n
                    pw[ch_idx1, ch_idx2, cps_idx, ep_idx] = @views pw[ch_idx2, ch_idx1, cps_idx, ep_idx]
                    ph[ch_idx1, ch_idx2, cps_idx, ep_idx] = @views ph[ch_idx2, ch_idx1, cps_idx, ep_idx]
                end
            end
        end
    end

    return (pw=pw, ph=ph, f=f)

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
- `pw::Array{Float64, 4}`: cross power spectrum power
- `ph::Array{Float64, 4}`: cross power spectrum phase (in radians)
- `f::Vector{Float64}`: cross power spectrum frequencies
"""
function cps(s1::AbstractArray, s2::AbstractArray; fs::Int64)

    @assert size(s1) == size(s2) "s1 and s2 must have the same size."

    ch_n = size(s1, 1)
    ep_n = size(s1, 3)

    pw, ph, f = @views cps(s1[1, :, 1], s2[1, :, 1], fs=fs)

    pw = zeros(ch_n, length(pw), ep_n)
    ph = zeros(ch_n, length(ph), ep_n)
    f = zeros(ch_n, length(f), ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            pw[ch_idx, :, ep_idx], ph[ch_idx, :, ep_idx], f[ch_idx, :, ep_idx] = @views cps(s1[ch_idx, :, ep_idx], s2[ch_idx, :, ep_idx], fs=fs)
        end
    end

    return (pw=pw, ph=ph, f=f)

end

"""
    cps(obj; ch)

Calculate cross power spectrum.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

# Returns

Named tuple containing:
- `pw::Array{Float64, 4}`: cross power spectrum power
- `ph::Array{Float64, 4}`: cross power spectrum phase (in radians)
- `f::Vector{Float64, 4}`: cross power spectrum frequencies
"""
function cps(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj))

    _check_channels(obj, ch)

    pw, ph, f = cps(obj.data[ch, :, :], fs=sr(obj))

    return (pw=pw, ph=ph, f=f)

end

"""
    cps(obj1, obj2; ch1, ch2, ep1, ep2, norm)

Calculate cross power spectrum.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch1::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj1)`: index of channels, default is all signal channels
- `ch2::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj2)`: index of channels, default is all signal channels
- `ep1::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj1))`: default use all epochs
- `ep2::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj2))`: default use all epochs

# Returns

Named tuple containing:
- `pw::Array{Float64, 3}`: cross power spectrum power
- `ph::Array{Float64, 3}`: cross power spectrum phase (in radians)
- `f::Vector{Float64, 3}`: cross power spectrum frequencies
"""
function cps(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch1::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj1), ch2::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj2), ep1::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj1)), ep2::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj2)))

    @assert sr(obj1) == sr(obj2) "OBJ1 and OBJ2 must have the same sampling rate."

    _check_channels(obj1, ch1)
    _check_channels(obj2, ch2)
    @assert length(ch1) == length(ch2) "ch1 and ch2 must have the same length."
    
    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    @assert length(ep1) == length(ep2) "ep1 and ep2 must have the same length."
    @assert epoch_len(obj1) == epoch_len(obj2) "OBJ1 and OBJ2 must have the same epoch lengths."

    pw, ph, f = @views cps(reshape(obj1.data[ch1, :, ep1], length(ch1), :, length(ep1)), reshape(obj2.data[ch2, :, ep2], length(ch2), :, length(ep2)), fs=sr(obj1))

    return (pw=pw, ph=ph, f=f)

end
