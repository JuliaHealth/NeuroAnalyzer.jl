export hz2rads
export rads2hz
export s2t
export t2s
export freqs
export frqinst

"""
    hz2rads(f)

Convert frequency `f` in Hz to rad/s.

# Arguments

- `f::Real`

# Returns

- `f_rads::Float64`
"""
function hz2rads(f::Real)
    return 2 * pi * f
end

"""
    rads2hz(f)

Convert frequency `f` in rad/s to Hz.

# Arguments

- `f::Real`

# Returns

- `f_rads::Float64`
"""
function rads2hz(f::Real)
    return f / 2 * pi
end

"""
    t2f(t)

Convert cycle length in ms `t` to frequency.

# Arguments

- `t::Real`: cycle length in ms

# Returns

- `f::Float64`: frequency in Hz
"""
function t2f(t::Real)
    t <= 0 && throw(ArgumentError("t must be > 0."))
    return round(1000 / t, digits=2)
end

"""
    f2t(f)

Convert frequency `f` to cycle length in ms.

# Arguments

- `f::Real`: frequency in Hz

# Returns

- `f::Float64`: cycle length in ms
"""
function f2t(f::Real)
    f <= 0 && throw(ArgumentError("f must be > 0."))
    return round(1000 / f, digits=2)
end

"""
    freqs(t)

Return vector of frequencies and Nyquist frequency for time vector.

# Arguments

- `t::AbstractVector, AbstractRange}`

# Returns

- `hz::Vector{Float64}`
- `nyquist_freq::Float64`
"""
function freqs(t::Union{AbstractVector, AbstractRange})

    typeof(t) <: AbstractRange && (t = collect(t))

    # sampling interval
    dt = t[2] - t[1]
    # sampling rate
    fs = round(Int64, 1 / dt)
    # Nyquist frequency
    nyquist_freq = fs / 2
    # frequency array
    hz = linspace(0, nyquist_freq, floor(Int64, length(t) / 2))

    return hz, nyquist_freq
end

"""
    freqs(signal, fs)

Return vector of frequencies and Nyquist frequency for signal.

# Arguments

- `signal::Vector{Float64}`
- `fs::Int64`

# Returns

- `hz::Vector{Float64`
- `nyquist_freq::Float64`
"""
function freqs(signal::Vector{Float64}, fs::Int64)

    fs < 0 && throw(ArgumentError("Sampling rate must be >0 Hz."))

    # Nyquist frequency
    nyquist_freq = fs / 2
    # frequency array
    hz = linspace(0, nyquist_freq, floor(Int64, length(signal) / 2))

    return hz, nyquist_freq
end

"""
    freqs(obj)

Return vector of frequencies and Nyquist frequency.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

Named tuple containing:
- `hz::Vector{Float64}`
- `nyquist::Float64`
"""
function freqs(obj::NeuroAnalyzer.NEURO)
    hz, nyq = freqs(obj.data[1, :, 1], sr(obj))
    return (hz=hz, nyquist=nyq)
end

"""
    frqinst(signal; fs)

Calculate instantaneous frequency.

# Arguments

- `signal::AbstractVector`
- `fs::Int64`

# Returns

- `frqinst::Vector{Float64}`
"""
function frqinst(signal::AbstractVector; fs::Int64)
    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))
    _, _, _, h_phases = hspectrum(signal)
    return 256 * derivative(h_phases) / (2*pi)
end

"""
    frqinst(obj; channel)

Calculate instantaneous frequency.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

# Returns

- `frqinst::Array{Float64, 3}`
"""
function frqinst(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj))

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    frqinst = zeros(ch_n, epoch_len(obj), ep_n)
    fs = sr(obj)

    _info("frqinst() uses Hilbert transform, the signal should be narrowband for best results.")

    # initialize progress bar
    progress_bar == true && (p = Progress(ep_n * ch_n, 1))
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            frqinst[ch_idx, :, ep_idx] = @views frqinst(obj.data[channel[ch_idx], :, ep_idx], fs=fs)
        end
        # update progress bar
        progress_bar == true && next!(p)
    end
    return frqinst
end
