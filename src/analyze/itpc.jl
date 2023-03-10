export itpc
export itpc_s

"""
    itpc(signal; t)

Calculate ITPC (Inter-Trial-Phase Clustering) at sample number `t` over epochs/trials.

# Arguments

- `signal::AbstractArray`: one channel over epochs/trials
- `t::Int64`: time point (sample number) at which ITPC is calculated
- `w::Union{AbstractVector, Nothing}`: optional vector of epochs/trials weights for wITPC calculation

# Returns

Named tuple containing:
- `itpc_value::Float64`: ITPC value
- `itpcz::Float64`: Rayleigh's ITPC Z value
- `itpc_angle::Float64`: ITPC angle
- `itpc_phases::Vector{Float64}`: phases at time `t` averaged across trials/epochs
"""
function itpc(signal::AbstractArray; t::Int64, w::Union{AbstractVector, Nothing}=nothing)

    t < 1 && throw(ArgumentError("t must be ≥ 1."))
    t > size(signal, 2) && throw(ArgumentError("t must be ≤ $(size(signal, 2))."))
    size(signal, 1) == 1 || throw(ArgumentError("signal must have 1 channel."))
    ep_n = size(signal, 3)

    w === nothing && (w = ones(ep_n))
    # scale w if w contains negative values
    any(i -> i < 0, w) && (w .+= abs(minimum(w)))
    length(w) == ep_n || throw(ArgumentError("Length of w should be equal to number of epochs ($ep_n)."))
    
    s_phase = zeros(size(signal, 2), ep_n)
    @inbounds @simd for ep_idx in 1:ep_n
        _, _, _, s_phase[:, ep_idx] = @views hspectrum(signal[1, :, ep_idx])
    end
 
    itpc_phases = @view s_phase[t, :]
    itpc_value = abs.(mean(exp.(1im .* itpc_phases .* w)))
    itpc_angle = angle.(mean(exp.(1im .* itpc_phases .* w)))
    itpcz = ep_n * itpc_value^2

    return (itpc_value=itpc_value, itpcz=itpcz, itpc_angle=itpc_angle, itpc_phases=itpc_phases)
end

"""
    itpc(obj; channel, t, w)

Calculate ITPC (Inter-Trial-Phase Clustering) at sample number `t` over epochs/trials.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `t::Int64`: time point (sample number) at which ITPC is calculated
- `w::Union{Vector{<:Real}, Nothing}=nothing`: optional vector of epochs/trials weights for wITPC calculation

# Returns

Named tuple containing:
- `itpc_value::Vector{Float64}`: ITPC or wITPC value
- `itpcz::Vector{Float64}`: Rayleigh's ITPC Z value
- `itpc_angle::Vector{Float64}`: ITPC angle
- `itpc_phases::Array{Float64, 2}`: phase difference (channel2 - channel1)
"""
function itpc(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), t::Int64, w::Union{Vector{<:Real}, Nothing}=nothing)

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)
    t < 1 && throw(ArgumentError("t must be ≥ 1."))
    t > epoch_len(obj) && throw(ArgumentError("t must be ≤ $(epoch_len(obj))."))
    ep_n < 2 && throw(ArgumentError("OBJ must contain ≥ 2 epochs."))

    itpc_value = zeros(ch_n)
    itpcz = zeros(ch_n)
    itpc_angle = zeros(ch_n)
    itpc_phases = zeros(ch_n, ep_n)

    Threads.@threads for ch_idx in 1:ch_n
        @inbounds itpc_value[ch_idx], itpcz[ch_idx], itpc_angle[ch_idx], itpc_phases[ch_idx, :] = @views itpc(reshape(obj.data[channel[ch_idx], :, :], 1, :, ep_n), t=t, w=w)
    end
    return (itpc_value=itpc_value, itpcz=itpcz, itpc_angle=itpc_angle, itpc_phases=itpc_phases)
end

"""
    itpc_s(obj; <keyword arguments>)

Calculate spectrogram of ITPC (Inter-Trial-Phase Clustering).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Int64`
- `frq_lim::Tuple{Real, Real}`: frequency bounds for the spectrogram
- `frq_n::Int64`: number of frequencies
- `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies
- `w::Union{Vector{<:Real}, Nothing}=nothing`: optional vector of epochs/trials weights for wITPC calculation

# Returns

Named tuple containing:
- `itpc_spec::Array{Float64, 3}`: spectrogram of ITPC values
- `itpc_z_s::Array{Float64, 3}`: spectrogram ITPCz values
- `itpc_frq::Vector{Float64}`: frequencies list
"""
function itpc_s(obj::NeuroAnalyzer.NEURO; channel::Int64, frq_lim::Tuple{Real, Real}, frq_n::Int64, frq::Symbol=:log, w::Union{Vector{<:Real}, Nothing}=nothing)

    _check_var(frq, [:log, :lin], "frq")
    frq_lim = tuple_order(frq_lim)
    frq_lim[1] < 0 && throw(ArgumentError("Lower frequency bound must be ≥ 0."))
    frq_lim[2] > sr(obj) ÷ 2 && throw(ArgumentError("Upper frequency bound must be ≤ $(sr(obj) ÷ 2)."))
    frq_n < 2 && throw(ArgumentError("frq_n must be ≥ 2."))
    if frq === :log
        frq_lim[1] == 0 && (frq_lim = (0.01, frq_lim[2]))
        frq_lim = (frq_lim[1], frq_lim[2])
        frq_list = logspace(log10(frq_lim[1]), log10(frq_lim[2]), frq_n)
    else
        frq_list = linspace(frq_lim[1], frq_lim[2], frq_n)
    end

    _check_channels(obj, channel)
    ep_n = epoch_n(obj)
    ep_len = epoch_len(obj)
    ep_n < 2 && throw(ArgumentError("OBJ must contain ≥ 2 epochs."))

    itpc_spec = zeros(frq_n, ep_len)
    itpc_z_s = zeros(frq_n, ep_len)

    # initialize progress bar
    progress_bar == true && (p = Progress(frq_n, 1))

    Threads.@threads for frq_idx in 1:frq_n
        # create Morlet wavelet
        kernel = generate_morlet(sr(obj), frq_list[frq_idx], 1, ncyc=10)
        half_kernel = floor(Int64, length(kernel) / 2) + 1
        s_conv = zeros(Float64, 1, ep_len, ep_n)
        # convolute with Morlet wavelet
        @inbounds @simd for ep_idx in 1:ep_n
            s_conv[1, :, ep_idx] = @views DSP.conv(obj.data[channel, :, ep_idx], kernel)[(half_kernel - 1):(end - half_kernel)]
        end
        # calculate ITPC of the convoluted signals
        @inbounds @simd for t_idx in 1:ep_len
            itpc_value, itpc_z, _, _ = itpc(s_conv, t=t_idx, w=w)
            itpc_spec[frq_idx, t_idx] = itpc_value
            itpc_z_s[frq_idx, t_idx] = itpc_z
        end

        # update progress bar
        progress_bar == true && next!(p)
    end

    return (itpc_spec=itpc_spec, itpc_z_s=itpc_z_s, itpc_f=frq_list)
end
