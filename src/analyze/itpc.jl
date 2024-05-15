export itpc
export itpc_spec

"""
    itpc(s; t, w)

Calculate ITPC (Inter-Trial-Phase Clustering) at sample number `t` over epochs.

# Arguments

- `s::AbstractArray`: one channel over epochs
- `t::Int64`: time point (sample number) at which ITPC is calculated
- `w::Union{AbstractVector, Nothing}`: optional vector of epochs/trials weights for wITPC calculation

# Returns

Named tuple containing:
- `itpc_value::Float64`: ITPC value
- `itpcz_value::Float64`: Rayleigh's ITPC Z value
- `itpc_angle::Float64`: ITPC angle
- `itpc_phases::Vector{Float64}`: phases at time `t` averaged across trials/epochs
"""
function itpc(s::AbstractArray; t::Int64, w::Union{AbstractVector, Nothing}=nothing)

    @assert t >= 1 "t must be ≥ 1."
    @assert t <= size(s, 2) "t must be ≤ $(size(s, 2))."
    @assert size(s, 1) == 1 "s must have 1 channel."

    ep_n = size(s, 3)

    w === nothing && (w = ones(ep_n))
    # scale w if w contains negative values
    any(i -> i < 0, w) && (w .+= abs(minimum(w)))
    @assert length(w) == ep_n "Length of w should be equal to number of epochs ($ep_n)."
    
    s_phase = zeros(size(s, 2), ep_n)
    @inbounds for ep_idx in 1:ep_n
        _, _, _, s_phase[:, ep_idx] = @views hspectrum(s[1, :, ep_idx])
    end
 
    itpc_phases = @view s_phase[t, :]
    itpc_value = abs.(mean(exp.(1im .* itpc_phases .* w)))
    itpc_angle = angle.(mean(exp.(1im .* itpc_phases .* w)))
    itpcz_value = ep_n * itpc_value^2

    return (itpc_value=itpc_value, itpcz_value=itpcz_value, itpc_angle=itpc_angle, itpc_phases=itpc_phases)

end

"""
    itpc(obj; <keyword arguments>)

Calculate ITPC (Inter-Trial-Phase Clustering) at sample number `t` over epochs.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `t::Int64`: time point (sample number) at which ITPC is calculated
- `w::Union{Vector{<:Real}, Nothing}=nothing`: optional vector of epochs/trials weights for wITPC calculation

# Returns

Named tuple containing:
- `itpc_value::Vector{Float64}`: ITPC or wITPC value
- `itpcz_value::Vector{Float64}`: Rayleigh's ITPC Z value
- `itpc_angle::Vector{Float64}`: ITPC angle
- `itpc_phases::Array{Float64, 2}`: phase difference (channel2 - channel1)
"""
function itpc(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), t::Int64, w::Union{Vector{<:Real}, Nothing}=nothing)

    _check_channels(obj, ch)
    ch_n = length(ch)
    ep_n = nepochs(obj)
    @assert t >= 1 "t must be ≥ 1."
    @assert t <= epoch_len(obj) "t must be ≤ $(epoch_len(obj))."
    @assert ep_n >= 2 "OBJ must contain ≥ 2 epochs."

    itpc_value = zeros(ch_n)
    itpcz_value = zeros(ch_n)
    itpc_angle = zeros(ch_n)
    itpc_phases = zeros(ch_n, ep_n)

    Threads.@threads for ch_idx in 1:ch_n
        @inbounds itpc_value[ch_idx], itpcz_value[ch_idx], itpc_angle[ch_idx], itpc_phases[ch_idx, :] = @views itpc(reshape(obj.data[ch[ch_idx], :, :], 1, :, ep_n), t=t, w=w)
    end

    return (itpc_value=itpc_value, itpcz_value=itpcz_value, itpc_angle=itpc_angle, itpc_phases=itpc_phases)

end

"""
    itpc_spec(s; w)

Calculate spectrogram of ITPC (Inter-Trial-Phase Clustering).

# Arguments

- `s::AbstractArray`: one channel over epochs
- `w::Union{AbstractVector, Nothing}`: optional vector of epochs/trials weights for wITPC calculation

# Returns

Named tuple containing:
- `itpc_values::Vector{Float64}`: ITPC values
- `itpcz_values::Vector{Float64}`: Rayleigh's ITPC Z values
- `itpc_angles::Vector{Float64}`: ITPC angles
- `itpc_phases::Matrix{Float64}`: phases at time `t` averaged across trials/epochs
"""
function itpc_spec(s::AbstractArray; w::Union{AbstractVector, Nothing}=nothing)

    @assert size(s, 1) == 1 "s must have 1 channel."

    ep_n = size(s, 3)

    w === nothing && (w = ones(ep_n))
    # scale w if w contains negative values
    any(i -> i < 0, w) && (w .+= abs(minimum(w)))
    @assert length(w) == ep_n "Length of w should be equal to number of epochs ($ep_n)."
    
    itpc_phases = zeros(size(s, 2), ep_n)
    itpc_values = zeros(size(s, 2))
    itpc_angles = zeros(size(s, 2))
    itpcz_values = zeros(size(s, 2))

    @inbounds for ep_idx in 1:ep_n
        _, _, _, itpc_phases[:, ep_idx] = @views hspectrum(s[1, :, ep_idx])
    end

    for idx in 1:size(itpc_phases, 1)
        itpc_values[idx] = @views abs.(mean(exp.(1im .* itpc_phases[idx, :] .* w)))
        itpc_angles[idx] = @views angle.(mean(exp.(1im .* itpc_phases[idx, :] .* w)))
        itpcz_values[idx] = @views ep_n * itpc_values[idx]^2
    end

    return (itpc_values=itpc_values, itpcz_values=itpcz_values, itpc_angles=itpc_angles, itpc_phases=itpc_phases)

end

"""
    itpc_spec(obj; <keyword arguments>)

Calculate spectrogram of ITPC (Inter-Trial-Phase Clustering).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Int64`
- `frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2)`: frequency bounds for the spectrogram
- `frq_n::Int64=_tlength(frq_lim)`: number of frequencies
- `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
- `w::Union{Vector{<:Real}, Nothing}=nothing`: optional vector of epochs/trials weights for wITPC calculation

# Returns

Named tuple containing:
- `itpc_s::Array{Float64, 3}`: spectrogram of ITPC values
- `itpcz_s::Array{Float64, 3}`: spectrogram itpcz_value values
- `itpc_f::Vector{Float64}`: frequencies list
"""
function itpc_spec(obj::NeuroAnalyzer.NEURO; ch::Int64, frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2), frq_n::Int64=_tlength(frq_lim), frq::Symbol=:log, w::Union{Vector{<:Real}, Nothing}=nothing)

    _check_var(frq, [:log, :lin], "frq")
    _check_tuple(frq_lim, "frq_lim", (0, sr(obj) / 2))
    @assert frq_n >= 2 "frq_n must be ≥ 2."
    if frq === :log
        frq_lim = frq_lim[1] == 0 ? (0.01, frq_lim[2]) : (frq_lim[1], frq_lim[2])
        frq_lim = (frq_lim[1], frq_lim[2])
        frq_list = round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), frq_n), digits=3)
    else
        frq_list = linspace(frq_lim[1], frq_lim[2], frq_n)
    end

    _check_channels(obj, ch)
    ep_n = nepochs(obj)
    ep_len = epoch_len(obj)
    @assert ep_n >= 2 "OBJ must contain ≥ 2 epochs."

    itpc_s = zeros(frq_n, ep_len)
    itpcz_s = zeros(frq_n, ep_len)

    # initialize progress bar
    progress_bar && (progbar = Progress(frq_n, dt=1, barlen=20, color=:white))

    Threads.@threads for frq_idx in 1:frq_n
        # create Morlet wavelet
        kernel = generate_morlet(sr(obj), frq_list[frq_idx], 1, ncyc=10)
        half_kernel = floor(Int64, length(kernel) / 2) + 1
        s_conv = zeros(Float64, 1, ep_len, ep_n)
        # convolute with Morlet wavelet
        @inbounds for ep_idx in 1:ep_n
            s_conv[1, :, ep_idx] = @views DSP.conv(obj.data[ch, :, ep_idx], kernel)[(half_kernel - 1):(end - half_kernel)]
        end
        # calculate ITPC of the convoluted signals
        itpc_s[frq_idx, :], itpcz_s[frq_idx, :], _, _ = itpc_spec(s_conv, w=w)

        # update progress bar
        progress_bar && next!(progbar)
    end

    return (itpc_s=itpc_s, itpcz_s=itpcz_s, itpc_f=frq_list)

end
