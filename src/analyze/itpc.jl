export itpc
export itpc_spec

"""
    itpc(s; <keyword arguments>)

Calculate ITPC (Inter-Trial-Phase Clustering) at sample number `t` over epochs.

# Arguments

- `s::AbstractArray`: one channel over epochs
- `t::Int64`: time point (sample number) at which ITPC is calculated
- `w::Union{AbstractVector, Nothing}`: optional vector of epochs/trials weights for wITPC calculation

# Returns

Named tuple containing:
- `itpc_val::Float64`: ITPC value
- `itpcz_val::Float64`: Rayleigh's ITPC z value
- `itpc_ang::Float64`: ITPC angle
- `itpc_ph::Vector{Float64}`: phases at time `t` averaged across trials/epochs

# Source

1. Cohen, M. X. (2014). Analyzing Neural Time Series Data: Theory and Practice.Cambridge: MIT Press
"""
function itpc(s::AbstractArray; t::Int64, w::Union{AbstractVector, Nothing}=nothing)::@NamedTuple{itpc_val::Float64, itpcz_val::Float64, itpc_ang::Float64, itpc_ph::Vector{Float64}}

    _chk3d(s)
    @assert t >= 1 "t must be ≥ 1."
    @assert t <= size(s, 2) "t must be ≤ $(size(s, 2))."
    @assert size(s, 1) == 1 "s must have 1 channel."

    ep_n = size(s, 3)

    w === nothing && (w = ones(ep_n))
    # scale w if w contains negative values
    any(i -> i < 0, w) && (w .+= abs(minimum(w)))
    @assert length(w) == ep_n "Length of w ($(length(w))) and number of epochs ($ep_n) must be equal."

    s_phase = zeros(size(s, 2), ep_n)
    @inbounds for ep_idx in 1:ep_n
        _, _, _, s_phase[:, ep_idx] = @views htransform(s[1, :, ep_idx])
    end

    itpc_ph = @view s_phase[t, :]
    itpc_val = abs.(mean(exp.(1im .* itpc_ph .* w)))
    itpc_ang = DSP.angle.(mean(exp.(1im .* itpc_ph .* w)))
    itpcz_val = ep_n * itpc_val^2

    return (itpc_val=itpc_val, itpcz_val=itpcz_val, itpc_ang=itpc_ang, itpc_ph=itpc_ph)

end

"""
    itpc(obj; <keyword arguments>)

Calculate ITPC (Inter-Trial-Phase Clustering) at sample number `t` over epochs.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `t::Int64`: time point (sample number) at which ITPC is calculated
- `w::Union{Vector{<:Real}, Nothing}=nothing`: optional vector of epochs/trials weights for wITPC calculation

# Returns

Named tuple containing:
- `itpc_val::Vector{Float64}`: ITPC or wITPC value
- `itpcz_val::Vector{Float64}`: Rayleigh's ITPC z value
- `itpc_ang::Vector{Float64}`: ITPC angle
- `itpc_ph::Matrix{Float64}`: phase difference (channel2 - channel1)
"""
function itpc(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, t::Int64, w::Union{Vector{<:Real}, Nothing}=nothing)::@NamedTuple{itpc_val::Vector{Float64}, itpcz_val::Vector{Float64}, itpc_ang::Vector{Float64}, itpc_ph::Matrix{Float64}}

    ch = exclude_bads ? get_channel(obj, ch=ch, exclude="bad") : get_channel(obj, ch=ch, exclude="")
    ch_n = length(ch)
    ep_n = nepochs(obj)
    @assert t >= 1 "t must be ≥ 1."
    @assert t <= epoch_len(obj) "t must be ≤ $(epoch_len(obj))."
    @assert ep_n >= 2 "OBJ must contain ≥ 2 epochs."

    itpc_val = zeros(ch_n)
    itpcz_val = zeros(ch_n)
    itpc_ang = zeros(ch_n)
    itpc_ph = zeros(ch_n, ep_n)

    Threads.@threads :greedy for ch_idx in 1:ch_n
        @inbounds itpc_val[ch_idx], itpcz_val[ch_idx], itpc_ang[ch_idx], itpc_ph[ch_idx, :] = @views itpc(reshape(obj.data[ch[ch_idx], :, :], 1, :, ep_n), t=t, w=w)
    end

    return (itpc_val=itpc_val, itpcz_val=itpcz_val, itpc_ang=itpc_ang, itpc_ph=itpc_ph)

end

"""
    itpc_spec(s; <keyword arguments>)

Calculate spectrogram of ITPC (Inter-Trial-Phase Clustering).

# Arguments

- `s::AbstractArray`: one channel over epochs
- `w::Union{AbstractVector, Nothing}`: optional vector of epochs/trials weights for wITPC calculation

# Returns

Named tuple containing:
- `itpc_val::Vector{Float64}`: ITPC values
- `itpcz_val::Vector{Float64}`: Rayleigh's ITPC z values
- `itpc_ang::Vector{Float64}`: ITPC angles
- `itpc_ph::Matrix{Float64}`: phases at time `t` averaged across trials/epochs
"""
function itpc_spec(s::AbstractArray; w::Union{AbstractVector, Nothing}=nothing)::@NamedTuple{itpc_val::Vector{Float64}, itpcz_val::Vector{Float64}, itpc_ang::Vector{Float64}, itpc_ph::Matrix{Float64}}

    _chk3d(s)
    @assert size(s, 1) == 1 "s must have 1 channel."

    ep_n = size(s, 3)

    w === nothing && (w = ones(ep_n))
    # scale w if w contains negative values
    any(i -> i < 0, w) && (w .+= abs(minimum(w)))
    @assert length(w) == ep_n "Length of w ($(length(w))) and number of epochs ($ep_n) must be equal."

    itpc_ph = zeros(size(s, 2), ep_n)
    itpc_val = zeros(size(s, 2))
    itpc_ang = zeros(size(s, 2))
    itpcz_val = zeros(size(s, 2))

    @inbounds for ep_idx in 1:ep_n
        _, _, _, itpc_ph[:, ep_idx] = @views htransform(s[1, :, ep_idx])
    end

    for idx in axes(itpc_ph, 1)
        itpc_val[idx] = @views abs.(mean(exp.(1im .* itpc_ph[idx, :] .* w)))
        itpc_ang[idx] = @views DSP.angle.(mean(exp.(1im .* itpc_ph[idx, :] .* w)))
        itpcz_val[idx] = @views ep_n * itpc_val[idx]^2
    end

    return (itpc_val=itpc_val, itpcz_val=itpcz_val, itpc_ang=itpc_ang, itpc_ph=itpc_ph)

end

"""
    itpc_spec(obj; <keyword arguments>)

Calculate spectrogram of ITPC (Inter-Trial-Phase Clustering).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::String`: channel to analyze
- `frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2)`: frequency bounds for the spectrogram
- `frq_n::Int64=_tlength(frq_lim)`: number of frequencies
- `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
- `w::Union{Vector{<:Real}, Nothing}=nothing`: optional vector of epochs/trials weights for wITPC calculation

# Returns

Named tuple containing:
- `itpc_s::Matrix{Float64}`: spectrogram of ITPC values
- `itpcz_s::Matrix{Float64}`: spectrogram of ITPCZ values
- `itpc_f::Vector{Float64}`: frequencies list
"""
function itpc_spec(obj::NeuroAnalyzer.NEURO; ch::String, frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2), frq_n::Int64=_tlength(frq_lim), frq::Symbol=:log, w::Union{Vector{<:Real}, Nothing}=nothing)::@NamedTuple{itpc_s::Matrix{Float64}, itpcz_s::Matrix{Float64}, itpc_f::Vector{Float64}}

    _check_var(frq, [:log, :lin], "frq")
    _check_tuple(frq_lim, "frq_lim", (0, sr(obj) / 2))
    @assert frq_n >= 2 "frq_n must be ≥ 2."
    if frq === :log
        frq_lim = frq_lim[1] == 0 ? (0.01, frq_lim[2]) : (frq_lim[1], frq_lim[2])
        frq_lim = (frq_lim[1], frq_lim[2])
        frq_list = round.(log10space(log10(frq_lim[1]), log10(frq_lim[2]), frq_n), digits=3)
    else
        frq_list = linspace(frq_lim[1], frq_lim[2], frq_n)
    end

    ch = exclude_bads ? get_channel(obj, ch=ch, exclude="bad") : get_channel(obj, ch=ch, exclude="")[1]
    ep_n = nepochs(obj)
    ep_len = epoch_len(obj)
    @assert ep_n >= 2 "OBJ must contain ≥ 2 epochs."

    itpc_s = zeros(frq_n, ep_len)
    itpcz_s = zeros(frq_n, ep_len)

    # initialize progress bar
    progress_bar && (progbar = Progress(frq_n, dt=1, barlen=20, color=:white))

    Threads.@threads :greedy for frq_idx in 1:frq_n
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
