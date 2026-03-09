export itpc
export itpc_spec

"""
    itpc(s; <keyword arguments>)

Calculate ITPC (Inter-Trial Phase Clustering) at sample `t` over epochs. ITPC measures the consistency of the instantaneous phase across epochs at a given time point (or across all time points for the spectrogram variant):

ITPC value = |mean( w .* exp(i·φ) )| (0 = random, 1 = perfect lock)
ITPC z = N · ITPC² (Rayleigh statistic)
ITPC angle = angle( mean( w .* exp(i·φ) ) ) (preferred phase)

The weighted variant (wITPC) allows per-epoch importance weights.

# Arguments

- `s::AbstractArray`: one channel over epochs
- `t::Int64`: time point (sample number) at which ITPC is calculated
- `w::Union{AbstractVector, Nothing}`: optional vector of epochs/trials weights for wITPC calculation

# Returns

Named tuple containing:

- `itpcv::Float64`: ITPC (or wITPC) value
- `itpcz::Float64`: Rayleigh's ITPC z statistic
- `itpca::Float64`: ITPC angle (preferred phase)
- `itpcph::Vector{Float64}`: instantaneous phases at time `t` across epochs

# Reference

Cohen, M. X. (2014). Analyzing Neural Time Series Data: Theory and Practice.Cambridge: MIT Press
"""
function itpc(
    s::AbstractArray;
    t::Int64,
    w::Union{AbstractVector, Nothing} = nothing
)::@NamedTuple{
    itpcv::Float64,
    itpcz::Float64,
    itpca::Float64,
    itpcph::Vector{Float64}
}

    # validate that the input is a proper 3-D array (channels × samples × epochs)
    _chk3d(s)

    @assert t >= 1 "t must be ≥ 1."
    @assert t <= size(s, 2) "t must be ≤ $(size(s, 2))."
    @assert size(s, 1) == 1 "s must have 1 channel."

    # number of epochs
    ep_n = size(s, 3)

    isnothing(w) && (w = ones(ep_n))
    # scale w if w contains negative values
    any(i -> i < 0, w) && (w .+= abs(minimum(w)))
    @assert length(w) == ep_n "Length of w ($(length(w))) and number of epochs ($ep_n) must be equal."

    # compute instantaneous phase for every epoch
    s_phase = zeros(size(s, 2), ep_n)
    @inbounds for ep_idx in 1:ep_n
        s_phase[:, ep_idx] = htransform(@view(s[1, :, ep_idx])).ph
    end

    itpcph = s_phase[t, :]

    # compute mean weighted complex phasor once; reuse for val, ang, and z_val
    mphasor = mean(Base.cis.(itpcph) .* w)

    itpcv  = abs(mphasor)
    itpca  = DSP.angle(mphasor)
    itpcz = ep_n * itpcv^2

    return (itpcv = itpcv, itpcz = itpcz, itpca = itpca, itpcph = itpcph)

end

"""
    itpc(obj; <keyword arguments>)

Calculate ITPC (Inter-Trial Phase Clustering) at time `t` over epochs. ITPC measures the consistency of the instantaneous phase across epochs at a given time point (or across all time points for the spectrogram variant):

ITPC value = |mean( w .* exp(i·φ) )| (0 = random, 1 = perfect lock)
ITPC z = N · ITPC² (Rayleigh statistic)
ITPC angle = angle( mean( w .* exp(i·φ) ) ) (preferred phase)

The weighted variant (wITPC) allows per-epoch importance weights.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `t::Real`: time point in seconds at which ITPC is calculated
- `w::Union{Vector{<:Real}, Nothing}=nothing`: optional per-epoch weights

# Returns

Named tuple containing:

- `itpcv::Vector{Float64}`: ITPC or wITPC value per channel
- `itpcz::Vector{Float64}`: Rayleigh's ITPC z statistic per channel
- `itpca::Vector{Float64}`: ITPC angle per channel
- `itpcph::Matrix{Float64}`: instantaneous phases at `t`, shape `(channels, epochs)`
"""
function itpc(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    t::Real,
    w::Union{Vector{<:Real}, Nothing} = nothing
)::@NamedTuple{
    itpcv::Vector{Float64},
    itpcz::Vector{Float64},
    itpca::Vector{Float64},
    itpcph::Matrix{Float64}
}

    # number of epochs
    ep_n = nepochs(obj)
    @assert ep_n >= 2 "OBJ must contain ≥ 2 epochs."

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")
    ch_n = length(ch)

    # get time point index
    t_idx = vsearch(t, obj.epoch_time)

    # pre-allocate outputs
    itpcv = zeros(ch_n)
    itpcz = zeros(ch_n)
    itpca = zeros(ch_n)
    itpcph = zeros(ch_n, ep_n)

    Threads.@threads :dynamic for ch_idx in 1:ch_n
        @inbounds begin
            itpc_data = itpc(
                reshape(
                    @view(obj.data[ch[ch_idx], :, :]), 1, size(obj.data, 2), ep_n
                ),
                t = t_idx,
                w = w
            )
            itpcv[ch_idx] = itpc_data.itpcv
            itpcz[ch_idx] = itpc_data.itpcz
            itpca[ch_idx] = itpc_data.itpca
            itpcph[ch_idx, :] = itpc_data.itpcph
        end
    end

    return (itpcv = itpcv, itpcz = itpcz, itpca = itpca, itpcph = itpcph)

end

"""
    itpc_spec(s; <keyword arguments>)

Calculate the ITPC spectrogram (ITPC at every time point). ITPC measures the consistency of the instantaneous phase across epochs at a given time point (or across all time points for the spectrogram variant):

ITPC value = |mean( w .* exp(i·φ) )| (0 = random, 1 = perfect lock)
ITPC z = N · ITPC² (Rayleigh statistic)
ITPC angle = angle( mean( w .* exp(i·φ) ) ) (preferred phase)

The weighted variant (wITPC) allows per-epoch importance weights.

# Arguments

- `s::AbstractArray`: single-channel array (1 × samples × epochs)
- `w::Union{AbstractVector, Nothing}`: optional per-epoch weights

# Returns

Named tuple containing:

- `itpcv::Vector{Float64}`: ITPC or wITPC value per channel
- `itpcz::Vector{Float64}`: Rayleigh's ITPC z statistic per channel
- `itpca::Vector{Float64}`: ITPC angle per channel
- `itpcph::Matrix{Float64}`: instantaneous phases at `t`, shape `(channels, epochs)`
"""
function itpc_spec(
    s::AbstractArray;
    w::Union{AbstractVector, Nothing} = nothing
)::@NamedTuple{
    itpcv::Vector{Float64},
    itpcz::Vector{Float64},
    itpca::Vector{Float64},
    itpcph::Matrix{Float64}
}

    # validate that the input is a proper 3-D array (channels × samples × epochs)
    _chk3d(s)
    @assert size(s, 1) == 1 "s must have 1 channel."

    # number of epochs
    ep_n = size(s, 3)

    w === nothing && (w = ones(ep_n))
    # scale w if w contains negative values
    any(i -> i < 0, w) && (w .+= abs(minimum(w)))
    @assert length(w) == ep_n "Length of w ($(length(w))) and number of epochs ($ep_n) must be equal."

    # pre-allocate outputs
    itpcph = zeros(size(s, 2), ep_n)
    itpcv = zeros(size(s, 2))
    itpca = zeros(size(s, 2))
    itpcz = zeros(size(s, 2))

    @inbounds for ep_idx in 1:ep_n
        itpcph[:, ep_idx] = htransform(@view(s[1, :, ep_idx])).ph
    end

    for idx in axes(itpcph, 1)
        # hoist the mean phasor so it is computed once per time point, not twice
        mphasor = mean(Base.cis.(@view(itpcph[idx, :])) .* w)
        itpcv[idx]  = abs(mphasor)
        itpca[idx]  = DSP.angle(mphasor)
        itpcz[idx] = ep_n * itpcv[idx]^2
    end

    return (itpcv = itpcv, itpcz = itpcz, itpca = itpca, itpcph = itpcph)

end

"""
    itpc_spec(obj; <keyword arguments>)

Calculate the ITPC spectrogram (ITPC at each frequency and time point). ITPC measures the consistency of the instantaneous phase across epochs at a given time point (or across all time points for the spectrogram variant):

ITPC value = |mean( w .* exp(i·φ) )| (0 = random, 1 = perfect lock)
ITPC z = N · ITPC² (Rayleigh statistic)
ITPC angle = angle( mean( w .* exp(i·φ) ) ) (preferred phase)

The weighted variant (wITPC) allows per-epoch importance weights.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::String`: channel to analyze
- `flim::Tuple{Real, Real}=(0, sr(obj) / 2)`: frequency bounds for the spectrogram
- `nfrq::Int64=_tlength(flim)`: number of frequencies
- `frq::Symbol=:log`: frequency scaling - `:lin` or `:log`
- `w::Union{Vector{<:Real}, Nothing}=nothing`: optional vector of epochs/trials weights for wITPC calculation

# Returns

Named tuple containing:

- `itpcs::Matrix{Float64}`: ITPC spectrogram, shape `(frequencies, samples)`
- `itpczs::Matrix{Float64}`: ITPCZ spectrogram, shape `(frequencies, samples)`
- `f::Vector{Float64}`: frequency vector
"""
function itpc_spec(
    obj::NeuroAnalyzer.NEURO;
    ch::String,
    flim::Tuple{Real, Real} = (0, sr(obj) / 2),
    nfrq::Int64 = _tlength(flim),
    frq::Symbol = :log,
    w::Union{Vector{<:Real}, Nothing} = nothing
)::@NamedTuple{
    itpcs::Matrix{Float64},
    itpczs::Matrix{Float64},
    f::Vector{Float64}
}

    _check_var(frq, [:log, :lin], "frq")
    _check_tuple(flim, (0, sr(obj) / 2), "flim")
    @assert nfrq >= 2 "nfrq must be ≥ 2."

    # build frequency vector; log scale requires a strictly positive lower bound
    if frq === :log
        @assert flim[1] > 0 "For :log scale, lower flim bound must be > 0 Hz."
        f = round.(logspace(flim[1], flim[2], nfrq), digits = 3)
    else
        f = linspace(flim[1], flim[2], nfrq)
    end

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")

    # number of channels
    ch_n = length(ch)
    # number of epochs
    ep_n = nepochs(obj)
    # epoch length
    ep_len = epoch_len(obj)
    @assert ep_n >= 2 "OBJ must contain ≥ 2 epochs."

    # pre-allocate outputs
    itpcs = zeros(nfrq, ep_len)
    itpczs = zeros(nfrq, ep_len)

    # initialize progress bar
    progbar = Progress(nfrq, dt = 1, barlen = 20, color = :white, enabled = progress_bar)

    Threads.@threads :dynamic for frq_idx in 1:nfrq

        # build Morlet wavelet and compute half-kernel offset for trimming
        kernel = generate_morlet(sr(obj), f[frq_idx], 1, ncyc = 10)
        half_kernel = floor(Int64, length(kernel) / 2) + 1

        # convolve each epoch with the Morlet kernel
        s_conv = zeros(Float64, 1, ep_len, ep_n)
        @inbounds for ep_idx in 1:ep_n
            s_conv[1, :, ep_idx] = DSP.conv(
                @view(obj.data[ch, :, ep_idx]),
                kernel
            )[(half_kernel - 1):(end - half_kernel)]
        end

        # compute the ITPC spectrogram for the convolved signal
        itpc_data = itpc_spec(s_conv, w = w)
        itpcs[frq_idx, :]  = itpc_data.itpcv
        itpczs[frq_idx, :] = itpc_data.itpcz

        progress_bar && next!(progbar)

    end

    return (itpcs = itpcs, itpczs = itpczs, f = f)

end
