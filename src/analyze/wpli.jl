export wpli

"""
    _wpli(s1, s2; <keyword arguments>)

Calculate weighted PLI (Phase Lag Index) for two 1-D signal vectors (Vinck et al. 2011).

The imaginary part of the cross-spectrum is calculated at every time point from the analytic signals (Hilbert transform): `X = Im(z₁ · conj(z₂))`; the wPLI is the asymmetry of its distribution around zero, weighted by its magnitude:

`wPLI = |Σₜ X| / Σₜ |X|`

and the debiased squared wPLI: `(( Σₜ X)² − Σₜ X²) / ((Σₜ |X|)² − Σₜ X²)`.

As PLI, it is insensitive to zero-lag (volume-conducted) coupling, and less sensitive to noise near zero lag, because small imaginary components get small weights. The signals should be narrowband (band-pass filtered).

# Arguments

- `s1::AbstractVector`: signal vector
- `s2::AbstractVector`: signal vector
- `debiased::Bool=false`: if `true`, calculate debiased wPLI

# Returns

Named tuple:

- `pv::Float64`: wPLI value
- `sd::Vector{Float64}`: signal difference (s1 - s2)
- `phd::Vector{Float64}`: phase difference (s1 - s2)
- `s1ph::Vector{Float64}`: signal 1 phase
- `s2ph::Vector{Float64}`: signal 2 phase
"""
function _wpli(
    s1::AbstractVector,
    s2::AbstractVector;
    debiased::Bool = false,
)::@NamedTuple{
    pv::Float64,
    sd::Vector{Float64},
    phd::Vector{Float64},
    s1ph::Vector{Float64},
    s2ph::Vector{Float64},
}
    # validate
    length(s1) == length(s2) ||
        throw(ArgumentError("Both signals must have the same length."))

    pli_data = _pli(s1, s2)
    sd = pli_data.sd
    phd = pli_data.phd
    s1ph = pli_data.s1ph
    s2ph = pli_data.s2ph

    # imaginary part of the cross-spectrum at every time point (analytic signals)
    im_xy = imag.(htransform(s1).c .* conj.(htransform(s2).c))
    # wPLI = |Σ Im| / Σ |Im|; debiased squared wPLI (Vinck et al. 2011)
    num = sum(im_xy)
    denom = sum(abs, im_xy)
    sum_sq = sum(abs2, im_xy)
    if denom == 0
        pv = 0.0
    elseif debiased
        pv = (num^2 - sum_sq) / (denom^2 - sum_sq)
    else
        pv = abs(num) / denom
    end

    return (; pv, sd, phd, s1ph, s2ph)
end

"""
    _wpli(s1, s2; <keyword arguments>)

Calculate weighted PLI (Phase Lag Index) across epochs (trials) for two signal matrices (Vinck et al. 2011).

For each epoch, the imaginary part of the cross-spectrum is averaged over time (analytic signals): `Xₑ = meanₜ Im(z₁ · conj(z₂))`; the wPLI is calculated across epochs: `|Σₑ Xₑ| / Σₑ |Xₑ|` (debiased: `((Σₑ Xₑ)² − Σₑ Xₑ²) / ((Σₑ |Xₑ|)² − Σₑ Xₑ²)`). This is the standard definition, unbiased (debiased version) for independent epochs.

# Arguments

- `s1::AbstractMatrix`: signal matrix, shape (samples, epochs)
- `s2::AbstractMatrix`: signal matrix, shape (samples, epochs)
- `debiased::Bool=false`: if `true`, calculate debiased wPLI

# Returns

- `Float64`: wPLI value
"""
function _wpli(s1::AbstractMatrix, s2::AbstractMatrix; debiased::Bool = false)::Float64
    # validate
    size(s1) == size(s2) || throw(ArgumentError("Both signals must have the same size."))
    size(s1, 2) >= 2 || throw(ArgumentError("At least 2 epochs are required."))

    # imaginary part of the cross-spectrum of each epoch (mean over time)
    x = [mean(imag.(htransform(@view(s1[:, ep])).c .* conj.(htransform(@view(s2[:, ep])).c)))
         for ep in axes(s1, 2)]
    num = sum(x)
    denom = sum(abs, x)
    sum_sq = sum(abs2, x)
    denom == 0 && return 0.0

    return debiased ? (num^2 - sum_sq) / (denom^2 - sum_sq) : abs(num) / denom
end

"""
    wpli(obj1, obj2; <keyword arguments>)

Calculate weighted PLI (Phase Lag Index) for two NEURO objects.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: input NEURO object
- `obj2::NeuroAnalyzer.NEURO`: input NEURO object
- `ch1::Union{String, Vector{String}, Regex}`: channel name(s) in `obj1`
- `ch2::Union{String, Vector{String}, Regex}`: channel name(s) in `obj2`
- `ep1::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}}=_c(nepochs(obj1))`: epoch number(s) in `obj1`
- `ep2::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}}=_c(nepochs(obj2))`: epoch number(s) in `obj2`
- `debiased::Bool=false`: if `true`, calculate debiased wPLI
- `across_epochs::Bool=false`: if `true`, calculate wPLI across the selected epochs (the standard definition, Vinck et al. 2011; `pv` has one column); otherwise for each epoch separately (across its time points)

# Returns

Named tuple:

- `pv::Matrix{Float64}`: wPLI values, shape (channels, epochs) or (channels, 1) if `across_epochs = true`

# Notes

wPLI calculated within single epochs is biased upwards for short epochs of narrowband signals (few independent time points), even for zero-lag coupling; with several epochs, `across_epochs = true` (optionally with `debiased = true`) is recommended.
- `sd::Array{Float64, 3}`: signal difference (s1 - s2)
- `phd::Array{Float64, 3}`: phase difference (s1 - s2)
- `s1ph::Array{Float64, 3}`: signal 1 phase
- `s2ph::Array{Float64, 3}`: signal 2 phase
"""
function wpli(
    obj1::NeuroAnalyzer.NEURO,
    obj2::NeuroAnalyzer.NEURO;
    ch1::Union{String, Vector{String}, Regex},
    ch2::Union{String, Vector{String}, Regex},
    ep1::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}} = _c(nepochs(obj1)),
    ep2::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}} = _c(nepochs(obj2)),
    debiased::Bool = false,
    across_epochs::Bool = false,
)::@NamedTuple{
    pv::Matrix{Float64},
    sd::Array{Float64, 3},
    phd::Array{Float64, 3},
    s1ph::Array{Float64, 3},
    s2ph::Array{Float64, 3},
}
    # resolve channel names to integer indices, optionally skipping bad channels
    ch1 =
        exclude_bads ? get_channel(obj1; ch = ch1, exclude = "bad") :
        get_channel(obj1; ch = ch1, exclude = "")
    ch2 =
        exclude_bads ? get_channel(obj2; ch = ch2, exclude = "bad") :
        get_channel(obj2; ch = ch2, exclude = "")
    isempty(ch1) && throw(ArgumentError("No channels selected."))
    isempty(ch2) && throw(ArgumentError("No channels selected."))
    length(ch1) == length(ch2) ||
        throw(
            ArgumentError(
                "Lengths of ch1 ($(length(ch1))) and ch2 ($(length(ch2))) must be equal.",
            ),
        )

    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    length(ep1) == length(ep2) ||
        throw(
            ArgumentError(
                "Lengths of ep1 ($(length(ep1))) and ep2 ($(length(ep2))) must be equal.",
            ),
        )
    epoch_len(obj1) == epoch_len(obj2) ||
        throw(ArgumentError("OBJ1 and OBJ2 must have the same epoch lengths."))
    ep1 = _n2v(ep1)
    ep2 = _n2v(ep2)

    ch_n = length(ch1)
    ep_n = length(ep1)

    pv = zeros(ch_n, ep_n)
    sd = zeros(ch_n, epoch_len(obj1), ep_n)
    phd = zeros(ch_n, epoch_len(obj1), ep_n)
    s1ph = zeros(ch_n, epoch_len(obj1), ep_n)
    s2ph = zeros(ch_n, epoch_len(obj1), ep_n)

    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        wpli_data = _wpli(
            @view(obj1.data[ch1[ch_idx], :, ep1[ep_idx]]),
            @view(obj2.data[ch2[ch_idx], :, ep2[ep_idx]]),
            debiased = debiased,
        )
        pv[ch_idx, ep_idx] = wpli_data.pv
        sd[ch_idx, :, ep_idx] = wpli_data.sd
        phd[ch_idx, :, ep_idx] = wpli_data.phd
        s1ph[ch_idx, :, ep_idx] = wpli_data.s1ph
        s2ph[ch_idx, :, ep_idx] = wpli_data.s2ph
    end

    # wPLI across epochs
    if across_epochs
        ep_n >= 2 || throw(ArgumentError("across_epochs requires at least 2 epochs."))
        pv = zeros(ch_n, 1)
        for ch_idx in 1:ch_n
            pv[ch_idx, 1] = _wpli(@view(obj1.data[ch1[ch_idx], :, ep1]),
                                  @view(obj2.data[ch2[ch_idx], :, ep2]);
                                  debiased = debiased)
        end
    end

    return (; pv, sd, phd, s1ph, s2ph)
end

"""
    wpli(obj; <keyword arguments>)

Calculate weighted PLI (Phase Lag Index) for a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `debiased::Bool=false`: if `true`, calculate debiased wPLI

# Returns

- `Array{Float64, 3}`: wPLI value
"""
function wpli(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    debiased::Bool = false,
)::Array{Float64, 3}
    # resolve channel names to integer indices, optionally skipping bad channels
    ch =
        exclude_bads ?
        get_channel(obj; ch = ch, exclude = "bad") :
        get_channel(obj; ch = ch, exclude = "")
    isempty(ch) && throw(ArgumentError("No channels selected."))

    # number of channels
    ch_n = length(ch)
    # number of epochs
    ep_n = nepochs(obj)

    # pre-allocate output
    pv = zeros(ch_n, ch_n, ep_n)

    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx1, ep_idx = idx[1], idx[2]
        for ch_idx2 = 1:ch_idx1
            pv[ch_idx1, ch_idx2, ep_idx] = _wpli(
                @view(obj.data[ch[ch_idx1], :, ep_idx]),
                @view(obj.data[ch[ch_idx2], :, ep_idx]),
                debiased = debiased,
            ).pv
        end
    end

    # mirror the lower triangle to the upper triangle to produce the full symmetric matrix
    return _copy_lt2ut(pv)
end
