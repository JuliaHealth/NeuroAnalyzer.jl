export denoise_wien
export denoise_wien!

"""
    _smooth_bins(p, h)

Moving average of a spectrum over `2h + 1` frequency bins (fewer at the edges).

# Arguments

- `p::AbstractVector`: power spectrum
- `h::Int64`: half-width in bins; `0` returns a copy of `p`

# Returns

- `Vector{Float64}`
"""
function _smooth_bins(p::AbstractVector, h::Int64)::Vector{Float64}
    h == 0 && return collect(Float64, p)
    n = length(p)
    c = vcat(0.0, cumsum(p))
    return [(c[min(n, i + h) + 1] - c[max(1, i - h)]) / (min(n, i + h) - max(1, i - h) + 1) for i = 1:n]
end

"""
    _wien_gain(S, N, k)

Wiener gain `Ŝ / (Ŝ + N)` from the power spectrum `S` of a reference (an average of `k` signals) and the noise power spectrum `N` of a signal. The noise left in the average (`N / k`) is subtracted from `S` to estimate the power of the signal common to all averaged signals (`Ŝ ≥ 0`).

# Arguments

- `S::AbstractVector`: power spectrum of the reference
- `N::AbstractVector`: power spectrum of the noise
- `k::Int64`: number of averaged signals

# Returns

- `Vector{Float64}`: gain in [0, 1] at each frequency
"""
function _wien_gain(S::AbstractVector, N::AbstractVector, k::Int64)::Vector{Float64}
    Sc = max.(S .- N ./ k, 0.0)
    return [Sc[i] + N[i] > 0 ? Sc[i] / (Sc[i] + N[i]) : 0.0 for i in eachindex(Sc)]
end

"""
    _denoise_wien(s; <keyword arguments>)

Denoise a 3-D signal array with a Wiener filter that keeps the signal common to a set of signals.

The reference is the average of the channels (in each epoch, `ref = :channels`) or of the epochs (for each channel, `ref = :epochs`); the noise of each signal is its difference from the reference. At each frequency, the gain is `Ŝ / (Ŝ + N)`, where `Ŝ` is the power of the common signal (the power of the reference minus the noise left in it) and `N` the noise power of the signal (averaged over epochs for `ref = :epochs`). The power spectra are smoothed over `fsmooth` Hz. The gain is real, so the filter does not shift the phase.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `fs::Int64`: sampling rate in Hz
- `ref::Symbol=:channels`: reference: `:channels` (mean of the channels; ≥ 2 channels) or `:epochs` (mean of the epochs of each channel; ≥ 2 epochs)
- `fsmooth::Real=2`: width (Hz) of the moving average smoothing the power spectra; `0` for no smoothing

# Returns

- `Array{Float64, 3}`: denoised array of the same shape as `s`
"""
function _denoise_wien(
    s::AbstractArray;
    fs::Int64,
    ref::Symbol = :channels,
    fsmooth::Real = 2,
)::Array{Float64, 3}
    # validate
    _chk3d(s)
    _check_var(ref, [:channels, :epochs], "ref")
    fs >= 1 || throw(ArgumentError("fs must be ≥ 1."))
    fsmooth >= 0 || throw(ArgumentError("fsmooth must be ≥ 0."))
    ch_n, n, ep_n = size(s)
    ref === :channels && ch_n < 2 &&
        throw(ArgumentError("ref = :channels requires ≥ 2 channels."))
    ref === :epochs && ep_n < 2 &&
        throw(ArgumentError("ref = :epochs requires ≥ 2 epochs."))

    # half-width of the spectral smoothing in frequency bins (bin width fs / n)
    h = round(Int64, fsmooth / 2 / (fs / n))
    pw(x) = _smooth_bins(abs2.(rfft(x)), h)
    flt(x, g) = irfft(rfft(x) .* g, n)

    s_new = similar(s, Float64)
    if ref === :channels
        Threads.@threads :static for ep_idx = 1:ep_n
            x = @view s[:, :, ep_idx]
            m = vec(mean(x; dims = 1))
            S = pw(m)
            for ch_idx = 1:ch_n
                N = pw(x[ch_idx, :] .- m)
                s_new[ch_idx, :, ep_idx] = flt(x[ch_idx, :], _wien_gain(S, N, ch_n))
            end
        end
    else
        Threads.@threads :static for ch_idx = 1:ch_n
            x = @view s[ch_idx, :, :]
            m = vec(mean(x; dims = 2))
            S = pw(m)
            N = mean(pw(x[:, ep_idx] .- m) for ep_idx = 1:ep_n)
            g = _wien_gain(S, N, ep_n)
            for ep_idx = 1:ep_n
                s_new[ch_idx, :, ep_idx] = flt(x[:, ep_idx], g)
            end
        end
    end

    return s_new
end

"""
    denoise_wien(obj; <keyword arguments>)

Denoise selected channels of a NEURO object with a Wiener filter that keeps the signal common to the channels (or epochs).

The reference is the mean of the selected channels in each epoch (`ref = :channels`) or the mean of the epochs of each channel (`ref = :epochs`, e.g. the ERP). The noise of each signal is its difference from the reference. At each frequency, the signal is multiplied by `Ŝ / (Ŝ + N)`, where `Ŝ` is the power of the common signal (power of the reference minus the noise left in it) and `N` the noise power of the signal: frequencies at which the signal follows the reference are kept, frequencies dominated by activity not shared with the reference are attenuated. The filter does not shift the phase.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `ref::Symbol=:channels`: reference:
    - `:channels`: mean of the selected channels (in each epoch); requires ≥ 2 channels
    - `:epochs`: mean of the epochs of each channel; requires ≥ 2 epochs
- `fsmooth::Real=2`: width (Hz) of the moving average smoothing the power spectra; `0` for no smoothing

# Returns

- `NeuroAnalyzer.NEURO`: new object with denoised channels
"""
function denoise_wien(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    ref::Symbol = :channels,
    fsmooth::Real = 2,
)::NeuroAnalyzer.NEURO
    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))

    # create new dataset
    obj_tmp = deepcopy(obj)

    obj_tmp.data[ch, :, :] = _denoise_wien(
        @view(obj.data[ch, :, :]); fs = sr(obj), ref = ref, fsmooth = fsmooth,
    )
    push!(obj_tmp.history, "denoise_wien(obj; ch=$ch, ref=$ref, fsmooth=$fsmooth)")

    return obj_tmp
end

"""
    denoise_wien!(obj; <keyword arguments>)

Denoise selected channels of a NEURO object in-place with a Wiener filter (see `denoise_wien()`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; modified in-place
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `ref::Symbol=:channels`: reference: `:channels` (mean of the selected channels) or `:epochs` (mean of the epochs of each channel)
- `fsmooth::Real=2`: width (Hz) of the moving average smoothing the power spectra

# Returns

- `Nothing`
"""
function denoise_wien!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    ref::Symbol = :channels,
    fsmooth::Real = 2,
)::Nothing
    obj_tmp = denoise_wien(obj; ch = ch, ref = ref, fsmooth = fsmooth)
    obj.data = obj_tmp.data
    obj.history = obj_tmp.history
    obj_tmp = nothing

    return nothing
end
