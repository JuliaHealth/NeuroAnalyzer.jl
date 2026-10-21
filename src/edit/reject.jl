export channel_reject
export channel_reject!
export epoch_reject

const _REJECT_METHODS = [:flat, :rmse, :rmsd, :euclid, :var, :p2p, :tkeo, :kurt, :z, :ransac, :amp]

# methods comparing a per-channel statistic between channels (channel_reject) or epochs
# (epoch_reject) with robust z-scores; :var is two-sided (high and low variance)
const _REJECT_ROBUST = (:rmse, :rmsd, :euclid, :var, :p2p, :tkeo, :kurt, :z)

# ------------------------------------------------------------------------------------------------
# Internal helpers
# ------------------------------------------------------------------------------------------------

"""
    _robust_z(x)

Robust z-scores of `x`: (x - median) / (1.4826 × MAD). If the MAD is 0 (more than half of the values are identical), the mean absolute deviation is used instead; constant input gives zeros.
"""
function _robust_z(x::AbstractVector{<:Real})::Vector{Float64}
    m = median(x)
    s = 1.4826 * median(abs.(x .- m))
    s == 0 && (s = 1.2533 * mean(abs.(x .- m)))
    s == 0 && return zeros(length(x))
    return (x .- m) ./ s
end

"""
    _reject_windows(n, w)

Non-overlapping windows of `w` samples covering a signal of `n` samples (at least one window; a signal shorter than `w` is one window).
"""
function _reject_windows(n::Int64, w::Int64)::Vector{UnitRange{Int64}}
    w >= n && return [1:n]
    return [i:(i + w - 1) for i = 1:w:(n - w + 1)]
end

"""
    _reject_stat(s, method; <keyword arguments>)

Statistic of each channel (row) of `s` (channels × samples, one epoch) for the robust detection methods; larger values are more abnormal (for `:var` both tails are abnormal).
"""
function _reject_stat(
    s::AbstractMatrix,
    method::Symbol;
    w::Int64,
    tkeo_method::Symbol = :pow,
    z::Real = 3.5,
)::Vector{Float64}
    ch_n, n = size(s)
    # small offset: the logarithm of 0 (flat signal) stays finite
    lg(x) = log(x + 1e-12)
    if method in (:rmse, :rmsd, :euclid)
        # distance from the median signal of all channels
        m = vec(median(s; dims = 1))
        d = s .- m'
        method === :rmse && return [lg(sqrt(mean(abs2, @view(d[i, :])))) for i = 1:ch_n]
        # :rmsd ignores a constant offset between the channel and the median signal
        method === :rmsd && return [lg(std(@view(d[i, :]))) for i = 1:ch_n]
        return [lg(sqrt(sum(abs2, @view(d[i, :])))) for i = 1:ch_n]
    elseif method === :var
        return [lg(var(@view(s[i, :]))) for i = 1:ch_n]
    elseif method === :p2p
        win = _reject_windows(n, w)
        return [lg(maximum(maximum(@view(s[i, r])) - minimum(@view(s[i, r])) for r in win))
                for i = 1:ch_n]
    elseif method === :kurt
        # excess kurtosis; a flat signal (undefined kurtosis) gets the median of the other
        # channels, so that it is not flagged
        k = [kurtosis(@view(s[i, :])) for i = 1:ch_n]
        k[isnan.(k)] .= all(isnan, k) ? 0.0 : median(filter(!isnan, k))
        return k
    elseif method === :tkeo
        return [lg(mean(abs.(_tkeo(s[i, :]; method = tkeo_method)))) for i = 1:ch_n]
    elseif method === :z
        # fraction of outlying samples: robust z-score within the channel > z
        return [count(abs.(_robust_z(@view(s[i, :]))) .> z) / n for i = 1:ch_n]
    end
    throw(ArgumentError("Unknown method :$method."))
end

"""
    _reject_absolute(s, method; <keyword arguments>)

Bad-channel flags of each channel (row) of `s` (channels × samples, one epoch) for the methods with absolute criteria (`:flat`, `:amp`, `:z`).
"""
function _reject_absolute(
    s::AbstractMatrix,
    method::Symbol;
    w::Int64,
    flat_tol::Real,
    flat_fr::Real,
    amp_t::Real,
    z::Real,
    p::Real,
)::Vector{Bool}
    ch_n, n = size(s)
    if method === :flat
        # fraction of windows in which the signal is (nearly) constant
        win = _reject_windows(n, w)
        return [count(std(@view(s[i, r])) < flat_tol for r in win) / length(win) > flat_fr
                for i = 1:ch_n]
    elseif method === :amp
        return [maximum(abs, @view(s[i, :])) > amp_t for i = 1:ch_n]
    end
    throw(ArgumentError("Unknown method :$method."))
end

"""
    _reject_neighbors(s, xyz; <keyword arguments>)

Statistic of the `:ransac` method: the signal of each channel is predicted from its 3 nearest neighbors (inverse-distance weighted mean); the statistic is the median correlation between the channel and the prediction in windows of `w` samples (a flat channel has correlation 0).
"""
function _reject_neighbors(s::AbstractMatrix, xyz::AbstractMatrix; w::Int64)::Vector{Float64}
    ch_n, n = size(s)
    k = min(3, ch_n - 1)
    win = _reject_windows(n, w)
    rc = zeros(ch_n)
    for i = 1:ch_n
        d = [i == j ? Inf : sqrt(sum(abs2, @view(xyz[i, :]) .- @view(xyz[j, :]))) for j = 1:ch_n]
        nb = sortperm(d)[1:k]
        wt = 1 ./ max.(d[nb], 1e-6)
        pred = vec(sum(s[nb, :] .* wt; dims = 1)) ./ sum(wt)
        r = [cor(@view(s[i, rw]), @view(pred[rw])) for rw in win]
        # a constant window (e.g. flat channel) has no correlation
        r[isnan.(r)] .= 0.0
        rc[i] = median(r)
    end
    return rc
end

"""
    _reject_flags(obj; <keyword arguments>)

Flags (channels × epochs) of one detection method; `across` sets how the robust statistics are compared: `:channels` (between the channels of each epoch, for `channel_reject()`) or `:epochs` (between the epochs of each channel, for `epoch_reject()`).
"""
function _reject_flags(
    obj::NeuroAnalyzer.NEURO,
    ch::Vector{Int64},
    method::Symbol;
    across::Symbol,
    w::Int64,
    flat_tol::Real,
    flat_fr::Real,
    p::Real,
    tkeo_method::Symbol,
    z::Real,
    ransac_r::Real,
    ransac_tr::Real,
    amp_t::Real,
)::Matrix{Bool}
    ch_n = length(ch)
    ep_n = nepochs(obj)
    flags = falses(ch_n, ep_n)

    if method === :z && across === :epochs
        # fraction of samples of each epoch outlying from the whole channel (robust z-score
        # computed from all epochs of the channel) > z; the epoch is flagged if the fraction
        # exceeds 1 - p
        Threads.@threads :static for ch_idx = 1:ch_n
            x = @view(obj.data[ch[ch_idx], :, :])
            xz = reshape(_robust_z(vec(x)), size(x))
            for ep_idx = 1:ep_n
                flags[ch_idx, ep_idx] = count(>(z), abs.(@view(xz[:, ep_idx]))) / size(x, 1) > 1 - p
            end
        end
    elseif method in _REJECT_ROBUST
        st = zeros(ch_n, ep_n)
        Threads.@threads :static for ep_idx = 1:ep_n
            st[:, ep_idx] = _reject_stat(
                @view(obj.data[ch, :, ep_idx]), method; w = w, tkeo_method = tkeo_method, z = z,
            )
        end
        rz = similar(st)
        if across === :channels
            for ep_idx = 1:ep_n
                rz[:, ep_idx] = _robust_z(st[:, ep_idx])
            end
        else
            for ch_idx = 1:ch_n
                rz[ch_idx, :] = _robust_z(st[ch_idx, :])
            end
        end
        flags .= method === :var ? abs.(rz) .> z : rz .> z
        # :z: in addition, more than 1 - p of the samples must be outliers
        method === :z && (flags .&= st .> 1 - p)
    elseif method === :ransac
        _check_datatype(obj, ["eeg", "seeg", "ecog", "meg"])
        ch_list = labels(obj)[ch]
        loc_labels = lowercase.(obj.locs[!, :label])
        loc_idx = [findfirst(==(lowercase(l)), loc_labels) for l in ch_list]
        any(isnothing, loc_idx) &&
            throw(ArgumentError(":ransac requires locations of all channels."))
        xyz = Matrix{Float64}(obj.locs[Vector{Int64}(loc_idx), [:loc_x, :loc_y, :loc_z]])
        rc = zeros(ch_n, ep_n)
        Threads.@threads :static for ep_idx = 1:ep_n
            rc[:, ep_idx] = _reject_neighbors(@view(obj.data[ch, :, ep_idx]), xyz; w = w)
        end
        # abnormally low correlation compared with the other channels (or epochs) and below ransac_r
        rz = similar(rc)
        if across === :channels
            for ep_idx = 1:ep_n
                rz[:, ep_idx] = _robust_z(rc[:, ep_idx])
            end
        else
            for ch_idx = 1:ch_n
                rz[ch_idx, :] = _robust_z(rc[ch_idx, :])
            end
        end
        flags .= (rz .< -z) .& (rc .< ransac_r)
    else
        Threads.@threads :static for ep_idx = 1:ep_n
            flags[:, ep_idx] = _reject_absolute(
                @view(obj.data[ch, :, ep_idx]), method;
                w = w, flat_tol = flat_tol, flat_fr = flat_fr, amp_t = amp_t, z = z, p = p,
            )
        end
    end

    return flags
end

# shared validation of channel_reject() and epoch_reject()
function _reject_check(obj, ch, method, p, tc, z, ransac_r, ransac_tr, w)
    _in(p, (0, 1), "p")
    _in(tc, (0, 1), "tc")
    _in(ransac_r, (0, 1), "ransac_r")
    _in(ransac_tr, (0, 1), "ransac_tr")
    z > 0 || throw(ArgumentError("z must be > 0."))
    w >= 2 || throw(ArgumentError("w must be ≥ 2."))
    method = method isa Symbol ? [method] : method
    for m in method
        _check_var(m, _REJECT_METHODS, "method")
    end
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))
    return ch, method
end

"""
    channel_reject(obj; <keyword arguments>)

Detect bad channels using one or more artifact-detection methods.

Each method flags channels in each epoch; a channel is marked bad if any of the selected methods flags it in at least a fraction `tc` of the epochs (for a continuous recording, in its only epoch). The returned mask has one entry per *selected* channel.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `method::Union{Symbol, Vector{Symbol}}=[:flat, :rmse, :rmsd, :euclid, :var, :p2p, :tkeo, :kurt, :z, :ransac, :amp]`: one or more detection methods (see below); default: all available methods
- `w::Int64=sr(obj)`: window length in samples (default = 1 s), used by `:flat`, `:p2p` and `:ransac`
- `flat_tol::Float64=0.1`: a window is flat if the standard deviation of the signal in it is < `flat_tol` (signal units)
- `flat_fr::Float64=0.3`: a channel is flat if more than this fraction of its windows is flat
- `p::Float64=0.99`: `:z` flags only channels with more than `1 - p` of outlying samples
- `tc::Float64=0.2`: fraction of epochs in which a channel must be flagged to be marked bad
- `tkeo_method::Symbol=:pow`: TKEO variant; see `tkeo()` for options
- `z::Real=3.5`: robust z-score threshold of the statistical methods (and of outlying samples in `:z`); must be > 0
- `ransac_r::Float64=0.6`: `:ransac` flags only channels whose correlation with the signal predicted from their neighbors is below `ransac_r`
- `ransac_tr::Float64=0.4`: not used (kept for compatibility)
- `ransac_t::Float64=100.0`: not used (kept for compatibility)
- `amp_t::Real=400.0`: amplitude threshold (signal units); channels with `|signal| > amp_t` are flagged

# Detection methods

The statistical methods compute a statistic for each channel and epoch and flag channels whose robust z-score (median/MAD, computed across the channels of each epoch) exceeds `z`; they (and `:ransac`) need at least 3 channels:

- `:var`: variance (log); flags both high (noise, artifacts) and low (flat, bridged, disconnected) variance
- `:p2p`: largest peak-to-peak amplitude in windows of `w` samples (log): jumps, pops, large artifacts
- `:kurt`: kurtosis: spiky signals (electrode pops, spikes)
- `:tkeo`: mean Teager-Kaiser energy (log): high-frequency activity (muscle, noise)
- `:rmse`: RMS deviation from the median signal of all channels (log): channels different from the others
- `:rmsd`: as `:rmse`, ignoring a constant offset
- `:euclid`: as `:rmse`, Euclidean distance
- `:z`: fraction of outlying samples (robust z-score within the channel > `z`), which must also exceed `1 - p`: spikes, pops, bursts of artifacts

The other methods use absolute criteria:

- `:flat`: more than `flat_fr` of windows with SD < `flat_tol`: disconnected, saturated channels
- `:amp`: `|signal| > amp_t` (the signal should be filtered, e.g. high-pass, so that DC offsets do not exceed the threshold)
- `:ransac`: correlation with the signal predicted from the 3 nearest channels (inverse-distance weighted mean; median over windows of `w` samples) is abnormally low compared with the other channels (robust z-score < `-z`) and below `ransac_r`; requires channel locations

# Returns

- `Vector{Bool}`: bad-channel mask (length = number of selected channels); `true` = bad
"""
function channel_reject(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    method::Union{Symbol, Vector{Symbol}} = [
        :flat,
        :rmse,
        :rmsd,
        :euclid,
        :var,
        :p2p,
        :tkeo,
        :kurt,
        :z,
        :ransac,
        :amp,
    ],
    w::Int64 = sr(obj),
    flat_tol::Float64 = 0.1,
    flat_fr::Float64 = 0.3,
    p::Float64 = 0.99,
    tc::Float64 = 0.2,
    tkeo_method::Symbol = :pow,
    z::Real = 3.5,
    ransac_r::Float64 = 0.6,
    ransac_tr::Float64 = 0.4,
    ransac_t::Float64 = 100.0,
    amp_t::Real = 400.0,
)::Vector{Bool}
    ch, method = _reject_check(obj, ch, method, p, tc, z, ransac_r, ransac_tr, w)
    ch_n = length(ch)
    for m in method
        m in (_REJECT_ROBUST..., :ransac) && ch_n < 3 &&
            throw(ArgumentError(":$m requires ≥ 3 channels."))
    end

    bc = falses(ch_n)
    for m in method
        _info("Using :$m method")
        flags = _reject_flags(
            obj, ch, m;
            across = :channels, w = w, flat_tol = flat_tol, flat_fr = flat_fr, p = p,
            tkeo_method = tkeo_method, z = z, ransac_r = ransac_r, ransac_tr = ransac_tr,
            amp_t = amp_t,
        )
        # fraction of epochs in which the channel was flagged
        bc .|= vec(sum(flags; dims = 2)) ./ size(flags, 2) .>= max(tc, eps())
    end

    return collect(bc)
end

"""
    channel_reject!(obj; <keyword arguments>)

Detect bad channels and mark them in the `:bad_channel` field of the object header in-place. See `channel_reject()` for the methods and arguments.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- other arguments as in `channel_reject()`

# Returns

- `Nothing`
"""
function channel_reject!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    method::Union{Symbol, Vector{Symbol}} = [
        :flat,
        :rmse,
        :rmsd,
        :euclid,
        :var,
        :p2p,
        :tkeo,
        :kurt,
        :z,
        :ransac,
        :amp,
    ],
    w::Int64 = sr(obj),
    flat_tol::Float64 = 0.1,
    flat_fr::Float64 = 0.3,
    p::Float64 = 0.99,
    tc::Float64 = 0.2,
    tkeo_method::Symbol = :pow,
    z::Real = 3.5,
    ransac_r::Float64 = 0.6,
    ransac_tr::Float64 = 0.4,
    ransac_t::Float64 = 100.0,
    amp_t::Real = 400.0,
)::Nothing
    bc = channel_reject(
        obj;
        ch = ch,
        method = method,
        w = w,
        flat_tol = flat_tol,
        flat_fr = flat_fr,
        p = p,
        tc = tc,
        tkeo_method = tkeo_method,
        z = z,
        ransac_r = ransac_r,
        ransac_tr = ransac_tr,
        ransac_t = ransac_t,
        amp_t = amp_t,
    )
    obj.header.recording[:bad_channel][get_channel(obj; ch = ch)] = bc
    push!(obj.history, "channel_reject!(obj; ch=$ch, method=$method)")

    return nothing
end

"""
    epoch_reject(obj; <keyword arguments>)

Detect bad epochs using one or more artifact-detection methods.

Each method flags channels in each epoch; an epoch is marked bad if at least `nbad` channels are flagged in it by any of the selected methods. The statistical methods compare each epoch with the other epochs of the same channel (robust z-score across epochs), so that channels with naturally larger or smaller signals are not penalized.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `method`, `w`, `flat_tol`, `flat_fr`, `p`, `tc`, `tkeo_method`, `z`, `ransac_r`, `ransac_tr`, `ransac_t`, `amp_t`: as in `channel_reject()`; `tc` is not used
- `nbad::Int64=1`: minimum number of flagged channels to mark an epoch as bad; must be in `[1, number of selected channels]`

# Detection methods

As in `channel_reject()`; the statistical methods need at least 3 epochs. `:z` flags epochs in which more than `1 - p` of the samples of a channel are outlying (robust z-score > `z`, computed from all epochs of the channel).

With many channels and `nbad = 1`, the statistical methods are sensitive: the robust z-scores of e.g. kurtosis have heavy tails, so spiky but artifact-free epochs may be flagged. Increase `z` or `nbad` to make the detection stricter.

# Returns

- `Vector{Int64}`: sorted indices of bad epochs
"""
function epoch_reject(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    method::Union{Symbol, Vector{Symbol}} = [
        :flat,
        :rmse,
        :rmsd,
        :euclid,
        :var,
        :p2p,
        :tkeo,
        :kurt,
        :z,
        :ransac,
        :amp,
    ],
    w::Int64 = sr(obj),
    flat_tol::Float64 = 0.1,
    flat_fr::Float64 = 0.3,
    p::Float64 = 0.99,
    tc::Float64 = 0.2,
    tkeo_method::Symbol = :pow,
    z::Real = 3.5,
    ransac_r::Float64 = 0.6,
    ransac_tr::Float64 = 0.4,
    ransac_t::Float64 = 100.0,
    amp_t::Real = 400.0,
    nbad::Int64 = 1,
)::Vector{Int64}
    ch, method = _reject_check(obj, ch, method, p, tc, z, ransac_r, ransac_tr, w)
    ch_n = length(ch)
    nbad >= 1 || throw(ArgumentError("nbad must be ≥ 1."))
    nbad <= ch_n || throw(ArgumentError("nbad must be ≤ $ch_n (number of channels)."))
    for m in method
        m in _REJECT_ROBUST && nepochs(obj) < 3 &&
            throw(ArgumentError(":$m requires ≥ 3 epochs."))
        m in (:rmse, :rmsd, :euclid, :ransac) && ch_n < 3 &&
            throw(ArgumentError(":$m requires ≥ 3 channels."))
    end

    be = falses(nepochs(obj))
    for m in method
        _info("Using :$m method")
        flags = _reject_flags(
            obj, ch, m;
            across = :epochs, w = w, flat_tol = flat_tol, flat_fr = flat_fr, p = p,
            tkeo_method = tkeo_method, z = z, ransac_r = ransac_r, ransac_tr = ransac_tr,
            amp_t = amp_t,
        )
        be .|= vec(sum(flags; dims = 1)) .>= nbad
    end

    return findall(be)
end
