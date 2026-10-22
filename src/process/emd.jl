export emd

# indices of the local maxima and minima of h (plateaus: the first sample)
function _emd_extrema(h::AbstractVector)::Tuple{Vector{Int64}, Vector{Int64}}
    mx = Int64[]
    mn = Int64[]
    for i = 2:(length(h) - 1)
        h[i - 1] < h[i] >= h[i + 1] && push!(mx, i)
        h[i - 1] > h[i] <= h[i + 1] && push!(mn, i)
    end
    return mx, mn
end

# cubic-spline envelope through the extrema idx of h; up to 2 extrema at each end are
# mirrored about the ends of the signal (Rilling et al. 2003), so that the spline is
# anchored beyond the ends (reduces end effects)
function _emd_envelope(h::AbstractVector, x::AbstractVector, idx::Vector{Int64})::Vector{Float64}
    nm = min(2, length(idx))
    l = idx[nm:-1:1]
    r = idx[end:-1:(end - nm + 1)]
    xs = vcat(2 * x[1] .- x[l], x[idx], 2 * x[end] .- x[r])
    ys = vcat(h[l], h[idx], h[r])
    return Spline1D(Float64.(xs), Float64.(ys); k = min(3, length(xs) - 1), bc = "extrapolate")(x)
end

"""
    _emd(s, x; <keyword arguments>)

Perform Empirical Mode Decomposition (EMD) for a 1-D signal vector (Huang et al. 1998).

The signal is decomposed into Intrinsic Mode Functions (IMFs), from the fastest to the slowest oscillation, and a residual: `s = IMF₁ + IMF₂ + … + residual`. Each IMF is extracted from the current residual by sifting: the mean of the upper and lower cubic-spline envelopes (through the local maxima and minima) is subtracted repeatedly until the relative change `SD = Σ m² / Σ h²` (`m`: mean envelope, `h`: current IMF candidate) drops below `epsilon` (or `max_sift` iterations); the IMF is then subtracted from the residual. The decomposition stops when the residual has fewer than 2 maxima or 2 minima (a trend), or after `max_imf` IMFs.

# Arguments

- `s::AbstractVector`: signal vector
- `x::AbstractVector`: x-axis points (e.g. time points)
- `epsilon::Real=0.3`: sifting stops when SD drops below `epsilon`; smaller values give more sifting (narrower-band IMFs)
- `max_imf::Int64=20`: maximum number of IMFs
- `max_sift::Int64=100`: maximum number of sifting iterations per IMF

# Returns

- `Matrix{Float64}`: IMFs by rows, with the residual as the last row; the rows sum to `s`
"""
function _emd(
    s::AbstractVector,
    x::AbstractVector;
    epsilon::Real = 0.3,
    max_imf::Int64 = 20,
    max_sift::Int64 = 100,
)::Matrix{Float64}
    # validate
    epsilon > 0 || throw(ArgumentError("epsilon must be > 0."))
    max_imf >= 1 || throw(ArgumentError("max_imf must be ≥ 1."))
    max_sift >= 1 || throw(ArgumentError("max_sift must be ≥ 1."))
    length(s) == length(x) || throw(ArgumentError("Lengths of s and x must be equal."))

    r = collect(Float64, s)
    imfs = Vector{Vector{Float64}}()
    while length(imfs) < max_imf
        mx, mn = _emd_extrema(r)
        # too few extrema for a reliable envelope: the residual is the trend
        (length(mx) < 3 || length(mn) < 3) && break
        # sifting
        h = copy(r)
        n_sift = 0
        for k = 1:max_sift
            mx, mn = _emd_extrema(h)
            (length(mx) < 2 || length(mn) < 2) && break
            m = (_emd_envelope(h, x, mx) .+ _emd_envelope(h, x, mn)) ./ 2
            sd = sum(abs2, m) / max(sum(abs2, h), eps())
            h .-= m
            n_sift = k
            sd < epsilon && break
        end
        push!(imfs, h)
        r .-= h
        _info("IMF $(length(imfs)): $n_sift sifting iterations")
    end

    # IMFs by rows, residual last
    imf = Matrix{Float64}(undef, length(imfs) + 1, length(r))
    for idx in eachindex(imfs)
        imf[idx, :] = imfs[idx]
    end
    imf[end, :] = r

    return imf
end

"""
    emd(obj; <keyword arguments>)

Perform Empirical Mode Decomposition (EMD) of one channel in one epoch of a NEURO object.

The signal is decomposed into Intrinsic Mode Functions (IMFs), from the fastest to the slowest oscillation, and a residual (trend), which sum to the signal. See `_emd()` for the algorithm.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: channel name; must resolve to exactly one channel
- `ep::Int64`: epoch index
- `epsilon::Real=0.3`: sifting of each IMF stops when the relative change (SD criterion) drops below `epsilon`; smaller values give narrower-band IMFs
- `max_imf::Int64=20`: maximum number of IMFs

# Returns

- `Matrix{Float64}`: intrinsic mode functions (IMF) by rows, with the residual as the last row; the rows sum to the signal
"""
function emd(
    obj::NeuroAnalyzer.NEURO;
    ch::String,
    ep::Int64,
    epsilon::Real = 0.3,
    max_imf::Int64 = 20,
)::Matrix{Float64}
    # resolve channel name to a single integer index; [1] selects the first (and expected only) result from get_channel
    ch =
        exclude_bads ? get_channel(obj; ch = ch, exclude = "bad") :
        get_channel(obj; ch = ch, exclude = "")
    isempty(ch) && throw(ArgumentError("No channels selected."))
    length(ch) == 1 || throw(ArgumentError("ch must resolve to exactly one channel."))
    ch = ch[1]

    # validate
    _check_epochs(obj, ep)

    imf = _emd(@view(obj.data[ch, :, ep]), obj.epoch_time; epsilon = epsilon, max_imf = max_imf)
    _info("$(size(imf, 1) - 1) IMFs were calculated")

    return imf
end
