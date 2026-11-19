"""
    _tpt_peaks(acc, fs)

Pinches in a TPT recording: peaks of the movement envelope.

1. Single-sample glitches of the accelerometer are removed from each acceleration axis with a running median over 5 samples (± 0.04 s at 50 Hz).
2. The posture of the hand (and slow drifts) is removed from each axis with a running median over ± 2 seconds, which is not affected by the short pinches.
3. The movement magnitude is the Euclidean norm of the remaining 3-D acceleration; its envelope is the moving mean over ± 0.04 s. A pinch (closing, contact and opening of the fingers) is a burst of several short acceleration peaks; the envelope merges them into one peak.
4. Peaks of the envelope are pinches if their prominence is ≥ the larger of 6 × the robust standard deviation of the envelope (noise; recordings without pinches) and 20% of its 99th percentile above the median (smaller movements between pinches), and if they are ≥ 0.25 s apart.

The parameters were selected using 22 recordings with 0 to 10 pinches (cross-validated between two repeats of each number of pinches): the number of pinches was exact in 17/22 recordings, mean absolute error 0.7 pinch.

# Arguments

- `acc::AbstractMatrix`: acceleration, shape (3, samples)
- `fs::Int64`: sampling rate

# Returns

- `Vector{Int64}`: sample indices of the pinches
"""
function _tpt_peaks(acc::AbstractMatrix, fs::Int64)::Vector{Int64}
    size(acc, 1) == 3 || throw(ArgumentError("acc must have 3 rows (X, Y, Z)."))
    n = size(acc, 2)
    n < 3 && return Int64[]

    _runmed(x, h) = [median(view(x, max(1, j - h):min(n, j + h))) for j = 1:n]
    _runmean(x, h) = [mean(view(x, max(1, j - h):min(n, j + h))) for j = 1:n]

    # 1. remove glitches
    h_despike = max(1, round(Int64, 0.04 * fs))
    a = reduce(vcat, (_runmed(@view(acc[c, :]), h_despike)' for c = 1:3))

    # 2-3. remove posture, movement magnitude and its envelope
    h_post = 2 * fs
    mag = zeros(n)
    for c = 1:3
        mag .+= abs2.(@view(a[c, :]) .- _runmed(@view(a[c, :]), h_post))
    end
    env = _runmean(sqrt.(mag), max(1, round(Int64, 0.04 * fs)))

    # 4. peaks
    base = median(env)
    noise = 1.4826 * median(abs.(env .- base))
    prom = max(6 * noise, 0.2 * (quantile(env, 0.99) - base))
    prom > 0 || return Int64[]
    p_idx, _ = findpeaks1d(env; prominence = prom, distance = max(1, round(Int64, 0.25 * fs)))

    return sort(p_idx)
end
