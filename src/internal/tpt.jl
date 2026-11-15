"""
    _tpt_peaks(acc, fs)

Pinches in a TPT recording: peaks of the movement magnitude.

The posture of the hand (and slow drifts) is removed from each acceleration axis with a running median over 2 seconds, which is not affected by the short pinches; the movement magnitude is the Euclidean norm of the remaining 3-D acceleration, smoothed over 3 samples. Its peaks are pinches if their prominence is ≥ the larger of 6 × the robust standard deviation of the magnitude (noise) and 25% of its 99th percentile above the median (smaller movements between pinches), and if they are ≥ 0.12 s apart.

# Arguments

- `acc::AbstractMatrix`: acceleration, shape (3, samples)
- `fs::Int64`: sampling rate

# Returns

- `Vector{Int64}`: sample indices of the pinches
"""
function _tpt_peaks(acc::AbstractMatrix, fs::Int64)::Vector{Int64}
    _wip()
    size(acc, 1) == 3 || throw(ArgumentError("acc must have 3 rows (X, Y, Z)."))
    n = size(acc, 2)
    n < 3 && return Int64[]

    # posture: running median over 2 s
    h = fs
    mag = zeros(n)
    for j = 1:n
        w = max(1, j - h):min(n, j + h)
        mag[j] = sqrt(sum(abs2(acc[c, j] - median(view(acc, c, w))) for c = 1:3))
    end
    # merge jittery double peaks
    mag = [mean(view(mag, max(1, j - 1):min(n, j + 1))) for j = 1:n]

    base = median(mag)
    noise = 1.4826 * median(abs.(mag .- base))
    prom = max(6 * noise, 0.25 * (quantile(mag, 0.99) - base))
    prom > 0 || return Int64[]
    p_idx, _ = findpeaks1d(mag; prominence = prom, distance = max(1, round(Int64, 0.12 * fs)))

    return sort(p_idx)
end
