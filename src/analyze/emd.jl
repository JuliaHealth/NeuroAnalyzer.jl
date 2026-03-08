export emd

"""
    emd(s, x; <keyword arguments>)

Perform Empirical Mode Decomposition (EMD). Decomposes a signal into Intrinsic Mode Functions (IMFs) by iteratively "sieving" out oscillatory components from fastest to slowest.

Each sieving iteration:

1. Replace any exact zeros with eps() to avoid spline interpolation issues.
2. Compute upper/lower cubic-spline envelopes of the local extrema.
3. Subtract the mean envelope to get a candidate IMF.
4. Check the IMF stopping criterion (number of extrema ≈ number of roots).
5. If criterion met: accept the IMF, subtract it, reset the sieve counter.
6. If stopping condition SD < epsilon: append the residue and exit.

Returns a matrix whose rows are IMFs (1..end-1) and the final residue (end).

# Arguments

- `s::AbstractVector`: signal vector
- `x::AbstractVector`: x-axis points (e.g. time points)
- `epsilon::Real=0.3`: decomposition stops when the normalised sum of squared differences (SD criterion) drops below `epsilon`

# Returns

  - `imf::Matrix{Float64}`: intrinsic mode functions (IMF) by rows, with the residue as the last row; returns an empty `0×0` matrix if no IMFs found
"""
function emd(s::AbstractVector, x::AbstractVector; epsilon::Real = 0.3)::Matrix{Float64}

    @assert epsilon > 0 "epsilon must be > 0."

    # work on a mutable copy so the caller's signal is not modified
    s_tmp = copy(s)

    # accumulate accepted IMFs
    imf_v = Vector{Float64}[]
    res = zeros(length(s))
    sd = Inf
    n_sieves = 1

    while sd > epsilon

        # cubic spline interpolation can fail at exact zeros — nudge them.
        s_tmp[s_tmp .== 0] .= eps()

        # compute upper and lower envelopes from local maxima/minima
        e_max = env_up(s_tmp, x, d = 2)
        e_min = env_lo(s_tmp, x, d = 2)

        # mean envelope: the "trend" to subtract
        e_avg = @. (e_max + e_min) / 2

        # candidate IMF: signal minus its mean envelope
        imf_tmp = @. s_tmp - e_avg

        # count extrema and zero-crossings to test the IMF conditions
        maxs = findpeaks(imf_tmp, d = 2)
        mins = findpeaks(_flipx(imf_tmp), d = 2)
        n_extrema = length(maxs) + length(mins)
        n_roots = _zeros(imf_tmp)

        # residue after removing the candidate IMF
        res = @. s_tmp - imf_tmp
        sd = sum((s_tmp[i] - imf_tmp[i])^2 / s_tmp[i]^2 for i in eachindex(s_tmp))

        # IMF validity check 
        # a valid IMF must have the number of extrema and zero-crossings differ by at most one, and both must exceed 1 (non-trivial oscillation)
        if n_roots >= n_extrema - 1 &&
                n_roots <= n_extrema + 1 &&
                n_extrema >= n_roots - 1 &&
                n_extrema <= n_roots + 1 &&
                n_roots > 1 &&
                n_extrema > 1

            # accept this IMF; next iteration works on the residue
            push!(imf_v, imf_tmp)
            _info("IMF found: $(length(imf_v)), sieves: $n_sieves, SD: $(round(sd, digits = 2))")
            s_tmp = res
            n_sieves = 1
        else
            # IMF conditions not met; refine further using the candidate
            s_tmp = imf_tmp
            n_sieves += 1
        end

        # once SD drops below epsilon, append the residue and exit.
        sd < epsilon && push!(imf_v, res)

    end

    # assemble all accepted IMFs (and residue) into a rows-as-IMFs matrix
    if length(imf_v) > 0
        imf = Matrix{Float64}(undef, length(imf_v), length(imf_v[1]))
        for idx in eachindex(imf_v)
            imf[idx, :] .= imf_v[idx]
        end
    else
        imf = Matrix{Float64}(undef, 0, 0)
    end

    return imf

end

"""
    emd(obj; <keyword arguments>)

Perform Empirical Mode Decomposition (EMD). Decomposes a signal into Intrinsic Mode Functions (IMFs) by iteratively "sieving" out oscillatory components from fastest to slowest.

Each sieving iteration:

1. Replace any exact zeros with eps() to avoid spline interpolation issues.
2. Compute upper/lower cubic-spline envelopes of the local extrema.
3. Subtract the mean envelope to get a candidate IMF.
4. Check the IMF stopping criterion (number of extrema ≈ number of roots).
5. If criterion met: accept the IMF, subtract it, reset the sieve counter.
6. If stopping condition SD < epsilon: append the residue and exit.

Returns a matrix whose rows are IMFs (1..end-1) and the final residue (end).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::String`: channel name
- `ep::Int64`: epoch number
- `epsilon::Real=0.3`: decomposition stops when the normalised sum of squared differences (SD criterion) drops below `epsilon`

# Returns

- `imf::Matrix{Float64}`: intrinsic mode functions (IMF) by rows, with the residue as the last row
"""
function emd(obj::NeuroAnalyzer.NEURO; ch::String, ep::Int64, epsilon::Real = 0.3)::Matrix{Float64}

    # resolve channel name to a single integer index; [1] selects the first (and expected only) result from get_channel
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad")[1] : get_channel(obj, ch = ch, exclude = "")[1]

    _check_epochs(obj, ep)

    imf = emd(@view(obj.data[ch, :, ep]), obj.epoch_time, epsilon = epsilon)

    # report how many IMFs (excluding the residue) were found
    size(imf, 1) > 0 && _info("$(size(imf, 1) - 1) IMFs were calculated")

    return imf

end
