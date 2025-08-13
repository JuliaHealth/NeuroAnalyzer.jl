export emd

"""
    emd(s, x; <keyword arguments>)

Perform Empirical Mode Decomposition (EMD).

# Arguments

- `s::AbstractVector`: signal
- `x::AbstractVector`: x-axis points (e.g. time points)
- `epsilon::Real=0.3`: decomposition stops when sum of the difference is lower than `epsilon`

# Returns

- `imf::Matrix{Float64}`: intrinsic mode functions (IMF) (by rows) and residue (last row in the matrix)
"""
function emd(s::AbstractVector, x::AbstractVector; epsilon::Real=0.3)::Matrix{Float64}

    @assert epsilon > 0 "epsilon must be > 0."

    # s must not contain 0s
    s_tmp = deepcopy(s)
    s_tmp[s_tmp .== 0] .= eps()

    imf_v = Vector{Float64}[]
    res = zeros(length(s))

    sd = Inf
    n_sieves = 1
    while sd > epsilon
        # cubic spline envelopes of all local extremas
        e_max = env_up(s_tmp, x, d=2)
        e_min = env_lo(s_tmp, x, d=2)
        e_avg = @. (e_max + e_min) / 2
        imf_tmp = @. s_tmp - e_avg

        maxs = findpeaks(imf_tmp, d=2)
        mins = findpeaks(NeuroAnalyzer._flipx(imf_tmp), d=2)
        n_extrema = length(maxs) + length(mins)

        n_roots = NeuroAnalyzer._zeros(imf_tmp)

        res = @. s_tmp - imf_tmp
        sd = sum( @. abs2.(s_tmp - imf_tmp) / s_tmp^2 )

        # check IMF basic conditions
        if n_roots >= n_extrema - 1 &&
            n_roots <= n_extrema + 1 &&
            n_extrema >= n_roots - 1 &&
            n_extrema <= n_roots + 1 &&
            n_roots > 1 &&
            n_extrema > 1
            # also: e_avg should be near zero at any point

            # calculate stopping criterion
            push!(imf_v, imf_tmp)
            _info("IMF found: $(length(imf_v)), # sieves: $n_sieves, SD: $(round(sd, digits=2))")
            s_tmp = res
            n_sieves = 1
        else
            s_tmp = imf_tmp
            n_sieves += 1
        end

        # check stopping criterion
        if sd < epsilon
            # add the residue to the set
            push!(imf_v, res)
        end
    end

    if length(imf_v) > 0
        imf = zeros(length(imf_v), length(imf_v[1]))
        for idx in eachindex(imf_v)
            imf[idx, :] = imf_v[idx]
        end
    else
        imf = Matrix{Float64}(undef, 0, 0)
    end

    return imf

end

"""
    emd(obj; <keyword arguments>)

Perform Empirical Mode Decomposition (EMD).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::String`: channel name
- `ep::Int64`: epoch number
- `epsilon::Real=0.3`: decomposition stops when sum of the difference is lower than `epsilon`

# Returns

- `imf::Matrix{Float64}`: intrinsic mode functions (IMF) (by rows) and residue (last row in the matrix)
"""
function emd(obj::NeuroAnalyzer.NEURO; ch::String, ep::Int64, epsilon::Real=0.3)::Matrix{Float64}

    ch = exclude_bads ? get_channel(obj, ch=ch, exclude="bad")[1] : get_channel(obj, ch=ch, exclude="")[1]
    _check_epochs(obj, ep)
    imf = @views emd(obj.data[ch, :, ep], obj.epoch_time, epsilon=epsilon)
    _info("$(size(imf, 1) - 1) IMFs were calculated")

    return imf

end
