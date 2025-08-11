export emd

"""
    emd(s, x; <keyword arguments>)

Perform Empirical Mode Decomposition (EMD).

# Arguments

- `s::AbstractVector`: signal
- `x::AbstractVector`: x-axis points
- `epsilon::Union{Float64, Nothing}=nothing`: decomposition stops when sum of the difference is lower than `epsilon`
- `n::Union{Int64, Nothing}=nothing`: decomposition stops when `n` IMFs are created

# Returns

- `imf::Vector{Vector{Float64}}`: intrinsic mode functions (IMF) and residue (last one)
"""
function emd(s::AbstractVector, x::AbstractVector; epsilon::Union{Float64, Nothing}=nothing, n::Union{Int64, Nothing}=nothing)::Vector{Vector{Float64}}

    @assert !(epsilon isa Nothing && n isa Nothing) "Use either epsilon or n."
    epsilon isa Float64 && @assert n isa Nothing "Use either epsilon or n."
    n isa Int64 && @assert epsilon isa Nothing "Use either epsilon or n."

    # s must not contain 0s
    s_tmp = deepcopy(s)
    s_tmp[s_tmp .== 0] .= eps()

    imf = Vector{Float64}[]

    if epsilon isa Float64
        sd = Inf
        while sd > epsilon
            # cubic spline envelopes of local extremas
            e_up = env_up(s_tmp, x, d=2)
            e_lo = env_lo(s_tmp, x, d=2)
            e_m = @. (e_up + e_lo) / 2
            imf_tmp = @. s_tmp - e_m
            res = @. s_tmp - imf_tmp
            # calculate stopping criterion
            sd = sum( @. abs2(res - s_tmp) / s_tmp^2 )
            if sd > epsilon
                s_tmp = res
                push!(imf, imf_tmp)
            else
                # last IMF is the residue
                push!(imf, res)
            end
        end
    else
        # use n sifting
        res = nothing
        for _ in 1:n
            # cubic spline envelopes of local extremas
            e_up = env_up(s_tmp, x, d=2)
            e_lo = env_lo(s_tmp, x, d=2)
            e_m = @. (e_up + e_lo) / 2
            imf_tmp = @. s_tmp - e_m
            res = @. s_tmp - imf_tmp
            # calculate stopping criterion
            sd = sum( @. (res - s_tmp)^2 / s_tmp^2 )
            s_tmp = res
            push!(imf, imf_tmp)
        end
        # last IMF is the residue
        push!(imf, res)
    end

    return imf

end
