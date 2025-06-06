export linspace
export logspace
export cmax
export cmin
export tuple_order
export cums
export f_nearest
export ntapers

"""
    linspace(start, stop, length)

Generates a sequence of evenly spaced numbers between `start` and `stop`.

# Arguments

- `start::Number`
- `stop::Number`
- `n::Int64`: sequence length

# Returns

- `range::Vector{Float64}`
"""
function linspace(start::Number, stop::Number, n::Int64)::Vector{Float64}

    @assert n >= 2 "n must be ≥ 2."

    return collect(range(start, stop, n))

end

"""
    logspace(start, stop, n)

Generates a sequence of log10-spaced numbers between `start` and `stop`.

# Arguments

- `start::Number`
- `stop::Number`
- `n::Int64`: sequence length

# Returns

- `range::Vector{Float64}`
"""
function logspace(start::Number, stop::Number, n::Int64)::Vector{Float64}

    @assert n >= 2 "n must be ≥ 2."

    return collect(exp10.(range(start, stop, n)))

end

"""
    cmax(x)

Return maximum value of the complex vector.

# Arguments

- `x::Vector{ComplexF64}`

# Returns

- `cmax::ComplexF64`
"""
function cmax(x::Vector{<:Complex})::ComplexF64

    return argmax(abs2, x)

end

"""
    cmin(x)

Return minimum value of the complex vector.

# Arguments

- `x::Vector{ComplexF64}`

# Returns

- `cmin::ComplexF64`
"""
function cmin(x::Vector{<:Complex})::ComplexF64

    return argmin(abs2, x)

end

"""
    cextrema(x)

Return extreme values of the complex vector.

# Arguments

- `x::Vector{ComplexF64}`

# Returns

Tuple containing:
- `cmax::ComplexF64`
- `cmin::ComplexF64`
"""
function cextrema(x::Vector{<:Complex})::Tuple{ComplexF64, ComplexF64}

    return (cmax(x), cmin(x))

end

"""
    tuple_order(t, rev)

Order tuple elements in ascending or descending (`rev=true`) order.

# Arguments

- `t::Tuple{Real, Real}`
- `rev::Bool=false`

# Returns

- `t::Tuple{Real, Real}`
"""
function tuple_order(t::Tuple{Real, Real}, rev::Bool=false)::Tuple{Real, Real}

    (!rev && t[1] > t[2]) && (t = (t[2], t[1]))
    (rev && t[1] < t[2]) && (t = (t[2], t[1]))

    return t

end

"""
    cums(signal)

Calculate cumulative sum of a 3-dimensional array.

# Arguments

- `signal::Array{<:Real, 3}`

# Returns

- `signal_cs::Array{Float64, 3}`
"""
function cums(signal::Array{<:Real, 3})::Array{Float64, 3}

    ch_n, _, ep_n = size(signal)
    signal_cs = similar(signal)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads :greedy for ch_idx in 1:ch_n
            signal_cs[ch_idx, :, ep_idx] = @views cumsum(signal[ch_idx, :, ep_idx])
        end
    end

    return signal_cs

end

"""
    f_nearest(m, pos)

Find nearest position tuple in a matrix of positions.

# Arguments

- `m::Matrix{Tuple{Float64, Float64}}`: matrix of positions
- `p::Tuple{Float64, Float64}`: position tuple

# Returns

- `pos::Tuple{Int64, Int64}`: row and column in m
"""
function f_nearest(m::Matrix{Tuple{Float64, Float64}}, p::Tuple{Float64, Float64})::Tuple{Int64, Int64}

    d = zeros(size(m))

    @inbounds for idx1 in axes(m, 1)
        for idx2 in axes(m, 2)
            d[idx1, idx2] = euclidean(m[idx1, idx2], p)
        end
    end

    return (findmin(d)[2][1], findmin(d)[2][2])

end

"""
    ntapers(obj; df)

Return recommended number of Slepian tapers for multi-taper power spectrum analysis.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `df::Real`: frequency resolution (bandwidth); smallest distance between frequency peaks that we want to observe (e.g. 1 Hz)

# Returns

- `nt::Int64`
"""
function ntapers(obj::NeuroAnalyzer.NEURO; df::Real)::Int64

    _bin(df, (0, sr(obj) / 2))

    n = epoch_len(obj) / sr(obj)
    nt = round(Int64, df * n) - 1

    return nt

end
