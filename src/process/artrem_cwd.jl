export artrem_cwd
export artrem_cwd!

"""
    artrem_cwd(s, t; <keyword arguments>)

Remove an artifact from a signal using continuous wavelet decomposition (CWD).

The signal is transformed into the time-frequency domain via CWD, the coefficients within the specified time (`tseg`) and frequency (`fseg`) windows are zeroed, and the signal is reconstructed via the inverse CWD (iCWD).

# Arguments

- `s::AbstractVector`: signal vector
- `t::AbstractVector`: time-point vector corresponding to `s`; must not be empty
- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=2)`: wavelet to use; see ContinuousWavelets.jl for available wavelets
- `tseg::Tuple{Real, Real}`: time window of the artifact `(t_start, t_end)` in seconds; must be a strictly ascending pair within `[t[1], t[end]]`
- `fseg::Tuple{Real, Real}`: frequency window of the artifact `(f_low, f_high)` in Hz; must be a strictly ascending pair within the wavelet frequency range
- `type::Symbol=:nd`: inverse CWD reconstruction method:
    - `:pd`: PenroseDelta
    - `:nd`: NaiveDelta (default)
    - `:df`: DualFrames

# Returns

- `Vector{Float64}`: artifact-corrected signal of the same length as `s`

# Throws
- `ArgumentError`: if `fs < 1`, or if `tseg`/`fseg` are outside valid ranges

# See also

[`artrem_cwd(::NeuroAnalyzer.NEURO)`](@ref), [`artrem_cwd!`](@ref)
"""
function artrem_cwd(
    s::AbstractVector,
    t::AbstractVector;
    fs::Int64,
    wt::T = wavelet(Morlet(2π), β = 2),
    tseg::Tuple{Real, Real},
    fseg::Tuple{Real, Real},
    type::Symbol = :nd,
)::Vector{Float64} where {T <: CWT}

    fs >= 1 || throw(ArgumentError("fs must be ≥ 1."))   # was: !(fs >= 1) && throw(...)

    # compute the wavelet frequency axis, then validate the segment bounds
    f = cwtfrq(s; fs=fs, wt=wt)
    _check_tuple(tseg, (t[1], t[end]), "tseg")
    _check_tuple(fseg, (f[1], f[end]), "fseg")

    # forward CWD → complex coefficient matrix (frequencies × time)
    coef = cwd(s; wt=wt)

    # find the nearest frequency and time indices for the artifact window
    f_idx1 = vsearch(fseg[1], f)
    f_idx2 = vsearch(fseg[2], f)
    t_idx1 = vsearch(tseg[1], t)
    t_idx2 = vsearch(tseg[2], t)

    # zero the artifact region in the coefficient matrix
    coef[f_idx1:f_idx2, t_idx1:t_idx2] .= 0

    # reconstruct the clean signal via inverse CWD
    return vec(icwd(coef; wt=wt, type=type))

end

"""
    artrem_cwd(obj; <keyword arguments>)

Remove an artifact from one channel and one epoch of a NEURO object using continuous wavelet decomposition (CWD).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: channel name; must resolve to exactly one channel
- `ep::Int64`: epoch index; must be a valid epoch number
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=2)`: wavelet to use; see ContinuousWavelets.jl for available wavelets
- `tseg::Tuple{Real, Real}`: time window of the artifact `(t_start, t_end)` in seconds; must be a strictly ascending pair within `[t[1], t[end]]`
- `fseg::Tuple{Real, Real}`: frequency window of the artifact `(f_low, f_high)` in Hz; must be a strictly ascending pair within the wavelet frequency range
- `type::Symbol=:nd`: inverse CWD reconstruction method:
    - `:pd`: PenroseDelta
    - `:nd`: NaiveDelta (default)
    - `:df`: DualFrames

# Returns

- `NeuroAnalyzer.NEURO`: new object with the artifact removed from the specified channel and epoch

# Throws

- `ArgumentError`: if `ch` does not resolve to exactly one channel, or `ep` is out of range

# See also

[`artrem_cwd!`](@ref), [`artrem_cwd(::AbstractVector, ::AbstractVector)`](@ref)
"""
function artrem_cwd(
    obj::NeuroAnalyzer.NEURO;
    ch::String,
    ep::Int64,
    wt::T = wavelet(Morlet(2π), β = 2),
    tseg::Tuple{Real, Real},
    fseg::Tuple{Real, Real},
    type::Symbol = :nd,
)::NeuroAnalyzer.NEURO where {T <: CWT}

    # resolve channel names to integer indices
    ch = get_channel(obj, ch = ch)
    !(length(ch) == 1) && throw(ArgumentError("ch must resolve to exactly one channel."))
    ch = ch[1]

    _check_epochs(obj, ep)

    obj_new = deepcopy(obj)
    _log_off()
    obj_new.data[ch, :, ep] = artrem_cwd(
        @view(obj.data[ch, :, ep]),
        obj.epoch_time,
        fs = sr(obj),
        wt = wt,
        tseg = tseg,
        fseg = fseg,
        type = type
    )
    _log_on()
    push!(obj_new.history, "artrem_cwd(OBJ, ch=$ch, ep=$ep, wt=$wt, tseg=$tseg, fseg=$fseg, type=$type)")

    return obj_new

end

"""
    artrem_cwd!(obj; <keyword arguments>)

Remove an artifact from one channel and one epoch of a NEURO object in-place using continuous wavelet decomposition (CWD).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; modified in-place
- `ch::String`: channel name; must resolve to exactly one channel
- `ep::Int64`: epoch index; must be a valid epoch number
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=2)`: wavelet to use; see ContinuousWavelets.jl for available wavelets
- `tseg::Tuple{Real, Real}`: time window of the artifact `(t_start, t_end)` in seconds; must be a strictly ascending pair within `[t[1], t[end]]`
- `fseg::Tuple{Real, Real}`: frequency window of the artifact `(f_low, f_high)` in Hz; must be a strictly ascending pair within the wavelet frequency range
- `type::Symbol=:nd`: inverse CWD reconstruction method:
    - `:pd`: PenroseDelta
    - `:nd`: NaiveDelta (default)
    - `:df`: DualFrames

# Returns

- `Nothing`

# Throws
- `ArgumentError`: if `ch` does not resolve to exactly one channel, or `ep` is out of range

# See also

[`artrem_cwd`](@ref)
"""
function artrem_cwd!(
    obj::NeuroAnalyzer.NEURO;
    ch::String,
    ep::Int64,
    wt::T = wavelet(Morlet(2π), β = 2),
    tseg::Tuple{Real, Real},
    fseg::Tuple{Real, Real},
    type::Symbol = :nd,
)::Nothing where {T <: CWT}

    obj_new = artrem_cwd(obj, ch = ch, ep = ep, wt = wt, tseg = tseg, fseg = fseg, type = type)
    obj.data = obj_new.data
    obj.history = obj_new.history

    return nothing

end
