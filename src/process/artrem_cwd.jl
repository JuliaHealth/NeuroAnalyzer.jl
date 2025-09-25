export artrem_cwd
export artrem_cwd!

"""
    artrem_cwd(s, t; <keyword arguments>)

Remove artifacts using continuous wavelet decomposition (CWD).

# Arguments

- `s::AbstractVector`
- `t::AbstractVector`: time points
- `fs::Int64`: sampling rate
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets
- `tseg::Tuple{Real, Real}`: artifact time location
- `fseg::Tuple{Real, Real}`: artifact frequency location
- `type::Symbol=:nd`: inverse style type:
    - `:pd`: PenroseDelta
    - `:nd`: NaiveDelta
    - `:df`: DualFrames

# Returns

- `s_new::Vector{Float64}`
"""
function artrem_cwd(s::AbstractVector, t::AbstractVector; fs::Int64, wt::T=wavelet(Morlet(2π), β=2), tseg::Tuple{Real, Real}, fseg::Tuple{Real, Real}, type::Symbol=:nd) where {T <: CWT}

    @assert fs >= 1 "fs must be ≥ 1."

    f = cwtfrq(s, fs=fs, wt=wt)
    _check_tuple(tseg, "tseg", (t[1], t[end]))
    _check_tuple(fseg, "fseg", (f[1], f[end]))

    # perform continuous wavelet transformation
    s_new = cwd(s, wt=wt)

    # locate artifact
    f_idx1 = vsearch(fseg[1], f)
    f_idx2 = vsearch(fseg[2], f)
    t_idx1 = vsearch(tseg[1], t)
    t_idx2 = vsearch(tseg[2], t)

    # zero components
    s_new[f_idx1:f_idx2, t_idx1:t_idx2] .= 0

    # reconstruct
    s_new = vec(icwd(s_new, wt=wt, type=type))

    return s_new

end

"""
    artrem_cwd(obj; <keyword arguments>)

Remove artifacts using continuous wavelet decomposition (CWD).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::String`: channel name
- `ep::Int64`
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets
- `tseg::Tuple{Real, Real}`: artifact time location
- `fseg::Tuple{Real, Real}`: artifact frequency location
- `type::Symbol=:nd`: inverse style type:
    - `:pd`: PenroseDelta
    - `:nd`: NaiveDelta
    - `:df`: DualFrames

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function artrem_cwd(obj::NeuroAnalyzer.NEURO; ch::String, ep::Int64, wt::T=wavelet(Morlet(2π), β=2), tseg::Tuple{Real, Real}, fseg::Tuple{Real, Real}, type::Symbol=:nd) where {T <: CWT}

    ch = get_channel(obj, ch=ch)[1]
    _check_epochs(obj, ep)

    obj_new = deepcopy(obj)
    _log_off()
    obj_new.data[ch, :, ep] = @views artrem_cwd(obj.data[ch, :, ep], obj.epoch_time, fs=sr(obj), wt=wt, tseg=tseg, fseg=fseg, type=type)
    _log_on()
    reset_components!(obj_new)
    push!(obj_new.history, "artrem_cwd(OBJ, ch=$ch, ep=$ep, wt=$wt, tseg=$tseg, fseg=$fseg, type=$type)")

    return obj_new

end

"""
    artrem_cwd!(obj; <keyword arguments>)

Remove artifacts using continuous wavelet decomposition (CWD).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::String`: channel name
- `ep::Int64`
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets
- `tseg::Tuple{Real, Real}`: artifact time location
- `fseg::Tuple{Real, Real}`: artifact frequency location
- `type::Symbol=:nd`: inverse style type:
    - `:pd`: PenroseDelta
    - `:nd`: NaiveDelta
    - `:df`: DualFrames
"""
function artrem_cwd!(obj::NeuroAnalyzer.NEURO; ch::String, ep::Int64, wt::T=wavelet(Morlet(2π), β=2), tseg::Tuple{Real, Real}, fseg::Tuple{Real, Real}, type::Symbol=:nd) where {T <: CWT}

    obj_new = artrem_cwd(obj, ch=ch, ep=ep, wt=wt, tseg=tseg, fseg=fseg, type=type)
    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history

    return nothing

end

