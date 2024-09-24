export cph

"""
    cph(s1, s2; <keyword arguments>)

Calculate cross-phases.

# Arguments

- `s1::AbstractVector`
- `s2::AbstractVector`
- `fs::Int64`: sampling rate

# Returns

Named tuple containing:
- `ph::Vector{Float64}`: cross-power spectrum phase (in radians)
- `f::Vector{Float64}`: cross-power spectrum frequencies
"""
function cph(s1::AbstractVector, s2::AbstractVector; fs::Int64)::NamedTuple{(:ph, :f), Tuple{Vector{Float64}, Vector{Float64}}}

    @assert fs >= 1 "fs must be ≥ 1."
    @assert length(s1) == length(s2) "s1 and s2 must have the same length."

    p = mt_cross_power_spectra(hcat(s1, s2)', fs=fs)
    ph = angle.(imag.(p.power))[1, 2, :]
    f = Vector(p.freq)

    return (ph=ph, f=f)

end

"""
    cph(s; <keyword arguments>)

Calculate cross-phases.

# Arguments

- `s::AbstractArray`
- `fs::Int64`: sampling rate

# Returns

Named tuple containing:
- `ph::Array{Float64, 4}`: cross-power spectrum phase (in radians)
- `f::Vector{Float64}`: cross-power spectrum frequencies
"""
function cph(s::AbstractArray; fs::Int64)::NamedTuple{(:ph, :f), Tuple{Array{Float64, 4}, Vector{Float64}}}

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    _, f = @views cph(s[1, :, 1], s[1, :, 1], fs=fs)
    ph = zeros(ch_n, ch_n, length(f), ep_n)

    # initialize progress bar
    progress_bar && (progbar = Progress(ep_n * ch_n, dt=1, barlen=20, color=:white))

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx1 in 1:ch_n
           for ch_idx2 in 1:ch_idx1
                ph[ch_idx1, ch_idx2, :, ep_idx], _ = @views cph(s[ch_idx1, :, ep_idx], s[ch_idx2, :, ep_idx], fs=fs)
            end

        # update progress bar
        progress_bar && next!(progbar)
        end
    end

    @inbounds for cph_idx in axes(ph, 3)
        Threads.@threads for ep_idx in 1:ep_n
            for ch_idx1 in 1:(ch_n - 1)
                for ch_idx2 in (ch_idx1 + 1):ch_n
                    ph[ch_idx1, ch_idx2, cph_idx, ep_idx] = @views ph[ch_idx2, ch_idx1, cph_idx, ep_idx]
                end
            end
        end
    end

    return (ph=ph, f=f)

end

"""
    cph(s1, s2; <keyword arguments>)

Calculate cross-phases.

# Arguments

- `s1::AbstractArray`
- `s2::AbstractArray`
- `fs::Int64`: sampling rate

# Returns

Named tuple containing:
- `ph::Array{Float64, 4}`: cross-power spectrum phase (in radians)
- `f::Vector{Float64}`: cross-power spectrum frequencies
"""
function cph(s1::AbstractArray, s2::AbstractArray; fs::Int64)::NamedTuple{(:ph, :f), Tuple{Array{Float64, 4}, Vector{Float64}}}

    @assert size(s1) == size(s2) "s1 and s2 must have the same size."

    ch_n = size(s1, 1)
    ep_n = size(s1, 3)

    _, f = @views cph(s1[1, :, 1], s2[1, :, 1], fs=fs)
    ph = zeros(ch_n, length(f), ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            ph[ch_idx, :, ep_idx], _ = @views cph(s1[ch_idx, :, ep_idx], s2[ch_idx, :, ep_idx], fs=fs)
        end
    end

    return (ph=ph, f=f)

end

"""
    cph(obj; <keyword arguments>)

Calculate cross-phases.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}}`: channel name or list of channel names

# Returns

Named tuple containing:
- `ph::Array{Float64, 4}`: cross-power spectrum phase (in radians)
- `f::Vector{Float64, 4}`: cross-power spectrum frequencies
"""
function cph(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}})::NamedTuple{(:ph, :f), Tuple{Array{Float64, 4}, Vector{Float64}}}

    ch = get_channel(obj, ch=ch)

    ph, f = cph(obj.data[ch, :, :], fs=sr(obj))

    return (ph=ph, f=f)

end

"""
    cph(obj1, obj2; <keyword arguments>)

Calculate cross-phases.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch1::Union{String, Vector{String}}: list of channels
- `ch2::Union{String, Vector{String}}: list of channels
- `ep1::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj1))`: default use all epochs
- `ep2::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj2))`: default use all epochs

# Returns

Named tuple containing:
- `ph::Array{Float64, 3}`: cross-power spectrum phase (in radians)
- `f::Vector{Float64}`: cross-power spectrum frequencies
"""
function cph(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch1::Union{String, Vector{String}}, ch2::Union{String, Vector{String}}, ep1::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj1)), ep2::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj2)))::NamedTuple{(:ph, :f), Tuple{Array{Float64, 4}, Vector{Float64}}}

    @assert sr(obj1) == sr(obj2) "OBJ1 and OBJ2 must have the same sampling rate."
    @assert length(ch1) == length(ch2) "ch1 and ch2 must have the same length."
    @assert length(ep1) == length(ep2) "ep1 and ep2 must have the same length."
    @assert epoch_len(obj1) == epoch_len(obj2) "OBJ1 and OBJ2 must have the same epoch lengths."

    ch1 = get_channel(obj1, ch=ch1)
    ch2 = get_channel(obj2, ch=ch2)
    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    length(ep1) == 1 && (ep1 = [ep1])
    length(ep2) == 1 && (ep2 = [ep2])

    ph, f = @views cph(obj1.data[ch1, :, ep1], obj2.data[ch2, :, ep2], fs=sr(obj1))

    return (ph=ph, f=f)

end
