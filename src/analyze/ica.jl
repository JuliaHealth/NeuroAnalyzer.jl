export ica
export ica_reconstruct
export ica_reconstruct!

"""
    ica(signal, n, tol, iter, f)

Calculate `n` first Independent Components.

# Arguments

- `signal::AbstractArray`
- `n::Int64`: number of PCs
- `tol::Float64=1.0e-6`: tolerance for ICA
- `iter::Int64=100`: maximum number of iterations
- `f::Symbol=:tanh`: neg-entropy functor (:tanh or :gaus)

# Returns

Named tuple containing:
- `ic::Array{Float64, 3}:`: IC(1)..IC(n) × epoch
- `ic_mw::Array{Float64, 3}:`: IC(1)..IC(n) × epoch
"""
function ica(signal::AbstractArray; n::Int64, tol::Float64=1.0e-6, iter::Int64=100, f::Symbol=:tanh)

    _check_var(f, [:tanh, :gaus], "f")
    n < 0 && throw(ArgumentError("n must be ≥ 1."))
    n > size(signal, 1) && throw(ArgumentError("n must be ≤ $(size(signal, 1))."))
    ch_n, _, ep_n = size(signal)
    ic = zeros(n, size(signal, 2), ep_n)
    ic_mw = zeros(ch_n, n, ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        f === :tanh && (M = @views MultivariateStats.fit(ICA, signal[:, :, ep_idx], n, tol=tol, maxiter=iter, fun=MultivariateStats.Tanh(1.0)))
        f === :gaus && (M = @views MultivariateStats.fit(ICA, signal[:, :, ep_idx], n, tol=tol, maxiter=iter, fun=MultivariateStats.Gaus()))

        n == size(signal, 1) && (mw = inv(M.W)')
        n < size(signal, 1) && (mw = pinv(M.W)')

        for idx in 1:n
            ic[idx, :, ep_idx] = @views MultivariateStats.predict(M, signal[:, :, ep_idx])[idx, :]
        end

        ic_mw[:, :, ep_idx] = mw
    end

    return (ic=ic, ic_mw=ic_mw)
end

"""
    ica_reconstruct(signal, ic, ic_mw, channel)

Reconstructs `signal` via removal of `ic` ICA components.

# Arguments

- `signal::AbstractArray`
- `ic::AbstractArray:`: IC(1)..IC(n) × epoch
- `ic_mw::AbstractArray:`: IC(1)..IC(n) × epoch
- `ic_idx::Union{Int64, Vector{Int64}, AbstractRange} - list of ICs to remove

# Returns

- `s_reconstructed::Array{Float64, 3}`
"""
function ica_reconstruct(signal::AbstractArray; ic::AbstractArray, ic_mw::AbstractArray, ic_idx::Union{Int64, Vector{Int64}, AbstractRange})

    typeof(ic_idx) <: AbstractRange && (ic_idx = collect(ic_idx))
    if typeof(ic_idx) == Vector{Int64}
        sort!(ic_idx)
        for idx in ic_idx
            (idx < 1 || idx > size(ic_mw, 2)) && throw(ArgumentError("ic_idx must be ≥ 1 and ≤ $(size(ic_mw, 2))"))
        end
    else
        (ic_idx < 1 || ic_idx > size(ic_mw, 2)) && throw(ArgumentError("ic_idx must be ≥ 1 and ≤ $(size(ic_mw, 2))"))
    end
    ic_removal = setdiff(1:size(ic_mw, 2), ic_idx)

    s_reconstructed = similar(signal)

    ep_n = size(signal, 3)

    @inbounds @simd for ep_idx in 1:ep_n
        s_reconstructed[:, :, ep_idx] = ic_mw[:, ic_removal, ep_idx] * ic[ic_removal, :, ep_idx]
    end

    return s_reconstructed
end

"""
    ica(obj; <keyword arguments>)

Perform independent component analysis (ICA).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `n::Int64`: number of ICs
- `tol::Float64=1.0e-6`: tolerance for ICA
- `iter::Int64=100`: maximum number of iterations
- `f::Symbol=:tanh`: neg-entropy functor: :tanh, :gaus

# Returns

Named tuple containing:
- `ic::Array{Float64, 3}`: IC(1)..IC(n) × epoch (W * data)
- `ic_mw::Array{Float64, 3}`: IC(1)..IC(n) × epoch inv(W)
"""
function ica(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj)), n::Int64, tol::Float64=1.0e-6, iter::Int64=100, f::Symbol=:tanh)

    _check_channels(obj, channel)

    ic, ic_mw = @views ica(obj.data[channel, :, :], n=n, tol=tol, iter=iter, f=f)

    return (ic=ic, ic_mw=ic_mw)
end

"""
    ica_reconstruct(obj; channel, ic_idx)

Reconstruct signals using embedded ICA components (`:ic` and `:ic_mw`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `ic_idx::Union{Int64, Vector{Int64}, AbstractRange}`: list of ICs to remove

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function ica_reconstruct(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj)), ic_idx::Union{Int64, Vector{Int64}, AbstractRange})

    :ic in obj.header.component_names || throw(ArgumentError("OBJ does not contain :ic component. Perform ica() first."))
    :ic_mw in obj.header.component_names || throw(ArgumentError("OBJ does not contain :ic_mw component. Perform ica() first."))

    _check_channels(obj, channel)

    obj_new = deepcopy(obj)
    ic_idx = findfirst(isequal(:ic), obj.header.component_names)
    ic_mw_idx = findfirst(isequal(:ic_mw), obj.header.component_names)
    obj_new.data[channel, :, :] = @views ica_reconstruct(obj_new.data[channel, :, :], ic=obj_new.components[ic_idx], ic_mw=obj_new.components[ic_mw_idx], ic_idx=ic_idx)
    reset_components!(obj_new)
    push!(obj_new.header.history, "ica_reconstruct(OBJ, channel=$channel, ic_idx=$ic_idx)")

    return obj_new
end

"""
    ica_reconstruct!(obj; channel, ic_idx)

Reconstruct signals using embedded ICA components (`:ic` and `:ic_mw`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `ic_idx::Union{Int64, Vector{Int64}, AbstractRange} - list of ICs to remove
"""
function ica_reconstruct!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj)), ic_idx::Union{Int64, Vector{Int64}, AbstractRange})

    obj_tmp = ica_reconstruct(obj, channel=channel, ic_idx=ic_idx)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    reset_components!(obj)

    return nothing
end

"""
    ica_reconstruct(obj, ic, ic_mw; channel, ic_idx)

Reconstruct signals using external ICA components (`ic` and `ic_mw`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ic::Array{Float64, 3}`: IC(1)..IC(n) × epoch (W * data)
- `ic_mw::Array{Float64, 3}`: IC(1)..IC(n) × epoch inv(W)
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `ic_idx::Union{Int64, Vector{Int64}, AbstractRange} - list of ICs to remove

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function ica_reconstruct(obj::NeuroAnalyzer.NEURO, ic::Array{Float64, 3}, ic_mw::Array{Float64, 3}; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj)), ic_idx::Union{Int64, Vector{Int64}, AbstractRange})

    _check_channels(obj, channel)

    obj_new = deepcopy(obj)
    obj_new.data[channel, :, :] = @views ica_reconstruct(obj_new.data[channel, :, :], ic=ic, ic_mw=ic_mw, ic_idx=ic_idx)
    
    reset_components!(obj_new)
    push!(obj_new.header.history, "ica_reconstruct(OBJ, ic=$ic, ic_mw=$ic_mw, channel=$channel, ic_idx=$ic_idx)")

    return obj_new
end

"""
    ica_reconstruct!(obj, ic, ic_mw; channel, ic_idx)

Reconstruct signals using external ICA components (`ic` and `ic_mw`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ic::Array{Float64, 3}`: IC(1)..IC(n) × epoch (W * data)
- `ic_mw::Array{Float64, 3}`: IC(1)..IC(n) × epoch inv(W)
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `ic_idx::Union{Int64, Vector{Int64}, AbstractRange} - list of ICs to remove
"""
function ica_reconstruct!(obj::NeuroAnalyzer.NEURO, ic::Array{Float64, 3}, ic_mw::Array{Float64, 3}; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj)), ic_idx::Union{Int64, Vector{Int64}, AbstractRange})

    obj_tmp = ica_reconstruct(obj, ic, ic_mw, channel=channel, ic_idx=ic_idx)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    reset_components!(obj)

    return nothing
end
