export ica
export ica_reconstruct
export ica_reconstruct!

"""
    ica(s; <keyword arguments>)

Calculate `n` first Independent Components.

# Arguments

- `s::AbstractMatrix`
- `n::Int64`: number of ICs
- `tol::Float64=1.0e-6`: tolerance for ICA
- `iter::Int64=100`: maximum number of iterations
- `f::Symbol=:tanh`: neg-entropy functor:
    - `:tanh`
    - `:gaus`

# Returns

Named tuple containing:
- `ic::Array{Float64, 3}:`: IC(1)..IC(n) × epoch
- `ic_mw::Array{Float64, 3}:`: IC(1)..IC(n) × epoch
"""
function ica(s::AbstractMatrix; n::Int64, tol::Float64=1.0e-6, iter::Int64=100, f::Symbol=:tanh)

    _check_var(f, [:tanh, :gaus], "f")
    n < 0 && throw(ArgumentError("n must be ≥ 1."))
    n > size(s, 1) && throw(ArgumentError("n must be ≤ $(size(s, 1))."))

    ic = zeros(n, size(s, 2))

    f === :tanh && (M = @views MultivariateStats.fit(ICA, s[:, :], n, tol=tol, maxiter=iter, fun=MultivariateStats.Tanh(1.0)))
    f === :gaus && (M = @views MultivariateStats.fit(ICA, s[:, :], n, tol=tol, maxiter=iter, fun=MultivariateStats.Gaus()))

    if n == size(s, 1)
        ic_mw = inv(M.W)'
    else
        ic_mw = pinv(M.W)'
    end

    @inbounds for idx in 1:n
        ic[idx, :] = @views MultivariateStats.predict(M, s[:, :])[idx, :]
    end

    return (ic=ic, ic_mw=ic_mw)

end

"""
    ica(s; <keyword arguments>)

Calculate `n` first Independent Components.

# Arguments

- `s::AbstractArray`
- `n::Int64`: number of ICs
- `tol::Float64=1.0e-6`: tolerance for ICA
- `iter::Int64=100`: maximum number of iterations
- `f::Symbol=:tanh`: neg-entropy functor:
    - `:tanh`
    - `:gaus`

# Returns

Named tuple containing:
- `ic::Array{Float64, 3}:`: IC(1)..IC(n) × epoch
- `ic_mw::Array{Float64, 3}:`: IC(1)..IC(n) × epoch
"""
function ica(s::AbstractArray; n::Int64, tol::Float64=1.0e-6, iter::Int64=100, f::Symbol=:tanh)

    _check_var(f, [:tanh, :gaus], "f")
    n < 0 && throw(ArgumentError("n must be ≥ 1."))
    n > size(s, 1) && throw(ArgumentError("n must be ≤ $(size(s, 1))."))

    ch_n, _, ep_n = size(s)
    ic = zeros(n, size(s, 2), ep_n)
    
    ic_mw = zeros(ch_n, n, ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        f === :tanh && (M = @views MultivariateStats.fit(ICA, s[:, :, ep_idx], n, tol=tol, maxiter=iter, fun=MultivariateStats.Tanh(1.0)))
        f === :gaus && (M = @views MultivariateStats.fit(ICA, s[:, :, ep_idx], n, tol=tol, maxiter=iter, fun=MultivariateStats.Gaus()))

        if n == size(s, 1)
            mw = inv(M.W)'
        else
            mw = pinv(M.W)'
        end

        for idx in 1:n
            ic[idx, :, ep_idx] = @views MultivariateStats.predict(M, s[:, :, ep_idx])[idx, :]
        end

        ic_mw[:, :, ep_idx] = mw
    end

    return (ic=ic, ic_mw=ic_mw)

end

"""
    ica_reconstruct(s; ic, ic_mw, ic_idx)

Reconstructs `s` via removal of `ic` ICA components.

# Arguments

- `s::AbstractArray`
- `ic::AbstractArray:`: IC(1)..IC(n) × epoch
- `ic_mw::AbstractArray:`: IC(1)..IC(n) × epoch
- `ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange} - list of ICs to remove

# Returns

- `s_reconstructed::Array{Float64, 3}`
"""
function ica_reconstruct(s::AbstractArray; ic::AbstractArray, ic_mw::AbstractArray, ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange})

    typeof(ic_idx) <: AbstractRange && (ic_idx = collect(ic_idx))

    if typeof(ic_idx) == Vector{Int64}
        sort!(ic_idx)
        for idx in ic_idx
            (idx < 1 || idx > size(ic_mw, 2)) && throw(ArgumentError("ic_idx must be in [1, $(size(ic_mw, 2))]."))
        end
    else
        (ic_idx < 1 || ic_idx > size(ic_mw, 2)) && throw(ArgumentError("ic_idx must be in [1, $(size(ic_mw, 2))]."))
    end

    ic_removal = setdiff(1:size(ic_mw, 2), ic_idx)

    s_reconstructed = similar(s)

    ep_n = size(s, 3)

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
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `n::Int64`: number of ICs
- `tol::Float64=1.0e-6`: tolerance for ICA
- `iter::Int64=100`: maximum number of iterations
- `f::Symbol=:tanh`: neg-entropy functor:
    - `:tanh`
    - `:gaus`

# Returns

Named tuple containing:
- `ic::Array{Float64, 3}`: IC(1)..IC(n) × epoch (W * data)
- `ic_mw::Array{Float64, 3}`: IC(1)..IC(n) × epoch inv(W)
"""
function ica(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)), n::Int64, tol::Float64=1.0e-6, iter::Int64=100, f::Symbol=:tanh)

    _check_channels(obj, ch)

    ic, ic_mw = @views ica(obj.data[ch, :, :], n=n, tol=tol, iter=iter, f=f)

    return (ic=ic, ic_mw=ic_mw)
end

"""
    ica_reconstruct(obj; ch, ic_idx)

Reconstruct signals using embedded ICA components (`:ic` and `:ic_mw`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange}`: list of ICs to remove

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function ica_reconstruct(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)), ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange})

    :ic in keys(obj.components) || throw(ArgumentError("OBJ does not contain :ic component. Perform ica() first."))
    :ic_mw in keys(obj.components) || throw(ArgumentError("OBJ does not contain :ic_mw component. Perform ica() first."))

    _check_channels(obj, ch)

    obj_new = deepcopy(obj)

    ic_idx = component_idx(obj, c=:ic)
    ic_mw_idx = component_idx(obj, c=:ic_mw)
    obj_new.data[ch, :, :] = @views ica_reconstruct(obj_new.data[ch, :, :], ic=obj_new.components[ic_idx], ic_mw=obj_new.components[ic_mw_idx], ic_idx=ic_idx)

    reset_components!(obj_new)
    push!(obj_new.header.history, "ica_reconstruct(OBJ, ch=$ch, ic_idx=$ic_idx)")

    return obj_new

end

"""
    ica_reconstruct!(obj; ch, ic_idx)

Reconstruct signals using embedded ICA components (`:ic` and `:ic_mw`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange} - list of ICs to remove
"""
function ica_reconstruct!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)), ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange})

    obj_tmp = ica_reconstruct(obj, ch=ch, ic_idx=ic_idx)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    obj.components = obj_tmp.components

    return nothing

end

"""
    ica_reconstruct(obj; ic, ic_mw; ch, ic_idx)

Reconstruct signals using external ICA components (`ic` and `ic_mw`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ic::Array{Float64, 3}`: IC(1)..IC(n) × epoch (W * data)
- `ic_mw::Array{Float64, 3}`: IC(1)..IC(n) × epoch inv(W)
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange} - list of ICs to remove

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function ica_reconstruct(obj::NeuroAnalyzer.NEURO, ic::Array{Float64, 3}, ic_mw::Array{Float64, 3}; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)), ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange})

    _check_channels(obj, ch)

    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = @views ica_reconstruct(obj_new.data[ch, :, :], ic=ic, ic_mw=ic_mw, ic_idx=ic_idx)
    
    reset_components!(obj_new)
    push!(obj_new.header.history, "ica_reconstruct(OBJ, ic=$ic, ic_mw=$ic_mw, ch=$ch, ic_idx=$ic_idx)")

    return obj_new

end

"""
    ica_reconstruct!(obj, ic, ic_mw; ch, ic_idx)

Reconstruct signals using external ICA components (`ic` and `ic_mw`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ic::Array{Float64, 3}`: IC(1)..IC(n) × epoch (W * data)
- `ic_mw::Array{Float64, 3}`: IC(1)..IC(n) × epoch inv(W)
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange} - list of ICs to remove
"""
function ica_reconstruct!(obj::NeuroAnalyzer.NEURO, ic::Array{Float64, 3}, ic_mw::Array{Float64, 3}; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)), ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange})

    obj_tmp = ica_reconstruct(obj, ic, ic_mw, ch=ch, ic_idx=ic_idx)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    obj.components = obj_tmp.components

    return nothing

end
