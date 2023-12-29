export ica_decompose
export ica_reconstruct
export ica_reconstruct!
export ica_remove
export ica_remove!

"""
    ica_decompose(s; <keyword arguments>)

Calculate `n` first Independent Components using FastICA algorithm.

# Arguments

- `s::AbstractMatrix`
- `n::Int64`: number of ICs
- `iter::Int64=100`: maximum number of iterations
- `f::Symbol=:tanh`: neg-entropy functor:
    - `:tanh`
    - `:gaus`

# Returns

Named tuple containing:
- `ic::Matrix{Float64}`: components IC(1)..IC(n) (W * data), components are sorted by decreasing variance
- `ic_mw::Matrix{Float64}`: weighting matrix IC(1)..IC(n) (inv(W))
"""
function ica_decompose(s::AbstractMatrix; n::Int64, iter::Int64=100, f::Symbol=:tanh)

    _check_var(f, [:tanh, :gaus], "f")
    @assert n >= 1 "n must be ≥ 1."
    @assert n <= size(s, 1) "n must be ≤ $(size(s, 1))."

    f === :tanh && (f = MultivariateStats.Tanh(1.0))
    f === :gaus && (f = MultivariateStats.Gaus())

    Random.seed!(1234)

    tol = [0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0.5, 0.9, 0.99]
    M = nothing

    _warn("The input signal should be cleaned from major artifacts and HP filtered at 1-2 Hz prior to ICA decomposition.")
    _info("Attempting to calculate $n components")
    _info("Training will end when W change = $(tol[end]) or after $(iter * length(tol)) steps")
    _info("Data will be demeaned and pre-whitened")

    final_tol = nothing

    # initialize progress bar
    progress_bar == true && (progbar = Progress(iter * length(tol), dt=1, barlen=20, color=:white))

    @inbounds for tol_idx in eachindex(tol)
        for _ in 1:iter
            err = nothing
            try
                M = MultivariateStats.fit(ICA, s, n, maxiter=iter, tol=tol[tol_idx], fun=f)
            catch err
            end
            # if typeof(err) != MultivariateStats.ConvergenceException{Float64}
            if err === nothing
                # @info "Iteration: $iter_idx convergence error: $(err.lastchange)."
                final_tol = tol[tol_idx]
                break
            end

            # update progress bar
            progress_bar == true && next!(progbar)

        end
        final_tol !== nothing && break
    end

    if M === nothing
        _warn("The target sources could not be find.")
        return nothing
    end
    
    _info("Converged at: $final_tol")
    
    # inverse or pseudoinverse the weighting matrix
    ic_mw = n == size(s, 1) ? inv(M.W)' : pinv(M.W)'

    ic = MultivariateStats.predict(M, s)

    return (ic=ic, ic_mw=ic_mw[:, :])

end

"""
    ica_decompose(obj; <keyword arguments>)

Perform independent component analysis (ICA) using FastICA algorithm.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `n::Int64=length(ch)`: number of ICs, default is the number of channels
- `iter::Int64=100`: maximum number of iterations
- `f::Symbol=:tanh`: neg-entropy functor:
    - `:tanh`
    - `:gaus`

# Returns

Named tuple containing:
- `ic::Matrix{Float64}`: components IC(1)..IC(n) (W * data), components are sorted by decreasing variance
- `ic_mw::Matrix{Float64}`: weighting matrix IC(1)..IC(n) (inv(W))
- `ic_var::Vector{Float64}`: variance of components
"""
function ica_decompose(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), n::Int64=length(ch), iter::Int64=100, f::Symbol=:tanh)

    _check_channels(obj, ch)
    @assert nepochs(obj) == 1 "ica_decompose() should be applied to a continuous signal."

    signal_len(obj) / sr(obj) <= 10 && _warn("For ICA decomposition the signal length should be >10 seconds.")

    ic, ic_mw = @views ica_decompose(obj.data[ch, :, 1], n=n, iter=iter, f=f)

    v = var(obj.data[ch, :, 1])
    ic_var = ones(n)
    @inbounds for ic_idx in 1:n
        ic_back = @views ic_mw[:, ic_idx] * ic[ic_idx, :][:, :]'
        ic_var[ic_idx] = @views 100 * (1 - var(obj.data[ch, :, 1] - ic_back) / v)
    end

    # sort components by decreasing variance
    ic_var_idx = reverse(sortperm(ic_var))
    ic = ic[ic_var_idx, :]
    ic_var = ic_var[ic_var_idx]

    for ic_idx in 1:n
        _info("Component $(lpad(ic_idx, 2)): percent variance accounted for: $(round(ic_var[ic_idx], digits=2))")
    end

    return (ic=ic, ic_mw=ic_mw, ic_var=ic_var)

end

"""
    ica_reconstruct(; ic, ic_mw, ic_idx)

Reconstruct signal using ICA components.

# Arguments

- `ic::Matrix{Float64}`: components IC(1)..IC(n)
- `ic_mw::Matrix{Float64}`: weighting matrix IC(1)..IC(n)
- `ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange} - list of ICs to remove or keep
- `keep::Bool=false`: if `true`, then the ICs are kept instead of removed

# Returns

- `s_new::Matrix{Float64}`: reconstructed signal
"""
function ica_reconstruct(; ic::Matrix{Float64}, ic_mw::Matrix{Float64}, ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange}, keep::Bool=false)

    typeof(ic_idx) <: AbstractRange && (ic_idx = collect(ic_idx))
    @assert size(ic, 1) == size(ic_mw, 2) "ic and ic_mw dimensions do not match (ic: $(size(ic)), ic_mw: $(size(ic_mw)))."

    if typeof(ic_idx) == Vector{Int64}
        sort!(ic_idx)
        for idx in ic_idx
            @assert !(idx < 1 || idx > size(ic_mw, 2)) "ic_idx must be in [1, $(size(ic_mw, 2))]."
        end
    else
        @assert !(ic_idx < 1 || ic_idx > size(ic_mw, 2)) "ic_idx must be in [1, $(size(ic_mw, 2))]."
    end

    if keep == false
        ic_idx = setdiff(1:size(ic_mw, 2), ic_idx)
        s_new = @views ic_mw[:, ic_idx] * ic[ic_idx, :]
    else
        s_new = @views ic_mw[:, ic_idx] * ic[ic_idx, :]'
    end

    return s_new

end

"""
    ica_reconstruct(obj; ch, ic_idx, keep)

Reconstruct signals using embedded ICA components (`:ic` and `:ic_mw`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange}`: list of ICs to remove or keep or keep
- `keep::Bool=false`: if `true`, then the ICs are kept instead of removed

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function ica_reconstruct(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange}, keep::Bool=false)

    _check_channels(obj, ch)
    @assert :ic in keys(obj.components) "OBJ does not contain :ic component. Perform ica_decompose() first."
    @assert :ic_mw in keys(obj.components) "OBJ does not contain :ic_mw component. Perform ica_decompose() first."

    return ica_reconstruct(obj, obj.components[:ic], obj.components[:ic_mw], ch=ch, ic_idx=ic_idx, keep=keep)

end

"""
    ica_reconstruct!(obj; ch, ic_idx, keep)

Reconstruct signals using embedded ICA components (`:ic` and `:ic_mw`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all channels
- `ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange} - list of ICs to remove or keep
- `keep::Bool=false`: if `true`, then the ICs are kept instead of removed
"""
function ica_reconstruct!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange}, keep::Bool=false)

    obj_new = ica_reconstruct(obj, ch=ch, ic_idx=ic_idx, keep=keep)
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components

    return nothing

end

"""
    ica_reconstruct(obj, ic, ic_mw; ch, ic_idx, keep)

Reconstruct signals using external ICA components.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ic::Matrix{Float64}`: components IC(1)..IC(n)
- `ic_mw::Matrix{Float64}`: weighting matrix IC(1)..IC(n)
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange} - list of ICs to remove or keep
- `keep::Bool=false`: if `true`, then the ICs are kept instead of removed

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function ica_reconstruct(obj::NeuroAnalyzer.NEURO, ic::Matrix{Float64}, ic_mw::Matrix{Float64}; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange}, keep::Bool=false)

    _check_channels(obj, ch)

    @assert nepochs(obj) == 1 "ica_reconstruct() should be applied to a continuous signal."

    obj_new = deepcopy(obj)
    obj_new.data[ch, :, 1] = @views ica_reconstruct(ic=ic, ic_mw=ic_mw, ic_idx=ic_idx, keep=keep)[ch, :, :]
    
    reset_components!(obj_new)
    push!(obj_new.history, "ica_reconstruct(OBJ, ch=$ch, ic_idx=$ic_idx, keep=$keep)")

    return obj_new

end

"""
    ica_reconstruct!(obj, ic, ic_mw; ch, ic_idx, keep)

Reconstruct signals using external ICA components.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ic::Matrix{Float64}`: components IC(1)..IC(n)
- `ic_mw::Matrix{Float64}`: weighting matrix IC(1)..IC(n)
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all channels
- `ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange} - list of ICs to remove or keep
- `keep::Bool=false`: if `true`, then the ICs are kept instead of removed
"""
function ica_reconstruct!(obj::NeuroAnalyzer.NEURO, ic::Matrix{Float64}, ic_mw::Matrix{Float64}; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange}, keep::Bool=false)

    obj_new = ica_reconstruct(obj, ic, ic_mw, ch=ch, ic_idx=ic_idx, keep=keep)
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components

    return nothing

end

"""
    ica_remove(obj, ic, ic_mw; ch, ic_idx)

Remove external ICA components from the signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ic::Matrix{Float64}`: components IC(1)..IC(n)
- `ic_mw::Matrix{Float64}`: weighting matrix IC(1)..IC(n)
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange} - list of ICs to remove

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function ica_remove(obj::NeuroAnalyzer.NEURO, ic::Matrix{Float64}, ic_mw::Matrix{Float64}; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange})

    _check_channels(obj, ch)

    @assert nepochs(obj) == 1 "ica_remove() should be applied to a continuous signal."

    obj_new = deepcopy(obj)

    @inbounds for ica_idx in eachindex(ic_idx)
        Threads.@threads for ch_idx in eachindex(ch)
            obj_tmp = ica_reconstruct(obj, ic, ic_mw, ch=ch[ch_idx], ic_idx=ic_idx[ica_idx], keep=true)
            obj_new.data[ch[ch_idx], :, 1] = @views obj_new.data[ch[ch_idx], :, 1] - obj_tmp.data[ch[ch_idx], :, 1]
        end
    end
    
    reset_components!(obj_new)
    push!(obj_new.history, "ica_remove(OBJ, ch=$ch, ic_idx=$ic_idx)")

    return obj_new

end

"""
    ica_remove!(obj, ic, ic_mw; ch, ic_idx)

Remove external ICA components from the signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ic::Matrix{Float64}`: components IC(1)..IC(n)
- `ic_mw::Matrix{Float64}`: weighting matrix IC(1)..IC(n)
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all channels
- `ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange} - list of ICs to remove or keep
"""
function ica_remove!(obj::NeuroAnalyzer.NEURO, ic::Matrix{Float64}, ic_mw::Matrix{Float64}; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange}, keep::Bool=false)

    obj_new = ica_remove(obj, ic, ic_mw, ch=ch, ic_idx=ic_idx)
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components

    return nothing

end

"""
    ica_remove(obj, ic; ch, ic_idx)

Remove embedded ICA components (`:ic` and `:ic_mw`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange} - list of ICs to remove

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function ica_remove(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange})

    _check_channels(obj, ch)
    @assert :ic in keys(obj.components) "OBJ does not contain :ic component. Perform ica_decompose() first."
    @assert :ic_mw in keys(obj.components) "OBJ does not contain :ic_mw component. Perform ica_decompose() first."

    return ica_remove(obj, obj.components[:ic], obj.components[:ic_mw], ch=ch, ic_idx=ic_idx)

end

"""
    ica_remove!(obj; ch, ic_idx)

Remove embedded ICA components (`:ic` and `:ic_mw`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all channels
- `ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange} - list of ICs to remove or keep
"""
function ica_remove!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange})

    obj_new = ica_remove(obj, ch=ch, ic_idx=ic_idx)
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components

    return nothing

end
