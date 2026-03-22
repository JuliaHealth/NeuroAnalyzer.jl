export ica_decompose
export ica_reconstruct
export ica_reconstruct!
export ica_remove
export ica_remove!

"""
    ica_decompose(s; <keyword arguments>)

Decompose a signal into Independent Components (ICs) using the FastICA algorithm.

# Arguments

- `s::AbstractMatrix`: input signal, shape `(channels,  samples)`
- `n::Int`: number of independent components to extract
- `iter::Int=100`: max iterations per tolerance level
- `f::Symbol=:tanh`: nonlinear function for neg-entropy approximation (`:tanh` or `:gaus`)

# Returns

Named tuple:

- `ic::Matrix{Float64}`: independent components, shape `(n, samples)`
- `ic_mw::Matrix{Float64}`: weighting matrix, shape `(channels, n)`

# Throws

- `ArgumentError`: if `n` is not in `[1, size(s, 1)]` or if `f` is invalid

# See also

[`ica_decompose(::NeuroAnalyzer.NEURO)`](@ref)
"""
function ica_decompose(
    s::AbstractMatrix;
    n::Int64,
    iter::Int64 = 100,
    f::Symbol = :tanh
)::@NamedTuple{
    ic::Matrix{Float64},
    ic_mw::Matrix{Float64}
}

    # validate
    _check_var(f, [:tanh, :gaus], "f")
    n >= 1 || throw(ArgumentError("n must be ≥ 1."))
    n <= size(s, 1) || throw(ArgumentError("n must be ≤ number of channels."))

    # map symbols to MultivariateStats functors
    functor = f === :tanh ? MultivariateStats.Tanh(1.0) : MultivariateStats.Gaus()

    # ensure reproducibility for the random initialization in FastICA
    Random.seed!(1234)

    # tolerance schedule: try strictest first, fallback to looser if convergence fails
    tols = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.5, 0.9, 0.99]
    model = nothing
    final_tol = nothing

    _warn("Signal should be artifact-cleaned and HP filtered (1-2 Hz) before ICA.")
    _info("Attempting to calculate $n components across $(length(tols)) tolerance levels")
    _info("Training will end when W change = $(tol[end]) or after $(iter * length(tol)) steps")
    _info("Data will be demeaned and pre-whitened")


    M = nothing

    final_tol = nothing

    # initialize progress bar
    progbar = Progress(iter * length(tol), dt = 1, barlen = 20, color = :white, enabled = progress_bar)

    for tol in tols
        try
            # attempt fit with current tolerance
            model = MultivariateStats.fit(ICA, s, n; maxiter=iter, tol=tol, fun=functor)
            final_tol = tol
            println()
            break # exit loop if converged
        catch err
            # if it's not a convergence error, rethrow it; otherwise, update progress and continue
            !(err isa MultivariateStats.ConvergenceException) && rethrow(err)
            # skip progress for failed tolerance bracket
            update!(progbar, iter)
        end
    end

    model === nothing && throw(ErrorException("ICA failed to converge even at highest tolerance."))
    _info("Converged at tolerance: $final_tol")

    # W is the unmixing matrix; ic_mw is the mixing matrix (W⁻¹ or W⁺)
    # transpose used to align with signal reconstruction logic (channels × components)
    ic_mw = n == size(s, 1) ? inv(model.W)' : pinv(model.W)'
    ic_mw = Matrix(ic_mw)
    ic = MultivariateStats.predict(model, s)

    return (; ic, ic_mw)

end

"""
    ica_decompose(obj; <keyword arguments>)

Decompose selected channels of a NEURO object into Independent Components (ICs) using the FastICA algorithm. Sorts components by variance explained.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `n::Int`: number of independent components to extract
- `iter::Int=100`: max iterations per tolerance level
- `f::Symbol=:tanh`: nonlinear function for neg-entropy approximation (`:tanh` or `:gaus`)

# Returns

Named tuple:

- `ic::Matrix{Float64}`: independent components, shape `(n, samples)`
- `ic_mw::Matrix{Float64}`: weighting matrix, shape `(channels, n)`
- `ic_var::Vector{Float64}`: variance explained by each component

# Throws

- `ArgumentError`: if `obj` is not continuous or if `n` is invalid

# See also

[`ica_decompose(::AbstractMatrix)`](@ref)
"""
function ica_decompose(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    n::Int64 = length(ch),
    iter::Int64 = 100,
    f::Symbol = :tanh
)::@NamedTuple{
    ic::Matrix{Float64},
    ic_mw::Matrix{Float64},
    ic_var::Vector{Float64}
}

    # validate
    nepochs(obj) == 1 || throw(ArgumentError("ica_decompose() must be applied to continuous object."))
    signal_len(obj) / sr(obj) <= 10 && _warn("For ICA decomposition the signal length should be >10 seconds.")

    # resolve channel names to integer indices
    ch = get_channel(obj, ch = ch)
    length(ch) == 1 && (ch = ch[1])

    # perform decomposition on the selected slice
    ica_data = ica_decompose(@view(obj.data[ch, :, 1]), n = n, iter = iter, f = f)
    ic = ica_data.ic
    ic_mw = ica_data.ic_mw

    # calculate Variance Accounted For (VAF) per component
    total_var = var(@view(obj.data[ch, :, 1]))
    ic_var = Vector{Float64}(undef, n)
    
    for idx in 1:n
        # reconstruct signal using only the i-th component
        ic_back = @views ic_mw[:, idx] * ic[idx, :]'
        # VAF formula: 100 * (1 - var(residual) / var(original))
        ic_var[idx] = 100.0 * (1.0 - var(@view(obj.data[ch, :, 1]) .- ic_back) / total_var)
    end

    # sort components by descending variance
    p = sortperm(ic_var, rev=true)
    ic = ic[p, :]
    ic_var = ic_var[p]
    ic_mw = ic_mw[:, p]

    for i in 1:n
        _info("Component $(lpad(i, 2)): VAF = $(round(ic_var[i], digits=2))%")
    end

    return (; ic, ic_mw, ic_var)

end

"""
    ica_reconstruct(; <keyword arguments>)

Reconstruct a signal from independent components.

# Arguments

- `ic::Matrix{Float64}`: independent components, shape (n, samples)
- `ic_mw::Matrix{Float64}`: weighting matrix, shape (channels, n)
- `ic_idx::Union{Int64, Vector{Int64}, AbstractRange}`: indices of components to keep or remove
- `keep::Bool=false`: if `true`, keep specified components; otherwise, remove them

# Returns
- `Matrix{Float64}`: reconstructed signal, shape (channels, samples)

# Throws

- `ArgumentError`: if `ic_idx` is out of bounds or if dimensions of `ic` and `ic_mw` do not match

# See also

[`ica_reconstruct(::NeuroAnalyzer.NEURO)`](@ref)
"""
function ica_reconstruct(;
    ic::Matrix{Float64},
    ic_mw::Matrix{Float64},
    ic_idx::Union{Int64, Vector{Int64}, AbstractRange},
    keep::Bool = false
)::Matrix{Float64}

    # validate
    typeof(ic_idx) <: AbstractRange && (ic_idx = collect(ic_idx))
    size(ic, 1) == size(ic_mw, 2) || throw(ArgumentError("Dimension mismatch between ic ($(size(ic)))and ic_mw ($size(ic_mw)))."))

    # bounds check
    all(1 .<= idx_vec .<= size(ic_mw, 2)) || throw(ArgumentError("ic_idx must be in [1, $(size(ic_mw, 2))]."))

    # determine which indices to actually use for reconstruction
    target_idx = keep ? idx_vec : setdiff(1:size(ic_mw, 2), idx_vec)

    # reconstruction: Signal = MixingMatrix[:, target] * Components[target, :]
    return @views ic_mw[:, target_idx] * ic[target_idx, :]

end

"""
    ica_reconstruct(obj, ic, ic_mw; <keyword arguments>)

Reconstruct selected channels of a NEURO object from independent components.

# Arguments
- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `ic_idx::Union{Int64, Vector{Int64}, AbstractRange}`: indices of components to keep or remove
- `ic::Matrix{Float64}`: independent components, shape (n, samples)
- `ic_mw::Matrix{Float64}`: weighting matrix, shape (channels, n)
- `keep::Bool=false`: if `true`, keep specified components; otherwise, remove them

# Returns

- `NeuroAnalyzer.NEURO`: reconstructed NEURO object

# Throws

- `ArgumentError`: if `obj` is not continuous or if `ic_idx` is invalid

# See also

[`ica_reconstruct(::Matrix{Float64})`](@ref), [`ica_reconstruct!(::NeuroAnalyzer.NEURO)`](@ref)
"""
function ica_reconstruct(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    ic_idx::Union{Int64, Vector{Int64}, AbstractRange},
    ic::Matrix{Float64},
    ic_mw::Matrix{Float64},
    keep::Bool = false
)::NeuroAnalyzer.NEURO

    # validate
    nepochs(obj) == 1 || throw(ArgumentError("ica_reconstruct() must be applied to continuous object."))

    # resolve channel names to integer indices
    ch = get_channel(obj, ch = ch)
    length(ch) == 1 && (ch = ch[1])

    # reconstruction
    obj_new = deepcopy(obj)
    obj_new.data[ch, :, 1] = ica_reconstruct(ic = ic, ic_mw = ic_mw, ic_idx = ic_idx, keep = keep)[ch, :]

    push!(obj_new.history, "ica_reconstruct(OBJ, ch=$ch, ic_idx=$ic_idx, keep=$keep)")

    return obj_new

end

"""
    ica_reconstruct!(obj, ic, ic_mw; <keyword arguments>)

Reconstruct selected channels of a NEURO object in-place from independent components

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; modified in-place
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `ic_idx::Union{Int64, Vector{Int64}, AbstractRange}`: indices of components to keep or remove
- `ic::Matrix{Float64}`: independent components, shape (n, samples)
- `ic_mw::Matrix{Float64}`: weighting matrix, shape (channels, n)
- `keep::Bool=false`: if `true`, keep specified components; otherwise, remove them

# Returns

- `Nothing`

# See also

[`ica_reconstruct(::Matrix{Float64})`](@ref), [`ica_reconstruct(::NeuroAnalyzer.NEURO)`](@ref)
"""
function ica_reconstruct!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    ic_idx::Union{Int64, Vector{Int64}, AbstractRange},
    ic::Matrix{Float64},
    ic_mw::Matrix{Float64},
    keep::Bool = false
)::Nothing

    obj_new = ica_reconstruct(obj, ch = ch, ic_idx = ic_idx, ic = ic, ic_mw = ic_mw, keep = keep)
    obj.data = obj_new.data
    obj.history = obj_new.history

    return nothing

end

"""
    ica_remove(obj, ic, ic_mw; <keyword arguments>)

Remove independent components from a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s), default is all channels
- `ic_idx::Union{Int64, Vector{Int64}, AbstractRange}`: indices of components to keep or remove
- `ic::Matrix{Float64}`: independent components, shape (n, samples)
- `ic_mw::Matrix{Float64}`: weighting matrix, shape (channels, n)

# Returns

- `NeuroAnalyzer.NEURO`: reconstructed NEURO object

# Throws

- `ArgumentError`: if `obj` is not continuous or if `ic_idx` is invalid

# See also

[`ica_remove!`](@ref)
"""
function ica_remove(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    ic_idx::Union{Int64, Vector{Int64}, AbstractRange},
    ic::Matrix{Float64},
    ic_mw::Matrix{Float64}
)::NeuroAnalyzer.NEURO

    # validate
    nepochs(obj) == 1 || throw(ArgumentError("ica_remove() must be applied to continuous object."))

    # resolve channel names to integer indices
    ch = get_channel(obj, ch = ch)
    length(ch) == 1 && (ch = ch[1])
    ch_n = length(ch)

    # number of IC components
    ic_n = length(ic_idx)

    obj_new = deepcopy(obj)

    # calculate over components and channels
    @inbounds Threads.@threads :static for idx in CartesianIndices((ic_n, ch_n))
        ic_idx, ch_idx = idx[1], idx[2]
        obj_tmp = ica_reconstruct(
            obj,
            ch = labels(obj)[ch[ch_idx]],
            ic_idx = ic_idx[ica_idx],
            ic = ic,
            ic_mw = ic_mw,
            keep = true
        )
        obj_new.data[ch[ch_idx], :, 1] = @views obj_new.data[ch[ch_idx], :, 1] - obj_tmp.data[ch[ch_idx], :, 1]
    end

    push!(obj_new.history, "ica_remove(OBJ, ch=$ch, ic_idx=$ic_idx)")

    return obj_new

end

"""
    ica_remove!(obj, ic, ic_mw; <keyword arguments>)

Remove independent components from a NEURO object in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; modified in-place
- `ch::Union{String, Vector{String}, Regex}`: channel name(s), default is all channels
- `ic_idx::Union{Int64, Vector{Int64}, AbstractRange}`: indices of components to keep or remove
- `ic::Matrix{Float64}`: independent components, shape (n, samples)
- `ic_mw::Matrix{Float64}`: weighting matrix, shape (channels, n)

# Returns

- `Nothing`

# Throws

- `ArgumentError`: if `obj` is not continuous or if `ic_idx` is invalid

# See also

[`ica_remove`](@ref)
"""
function ica_remove!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    ic_idx::Union{Int64, Vector{Int64}, AbstractRange},
    ic::Matrix{Float64},
    ic_mw::Matrix{Float64}
)::Nothing

    obj_new = ica_remove(
        obj,
        ic,
        ic_mw,
        ch = ch,
        ic_idx = ic_idx,
        ic = ic,
        ic_mw = ic_mw
    )
    obj.data = obj_new.data
    obj.history = obj_new.history

    return nothing

end
