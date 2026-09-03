export ica_decompose
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

- `ic::Matrix{Float64}`: independent components, shape (n, samples)
- `ic_mw::Matrix{Float64}`: weighting (unmixing) matrix, shape (channels, n)
"""
function ica_decompose(
    s::AbstractMatrix;
    n::Int64,
    iter::Int64 = 100,
    f::Symbol = :tanh,
)::@NamedTuple{
    ic::Matrix{Float64},
    ic_mw::Matrix{Float64},
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
    tols = [1.0e-6, 1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2, 0.1, 0.5, 0.9, 0.99]
    model = nothing
    final_tol = nothing

    _warn("Signal should be artifact-cleaned and HP filtered (1-2 Hz) before ICA.")
    _info("Attempting to calculate $n components across $(length(tols)) tolerance levels")
    _info(
        "Training will end when W change = $(tols[end]) or after $(iter * length(tols)) steps",
    )
    _info("Data will be demeaned and pre-whitened")

    final_tol = nothing

    # initialize progress bar
    progbar = Progress(
        iter * length(tols);
        dt = 1,
        barlen = 20,
        color = :white,
        enabled = progress_bar,
    )

    for tol in tols
        try
            # attempt fit with current tolerance
            model =
                MultivariateStats.fit(ICA, s, n; maxiter = iter, tol = tol, fun = functor)
            final_tol = tol
            println()
            break # exit loop if converged
        catch err
            # if it's not a convergence error, rethrow it; otherwise, update progress and continue
            !(err isa MultivariateStats.ConvergenceException) && rethrow(err)
            # skip progress for failed tolerance bracket
            # update progress bar
            progress_bar && next!(progbar)
        end
    end

    model === nothing &&
        throw(ErrorException("ICA failed to converge even at highest tolerance."))
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

- `ic::Matrix{Float64}`: independent components, shape (n, samples)
- `ic_mw::Matrix{Float64}`: weighting (unmixing) matrix, shape (channels, n)
- `ic_var::Vector{Float64}`: variance explained by each component
"""
function ica_decompose(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    n::Int64 = length(ch),
    iter::Int64 = 100,
    f::Symbol = :tanh,
)::@NamedTuple{
    ic::Matrix{Float64},
    ic_mw::Matrix{Float64},
    ic_var::Vector{Float64},
}
    # validate
    nepochs(obj) == 1 ||
        throw(ArgumentError("ica_decompose() must be applied to continuous object."))
    signal_len(obj) / sr(obj) <= 10 &&
        _warn("For ICA decomposition the signal length should be >10 seconds.")

    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))
    length(ch) == 1 && (ch = ch[1])

    # perform decomposition on the selected slice
    ica_data = ica_decompose(@view(obj.data[ch, :, 1]); n = n, iter = iter, f = f)
    ic = ica_data.ic
    ic_mw = ica_data.ic_mw

    # calculate Variance Accounted For (VAF) per component
    total_var = var(@view(obj.data[ch, :, 1]))
    ic_var = Vector{Float64}(undef, n)

    for idx = 1:n
        # reconstruct signal using only the i-th component
        ic_back = @views ic_mw[:, idx] * ic[idx, :]'
        # VAF formula: 100 * (1 - var(residual) / var(original))
        ic_var[idx] = 100.0 * (1.0 - var(@view(obj.data[ch, :, 1]) .- ic_back) / total_var)
    end

    # sort components by descending variance
    p = sortperm(ic_var; rev = true)
    ic = ic[p, :]
    ic_var = ic_var[p]
    ic_mw = ic_mw[:, p]

    for i = 1:n
        _info("Component $(lpad(i, 2)): VAF = $(round(ic_var[i], digits = 2))%")
    end

    return (; ic, ic_mw, ic_var)
end

"""
    ica_remove(; <keyword arguments>)

Reconstruct a signal by removing independent component(s).

# Arguments

- `ic::Matrix{Float64}`: independent components, shape (n, samples)
- `ic_mw::Matrix{Float64}`: weighting (unmixing) matrix, shape (channels, n)
- `ic_idx::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}}`: indices of components to remove
- `keep::Bool=false`: if `true`, keep specified components and remove all other components

# Returns
- `Matrix{Float64}`: reconstructed signal, shape (channels, samples)
"""
function ica_remove(;
    ic::Matrix{Float64},
    ic_mw::Matrix{Float64},
    ic_idx::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}},
    keep::Bool = false,
)::Matrix{Float64}
    # validate
    size(ic, 1) == size(ic_mw, 2) || throw(
        ArgumentError(
            "Dimension mismatch between ic ($(size(ic)))and ic_mw ($size(ic_mw))).",
        ),
    )
    # bounds check
    all(1 .<= ic_idx .<= size(ic_mw, 2)) ||
        throw(ArgumentError("ic_idx must be in [1, $(size(ic_mw, 2))]."))

    # determine which indices to actually use for reconstruction
    target_idx = keep ? setdiff(1:size(ic_mw, 2), ic_idx) : ic_idx
    # in case a single target_idx
    target_idx = _n2v(target_idx)

    # zero target components
    ic_modified = deepcopy(ic)
    ic_modified[target_idx, :] .= 0

    # reconstruction: signal = MixingMatrix * Components
    return ic_mw * ic_modified
end

"""
    ica_remove(obj, ic, ic_mw; <keyword arguments>)

Reconstruct selected channels of a NEURO object by removing independent component(s).

# Arguments
- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `ic::Matrix{Float64}`: independent components, shape (n, samples)
- `ic_mw::Matrix{Float64}`: weighting (unmixing) matrix, shape (channels, n)
- `ic_idx::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}}`: indices of components to remove
- `keep::Bool=false`: if `true`, keep specified components and remove all other components

# Returns

- `NeuroAnalyzer.NEURO`: reconstructed NEURO object
"""
function ica_remove(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    ic::Matrix{Float64},
    ic_mw::Matrix{Float64},
    ic_idx::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}},
    keep::Bool = false,
)::NeuroAnalyzer.NEURO
    # validate
    nepochs(obj) == 1 ||
        throw(ArgumentError("ica_remove() must be applied to continuous object."))

    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))
    length(ch) == 1 && (ch = ch[1])

    # reconstruction

    # create new dataset
    obj_tmp = deepcopy(obj)

    obj_tmp.data[ch, :, 1] =
        ica_remove(; ic = ic, ic_mw = ic_mw, ic_idx = ic_idx, keep = keep)

    push!(obj_tmp.history, "ica_remove(obj; ch=$ch, ic_idx=$ic_idx, keep=$keep)")

    return obj_tmp
end

"""
    ica_remove!(obj, ic, ic_mw; <keyword arguments>)

Reconstruct selected channels of a NEURO object by removing independent component(s) in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s), default is all channels
- `ic_idx::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}}`: indices of components to keep or remove
- `ic::Matrix{Float64}`: independent components, shape (n, samples)
- `ic_mw::Matrix{Float64}`: weighting (unmixing) matrix, shape (channels, n)
- `keep::Bool=false`: if `true`, keep specified components; otherwise, remove them

# Returns

- `Nothing`
"""
function ica_remove!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    ic::Matrix{Float64},
    ic_mw::Matrix{Float64},
    ic_idx::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}},
    keep::Bool = false,
)::Nothing
    obj_tmp = ica_remove(obj; ch = ch, ic_idx = ic_idx, ic = ic, ic_mw = ic_mw, keep = keep)
    obj.data = obj_tmp.data
    obj.history = obj_tmp.history
    obj_tmp = nothing

    return nothing
end
