export ica_decompose
export ica_remove
export ica_remove!

"""
    _ica_decompose(s; <keyword arguments>)

Decompose a signal into Independent Components (ICs) using the FastICA algorithm.

# Arguments

- `s::AbstractMatrix`: input signal, shape `(channels,  samples)`
- `n::Int`: number of independent components to extract
- `iter::Int=100`: max iterations per tolerance level
- `f::Symbol=:tanh`: nonlinear function for neg-entropy approximation (`:tanh` or `:gaus`)
- `rng::Int64=42`: seed for random initialization in FastICA

# Returns

Named tuple:

- `ic::Matrix{Float64}`: independent components, shape (n, samples)
- `ic_mw::Matrix{Float64}`: mixing matrix (component weights at the channels), shape (channels, n)
"""
function _ica_decompose(
    s::AbstractMatrix;
    n::Int64,
    iter::Int64 = 100,
    f::Symbol = :tanh,
    rng::Int64 = 42,
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
    Random.seed!(rng)

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
            break # exit loop if converged
        catch err
            # if it's not a convergence error, rethrow it; otherwise, update progress and continue
            !(err isa MultivariateStats.ConvergenceException) && rethrow(err)
            # skip progress for failed tolerance bracket
            # update progress bar
            progress_bar && next!(progbar)
        end
    end
    # the loop stops at the first tolerance that converges
    progress_bar && finish!(progbar)

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

Decompose selected channels of a NEURO object into Independent Components (ICs). Sorts components by variance explained.

Methods:

- `:fastica` (default): FastICA (MultivariateStats), neg-entropy approximation with the nonlinearity `f`
- `:amica`: AMICA (Adaptive Mixture ICA): the AMICA program installed by the user (e.g. `amica15ub`, from the AMICA plugin of EEGLAB or https://github.com/sccn/amica) is run on the signal (written to a temporary file, with the default parameters of `runamica15()` of EEGLAB) and its output is read (see `import_amica()`); the components of model 1 are returned
- `:infomax`: extended Infomax (as `runica()` of EEGLAB, the default of EEGLAB and the decomposition ICLabel was trained on): natural-gradient learning on random blocks of samples after PCA reduction (if `n` < number of channels) and sphering; with `extended = true`, components with sub-Gaussian (e.g. line noise) and super-Gaussian (e.g. eye blinks) distributions are separated; `extended = false`: logistic Infomax (super-Gaussian sources only)

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `exclude_bad_channels::Bool=NeuroAnalyzer.exclude_bad_channels`: if `true`, channels marked as bad are not used (default: the `exclude_bad_channels` preference, see `set_exclude_bad_channels()`)
- `n::Union{Nothing, Int64}=nothing`: number of independent components to extract; default: the number of selected channels
- `method::Symbol=:fastica`: ICA algorithm: `:fastica`, `:infomax` or `:amica`
- `iter::Union{Nothing, Int64}=nothing`: `:fastica`: maximum number of iterations per tolerance level (default: 100); `:infomax`: maximum number of training steps (passes through the data; default: 512); `:amica`: maximum number of iterations (default: 2000)
- `f::Symbol=:tanh`: nonlinear function for neg-entropy approximation (`:tanh` or `:gaus`) (`:fastica`)
- `extended::Bool=true`: extended Infomax (sub- and super-Gaussian sources) (`:infomax`)
- `rng::Int64=42`: seed for the random initialization (`:fastica`) or the random order of the samples (`:infomax`)
- `amica::String="amica15ub"`: AMICA program: name in the `PATH` or path (`:amica`)
- `amica_models::Int64=1`: number of ICA models (`:amica`)
- `amica_threads::Int64=Sys.CPU_THREADS`: number of threads (`:amica`)
- `amica_outdir::Union{Nothing, String}=nothing`: folder for the AMICA output, kept (e.g. for `import_amica()` of other models); default: a temporary folder, deleted (`:amica`)

# Returns

Named tuple:

- `ic::Matrix{Float64}`: independent components, shape (n, samples)
- `ic_mw::Matrix{Float64}`: mixing matrix (component weights at the channels), shape (channels, n)
- `ic_var::Vector{Float64}`: variance explained by each component

# References

1. Hyvärinen A, Oja E. Independent component analysis: algorithms and applications. Neural Netw. 2000;13(4-5):411-430.
2. Lee TW, Girolami M, Sejnowski TJ. Independent component analysis using an extended infomax algorithm for mixed subgaussian and supergaussian sources. Neural Comput. 1999;11(2):417-441.
3. Palmer JA, Kreutz-Delgado K, Makeig S. AMICA: An Adaptive Mixture of Independent Component Analyzers with Shared Components. Technical report, Swartz Center for Computational Neuroscience, UCSD; 2011.
"""
function ica_decompose(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    n::Union{Nothing, Int64} = nothing,
    method::Symbol = :fastica,
    iter::Union{Nothing, Int64} = nothing,
    f::Symbol = :tanh,
    extended::Bool = true,
    rng::Int64 = 42,
    amica::String = "amica15ub",
    amica_models::Int64 = 1,
    amica_threads::Int64 = Sys.CPU_THREADS,
    amica_outdir::Union{Nothing, String} = nothing,
    exclude_bad_channels::Bool = NeuroAnalyzer.exclude_bad_channels,
)::@NamedTuple{
    ic::Matrix{Float64},
    ic_mw::Matrix{Float64},
    ic_var::Vector{Float64},
}
    # validate
    _check_var(method, [:fastica, :infomax, :amica], "method")
    nepochs(obj) == 1 ||
        throw(ArgumentError("ica_decompose() must be applied to continuous object."))
    signal_len(obj) / sr(obj) <= 10 &&
        _warn("For ICA decomposition the signal length should be >10 seconds.")

    # resolve channel names to integer indices
    ch = _get_ch(obj, ch, exclude_bad_channels)
    length(ch) >= 2 || throw(ArgumentError("At least 2 channels must be selected."))
    isnothing(n) && (n = length(ch))

    # Infomax and AMICA need full-rank data (e.g. average reference: rank = number of channels - 1)
    if method in (:infomax, :amica)
        r = _data_rank(@view(obj.data[ch, :, 1]))
        if n > r
            _warn("Data rank ($r) < number of components ($n): n reduced to $r.")
            n = r
        end
    end

    # perform decomposition on the selected slice
    ica_data = if method === :fastica
        _ica_decompose(
            @view(obj.data[ch, :, 1]);
            n = n,
            iter = something(iter, 100),
            f = f,
            rng = rng,
        )
    elseif method === :infomax
        _ica_infomax(
            @view(obj.data[ch, :, 1]);
            n = n,
            iter = something(iter, 512),
            extended = extended,
            rng = rng,
        )
    else
        _ica_amica(
            @view(obj.data[ch, :, 1]);
            n = n,
            iter = something(iter, 2000),
            amica = amica,
            models = amica_models,
            threads = amica_threads,
            outdir = amica_outdir,
        )
    end
    ic = ica_data.ic
    ic_mw = ica_data.ic_mw

    return _ica_vaf_sort(obj.data[ch, :, 1], ic, ic_mw)
end

# variance accounted for (VAF, %) by each component of a decomposition of the signal `s`
# (channels × samples): 100 × (1 - var(s - component back-projection) / var(s)) of the
# demeaned signal; the components are sorted by decreasing VAF
function _ica_vaf_sort(
    s::AbstractMatrix,
    ic::AbstractMatrix,
    ic_mw::AbstractMatrix,
)::@NamedTuple{
    ic::Matrix{Float64},
    ic_mw::Matrix{Float64},
    ic_var::Vector{Float64},
}
    d = Matrix{Float64}(s)
    d .-= mean(d; dims = 2)
    total_var = var(d)
    n = size(ic, 1)
    ic_var = Vector{Float64}(undef, n)
    for idx = 1:n
        # reconstruct signal using only the i-th component
        ic_back = @views ic_mw[:, idx] * ic[idx, :]'
        ic_var[idx] = 100.0 * (1.0 - var(d .- ic_back) / total_var)
    end

    # sort components by descending variance
    p = sortperm(ic_var; rev = true)
    for i = 1:n
        _info("Component $(lpad(i, 2)): VAF = $(round(ic_var[p[i]], digits = 2))%")
    end

    return (; ic = Matrix{Float64}(ic[p, :]), ic_mw = Matrix{Float64}(ic_mw[:, p]), ic_var = ic_var[p])
end

"""
    ica_remove(; <keyword arguments>)

Reconstruct a signal from independent components, removing the selected component(s). The signal is the mixture of the remaining components only: channel means and the part of the signal outside the ICA subspace (if fewer components than channels) are not included; `ica_remove(obj; ...)` keeps them.

# Arguments

- `ic::Matrix{Float64}`: independent components, shape (n, samples)
- `ic_mw::Matrix{Float64}`: mixing matrix (component weights at the channels), shape (channels, n)
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
            "Dimension mismatch between ic $(size(ic)) and ic_mw $(size(ic_mw)).",
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

Remove independent component(s) from selected channels of a NEURO object: the back-projection of the removed components (`ic_mw[:, i] * ic[i, :]`) is subtracted from the signal, so channel means and the signal outside the ICA subspace are kept.

# Arguments
- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `exclude_bad_channels::Bool=NeuroAnalyzer.exclude_bad_channels`: if `true`, channels marked as bad are not used (default: the `exclude_bad_channels` preference, see `set_exclude_bad_channels()`)
- `ic::Matrix{Float64}`: independent components, shape (n, samples)
- `ic_mw::Matrix{Float64}`: mixing matrix (component weights at the channels), shape (channels, n)
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
    exclude_bad_channels::Bool = NeuroAnalyzer.exclude_bad_channels,
)::NeuroAnalyzer.NEURO
    obj_new = deepcopy(obj)
    ica_remove!(obj_new; ch = ch, exclude_bad_channels = exclude_bad_channels, ic_idx = ic_idx, ic = ic, ic_mw = ic_mw, keep = keep)

    return obj_new
end

"""
    ica_remove!(obj, ic, ic_mw; <keyword arguments>)

Reconstruct selected channels of a NEURO object by removing independent component(s) in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s), default is all channels
- `exclude_bad_channels::Bool=NeuroAnalyzer.exclude_bad_channels`: if `true`, channels marked as bad are not used (default: the `exclude_bad_channels` preference, see `set_exclude_bad_channels()`)
- `ic_idx::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}}`: indices of components to keep or remove
- `ic::Matrix{Float64}`: independent components, shape (n, samples)
- `ic_mw::Matrix{Float64}`: mixing matrix (component weights at the channels), shape (channels, n)
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
    exclude_bad_channels::Bool = NeuroAnalyzer.exclude_bad_channels,
)::Nothing
    # validate
    nepochs(obj) == 1 ||
        throw(ArgumentError("ica_remove() must be applied to continuous object."))

    # resolve channel names to integer indices
    ch = _get_ch(obj, ch, exclude_bad_channels)
    isempty(ch) && throw(ArgumentError("No channels selected."))
    size(ic_mw, 1) == length(ch) ||
        throw(
            ArgumentError(
                "ic_mw has $(size(ic_mw, 1)) rows, $(length(ch)) channels are selected.",
            ),
        )
    size(ic, 2) == epoch_len(obj) ||
        throw(ArgumentError("ic has $(size(ic, 2)) samples, OBJ has $(epoch_len(obj))."))
    size(ic, 1) == size(ic_mw, 2) ||
        throw(
            ArgumentError(
                "Dimension mismatch between ic $(size(ic)) and ic_mw $(size(ic_mw)).",
            ),
        )
    all(1 .<= ic_idx .<= size(ic_mw, 2)) ||
        throw(ArgumentError("ic_idx must be in [1, $(size(ic_mw, 2))]."))

    # components to remove
    rm_idx = keep ? setdiff(1:size(ic_mw, 2), ic_idx) : unique(collect(ic_idx))

    # subtract the back-projection of the removed components
    obj.data[ch, :, 1] = obj.data[ch, :, 1] .- ic_mw[:, rm_idx] * ic[rm_idx, :]

    _history!(obj, "ica_remove", (; ch, ic_idx, keep, ic, ic_mw, exclude_bad_channels))

    return nothing
end
