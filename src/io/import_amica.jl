export import_amica

"""
    import_amica(obj; <keyword arguments>)

Import an AMICA (Adaptive Mixture ICA) decomposition: the output folder of the AMICA program (`outdir`, e.g. `amicaout/`; AMICA is run separately, e.g. from EEGLAB or as a standalone program) is read and the independent components of the selected channels of `obj` are calculated. The result has the same form as that of `ica_decompose()` and is used in the same way (`ica_remove()`, `plot_ica()`).

The files `W` (unmixing weights of each model), `S` (sphering matrix), `mean` (channel means) and `gm` (model weights) are read (raw little-endian Float64, as written by AMICA and read by `loadmodout15()` of EEGLAB). The unmixing matrix of model `h` is `W[:, :, h] × S[1:nw, :]` (`nw`: number of components, < number of channels if AMICA reduced the data with PCA), applied to the signal minus the channel means.

The signal (selected channels, in the same order) must be the signal AMICA was run on (the same channels and preprocessing).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; must be continuous (1 epoch)
- `ch::Union{String, Vector{String}, Regex}`: channel name(s), as in the AMICA decomposition
- `path::String`: AMICA output folder
- `model::Int64=1`: model of a multi-model decomposition
- `exclude_bad_channels::Bool=NeuroAnalyzer.exclude_bad_channels`: if `true`, channels marked as bad are not used (default: the `exclude_bad_channels` preference, see `set_exclude_bad_channels()`)

# Returns

Named tuple:

- `ic::Matrix{Float64}`: independent components, shape (components, samples)
- `ic_mw::Matrix{Float64}`: mixing matrix (component weights at the channels), shape (channels, components)
- `ic_var::Vector{Float64}`: variance explained by each component (%); the components are sorted by it

# References

1. Palmer JA, Kreutz-Delgado K, Makeig S. AMICA: An Adaptive Mixture of Independent Component Analyzers with Shared Components. Technical report, Swartz Center for Computational Neuroscience, UCSD; 2011.
2. Delorme A, Palmer J, Onton J, Oostenveld R, Makeig S. Independent EEG sources are dipolar. PLoS One. 2012;7(2):e30135.
"""
function import_amica(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    path::String,
    model::Int64 = 1,
    exclude_bad_channels::Bool = NeuroAnalyzer.exclude_bad_channels,
)::@NamedTuple{
    ic::Matrix{Float64},
    ic_mw::Matrix{Float64},
    ic_var::Vector{Float64},
}
    # validate
    nepochs(obj) == 1 ||
        throw(ArgumentError("import_amica() must be applied to a continuous object."))
    um = _amica_unmix(path; model = model)
    nx = um.nx

    # channels
    ch = _get_ch(obj, ch, exclude_bad_channels)
    length(ch) == nx || throw(
        ArgumentError(
            "The AMICA decomposition has $nx channels, $(length(ch)) channels are selected.",
        ),
    )

    _info("AMICA: model $model of $(um.n_models) (weight $(round(um.gm; digits = 3))), $(um.nw) components, $nx channels")

    # components of the signal minus the channel means
    x = obj.data[ch, :, 1] .- um.mean
    ic = um.unmix * x
    ic_mw = um.nw == nx ? inv(um.unmix) : LinearAlgebra.pinv(um.unmix)

    return _ica_vaf_sort(obj.data[ch, :, 1], ic, ic_mw)
end

# unmixing matrix of model `model` of an AMICA output folder: `W[:, :, model] × S[1:nw, :]`
# (components × channels), applied to the signal minus `mean`
function _amica_unmix(
    path::String;
    model::Int64 = 1,
)::@NamedTuple{
    unmix::Matrix{Float64},
    mean::Vector{Float64},
    nx::Int64,
    nw::Int64,
    n_models::Int64,
    gm::Float64,
}
    isdir(path) || throw(ArgumentError("AMICA output folder $path not found."))
    for f in ("W", "S", "mean", "gm")
        isfile(joinpath(path, f)) ||
            throw(ArgumentError("File $f not found in $path (not an AMICA output folder?)."))
    end

    # number of models
    gm = _amica_read(joinpath(path, "gm"))
    n_models = length(gm)
    n_models >= 1 || throw(ArgumentError("$(joinpath(path, "gm")) is empty."))
    1 <= model <= n_models ||
        throw(ArgumentError("model must be in [1, $n_models] (number of models)."))

    # sphering matrix: nx × nx (nx: number of channels)
    s_flat = _amica_read(joinpath(path, "S"))
    nx = isqrt(length(s_flat))
    nx^2 == length(s_flat) ||
        throw(ArgumentError("S is not a square matrix ($(length(s_flat)) values)."))
    sph = reshape(s_flat, nx, nx)

    # unmixing weights: nw × nw × models (nw: number of components)
    w_flat = _amica_read(joinpath(path, "W"))
    length(w_flat) % n_models == 0 ||
        throw(ArgumentError("W does not match the number of models ($n_models)."))
    nw = isqrt(length(w_flat) ÷ n_models)
    nw^2 * n_models == length(w_flat) ||
        throw(ArgumentError("W is not a set of square matrices ($(length(w_flat)) values)."))
    nw <= nx || throw(ArgumentError("Number of components ($nw) > number of channels ($nx)."))
    w = reshape(w_flat, nw, nw, n_models)

    mn = _amica_read(joinpath(path, "mean"))
    length(mn) == nx ||
        throw(ArgumentError("mean has $(length(mn)) values, $nx channels expected."))

    return (; unmix = w[:, :, model] * sph[1:nw, :], mean = mn, nx, nw, n_models, gm = gm[model])
end

# values of an AMICA output file (raw little-endian Float64)
function _amica_read(file::String)::Vector{Float64}
    b = read(file)
    length(b) % 8 == 0 ||
        throw(ArgumentError("$file is not a Float64 file ($(length(b)) bytes)."))
    return ltoh.(collect(reinterpret(Float64, b)))
end
