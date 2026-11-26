export study_measure!

"""
    study_measure!(study, kind; <keyword arguments>)

Compute a measure for every included recording and condition and cache it in the study (`study.measures[kind]`, saved with the study by `save_study()`). A cached measure is recomputed only when its parameters, the included recordings or their files (e.g. after `study_apply!()`) change.

The recordings must be compatible: the same data type, sampling rate and channel labels (in the same order) and, for `:erp` and `:tf`, the same epoch time axis (see `study_check()`). Channels marked as bad in a recording are `NaN` in its measure, so that they are left out of the grand averages. Rejected epochs (see `epoch_reject!()`) are not used.

# Arguments

- `study::NeuroAnalyzer.STUDY`: modified in-place
- `kind::Symbol`: measure:
    - `:erp`: average of the epochs (ERP/ERF), shape (recording, condition, channel, time)
    - `:psd`: power spectral density (`psd()`) averaged over epochs, shape (recording, condition, channel, frequency)
    - `:tf`: time-frequency power (`spectrogram()`) averaged over epochs, shape (recording, condition, channel, frequency, time)
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `condition::Union{Nothing, String, Vector{String}}=nothing`: epoch condition(s) (epoch metadata, see `get_epochs()`); `nothing`: the conditions of all recordings (all epochs as condition `""` for recordings without conditions)
- `bl::Tuple{Real, Real}=(0, 0)`: baseline window in seconds (`:erp`), its mean is subtracted from the average
- `force::Bool=false`: if `true`, recompute a cached measure
- `kwargs...`: arguments of `psd()` (`:psd`) or `spectrogram()` (`:tf`), e.g. `method`, `flim`, `wlen`; powers are always computed in linear units (dB: see `study_average()`)

# Returns

- `Nothing`

# Notes

The measure is a named tuple with the fields:

- `kind::Symbol`
- `data::Array{Float64}`: values (see `kind`); `NaN` for bad channels and for conditions missing in a recording
- `labels::Vector{String}`: channel labels
- `unit::Vector{String}`: channel units
- `data_type::String`: data type of the recordings
- `conditions::Vector{String}`: conditions
- `n_epochs::Matrix{Int64}`: number of epochs, shape (recording, condition)
- `files::Vector{String}`: files of the recordings (as in `study.design.file`)
- `mtime::Vector{Float64}`: modification times of the files
- `t::Vector{Float64}`: time points in seconds (`:erp`, `:tf`)
- `f::Vector{Float64}`: frequencies in Hz (`:psd`, `:tf`)
- `params::NamedTuple`: parameters
"""
function study_measure!(
    study::NeuroAnalyzer.STUDY,
    kind::Symbol;
    ch::Union{String, Vector{String}, Regex},
    condition::Union{Nothing, String, Vector{String}} = nothing,
    bl::Tuple{Real, Real} = (0, 0),
    force::Bool = false,
    kwargs...,
)::Nothing
    # validate
    _check_var(kind, [:erp, :psd, :tf], "kind")
    haskey(kwargs, :db) && throw(
        ArgumentError(
            "Powers are computed in linear units; dB conversion: study_average(; db = true).",
        ),
    )
    kind === :erp && !isempty(kwargs) && throw(
        ArgumentError("Unknown arguments for :erp: $(join(string.(keys(kwargs)), ", "))."),
    )
    bl != (0, 0) && kind !== :erp && throw(ArgumentError("bl is used only by :erp."))
    rows = findall(study.design.include)
    isempty(rows) && throw(ArgumentError("No recordings included."))

    files = study.design.file[rows]
    paths = [_study_path(study, f) for f in files]
    missing_files = Base.filter(!isfile, paths)
    isempty(missing_files) ||
        throw(ArgumentError("Files not found: $(join(missing_files, ", "))."))
    mt = mtime.(paths)
    params = kind === :erp ? (; ch, condition, bl) : (; ch, condition, kwargs...)

    # cached measure
    if !force && haskey(study.measures, kind)
        m = study.measures[kind]
        if isequal(m.params, params) && m.files == files && m.mtime == mt
            _info("Measure :$kind is up to date.")
            return nothing
        end
    end

    ref = nothing
    ref_file = ""
    lbl = String[]
    unit = String[]
    t = Float64[]
    f = Float64[]
    # per recording: condition => values (channel × axis...) and number of epochs
    res = Vector{Dict{String, Tuple{Array{Float64}, Int64}}}(undef, length(rows))
    conds = String[]

    p = Progress(length(rows); dt = 1, barlen = 20, color = :white, enabled = progress_bar)
    for (k, path) in enumerate(paths)
        obj = _silent(() -> load_neuro(path))
        chi = _get_ch(obj, ch, false)
        isempty(chi) && throw(ArgumentError("No channels selected in $(files[k])."))
        sig = _study_sig(obj, chi)
        if isnothing(ref)
            ref = sig
            ref_file = files[k]
            lbl = sig.labels
            unit = obj.header.recording[:unit][chi]
            kind === :erp && (t = copy(obj.epoch_time))
        else
            iss = _study_compat(ref, sig; time = kind !== :psd)
            isempty(iss) || throw(
                ArgumentError(
                    "Recording $(files[k]) is not compatible with $ref_file: $iss (see study_check()).",
                ),
            )
        end
        bad = obj.header.recording[:bad_channel][chi]
        cs = if isnothing(condition)
            c = epoch_conditions(obj)
            isempty(c) ? [""] : c
        else
            condition isa String ? [condition] : condition
        end

        res[k] = Dict{String, Tuple{Array{Float64}, Int64}}()
        for c in cs
            ep = get_epochs(obj; condition = isempty(c) ? nothing : c)
            isempty(ep) && continue
            o = length(ep) == nepochs(obj) ? obj : keep_epoch(obj; ep = ep)
            x = if kind === :erp
                x = dropdims(mean(@view(obj.data[chi, :, ep]); dims = 3); dims = 3)
                if bl != (0, 0)
                    _check_tuple(bl, (obj.epoch_time[1], obj.epoch_time[end]), "bl")
                    b1, b2 = vsearch(bl[1], obj.epoch_time), vsearch(bl[2], obj.epoch_time)
                    x .-= mean(x[:, b1:b2]; dims = 2)
                end
                x
            elseif kind === :psd
                r = _silent(
                    () -> psd(o; ch = labels(obj)[chi], db = false,
                              exclude_bad_channels = false, kwargs...),
                )
                isempty(f) && (f = r.f)
                r.f ≈ f || throw(
                    ArgumentError(
                        "Frequencies of $(files[k]) differ from those of $ref_file.",
                    ),
                )
                dropdims(mean(r.p; dims = 3); dims = 3)
            else
                r = _silent(
                    () -> spectrogram(o; ch = labels(obj)[chi], db = false,
                                      exclude_bad_channels = false, kwargs...),
                )
                if isempty(f)
                    f = r.f
                    t = r.t
                end
                (r.f ≈ f && length(r.t) == length(t) && r.t ≈ t) || throw(
                    ArgumentError(
                        "Frequencies or time points of $(files[k]) differ from those of $ref_file.",
                    ),
                )
                # (freq, time, channel) → (channel, freq, time)
                permutedims(dropdims(mean(r.p; dims = 4); dims = 4), (3, 1, 2))
            end
            # bad channels are left out of the averages
            selectdim(x, 1, bad) .= NaN
            res[k][c] = (x, length(ep))
            c in conds || push!(conds, c)
        end
        next!(p)
    end
    finish!(p)

    # conditions: as requested, otherwise in the order of their first appearance
    isnothing(condition) ||
        (conds = [c for c in (condition isa String ? [condition] : condition) if c in conds])
    isempty(conds) && throw(ArgumentError("No epochs of the condition(s) $(repr(condition))."))

    sz = size(first(v[1] for r in res for v in values(r)))
    data = fill(NaN, length(rows), length(conds), sz...)
    n_ep = zeros(Int64, length(rows), length(conds))
    for (k, r) in enumerate(res), (ci, c) in enumerate(conds)
        haskey(r, c) || continue
        selectdim(selectdim(data, 1, k), 1, ci) .= r[c][1]
        n_ep[k, ci] = r[c][2]
    end
    nm = count(==(0), n_ep)
    nm > 0 && _warn("$nm recording × condition combination(s) without epochs (NaN).")

    study.measures[kind] = (
        kind = kind,
        data = data,
        labels = lbl,
        unit = unit,
        data_type = ref.data_type,
        conditions = conds,
        n_epochs = n_ep,
        files = files,
        mtime = mt,
        t = t,
        f = f,
        params = params,
    )
    kw = join(["$k = $(repr(v))" for (k, v) in pairs(params)], ", ")
    push!(study.history, "study_measure!(study, :$kind; $kw)")

    return nothing
end

# measure `kind` of the study; an error if not computed, a warning if out of date
function _study_measure(study::NeuroAnalyzer.STUDY, kind::Symbol)
    haskey(study.measures, kind) ||
        throw(ArgumentError("Measure :$kind is not computed (see study_measure!())."))
    m = study.measures[kind]
    paths = [_study_path(study, f) for f in m.files]
    if !all(isfile, paths) || mtime.(paths) != m.mtime
        _warn("Files of measure :$kind have changed; run study_measure!() to update it.")
    end
    return m
end
