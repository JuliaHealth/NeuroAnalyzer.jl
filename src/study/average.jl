export study_average
export study_erp

"""
    study_average(study, kind; <keyword arguments>)

Grand average of a study measure (see `study_measure!()`) over subjects, for each level of the design factor(s) `by` (e.g. group) and each condition. The recordings of one subject (e.g. sessions, if `:session` is not in `by`) are averaged first, so each subject counts once. Bad channels (`NaN`) are left out, so the number of subjects may differ between channels.

# Arguments

- `study::NeuroAnalyzer.STUDY`
- `kind::Symbol`: measure (`:erp`, `:psd` or `:tf`)
- `by::Union{Symbol, Vector{Symbol}}=:group`: design column(s) defining the levels, e.g. `:group`, `[:group, :session]`; `Symbol[]`: all subjects together
- `condition::Union{Nothing, String, Vector{String}}=nothing`: condition(s); `nothing`: all conditions of the measure
- `db::Bool=true`: for `:psd` and `:tf`, convert the power of each subject to dB before averaging

# Returns

Named tuple:

- `kind::Symbol`
- `levels::DataFrame`: one row per level: the `by` columns, `condition` and `n` (number of subjects)
- `data::Array{Float64}`: mean over subjects, shape (level, channel, time or frequency) or (level, channel, frequency, time) for `:tf`
- `sd::Array{Float64}`: standard deviation over subjects (same shape)
- `n::Matrix{Int64}`: number of subjects, shape (level, channel)
- `labels::Vector{String}`: channel labels
- `unit::Vector{String}`: channel units
- `t::Vector{Float64}`: time points in seconds (`:erp`, `:tf`)
- `f::Vector{Float64}`: frequencies in Hz (`:psd`, `:tf`)
- `db::Bool`: powers in dB

# Notes

The averages can be plotted with the plots of arrays, e.g. `plot_psd(f, data[1, :, :]; clabels = labels)`; for ERPs `study_erp()` returns a NEURO object for `plot_erp()` and `plot_topo()`.
"""
function study_average(
    study::NeuroAnalyzer.STUDY,
    kind::Symbol;
    by::Union{Symbol, Vector{Symbol}} = :group,
    condition::Union{Nothing, String, Vector{String}} = nothing,
    db::Bool = true,
)::NamedTuple
    _check_var(kind, [:erp, :psd, :tf], "kind")
    m = _study_measure(study, kind)
    by = by isa Symbol ? [by] : by
    for b in by
        b in propertynames(study.design) ||
            throw(ArgumentError("Design table has no column $b."))
    end
    d, idx = _study_measure_rows(study, m)
    ci = _study_conditions(m, condition)
    X = selectdim(m.data, 1, idx)
    db = db && kind !== :erp
    if db
        X = map(x -> x > 0 ? 10 * log10(x) : NaN, X)
    end

    # levels: combinations of the `by` columns (in order of appearance) and conditions
    lv = isempty(by) ? DataFrame() : unique(d[!, by])
    nlv = max(1, DataFrames.nrow(lv))
    sz = size(m.data)[3:end]
    data = fill(NaN, nlv * length(ci), sz...)
    sd = fill(NaN, nlv * length(ci), sz...)
    n = zeros(Int64, nlv * length(ci), sz[1])
    levels = DataFrame()
    row = 0
    for l = 1:nlv
        sel = isempty(by) ? trues(DataFrames.nrow(d)) :
              [all(d[r, b] == lv[l, b] for b in by) for r = 1:DataFrames.nrow(d)]
        for c in ci
            row += 1
            S = _study_subject_avg(selectdim(X, 2, c), d, sel, :session in by)
            mu, s, nn = _nanstats(S)
            selectdim(data, 1, row) .= mu
            selectdim(sd, 1, row) .= s
            n[row, :] = nn[:, ntuple(_ -> 1, ndims(nn) - 1)...]
            r = isempty(by) ? DataFrame() : DataFrame(lv[l:l, :])
            r[!, :condition] = [m.conditions[c]]
            r[!, :n] = [maximum(n[row, :])]
            levels = vcat(levels, r)
        end
    end

    return (
        kind = kind,
        levels = levels,
        data = data,
        sd = sd,
        n = n,
        labels = m.labels,
        unit = m.unit,
        t = m.t,
        f = m.f,
        db = db,
    )
end

"""
    study_erp(study; <keyword arguments>)

Grand average ERP (ERF for MEG) of a study as a NEURO object, e.g. for `plot_erp()` and `plot_topo()`: epoch 1 is the grand average over subjects, the following epochs are the averages of the subjects (as `average_epochs()`). Requires the `:erp` measure (see `study_measure!()`).

The header and channel locations are those of the first recording. Channels marked as bad in a subject are left out of the grand average; in the epoch of that subject they are replaced by the grand average.

# Arguments

- `study::NeuroAnalyzer.STUDY`
- `condition::Union{Nothing, String}=nothing`: condition; may be omitted if the measure has one condition
- `group::Union{Nothing, String}=nothing`: group of subjects; `nothing`: all subjects

# Returns

- `NeuroAnalyzer.NEURO`
"""
function study_erp(
    study::NeuroAnalyzer.STUDY;
    condition::Union{Nothing, String} = nothing,
    group::Union{Nothing, String} = nothing,
)::NeuroAnalyzer.NEURO
    m = _study_measure(study, :erp)
    if isnothing(condition)
        length(m.conditions) == 1 || throw(
            ArgumentError(
                "condition must be given (conditions: $(join(repr.(m.conditions), ", "))).",
            ),
        )
        condition = m.conditions[1]
    end
    c = only(_study_conditions(m, condition))
    d, idx = _study_measure_rows(study, m)
    sel = isnothing(group) ? trues(DataFrames.nrow(d)) : d.group .== group
    any(sel) || throw(ArgumentError("No recordings of group $(repr(group))."))
    S = _study_subject_avg(selectdim(selectdim(m.data, 1, idx), 2, c), d, sel, false)
    subj = unique(d.subject[sel])
    # subjects without epochs of the condition
    has = [!all(isnan, selectdim(S, 1, k)) for k in axes(S, 1)]
    any(has) || throw(ArgumentError("No epochs of condition $(repr(condition))."))
    S = S[has, :, :]
    subj = subj[has]
    ga, _, _ = _nanstats(S)
    any(isnan, ga) && _warn("Channels bad in all subjects are NaN.")
    # bad channels of a subject: the grand average
    for k in axes(S, 1)
        v = @view S[k, :, :]
        v[isnan.(v)] .= ga[isnan.(v)]
    end

    # header of the first recording of the selected subjects
    obj = _silent(() -> load_neuro(_study_path(study, d.file[findfirst(sel)])))
    _silent(() -> keep_channel!(obj; ch = m.labels))
    labels(obj) == m.labels || _silent(() -> reorder_channel!(obj; order = m.labels))
    obj.data = cat(reshape(ga, size(ga)..., 1), permutedims(S, (2, 3, 1)); dims = 3)
    obj.time_pts, _ = _get_t(obj)
    obj.epoch_time = copy(m.t)
    obj.markers = obj.markers[1:0, :]
    obj.epochs = DataFrame(
        :condition => fill(condition, size(obj.data, 3)),
        :source    => ["grand average"; subj],
        :reject    => falses(size(obj.data, 3)),
        :onset     => fill(NaN, size(obj.data, 3)),
    )
    dt = m.data_type
    obj.header.recording[:data_type] = dt == "eeg" ? "erp" : dt == "meg" ? "erf" : dt
    obj.header.recording[:bad_channel] = falses(nchannels(obj))
    obj.header.recording[:file_name] = ""
    obj.history = [
        "study_erp(study; condition = $(repr(condition)), group = $(repr(group))) # $(length(subj)) subjects" *
        (isempty(study.name) ? "" : ", study $(study.name)"),
    ]

    return obj
end

# design rows of the included recordings and their rows in the measure
function _study_measure_rows(study::NeuroAnalyzer.STUDY, m)
    d = study.design[study.design.include, :]
    idx = indexin(d.file, m.files)
    any(isnothing, idx) && throw(
        ArgumentError(
            "Recordings not in measure :$(m.kind): $(join(d.file[isnothing.(idx)], ", ")); run study_measure!() to update it.",
        ),
    )
    return d, Int64.(idx)
end

# condition numbers in the measure
function _study_conditions(m, condition)::Vector{Int64}
    isnothing(condition) && return collect(eachindex(m.conditions))
    cs = condition isa String ? [condition] : condition
    ci = indexin(cs, m.conditions)
    any(isnothing, ci) && throw(
        ArgumentError(
            "Conditions not in measure :$(m.kind): $(join(repr.(cs[isnothing.(ci)]), ", ")) (conditions: $(join(repr.(m.conditions), ", "))).",
        ),
    )
    return Int64.(ci)
end

# averages of the subjects of the selected recordings (the recordings of one subject are
# averaged, unless `sessions`), shape (subject, values...)
function _study_subject_avg(
    X::AbstractArray, d::DataFrame, sel::AbstractVector{Bool}, sessions::Bool,
)::Array{Float64}
    rows = findall(sel)
    keys_ = sessions ? collect(eachindex(rows)) : indexin(d.subject[rows], unique(d.subject[rows]))
    nk = maximum(keys_)
    S = fill(NaN, nk, size(X)[2:end]...)
    for k = 1:nk
        rr = rows[keys_ .== k]
        mu, _, _ = _nanstats(selectdim(X, 1, rr))
        selectdim(S, 1, k) .= mu
    end
    return S
end

# mean, standard deviation and number of non-NaN values along dimension 1
function _nanstats(X::AbstractArray)
    n1 = size(X, 1)
    rest = size(X)[2:end]
    M = reshape(collect(Float64, X), n1, :)
    mu = fill(NaN, size(M, 2))
    sd = fill(NaN, size(M, 2))
    n = zeros(Int64, size(M, 2))
    @inbounds for j in axes(M, 2)
        v = Base.filter(!isnan, @view M[:, j])
        n[j] = length(v)
        n[j] > 0 && (mu[j] = mean(v))
        n[j] > 1 && (sd[j] = std(v))
    end
    return reshape(mu, rest), reshape(sd, rest), reshape(n, rest)
end
