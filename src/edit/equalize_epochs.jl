export equalize_epochs
export equalize_epochs!

"""
    equalize_epochs(obj; <keyword arguments>)

Return a copy of `obj` with the same number of epochs in each condition (epoch metadata, `obj.epochs`, see `epoch_conditions()`): epochs of the conditions with more epochs are deleted until each condition has as many epochs as the smallest one.

Only epochs that are not marked as rejected (see `epoch_reject!()`) are counted and selected; rejected epochs, epochs of other conditions and epochs without a condition are kept.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `condition::Union{Nothing, Vector{String}}=nothing`: conditions to equalize (≥ 2); `nothing`: all conditions
- `method::Symbol=:random`: selection of the epochs to keep:
    - `:random`: a random subset of the epochs of each condition
    - `:mintime`: the epochs nearest in time (onset in the source recording, `obj.epochs.onset`) to the epochs of the condition with the fewest epochs: an order-preserving matching that minimizes the sum of the onset differences; only epochs of the same source recording are matched
- `rng::AbstractRNG=Random.default_rng()`: random number generator for `method = :random` (e.g. `Random.Xoshiro(1)` for reproducible results)

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object

# Examples

```julia
e = epoch(eeg; marker_value = ["s70", "s71"], offset = 0.5, ep_len = 2.5)
e = equalize_epochs(e; method = :mintime)
```
"""
function equalize_epochs(
    obj::NeuroAnalyzer.NEURO;
    condition::Union{Nothing, Vector{String}} = nothing,
    method::Symbol = :random,
    rng::AbstractRNG = Random.default_rng(),
)::NeuroAnalyzer.NEURO
    obj_new = deepcopy(obj)
    equalize_epochs!(obj_new; condition = condition, method = method, rng = rng)

    return obj_new
end

"""
    equalize_epochs!(obj; <keyword arguments>)

Equalize the number of epochs across conditions in-place (see `equalize_epochs()`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `condition::Union{Nothing, Vector{String}}=nothing`: conditions to equalize (≥ 2); `nothing`: all conditions
- `method::Symbol=:random`: selection of the epochs to keep:
    - `:random`: a random subset of the epochs of each condition
    - `:mintime`: the epochs nearest in time to the epochs of the condition with the fewest epochs
- `rng::AbstractRNG=Random.default_rng()`: random number generator for `method = :random`

# Returns

- `Nothing`
"""
function equalize_epochs!(
    obj::NeuroAnalyzer.NEURO;
    condition::Union{Nothing, Vector{String}} = nothing,
    method::Symbol = :random,
    rng::AbstractRNG = Random.default_rng(),
)::Nothing
    # validate
    _check_var(method, [:random, :mintime], "method")
    datatype(obj) in ("erp", "erf") && throw(
        ArgumentError("equalize_epochs() requires epoched data, not $(datatype(obj)) (epoch 1 is the average)."),
    )
    cond = isnothing(condition) ? epoch_conditions(obj) : unique(condition)
    length(cond) >= 2 || throw(ArgumentError("At least 2 conditions are required (found: $(length(cond)))."))
    for c in cond
        c in epoch_conditions(obj) || throw(ArgumentError("Condition $c not found."))
    end

    # epochs of each condition (not rejected)
    ep_c = [get_epochs(obj; condition = c) for c in cond]
    n = minimum(length, ep_c)
    n > 0 || throw(
        ArgumentError("Condition $(cond[argmin(length.(ep_c))]) has no epochs that are not rejected."),
    )
    all(==(n), length.(ep_c)) && return nothing # already equal

    ep_tbl = _epochs_table(obj)
    # reference for :mintime: the (first) condition with the fewest epochs
    ref = ep_c[argmin(length.(ep_c))]

    ep_del = Int64[]
    for ep in ep_c
        length(ep) == n && continue
        keep = if method === :random
            sort(ep[randperm(rng, length(ep))[1:n]])
        else
            _mintime_match(ep_tbl, ref, ep)
        end
        append!(ep_del, setdiff(ep, keep))
    end

    n0 = length(obj.history)
    delete_epoch!(obj; ep = sort(ep_del))
    _history!(obj, "equalize_epochs", (; condition = cond, method); from = n0)

    _info("$(length(ep_del)) epoch(s) deleted, $n epoch(s) per condition")

    return nothing
end

# order-preserving matching of the reference epochs `ref` with the epochs `ep`
# (length(ep) ≥ length(ref)) that minimizes the sum of the absolute onset differences
# (dynamic programming); epochs of different source recordings are not matched (unless
# no other choice is possible); returns the matched epochs of `ep`
function _mintime_match(
    ep_tbl::DataFrame,
    ref::AbstractVector{<:Integer},
    ep::AbstractVector{<:Integer},
)::Vector{Int64}
    # both sequences in time order (by source recording, then onset)
    ord(x) = sort(collect(x); by = i -> (ep_tbl.source[i], ep_tbl.onset[i]))
    r = ord(ref)
    c = ord(ep)
    n = length(r)
    m = length(c)
    # epochs of different recordings: penalty larger than any onset difference
    pen = 1.0 + 2 * (maximum(ep_tbl.onset) - minimum(ep_tbl.onset))
    d(i, j) = ep_tbl.source[r[i]] == ep_tbl.source[c[j]] ?
              abs(ep_tbl.onset[r[i]] - ep_tbl.onset[c[j]]) : pen

    # cost[i + 1, j + 1]: first i references matched within the first j candidates
    cost = fill(Inf, n + 1, m + 1)
    cost[1, :] .= 0.0
    for i in 1:n, j in i:m
        skip = cost[i + 1, j]
        take = cost[i, j] + d(i, j)
        cost[i + 1, j + 1] = min(skip, take)
    end

    # backtrack
    keep = Int64[]
    i, j = n, m
    while i > 0
        if j > i && cost[i + 1, j + 1] == cost[i + 1, j]
            j -= 1
        else
            push!(keep, c[j])
            i -= 1
            j -= 1
        end
    end

    return sort(keep)
end
