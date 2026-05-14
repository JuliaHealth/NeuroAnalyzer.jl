"""
    _select_indices(input, n_total, default)

Core logic shared by `_select_channels` and `_select_epochs`.

Resolves `input` to a concrete, sorted index selection according to these rules:

|    `input` value    |                   Behavior                   |
| ------------------- | -------------------------------------------- |
| `0` (sentinel)      | Expands to `collect(1:eff_default)`          |
| `Int64` ≠ 0         | Returned as-is                               |
| `AbstractUnitRange` | Materialized to `Vector{Int64}`, then sorted |
| `Vector{Int64}`     | Sorted in-place (no copy) when length > 1    |

The effective default is derived from `default`:
- `0`          → `n_total` (select everything)
- `> n_total`  → clamped down to `n_total`
- `< 0`        → clamped up to `1`

# Arguments

- `input::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}}`: raw user-supplied selection
- `n_total::Int64`: total number of available indices (channels or epochs)
- `default ::Int64=0`: upper bound for the "select all" expansion; `0` means `n_total`

# Returns

- `Union{Int64, Vector{Int64}}`

# Notes

- The sentinel `0` is only meaningful for scalar `Int64` input. Passing a `Vector` or range whose first element is `0` is an upstream error and is not handled here.
- Out-of-range indices in a vector input are **not** validated; that is the caller's responsibility.
"""
function _select_indices(
    input::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}},
    n_total::Int64,
    default::Int64 = 0,
)::Union{Int64, Vector{Int64}}

    # resolve the effective upper bound for "select all" expansion
    eff_default = clamp(default == 0 ? n_total : default, 1, n_total)

    # scalar sentinel 0 → expand to the full default range
    input isa Int64 && input == 0 && return collect(1:eff_default)

    # materialize a range to a concrete vector
    result = input isa AbstractUnitRange{Int64} ? collect(input) : input

    # sort multi-element vectors in-place (single elements are trivially sorted)
    result isa Vector{Int64} && length(result) > 1 && sort!(result)

    return result
end

# ─────────────────────────────────────────────────────────────────────────────

"""
    _select_channels(obj, input, n_total, def_chn)

Resolve a channel selection for `obj`.

`channel = 0` selects all channels up to `def_chn` (or `nchannels(obj)` when `def_chn = 0`). See `_select_indices` for the full resolution rules.

# Returns

- `Union{Int64, Vector{Int64}}`: channel index/indices
"""
function _select_channels(
    obj::NeuroAnalyzer.NEURO,
    channel::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}},
    def_chn::Int64 = 0,
)::Union{Int64, Vector{Int64}}
    return _select_indices(channel, nchannels(obj), def_chn)
end

"""
    _select_epochs(obj, input, n_total, def_ep)

Resolve an epoch selection for `obj`.

`epoch = 0` selects all epochs up to `def_ep` (or `nepochs(obj)` when `def_ep = 0`).  See `_select_indices` for the full resolution rules.

# Returns

- `Union{Int64, Vector{Int64}}`: epoch index/indices
"""
function _select_epochs(
    obj::NeuroAnalyzer.NEURO,
    epoch::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}},
    def_ep::Int64 = 0,
)::Union{Int64, Vector{Int64}}
    return _select_indices(epoch, nepochs(obj), def_ep)
end
