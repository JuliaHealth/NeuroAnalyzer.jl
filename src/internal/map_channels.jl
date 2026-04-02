"""
    _map_channels(ch, chs)

Map channel indices from a full channel list into their positions within a subset `chs`, returning both the remapped indices and the original indices.

# Arguments

- `ch::Union{Int64, Vector{Int64}}`: channel index or indices to remap
- `chs::Vector{Int64}`: the channel subset to map into

# Returns

- `Tuple{Union{Int64, Vector{Int64}}, Union{Int64, Vector{Int64}}}`: `(mapped, original)` where `mapped` contains the positions of each element of `ch` within `chs`, and `original` is the unmodified input `ch`
"""
function _map_channels(
    ch::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}},
    chs = Vector{Int64},
)::Tuple{Union{Int64, Vector{Int64}}, Union{Int64, Vector{Int64}}}
    ch_orig = ch

    if ch isa Int64
        ch_mapped = vsearch(ch, chs)
    else
        ch_orig   = copy(ch)
        ch_mapped = [vsearch(ch[idx], chs) for idx in eachindex(ch)]
    end

    return ch_mapped, ch_orig
end
