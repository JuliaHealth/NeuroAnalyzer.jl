export corr

"""
    corr(obj1, obj2; <keyword arguments>)

Compute the Pearson correlation between paired channels across two NEURO objects.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `ch1::Union{String, Vector{String}}: channel name(s)
- `ch2::Union{String, Vector{String}}: channel name(s)
- `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1))`: epoch number(s)
- `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2))`: epoch number(s)

# Returns

- `cr::Matrix{Float64}`: correlation coefficients of shape `(channels, epochs)`
"""
function corr(
    obj1::NeuroAnalyzer.NEURO,
    obj2::NeuroAnalyzer.NEURO;
    ch1::Union{String, Vector{String}},
    ch2::Union{String, Vector{String}},
    ep1::Union{Int64, Vector{Int64}, AbstractRange} = _c(nepochs(obj1)),
    ep2::Union{Int64, Vector{Int64}, AbstractRange} = _c(nepochs(obj2)),
)::Matrix{Float64}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch1 = exclude_bads ? get_channel(obj1, ch = ch1, exclude = "bad") : get_channel(obj1, ch = ch1, exclude = "")
    ch2 = exclude_bads ? get_channel(obj2, ch = ch2, exclude = "bad") : get_channel(obj2, ch = ch2, exclude = "")
    @assert length(ch1) == length(ch2) "Lengths of ch1 ($(length(ch1))) and ch2 ($(length(ch2))) must be equal."

    # validate epoch indices and ensure both objects have matching epoch structure
    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    # normalize scalar epoch arguments to vectors so indexing is uniform
    isa(ep1, Int64) && (ep1 = [ep1])
    isa(ep2, Int64) && (ep2 = [ep2])
    @assert length(ep1) == length(ep2) "Lengths of ep1 ($(length(ep1))) and ep2 ($(length(ep2))) must be equal."
    @assert epoch_len(obj1) == epoch_len(obj2) "OBJ1 and OBJ2 must have the same epoch lengths."

    # number of channels
    ch_n = length(ch1)
    # number of epochs
    ep_n = length(ep1)

    # pre-allocate output
    cr = zeros(ch_n, ep_n)

    # calculate over channel and epochs
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        cr[ch_idx, ep_idx] = cor(
            @view(obj1.data[ch1[ch_idx], :, ep1[ep_idx]]),
            @view(obj2.data[ch2[ch_idx], :, ep2[ep_idx]]),
        )
    end

    return cr

end

"""
    corr(obj; <keyword arguments>)

Compute the Pearson correlation between all channel pairs within a single NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)

# Returns

- `cr::Array{Float64, 3}`: correlation coefficient of shape `(channels, channels, epochs)`
"""
function corr(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::Array{Float64, 3}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")

    # number of channels
    ch_n = length(ch)
    # number of epochs
    ep_n = nepochs(obj)

    # pre-allocate symmetric output: channels × channels × epochs
    # only the lower triangle is computed
    cr = zeros(ch_n, ch_n, ep_n)

    # calculate over channel and epochs
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx1, ep_idx = idx[1], idx[2]
        for ch_idx2 in 1:ch_idx1
            cr[ch_idx1, ch_idx2, ep_idx] = cor(
                @view(obj.data[ch[ch_idx1], :, ep_idx]),
                @view(obj.data[ch[ch_idx2], :, ep_idx]),
            )
        end
    end

    # mirror the lower triangle to the upper triangle to produce the full symmetric matrix
    cr = _copy_lt2ut(cr)

    return cr

end
