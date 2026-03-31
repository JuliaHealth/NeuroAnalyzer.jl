export normpower
export normpower!

"""
    normpower(s)

Return a signal with normalized power (amplitudes divided by the root-mean-squared value of the entire signal).

# Arguments

- `s::AbstractVector`: signal vector

# Returns

- `Vector{Float64}`
"""
function normpower(s::AbstractVector)::Vector{Float64}
    return s .* amp(s).rms_amp
end

"""
    normpower(s)

Return a signal with normalized power (amplitudes divided by the root-mean-squared value of the entire signal) for a 3-D signal array.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)

# Returns

- `Array{Float64, 3}`
"""
function normpower(s::AbstractArray)::Array{Float64, 3}

    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # pre-allocate output
    s_new = similar(s, Float64)

    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        s_new[ch_idx, :, ep_idx] = normpower(@view(s[ch_idx, :, ep_idx]))
    end

    return s_new
end

"""
    normpower(obj; <keyword arguments>)

Return a signal with normalized power (amplitudes divided by the root-mean-squared value of the entire signal).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function normpower(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
)::NeuroAnalyzer.NEURO

    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))

    # create new dataset
    obj_new = deepcopy(obj)

    obj_new.data[ch, :, :] = normpower(obj.data[ch, :, :])
    push!(obj_new.history, "normpower(obj; ch=$ch)")

    return obj_new
end

"""
    normpower!(obj; <keyword arguments>)

Return a signal with normalized power (amplitudes divided by the root-mean-squared value of the entire signal).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)

# Returns

- `Nothing`
"""
function normpower!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
)::Nothing
    obj_new = normpower(obj; ch = ch)
    obj.data = obj_new.data
    obj.history = obj_new.history

    return nothing
end
