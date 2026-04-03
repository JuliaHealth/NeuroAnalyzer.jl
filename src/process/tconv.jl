export tconv
export tconv!

"""
    tconv(s; <keyword arguments>)

Performs convolution in the time domain.

# Arguments

- `s::AbstractVector`: signal vector
- `kernel::AbstractVector`

# Returns

- `Union{Vector{Float64}, Vector{ComplexF64}}`: convoluted signal
"""
function tconv(
    s::AbstractVector;
    kernel::AbstractVector,
)::Union{Vector{Float64}, Vector{ComplexF64}}
    s_conv = DSP.conv(s, kernel)
    return _remove_kernel(s_conv, kernel)
end

"""
    tconv(s; <keyword arguments>)

Perform convolution in the time domain for a 3-D signal array.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `kernel::AbstractVector`: convolution kernel

# Returns

- `Union{Array{Float64, 3}, Array{ComplexF64, 3}}`: convoluted signal
"""
function tconv(
    s::AbstractArray;
    kernel::AbstractVector,
)::Union{Array{Float64, 3}, Array{ComplexF64, 3}}
    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # pre-allocate output
    s_new = zeros(eltype(kernel), size(s))

    # initialize progress bar
    progbar =
        Progress(ep_n * ch_n; dt = 1, barlen = 20, color = :white, enabled = progress_bar)

    # calculate over channels and epochs
    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        s_new[ch_idx, :, ep_idx] = tconv(@view(s[ch_idx, :, ep_idx]), kernel = kernel)
        # update progress bar
        progress_bar && next!(progbar)
    end

    return s_new
end

"""
    tconv(obj; <keyword arguments>)

Perform convolution in the time domain.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `kernel::AbstractVector`: convolution kernel

# Returns

- `Union{NeuroAnalyzer.NEURO, Array{ComplexF64, 3}}`: convoluted signal
"""
function tconv(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    kernel::AbstractVector,
)::Union{NeuroAnalyzer.NEURO, Array{ComplexF64, 3}}
    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))

    # create new dataset
    obj_new = deepcopy(obj)
    _info("Group delay: $(_group_delay(kernel)) samples")

    if eltype(kernel) == ComplexF64
        return tconv(obj.data[ch, :, :]; kernel = kernel)
    else
        obj_new.data[ch, :, :] = tconv(obj.data[ch, :, :]; kernel = kernel)
        push!(obj_new.history, "tconv(obj; ch=$ch, kernel=kernel)")
        return obj_new
    end
end

"""
    tconv!(obj; <keyword arguments>)

Perform convolution in the time domain.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `kernel::AbstractVector`: convolution kernel
"""
function tconv!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    kernel::AbstractVector,
)::Union{Nothing, Array{ComplexF64, 3}}
    if eltype(kernel) == ComplexF64
        return tconv(obj.data; ch = ch, kernel = kernel)
    else
        obj_new = tconv(obj; ch = ch, kernel = kernel)
        obj.data = obj_new.data
        obj.history = obj_new.history
        return nothing
    end
end
