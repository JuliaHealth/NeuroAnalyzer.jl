export denoise_wien
export denoise_wien!

"""
    denoise_wien(s)

Perform Wiener deconvolution denoising on a 3-D signal array.

For each epoch, the cross-channel mean signal is used as the reference, and a noise estimate is generated as white noise scaled to the mean signal power. The Wiener filter is then applied independently to each channel.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs).

# Returns

- `Array{Float64, 3}`: denoised array of the same shape as `s`

# Notes

- The noise estimate is random (`rand`); results are not reproducible unless a random seed is set by the caller
"""
function denoise_wien(s::AbstractArray)::AbstractArray
    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # pre-allocate output
    s_new = similar(s, Float64)

    # thread over epochs; the inner channel loop must remain sequential because
    # s_m and noise are shared across all channels of the same epoch.
    Threads.@threads :static for ep_idx in 1:ep_n
        # cross-channel mean for this epoch: shape (samples,)
        s_m = vec(mean(@view(s[:, :, ep_idx]), dims = 1))
        m = mean(s_m)
        # Noise estimate: white noise at mean signal power
        noise = rand(Float64, length(s_m)) .* m
        for ch_idx in 1:ch_n
            @inbounds s_new[ch_idx, :, ep_idx] =
                wiener(@view(s[ch_idx, :, ep_idx]), s_m, noise)
        end
    end

    return s_new
end

"""
    denoise_wien(obj; <keyword arguments>)

Perform Wiener deconvolution denoising on selected channels of a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)

# Returns

- `NeuroAnalyzer.NEURO`: new object with denoised channels
"""
function denoise_wien(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
)::NeuroAnalyzer.NEURO
    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))

    # create new dataset
    obj_new = deepcopy(obj)

    obj_new.data[ch, :, :] = denoise_wien(@view(obj.data[ch, :, :]))
    push!(obj_new.history, "denoise_wien(obj; ch=$ch)")

    return obj_new
end

"""
    denoise_wien!(obj; <keyword arguments>)

Perform Wiener deconvolution denoising in-place on selected channels of a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; modified in-place
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)

# Returns

- `Nothing`
"""
function denoise_wien!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
)::Nothing
    obj_new = denoise_wien(obj; ch = ch)
    obj.data = obj_new.data
    obj.history = obj_new.history

    return nothing
end
