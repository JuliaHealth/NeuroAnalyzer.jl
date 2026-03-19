export denoise_wien
export denoise_wien!

"""
    denoise_wien(s)

Perform Wiener deconvolution denoising.

# Arguments

- `s::AbstractArray`

# Returns

- `s_new::AbstractArray`
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

    # calculate over channel and epochs
    @inbounds Threads.@threads :dynamic for ep_idx in 1:ep_n
        s_m = @views mean(s[:, :, ep_idx], dims = 1)'[:, 1]
        m = mean(s_m)
        noise = rand(Float64, size(s_m)) .* m
        for ch_idx in 1:ch_n
            s_new[ch_idx, :, ep_idx] = wiener(@view(s[ch_idx, :, ep_idx]), s_m, noise)
        end
    end

    return s_new

end

"""
    denoise_wien(obj; <keyword arguments>)

Perform Wiener deconvolution denoising.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)

# Returns

- `obj_new::NeuroAnalyzer.NEURO`: output NEURO object
"""
function denoise_wien(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::NeuroAnalyzer.NEURO

    ch = get_channel(obj, ch = ch)
    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = @views denoise_wien(obj_new.data[ch, :, :])
    push!(obj_new.history, "denoise_wien(OBJ, ch=$ch)")

    return obj_new
end

"""
    denoise_wien!(obj; <keyword arguments>)

Perform Wiener deconvolution denoising.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)

# Returns

- `Nothing`
"""
function denoise_wien!(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::Nothing

    obj_new = denoise_wien(obj, ch = ch)
    obj.data = obj_new.data
    obj.history = obj_new.history

    return nothing

end
