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

    _chk3d(s)
    ch_n, _, ep_n = size(s)
    s_new = similar(s)

    @inbounds for ep_idx in 1:ep_n
        s_m = @views mean(s[:, :, ep_idx], dims=1)'[:, 1]
        m = mean(s_m)
        noise = rand(Float64, size(s_m)) .* m
        Threads.@threads for ch_idx in 1:ch_n
            s_new[ch_idx, :, ep_idx] = @views wiener(s[ch_idx, :, ep_idx], s_m, noise)
        end
    end

    return s_new

end

"""
    denoise_wien(obj; <keyword arguments>)

Perform Wiener deconvolution denoising.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names

# Returns
- `obj_new::NeuroAnalyzer.NEURO`
"""
function denoise_wien(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::NeuroAnalyzer.NEURO

    ch = get_channel(obj, ch=ch)
    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = @views denoise_wien(obj_new.data[ch, :, :])
    reset_components!(obj_new)
    push!(obj_new.history, "denoise_wien(OBJ, ch=$ch)")

    return obj_new
end

"""
    denoise_wien!(obj; <keyword arguments>)

Perform Wiener deconvolution denoising.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names

# Returns

Nothing
"""
function denoise_wien!(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::Nothing

    obj_new = denoise_wien(obj, ch=ch)
    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history

    return nothing

end
