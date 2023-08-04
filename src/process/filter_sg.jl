export filter_sg
export filter_sg!

"""
    filter_sg(s; order, window)

Filter using Savitzky-Golay filter.

# Arguments

- `s::AbstractVector`
- `order::Int64=6`: order of the polynomial used to fit the samples; must be less than `window`
- `window::Int64=11`: length of the filter window (i.e., the number of coefficients); must be an odd number

# Returns

- `s_filtered::Vector{Float64}`
"""
function filter_sg(s::AbstractVector; order::Int64=6, window::Int64=11)

    @assert !(window < 1 || window > length(s)) "window must be in [1, $(length(s))]."
    @assert isodd(window) "window must be an odd number."
    @assert order > 1 "order must be > 1."
    @assert order < window "order must be < $window."

    s_filtered = savitzky_golay(s, window, order).y

    return s_filtered

end

"""
    filter_sg(s; order, window)

Filter using Savitzky-Golay filter.

# Arguments

- `s::AbstractArray`
- `order::Int64=6`: order of the polynomial used to fit the samples; must be less than `window`
- `window::Int64=11`: length of the filter window (i.e., the number of coefficients); must be an odd number

# Returns

- `s_filtered::Array{Float64, 3}`: convoluted signal
"""
function filter_sg(s::AbstractArray; order::Int64=6, window::Int64=11)

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    s_filtered = similar(s)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            s_filtered[ch_idx, :, ep_idx] = @views filter_sg(s[ch_idx, :, ep_idx], order=order, window=window)
        end
    end

    return s_filtered

end

"""
    filter_sg(obj; ch, order, window)

Filter using Savitzky-Golay filter.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `order::Int64=6`: order of the polynomial used to fit the samples; must be less than `window`
- `window::Int64=11`: length of the filter window (i.e., the number of coefficients); must be an odd number

# Returns

- `obj_new::NeuroAnalyzer.NEURO`: convoluted signal
"""
function filter_sg(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), order::Int64=6, window::Int64=11)

    _check_channels(obj, ch)

    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = filter_sg(obj.data[ch, :, :], order=order, window=window)
    reset_components!(obj_new)
    push!(obj_new.history, "filter_sg(OBJ, ch=$ch, order=$order, window=$window")

    return obj_new

end

"""
    filter_sg!(obj; ch, order, window)

Filter using Savitzky-Golay filter.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `order::Int64=6`: order of the polynomial used to fit the samples; must be less than `window`
- `window::Int64=11`: length of the filter window (i.e., the number of coefficients); must be an odd number

"""
function filter_sg!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), order::Int64=6, window::Int64=11)

    obj_new = filter_sg(obj, ch=ch, order=order, window=window)
    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history

    return nothing

end
