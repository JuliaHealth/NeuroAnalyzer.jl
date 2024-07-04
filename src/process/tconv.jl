export tconv
export tconv!

"""
    tconv(s; kernel)

Performs convolution in the time domain.

# Arguments

- `s::AbstractVector`
- `kernel::AbstractVector`

# Returns

- `s_new::Vector{Float64}`: convoluted signal
"""
function tconv(s::AbstractVector; kernel::AbstractVector)

    s_new = DSP.conv(s, kernel)

    half_kernel = floor(Int, length(kernel) / 2)

    # remove in- and out- edges
    if mod(length(kernel), 2) == 0
        return s_new[half_kernel:(end - half_kernel)]
    else
        return s_new[half_kernel:(end - half_kernel - 1)]
    end

end

"""
    tconv(s; kernel)

Perform convolution in the time domain.

# Arguments

- `s::AbstractArray`
- `kernel::AbstractVector`: convolution kernel

# Returns

- `s_new::Array{Float64, 3}`: convoluted signal
"""
function tconv(s::AbstractArray; kernel::AbstractVector)

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    s_new = similar(s)

    # initialize progress bar
    progress_bar && (progbar = Progress(ep_n * ch_n, dt=1, barlen=20, color=:white))

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            s_new[ch_idx, :, ep_idx] = tconv(s[ch_idx, :, ep_idx], kernel=kernel)

            # update progress bar
            progress_bar && next!(progbar)
        end
    end

    return s_new

end

"""
    tconv(obj; ch, kernel)

Perform convolution in the time domain.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `kernel::AbstractVector`: convolution kernel

# Returns

- `obj_new::NeuroAnalyzer.NEURO`: convoluted signal
"""
function tconv(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), kernel::AbstractVector)

    _check_channels(obj, ch)
    isa(ch, Int64) && (ch = [ch])

    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = tconv(obj.data[ch, :, :], kernel=kernel)
    reset_components!(obj_new)
    push!(obj_new.history, "tconv(OBJ, ch=$ch, kernel=$kernel)")

    return obj_new

end

"""
    tconv!(obj; ch, kernel)

Perform convolution in the time domain.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nchannels(obj))`: index of channels, default is all channels
- `kernel::AbstractVector`: convolution kernel
"""
function tconv!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), kernel::AbstractVector)

    obj_new = tconv(obj, ch=ch, kernel=kernel)
    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history

    return nothing

end
