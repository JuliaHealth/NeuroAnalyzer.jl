export tconv

"""
    tconv(signal; kernel)

Performs convolution in the time domain.

# Arguments

- `s::AbstractVector`
- `kernel::AbstractVector`

# Returns

- `s_new::Vector{Float64}`
"""
function tconv(s::AbstractVector; kernel::AbstractVector)

    s_new = DSP.conv(s, kernel)

    half_kernel = floor(Int, length(kernel) / 2)

    # remove in- and out- edges
    if mod(length(kernel), 2) == 0 
        return real.(s_new)[half_kernel:(end - half_kernel)]
    else
        return real.(s_new)[(half_kernel + 1):(end - half_kernel)]
    end

end

"""
    tconv(s; kernel)

Perform convolution in the time domain.

# Arguments

- `s::AbstractArray`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `kernel::Union{Vector{<:Real}, Vector{ComplexF64}}`: kernel used for convolution

# Returns

- `s_new::Array{Float64, 3}`: convoluted signal
"""
function tconv(s::AbstractArray; kernel::Union{Vector{<:Real}, Vector{ComplexF64}})

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    s_new = similar(s)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            s_new[ch_idx, :, ep_idx] = tconv(s[ch_idx, :, ep_idx], kernel=kernel)
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
- `kernel::Union{Vector{<:Real}, Vector{ComplexF64}}`: kernel used for convolution

# Returns

- `s_new::Array{Float64, 3}`: convoluted signal
"""
function tconv(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), kernel::Union{Vector{<:Real}, Vector{ComplexF64}})

    _check_channels(obj, ch)
    s_new = tconv(obj.data[ch, :, :], kernel=kernel)

    return s_new

end
