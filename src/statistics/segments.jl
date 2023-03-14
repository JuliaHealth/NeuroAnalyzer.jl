export seg_mean

"""
    seg_mean(seg)

Calculate mean of a segment (e.g. spectrogram).

# Arguments

- `seg::AbstractArray`

# Returns

- `seg_mean::Vector{Float64}`: averaged segment
"""
function seg_mean(seg::AbstractArray)

    return reshape(mean(mean(seg, dims=1), dims=2), size(seg, 3))

end

"""
    seg2_mean(seg1, seg2)

Calculate mean of two segments (e.g. spectrograms).

# Arguments

- `seg1::AbstractArray`
- `seg2::AbstractArray`

# Returns

Named tuple containing:
- `seg1::Vector{Float64}`: averaged segment 1
- `seg2::Vector{Float64}`: averaged segment 2
"""
function seg_mean(seg1::AbstractArray, seg2::AbstractArray)

    return (seg1=seg_mean(seg1), seg2=seg_mean(seg2))
    
end
