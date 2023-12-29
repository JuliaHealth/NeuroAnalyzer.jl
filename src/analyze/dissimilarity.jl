export gfp
export gfp_norm
export diss

"""
    gfp(s)

Calculate GFP (Global Field Power).

# Arguments

- `s::AbstractVector`

# Returns

- `gfp::Float64`
"""
function gfp(s::AbstractVector)

    return sum(s.^2) / length(s)

end

"""
    gfp_norm(s)

Calculate signal normalized for GFP (Global Field Power).

# Arguments

- `s::AbstractVector`

# Returns

- `gfp_norm::Float64`
"""
function gfp_norm(s::AbstractVector)

    return s ./ gfp(s)

end

"""
    diss(s1, s2)

Calculate DISS (global dissimilarity) and spatial correlation between `s1` and `s2`.

# Arguments

- `s1::AbstractVector`
- `s2::AbstractVector`

# Returns

Named tuple containing:
- `gd::Float64`: global dissimilarity
- `sc::Float64`: spatial correlation
"""
function diss(s1::AbstractVector, s2::AbstractVector)
    
    @assert length(s1) == length(s2) "s1 and s2 must have the same length."

    gfp_norm1 = gfp_norm(s1)
    gfp_norm2 = gfp_norm(s2)
    gd = sqrt(sum((gfp_norm1 .- gfp_norm2).^2) / length(s1))
    sc = 0.5 * (2 - gd^2)

    return (gd=gd, sc=sc)

end

"""
    diss(s)

Calculate DISS (global dissimilarity) and spatial correlation (channels vs channels).

# Arguments

- `s::AbstractArray`

# Returns

Named tuple containing:
- `gd::Array{Float64, 3}`: global dissimilarity
- `sc::Array{Float64, 3}`: spatial correlation
"""
function diss(s::AbstractArray)

    ch_n = size(s, 1)
    ep_n = size(s, 3)
    
    # initialize progress bar
    progress_bar == true && (progbar = Progress(ep_n * ch_n, dt=1, barlen=20, color=:white))

    gd = zeros(ch_n, ch_n, ep_n)
    sc = zeros(ch_n, ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx1 in 1:ch_n
           for ch_idx2 in 1:ch_idx1
                gd[ch_idx1, ch_idx2, ep_idx], sc[ch_idx1, ch_idx2, ep_idx] = @views diss(s[ch_idx1, :, ep_idx], s[ch_idx2, :, ep_idx])
            end

        # update progress bar
        progress_bar == true && next!(progbar)
        end
    end

    # copy lower triangle to upper triangle
    @inbounds for ep_idx in 1:ep_n
        gd[:, :, ep_idx] = gd[:, :, ep_idx] + gd[:, :, ep_idx]' - diagm(diag(gd[:, :, ep_idx]))
        sc[:, :, ep_idx] = sc[:, :, ep_idx] + sc[:, :, ep_idx]' - diagm(diag(sc[:, :, ep_idx]))
    end

    return (gd=gd, sc=sc)

end

"""
    diss(s1, s2)

Calculate DISS (global dissimilarity) and spatial correlation (channels vs channels).

# Arguments

- `s1::AbstractArray`
- `s2::AbstractArray`

# Returns

Named tuple containing:
- `gd::Array{Float64, 3}`: global dissimilarity
- `sc::Array{Float64, 3}`: spatial correlation
"""
function diss(s1::AbstractArray, s2::AbstractArray)

    @assert size(s1) == size(s2) "s1 and s2 must have the same size."

    ch_n = size(s1, 1)
    ep_n = size(s1, 3)

    gd = zeros(ch_n, ep_n)
    sc = zeros(ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            gd[ch_idx, ep_idx], sc[ch_idx, ep_idx] = @views diss(s1[ch_idx, :, ep_idx], s2[ch_idx, :, ep_idx])
        end
    end

    return (gd=gd, sc=sc)

end

"""
    diss(obj; ch)

Calculate DISS (global dissimilarity) and spatial correlation (channels vs channels).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

# Returns

Named tuple containing:
- `gd::Array{Float64, 3}`: global dissimilarity
- `sc::Array{Float64, 3}`: spatial correlation
"""
function diss(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj))

    _check_channels(obj, ch)

    gd, sc = diss(obj.data[ch, :, :])

    return (gd=gd, sc=sc)

end

"""
    diss(obj1, obj2; ch1, ch2, ep1, ep2)

Calculate DISS (global dissimilarity) and spatial correlation (`ch1` of `obj1` vs `ch2` of `obj2`)..

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch1::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj1)`: index of channels, default is all signal channels
- `ch2::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj2)`: index of channels, default is all signal channels
- `ep1::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj1))`: default use all epochs
- `ep2::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj2))`: default use all epochs

# Returns

Named tuple containing:
- `gd::Array{Float64, 3}`: global dissimilarity
- `sc::Array{Float64, 3}`: spatial correlation
"""
function diss(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch1::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj1), ch2::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj2), ep1::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj1)), ep2::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj2)))

    @assert sr(obj1) == sr(obj2) "OBJ1 and OBJ2 must have the same sampling rate."

    _check_channels(obj1, ch1)
    _check_channels(obj2, ch2)
    @assert length(ch1) == length(ch2) "ch1 and ch2 must have the same length."
    
    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    @assert length(ep1) == length(ep2) "ep1 and ep2 must have the same length."
    @assert epoch_len(obj1) == epoch_len(obj2) "OBJ1 and OBJ2 must have the same epoch lengths."

    gd, sc = @views diss(reshape(obj1.data[ch1, :, ep1], length(ch1), :, length(ep1)), reshape(obj2.data[ch2, :, ep2], length(ch2), :, length(ep2)))

    return (gd=gd, sc=sc)
    
end
