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
function gfp(s::AbstractVector)::Float64

    return sum(s.^2) / length(s)

end

"""
    gfp_norm(s)

Calculate signal normalized for GFP (Global Field Power).

# Arguments

- `s::AbstractVector`

# Returns

- `gfp_norm::Vector{Float64}`
"""
function gfp_norm(s::AbstractVector)::Vector{Float64}

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
function diss(s1::AbstractVector, s2::AbstractVector)::@NamedTuple{gd::Float64, sc::Float64}

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
function diss(s::AbstractArray)::@NamedTuple{gd::Array{Float64, 3}, sc::Array{Float64, 3}}

    _chk3d(s)
    ch_n = size(s, 1)
    ep_n = size(s, 3)

    # initialize progress bar
    progress_bar && (progbar = Progress(ep_n * ch_n, dt=1, barlen=20, color=:white))

    gd = zeros(ch_n, ch_n, ep_n)
    sc = zeros(ch_n, ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx1 in 1:ch_n
           for ch_idx2 in 1:ch_idx1
                gd[ch_idx1, ch_idx2, ep_idx], sc[ch_idx1, ch_idx2, ep_idx] = @views diss(s[ch_idx1, :, ep_idx], s[ch_idx2, :, ep_idx])
            end

        # update progress bar
        progress_bar && next!(progbar)
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
- `gd::Matrix{Float64}`: global dissimilarity
- `sc::Matrix{Float64}`: spatial correlation
"""
function diss(s1::AbstractArray, s2::AbstractArray)::@NamedTuple{gd::Matrix{Float64}, sc::Matrix{Float64}}

    @assert size(s1) == size(s2) "s1 and s2 must have the same size."
    _chk3d(s1)
    _chk3d(s2)

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
    diss(obj; <keyword arguments>)

Calculate DISS (global dissimilarity) and spatial correlation (channels vs channels).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names

# Returns

Named tuple containing:
- `gd::Array{Float64, 3}`: global dissimilarity
- `sc::Array{Float64, 3}`: spatial correlation
"""
function diss(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::@NamedTuple{gd::Array{Float64, 3}, sc::Array{Float64, 3}}

    ch = get_channel(obj, ch=ch)
    gd, sc = diss(obj.data[ch, :, :])

    return (gd=gd, sc=sc)

end

"""
    diss(obj1, obj2; <keyword arguments>)

Calculate DISS (global dissimilarity) and spatial correlation (`ch1` of `obj1` vs `ch2` of `obj2`)..

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch1::Union{String, Vector{String}}: list of channels
- `ch2::Union{String, Vector{String}}: list of channels
- `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1))`: default use all epochs
- `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2))`: default use all epochs

# Returns

Named tuple containing:
- `gd::Array{Float64, 3}`: global dissimilarity
- `sc::Array{Float64, 3}`: spatial correlation
"""
function diss(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch1::Union{String, Vector{String}}, ch2::Union{String, Vector{String}}, ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1)), ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2)))::@NamedTuple{gd::Array{Float64, 3}, sc::Array{Float64, 3}}

    @assert sr(obj1) == sr(obj2) "OBJ1 and OBJ2 must have the same sampling rate."
    @assert length(ch1) == length(ch2) "ch1 and ch2 must have the same length."
    @assert length(ep1) == length(ep2) "ep1 and ep2 must have the same length."
    @assert epoch_len(obj1) == epoch_len(obj2) "OBJ1 and OBJ2 must have the same epoch lengths."

    ch1 = get_channel(obj1, ch=ch1)
    ch2 = get_channel(obj2, ch=ch2)
    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    length(ep1) == 1 && (ep1 = [ep1])
    length(ep2) == 1 && (ep2 = [ep2])

    gd, sc = @views diss(obj1.data[ch1, :, ep1], obj2.data[ch2, :, ep2])

    return (gd=gd, sc=sc)

end
