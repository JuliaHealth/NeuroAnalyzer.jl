export mutual_information

"""
    mutual_information(s1, s2)

Calculate mutual information.

# Arguments

- `s1::AbstractVector`
- `s2::AbstractVector`

# Returns

- `mutual_information::Float64`
"""
function mutual_information(s1::AbstractVector, s2::AbstractVector)::Float64

    return get_mutual_information(s1, s2)

end

"""
    mutual_information(s1, s2)

Calculate mutual information (channels of `s1` vs channels of `s2`).

# Arguments

- `s1::AbstractArray`
- `s2::AbstractArray`

# Returns

- `m::Matrix{Float64}`
"""
function mutual_information(s1::AbstractArray, s2::AbstractArray)::Matrix{Float64}

    @assert size(s1) == size(s2) "s1 and s2 must have the same size."
    _chk3d(s1)
    _chk3d(s2)

    ch_n = size(s1, 1)
    ep_n = size(s1, 3)

    m = zeros(ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads :greedy for ch_idx in 1:ch_n
            m[ch_idx, ep_idx] = @views mutual_information(s1[ch_idx, :, ep_idx], s2[ch_idx, :, ep_idx])
        end
    end

    return m

end

"""
    mutual_information(s)

Calculate mutual information (channels vs channels).

# Arguments

- `s::AbstractArray`

# Returns

- `m::Array{Float64, 3}`
"""
function mutual_information(s::AbstractArray)::Array{Float64, 3}

    _chk3d(s)
    ch_n = size(s, 1)
    ep_n = size(s, 3)

    # initialize progress bar
    progress_bar && (progbar = Progress(ep_n * ch_n, dt=1, barlen=20, color=:white))

    m = zeros(ch_n, ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads :greedy for ch_idx1 in 1:ch_n
           for ch_idx2 in 1:ch_idx1
                m[ch_idx1, ch_idx2, ep_idx] = @views mutual_information(s[ch_idx1, :, ep_idx], s[ch_idx2, :, ep_idx])
            end

        # update progress bar
        progress_bar && next!(progbar)
        end
    end

    # copy lower triangle to upper triangle
    m = _copy_lt2ut(m)

    return m

end

"""
    mutual_information(obj; <keyword arguments>)

Calculate mutual information between channels. Currently only one estimator (maximum likelihood) is available. Internally it uses `InformationMeasures.get_mutual_information()`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names

# Returns

- `m::Array{Float64, 3}`
"""
function mutual_information(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::Array{Float64, 3}

    ch = exclude_bads ? get_channel(obj, ch=ch, exclude="bad") : get_channel(obj, ch=ch, exclude="")
    m = @views mutual_information(obj.data[ch, :, :])

    return m

end

"""
    mutual_information(obj1, obj2; <keyword arguments>)

Calculate mutual information between two channels. Currently only one estimator (maximum likelihood) is available. Internally it uses `InformationMeasures.get_mutual_information()`.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `ch1::Union{String, Vector{String}}`: channel name or list of channel names
- `ch2::Union{String, Vector{String}}`: channel name or list of channel names
- `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1))`: default use all epochs
- `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2))`: default use all epochs

# Returns

- `m::Matrix{Float64}`
"""
function mutual_information(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch1::Union{String, Vector{String}}, ch2::Union{String, Vector{String}}, ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1)), ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2)))::Matrix{Float64}

    # check channels
    ch1 = exclude_bads ? get_channel(obj1, ch=ch1, exclude="bad") : get_channel(obj1, ch=ch1, exclude="")
    ch2 = exclude_bads ? get_channel(obj2, ch=ch2, exclude="bad") : get_channel(obj2, ch=ch2, exclude="")
    @assert length(ch1) == length(ch2) "ch1 and ch2 must have the same length."

    # check epochs
    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    @assert length(ep1) == length(ep2) "ep1 and ep2 must have the same length."
    @assert epoch_len(obj1) == epoch_len(obj2) "OBJ1 and OBJ2 must have the same epoch lengths."

    length(ep1) == 1 && (ep1 = [ep1])
    length(ep2) == 1 && (ep2 = [ep2])

    m = @views mutual_information(obj1.data[ch1, :, ep1], obj2.data[ch2, :, ep2])

    return m

end
