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
function mutual_information(s1::AbstractVector, s2::AbstractVector)

    return get_mutual_information(s1, s2)

end

"""
    mutual_information(s1, s2)

Calculate mutual information (channels of `s1` vs channels of `s2`).

# Arguments

- `s1::AbstractArray`
- `s2::AbstractArray`

# Returns

- `mutual_information::Array{Float64}`
"""
function mutual_information(s1::AbstractArray, s2::AbstractArray)

    @assert size(s1) == size(s2) "s1 and s2 must have the same size."

    ch_n = size(s1, 1)
    ep_n = size(s1, 3)

    m = zeros(ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
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

"""
function mutual_information(s::AbstractArray)

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    # initialize progress bar
    progress_bar && (progbar = Progress(ep_n * ch_n, dt=1, barlen=20, color=:white))

    m = zeros(ch_n, ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx1 in 1:ch_n
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

Calculate mutual information between channels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, AbstractRange}: list of channels

# Returns

- `mutual_information::Array{Float64, 3}`
"""
function mutual_information(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, AbstractRange})

    ch = get_channel(obj, ch=ch)
    ch_n = length(ch)
    ep_n = nepochs(obj)
    isa(ch, Int64) && (ch = [ch])

    m = zeros(ch_n, ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        # create half of the matrix
        Threads.@threads for ch_idx1 in 1:ch_n
            for ch_idx2 in 1:ch_idx1
                m[ch_idx1, ch_idx2, ep_idx] = @views mutual_information(obj.data[ch[ch_idx1], :, ep_idx], obj.data[ch[ch_idx2], :, ep_idx])
            end
        end

    end

    # copy to the other half
    m = _copy_lt2ut(m)

    return m

end

"""
    mutual_information(obj1, obj2; <keyword arguments>)

Calculate mutual information between two channels.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `ch1::Union{String, Vector{String}}: list of channels
- `ch2::Union{String, Vector{String}}: list of channels
- `ep1::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj1))`: default use all epochs
- `ep2::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj2))`: default use all epochs

# Returns

- `m::Array{Float64, 3}`
"""
function mutual_information(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch1::Union{String, Vector{String}}, ch2::Union{String, Vector{String}}, ep1::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj1)), ep2::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj2)))

    # check channels
    _check_channels(obj1, ch1)
    _check_channels(obj2, ch2)
    @assert length(ch1) == length(ch2) "ch1 and ch2 must have the same length."

    # check epochs
    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    @assert length(ep1) == length(ep2) "ep1 and ep2 must have the same length."
    @assert epoch_len(obj1) == epoch_len(obj2) "OBJ1 and OBJ2 must have the same epoch lengths."

    isa(ch1, Int64) && (ch1 = [ch1])
    isa(ch2, Int64) && (ch2 = [ch2])
    length(ep1) == 1 && (ep1 = [ep1])
    length(ep2) == 1 && (ep2 = [ep2])

    m = @views mutual_information(obj1.data[ch1, :, ep1], obj2.data[ch2, :, ep2])

    return m

end