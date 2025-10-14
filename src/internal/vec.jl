function _flipx(s::AbstractVector)::Vector{Float64}
    # center at Y=0
    m = mean(s)
    s_new = s .- m
    # flip the signal along the X axis
    s_new = .-s_new
    s_new .+= m
    return s_new
end

_zeros(s::AbstractVector)::Int64 = count(abs.(diff(sign.(s))) .!= 0)

function _split(s::AbstractVector; wlen::Int64, woverlap::Int64=round(Int64, wlen * 0.97))::Vector{Vector{Float64}}
    # get all segments available, the last one may be shorter
    s_new = Vector{Vector{Float64}}()
    n = length(s)
    idx = 1
    while (idx - 1 + wlen) <= n
        push!(s_new, s[idx:(idx - 1 + wlen)])
        idx += woverlap
    end
    push!(s_new, s[idx:n])
    return s_new
end

function _fsplit(s::AbstractVector; wlen::Int64, woverlap::Int64=round(Int64, wlen * 0.97))::Vector{Vector{Float64}}
    # get only complete segments
    s_new = Vector{Vector{Float64}}()
    n = length(s)
    idx = 1
    while (idx - 1 + wlen) <= n
        push!(s_new, s[idx:(idx - 1 + wlen)])
        idx += woverlap
    end
    return s_new
end

function _chunks(n::Int64; wlen::Int64, woverlap::Int64=round(Int64, wlen * 0.97))::Matrix{Int64}
    # get all chunks available, the last one may be shorter
    chunks_idx = Matrix{Int64}(undef, 0, 2)
    idx = 1
    while (idx - 1 + wlen) <= n
        chunks_idx = vcat(chunks_idx, [idx (idx - 1 + wlen)])
        idx += woverlap
    end
    chunks_idx = vcat(chunks_idx, [idx n])
    return chunks_idx
end

function _chunks(s::AbstractVector; wlen::Int64, woverlap::Int64=round(Int64, wlen * 0.97))::Matrix{Int64}
    # get all chunks available, the last one may be shorter
    chunks_idx = Matrix{Int64}(undef, 0, 2)
    n = length(s)
    idx = 1
    while (idx - 1 + wlen) <= n
        chunks_idx = vcat(chunks_idx, [idx (idx - 1 + wlen)])
        idx += woverlap
    end
    chunks_idx = vcat(chunks_idx, [idx n])
    return chunks_idx
end

function _fchunks(n::Int64; wlen::Int64, woverlap::Int64=round(Int64, wlen * 0.97))::Matrix{Int64}
    # get only full chunks
    chunks_idx = Matrix{Int64}(undef, 0, 2)
    idx = 1
    while (idx - 1 + wlen) <= n
        chunks_idx = vcat(chunks_idx, [idx (idx - 1 + wlen)])
        idx += woverlap
    end
    return chunks_idx
end

function _fchunks(s::AbstractVector; wlen::Int64, woverlap::Int64=round(Int64, wlen * 0.97))::Matrix{Int64}
    # get only full chunks
    chunks_idx = Matrix{Int64}(undef, 0, 2)
    n = length(s)
    idx = 1
    while (idx - 1 + wlen) <= n
        chunks_idx = vcat(chunks_idx, [idx (idx - 1 + wlen)])
        idx += woverlap
    end
    return chunks_idx
end

