function _dict2labeled_matrix(d::Dict; rev::Bool = true)::Tuple{Vector{String}, Vector{Vector{Float64}}}
    l = Vector{String}()
    v = Vector{Vector{Float64}}()
    for (kk, vv) in d
        push!(l, kk)
        push!(v, vv)
    end
    if rev
        return reverse!(l), reverse!(v)
    else
        return l, v
    end
end

function _labeled_matrix2dict(l::Vector{String}, v::Vector{Vector{Float64}})::Dict
    @assert length(l) == length(v) "Length of labels and values do not match."
    return Dict(zip(l, v))
end
