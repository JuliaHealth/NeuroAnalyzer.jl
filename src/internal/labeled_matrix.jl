function _dict2labeled_matrix(d::Dict; rev::Bool=true)
    l = Vector{String}()
    v = Vector{Vector{Float64}}()
    for (kk, vv) in d
        push!(l, kk)
        push!(v, vv)
    end
    rev == true && return reverse!(l), reverse!(v)
    rev == false && return l, v
end

function _labeled_matrix2dict(l::Vector{String}, v::Vector{Vector{Float64}})
    length(l) == length(v) || throw(ArgumentError("Length of labels and values do not match."))
    return Dict(zip(l, v))
end
