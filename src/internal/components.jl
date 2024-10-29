function _get_component(obj::NeuroAnalyzer.NEURO, c::Symbol)::Any
    # get component values
    @assert c in keys(obj.components) "Component $c not found."
    return obj.components[c]
end

function _select_cidx(obj::NeuroAnalyzer.NEURO, c::Symbol, c_idx::Union{Int64, Vector{Int64}, <:AbstractRange}, def_cidx::Int64=0)::Union{Int64, AbstractRange, Vector{Int64}}
    # select component channels, default is all or def_cidx
    c = _get_component(obj, c)
    if ndims(c) == 1
        c_idx == 0 && (c_idx = 1)
        _check_cidx(c, c_idx)
    else
        def_cidx > size(c, 1) && (def_cidx = size(c, 1))
        def_cidx == 0 && (def_cidx = size(c, 1))
        c_idx == 0 && (c_idx = 1:def_cidx)
        typeof(c_idx) <: AbstractRange && (c_idx = collect(c_idx))
        length(c_idx) > 1 && sort!(c_idx)
        for idx in c_idx
            @assert !(idx < 1 || idx > size(c, 1)) "c_idx must be in [1, $(size(c, 1))]."
        end
    end
    return c_idx
end

function _select_cidx(c::Union{AbstractVector, AbstractMatrix, AbstractArray}, c_idx::Union{Int64, Vector{Int64}, <:AbstractRange}, def_cidx::Int64=0)::Union{Int64, AbstractRange, Vector{Int64}}
    # select component channels, default is all or def_cidx
    if ndims(c) == 1
        c_idx == 0 && (c_idx = 1)
        _check_cidx(c, c_idx)
    else
        def_cidx > size(c, 1) && (def_cidx = size(c, 1))
        def_cidx == 0 && (def_cidx = size(c, 1))
        c_idx == 0 && (c_idx = 1:def_cidx)
        typeof(c_idx) <: AbstractRange && (c_idx = collect(c_idx))
        length(c_idx) > 1 && sort!(c_idx)
        for idx in c_idx
            @assert !(idx < 1 || idx > size(c, 1)) "c_idx must be in [1, $(size(c, 1))]."
        end
    end
    return c_idx
end