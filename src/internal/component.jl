function _get_component(obj::NeuroAnalyzer.NEURO, c::Symbol)
    # get component values
    c in obj.header.component_names || throw(ArgumentError("Component $c not found."))
    c_idx = findfirst(isequal(c), obj.header.component_names)
    c = obj.components[c_idx]
    return (c=c, c_idx=c_idx)
end

function _select_cidx(obj::NeuroAnalyzer.NEURO, c::Symbol, c_idx::Union{Int64, Vector{Int64}, AbstractRange}, def_cidx::Int64=0)
    # select component channels, default is all or def_cidx
    c, _ = _get_component(obj, c)
    def_cidx > size(c, 1) && (def_cidx = size(c, 1))
    def_cidx == 0 && (def_cidx = size(c, 1))
    c_idx == 0 && (c_idx = 1:def_cidx)
    typeof(c_idx) <: AbstractRange && (c_idx = collect(c_idx))
    length(c_idx) > 1 && sort!(c_idx)
    for idx in c_idx
        (idx < 1 || idx > size(c, 1)) && throw(ArgumentError("c_idx must be in [1, $(size(c, 1))]."))
    end
    return c_idx
end

function _select_cidx(c::AbstractArray, c_idx::Union{Int64, Vector{Int64}, AbstractRange}, def_cidx::Int64=0)
    # select component channels, default is all or def_cidx
    def_cidx > size(c, 1) && (def_cidx = size(c, 1))
    def_cidx == 0 && (def_cidx = size(c, 1))
    c_idx == 0 && (c_idx = 1:def_cidx)
    typeof(c_idx) <: AbstractRange && (c_idx = collect(c_idx))
    length(c_idx) > 1 && sort!(c_idx)
    for idx in c_idx
        (idx < 1 || idx > size(c, 1)) && throw(ArgumentError("c_idx must be in [1, $(size(c, 1))]."))
    end
    return c_idx
end