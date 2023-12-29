_info(s::String) = verbose == true && @info s

_warn(s::String) = verbose == true && @warn s

_deprecated(s::String) = verbose == true && @error "Function $s() is deprecated."
_deprecated(s1::String, s2::String) = verbose == true && @error "Function $s1() is deprecated, please use $s2() instead."

_wip() = allow_wip == true ? (@warn "This function has the WIP (Work In Progress) status and is not ready for production use.") : (@error "This function has the WIP (Work In Progress) status and is not ready for production use.")

_pl(x::Union{AbstractRange, AbstractVector}) = length(collect(x)) > 1 ? "s" : ""

_pl(x::Real) = x > 1 ? "s" : ""

# _get_range(signal::Union{AbstractVector, AbstractArray}) = round(abs(minimum(signal)) + abs(maximum(signal)), digits=0)
_get_range(signal::Union{AbstractVector, AbstractArray}) = round(rng(signal), digits=0)

_c(n) = collect(1:n)

_tuple_max(t::Tuple{Real, Real}) = abs(t[1]) > abs(t[2]) ? (-abs(t[1]), abs(t[1])) : (-abs(t[2]), abs(t[2]))

_s2v(s::Union{<:Number, Vector{<:Number}}) = typeof(s) <: Number ? [s] : s

function _v2s(x::Vector{<:Number})
    s_tmp = string.(x)
    s = ""
    if length(s_tmp) > 1
        for idx in 1:(length(s_tmp) - 1)
            s *= s_tmp[idx] * ", "
        end
    end
    s *= s_tmp[end]
    return s
end

function _copy_lt2ut(m::AbstractArray)
    if ndims(m) == 2
        return m + m' - diagm(diag(m))
    else
        Threads.@threads for ep_idx in 1:size(m, 3)
            @inbounds m[:, :, ep_idx] = m[:, :, ep_idx] + m[:, :, ep_idx]' - diagm(diag(m[:, :, ep_idx]))
        end
        return m
    end
end

_tlength(t::Tuple{Real, Real}) = length(t[1]:1:t[2])

function _s2i(s::String)
    s = replace(s, " "=>"")
    if occursin(":", s)
        return collect(parse(Int64, split(s, ":")[1]):parse(Int64, split(s, ":")[2]))
    elseif occursin(",", s)
        s = replace(s, "["=>"")
        s = replace(s, "]"=>"")
        return parse.(Int64, split(s, ","))
    elseif _check_sint(s)
        return parse.(Int64, s)
    end
end

function _i2s(s::Union{Int64, Vector{Int64}, AbstractRange})
    !isa(s, Int64) && (s = collect(s))
    s = string(s)
    s = replace(s, "["=>"")
    s = replace(s, "]"=>"")
    return s
end

function _s2tf(s::String)
    s = replace(s, " "=>"")
    s = replace(s, "("=>"")
    s = replace(s, ")"=>"")
    return (parse(Float64, split(s, ",")[1]), parse(Float64, split(s, ",")[2]))
end

function _s2ti(s::String)
    s = replace(s, " "=>"")
    s = replace(s, "("=>"")
    s = replace(s, ")"=>"")
    return (parse(Int64, split(s, ",")[1]), parse(Int64, split(s, ",")[2]))
end

function _detect_montage(clabels::Vector{String}, ch_type::Vector{String}, data_type::String)
    m = match.(r"(.+)\-(.+)", lowercase.(clabels[ch_type .== data_type]))
    if length(findall(!isnothing, m)) == length(clabels[ch_type .== data_type])
        r = String[]
        for idx in 1:length(m)
            push!(r, m[idx].captures[2])
        end
        if length(unique(r)) == 1
            occursin("a", lowercase(r[1])) && return "common (A)"
            occursin("m", lowercase(r[1])) && return "common (M)"
            return "common"
        else
            return "bipolar"
        end
    end
    m = match.(r"([a-z]+)([0-9]+[0-9]?)([a-z]+)([0-9]+)", lowercase.(clabels[ch_type .== data_type]))
    if length(findall(!isnothing, m)) == length(clabels[ch_type .== data_type])
        r = String[]
        for idx in 1:length(m)
            push!(r, m[idx].captures[3])
        end
        if length(unique(r)) == 1
            occursin("a", lowercase(r[1])) && return "common (A)"
            occursin("m", lowercase(r[1])) && return "common (M)"
            return "common"
        else
            return "bipolar"
        end
    else
        return "physical"
    end
end

function _fread(fid, n, t)
    (n > 1 && t === :c) && (t = :s)
    t === :s && (n = n)
    t === :c && (n = n)
    t === :l && (n *= 4)
    t === :ul && (n *= 4)
    t === :ui8 && (n = n)
    t === :ui16 && (n *= 2)
    t === :ui32 && (n *= 4)
    t === :ui64 && (n *= 8)
    t === :i && (n *= 4)
    t === :i8 && (n = n)
    t === :i16 && (n *= 2)
    t === :i32 && (n *= 4)
    t === :i64 && (n *= 8)
    t === :f16 && (n *= 2)
    t === :f32 && (n *= 4)
    t === :f64 && (n *= 8)
    header = zeros(UInt8, n)
    readbytes!(fid, header, n)
    t === :s && return Int64.(map(ltoh, reinterpret(UInt8, header)))
    t === :c && return Int64(map(ltoh, reinterpret(UInt8, header))[1])
    t === :l && return Int64(map(ltoh, reinterpret(Int32, header))[1])
    t === :ul && return Int64(map(ltoh, reinterpret(UInt32, header))[1])
    t === :ui8 && return Int64(map(ltoh, reinterpret(UInt8, header))[1])
    t === :ui16 && return Int64(map(ltoh, reinterpret(UInt16, header))[1])
    t === :ui32 && return Int64(map(ltoh, reinterpret(UInt32, header))[1])
    t === :ui64 && return Int64(map(ltoh, reinterpret(UInt64, header))[1])
    t === :i && return Int64(map(ltoh, reinterpret(Int32, header))[1])
    t === :i16 && return Int64(map(ltoh, reinterpret(Int16, header))[1])
    t === :i32 && return Int64(map(ltoh, reinterpret(Int32, header))[1])
    t === :i64 && return Int64(map(ltoh, reinterpret(Int64, header))[1])
    t === :f16 && return Float64(map(ltoh, reinterpret(Float16, header))[1])
    t === :f32 && return Float64(map(ltoh, reinterpret(Float32, header))[1])
    t === :f64 && return Float64(map(ltoh, reinterpret(Float64, header))[1])
end

function _vint2str(x::Vector{Int64})
    s = strip(String(Char.(x)))
    return replace(s, "\0"=>"")
end
