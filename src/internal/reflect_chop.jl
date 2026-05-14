"""Reflect the signal: 123 becomes 321123321."""
function _reflect(s::AbstractVector)::AbstractVector
    isempty(s) && throw(ArgumentError("s cannot be empty."))
    return vcat(s[end:-1:1], s, s[end:-1:1])
end

"""Reflect the signal: s2(reversed):s1:s3(reversed)."""
function _reflect(
    s1::AbstractVector,
    s2::AbstractVector,
    s3::AbstractVector,
)::AbstractVector
    isempty(s1) && throw(ArgumentError("s1 cannot be empty."))
    isempty(s2) && throw(ArgumentError("s2 cannot be empty."))
    isempty(s3) && throw(ArgumentError("s3 cannot be empty."))
    return vcat(s2[end:-1:1], s1, s3[end:-1:1])
end

"""Chop the reflected signal: 321123321 becomes 123."""
function _chop(s::AbstractVector)::AbstractVector
    isempty(s) && throw(ArgumentError("s cannot be empty."))
    return s[(length(s) ÷ 3 + 1):((length(s) ÷ 3) * 2)]
end

"""Chop `n` samples from the beginning and the end of a reflected signal."""
function _chop(s::AbstractVector, n::Int64)::AbstractVector
    isempty(s) && throw(ArgumentError("s cannot be empty."))
    n >= 0 || throw(ArgumentError("n must be ≥ 0."))
    if n == 0
        return s
    else
        s[(n + 1):(end - n)]
    end
end
