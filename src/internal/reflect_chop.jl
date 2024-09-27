_reflect(s::AbstractVector)::AbstractVector = vcat(s[end:-1:1], s, s[end:-1:1])
_reflect(s1::AbstractVector, s2::AbstractVector, s3::AbstractVector)::AbstractVector = vcat(s2[end:-1:1], s1, s3[end:-1:1])

_chop(s::AbstractVector)::AbstractVector = s[(length(s) ÷ 3 + 1):(length(s) ÷ 3) * 2]
_chop(s1::AbstractVector, n::Int64)::AbstractVector = s1[(n + 1):end - n]
