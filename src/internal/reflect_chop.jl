_reflect(signal::AbstractArray) = vcat(signal[end:-1:1], signal, signal[end:-1:1])
_reflect(s1::AbstractArray, s2::AbstractArray, s3::AbstractArray) = vcat(s2[end:-1:1], s1, s3[end:-1:1])

_chop(signal::AbstractArray) = signal[(length(signal) ÷ 3 + 1):(length(signal) ÷ 3) * 2]
_chop(s1::AbstractArray, n::Int64) = s1[(n + 1):end - n]
