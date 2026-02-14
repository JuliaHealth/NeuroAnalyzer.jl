function _fir_response(f::Vector{<:Real}, w = range(0, stop = π, length = 1024))::Array{ComplexF32}
    # code based on Matti Pastell "FIR filter design with Julia"
    n = length(w)
    h = Array{ComplexF32}(undef, n)
    sw = 0
    for i in 1:n
        for j in eachindex(f)
            sw += f[j] * exp(-im * w[i])^-j
        end
        h[i] = sw
        sw = 0
    end
    return h
end
