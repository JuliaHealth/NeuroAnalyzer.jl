to_mat(A::AbstractMatrix) = A
function to_mat(y::AbstractVector)
    ymat = zeros(eltype(y), length(y)-1, 2);
    ymat[:,1] .= @view y[1:end-1];
    ymat[:,2] .= @view y[2:end]
    return ymat
end
function to_mat(pairs::Vector{Pair{T1,T2}} where {T1<:Tuple{<:Real,<:Real}, T2<:Union{Real, Tuple{<:Real,<:Real}}})
    xmat, ymat = zeros(length(pairs), 2), zeros(length(pairs), 2)
    for i = 1:length(pairs)
        xmat[i,:] .= pairs[i].first
        ymat[i,:] .= pairs[i].second
    end
    return xmat, ymat
end

to_toeplitz(vals) = to_toeplitz(vals, vals)
function to_toeplitz(vals_left::AbstractVector{T1}, vals_top::AbstractVector{T2}) where {T1,T2}
    @assert vals_left[1] == vals_top[1]
    N_rows, N_cols = length(vals_left), length(vals_top)
    A = zeros(promote_type(T1,T2), N_rows, N_cols)
    vals_total = zeros(promote_type(T1,T2), N_rows+N_cols-1)
    vals_total[1:N_cols] .= @view(vals_top[end:-1:1])
    vals_total[N_rows+1:end] .= @view(vals_left[2:end])
    for j = 1:N_cols
        A[:,j] .= @view(vals_total[end-N_rows+2-j:end+1-j])
    end
    return A
end

to_hankel(vals) = to_hankel(vals, vals)
function to_hankel(vals_left::AbstractVector{T1}, vals_bottom::AbstractVector{T2}) where {T1,T2}
    @assert vals_left[end] == vals_bottom[1]
    N_rows, N_cols = length(vals_left), length(vals_bottom)
    A = zeros(promote_type(T1,T2), N_rows, N_cols)
    vals_total = zeros(promote_type(T1,T2), N_rows+N_cols-1)
    vals_total[1:N_rows] .= vals_left
    vals_total[N_rows+1:end] .= @view(vals_bottom[2:end])
    for j = 1:N_cols
        A[:,j] .= @view(vals_total[j:j+N_rows-1])
    end
    return A
end

function nan2zero!(x::Array{T}) where T
    for idx = 1:length(x)
        if isnan(x[idx])
            x[idx] = zero(T)
        end
    end
end

function freqz(h::Vector{T}; fs::Real = 1, N = 1000) where T
    ω = range(0, stop = π, length = N)
    H = zeros(Complex{T}, length(ω))
    for (idx,h_n) in enumerate(h)
        n = idx-1
        H .+= h_n * exp.(-n*1im*ω)
    end
    return H, ω * (fs/2π)
end