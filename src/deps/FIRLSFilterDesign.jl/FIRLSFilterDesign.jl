module FIRLSFilterDesign
"""
https://cnx.org/contents/6x7LNQOp@7/Linear-Phase-Fir-Filter-Design-By-Least-Squares
https://www.dsprelated.com/Plishowarticle/808.php
https://eeweb.engineering.nyu.edu/iselesni/EL713/zoom/linphase.pdf
"""
# using LinearAlgebra
import LinearAlgebra: I, Diagonal, UniformScaling

include("utils.jl")

export firls_design

abstract type FIR end
struct FIR_I <: FIR end
struct FIR_II <: FIR end
struct FIR_III <: FIR end
struct FIR_IV <: FIR end

@doc """
    infer_fir_type(is_odd, is_antisymmetric)

Determines the type of FIR filter to be designed, based on:
* Whether the number of filter coefficients is odd (`is_odd`)
* Whether the filter should be antisymmetric (`is_antisymmetric`)

The result is a `fir_type`, which can be:
1. `FIR_I`, when filter length is odd and the filter is not antisymmetric
2. `FIR_II`, when filter length is even and the filter is not antisymmetric
3. `FIR_III`, when filter length is odd and the filter is antisymmetric
4. `FIR_IV`, when filter length is even and the filter is antisymmetric

# Arguments
- `is_odd::Bool`
- `is_antisymmetric::Bool`

# Outputs
- `fir_type::FIR`
""" 
function infer_fir_type(is_odd, is_antisymmetric)
    if !is_antisymmetric
        fir_type = is_odd ? FIR_I : FIR_II
    else
        fir_type = is_odd ? FIR_III : FIR_IV
    end
    fir_type()
end

function validate_inputs(filter_order, bands_D, D, fs)
    @assert filter_order >= 0 "Filter order cannot be negative."
    @assert fs > 0 "Sampling frequency should be larger than 0 Hz"
    validate_bands_D(bands_D, D, fs)
end
function validate_inputs(filter_order, bands_DW, D, W, fs)
    @assert filter_order >= 0 "Filter order cannot be negative."
    @assert fs > 0 "Sampling frequency should be larger than 0 Hz"
    validate_bands_D(bands_DW, D, fs)
    validate_bands_W(bands_DW, W, fs)
end
function validate_inputs(filter_order, bands_D, D, bands_W, W, fs)
    @assert filter_order >= 0 "Filter order cannot be negative."
    @assert fs > 0 "Sampling frequency should be larger than 0 Hz"
    validate_bands_D(bands_D, D, fs)
    validate_bands_W(bands_W, W, fs)
end

function validate_bands_D(bands_D, D, fs)
    validate_bands(bands_D, fs)
    @assert size(bands_D,2) == size(D,2) == 2 "Frequency bands and desired response should be N x 2 matrices"
    @assert size(bands_D,1) == size(D,1) "Frequency bands and desired response should be N x 2 matrices"
end

function validate_bands_W(bands_W, W, fs)
    validate_bands(bands_W, fs)
    @assert size(bands_W,2) == size(W,2) == 2 "Frequency band and weight matrix should be N x 2."
    @assert size(bands_W,1) == size(W,1) "Number of frequency bands should match number of weights."
end

function validate_bands(fbands::Matrix{T}, fs) where T
    @assert fbands[1] == 0 "Frequency bands should start at 0."
    @assert fbands[end] == fs * (1//2) "Frequency bands should end at fs/2"
    @views @assert all((fbands[2:end,1] .- fbands[1:end-1,2]) .== zero(T)) "Frequency bands should cover entire [0,fs/2] interval, without gaps or overlaps."
end

@doc """
    get_flength_M(filter_order)

Determines the length of the FIR filter and the number of amplitude coefficients needed, based on its order.

# Arguments
- `filter_order::Integer` : the order of the FIR filter.

# Outputs
- `filter_length::Integer` : the number of filter coefficients
- `M::Integer` : the number of unique amplitude coefficients needed to form the filter is equal to ``M+1``
"""
function get_flength_M(filter_order)
    filter_length = filter_order + 1
    M = isodd(filter_length) ? (filter_length-1)÷2 : filter_length÷2-1
    return filter_length, M
end


to_matrix_simple(v::Vector) = hcat(v,v)
to_matrix_simple(m::Matrix) = m
knotpoints_to_matrix(v::Vector) = @views hcat(v[1:end-1], v[2:end])

@doc """
    firls_design(filter_order::Integer, bands_DW::Matrix, D::Matrix, W::Matrix, antisymmetric::Bool; fs::Real = 1, solver::Function = \\)

Designs a linear-phase FIR filter.

# Arguments
- `filter_order::Integer`   : the order of the FIR filter.
- `bands_DW::Matrix`        : a matrix of size (N,2) which contains rows of sequential frequency bands, spanning [0, fs/2].
- `D::Matrix`               : a matrix of size (N,2) which contains rows of amplitude responses for the frequency bands in `bands_DW`. The first and second columns indicate the amplitude response at the lower and upper bound of the frequency bands, interpolated linearly in between.
- `W::Matrix`               : a matrix of size (N,2) which contains rows of weighting coefficients for the frequency bands in `bands_DW`. The first and second columns indicate the weighting at the lower and upper bound of the frequency bands, interpolated linearly in between.
- `antisymmetric::Bool`     : a Boolean that signifies whether the filter coefficients will be anti-symmetric, as used in type III and IV FIR filters.
- `fs::Real`                : the sampling frequency.
- `solver::Function`        : the function that is called to solve the equation ``Qa = b``, with the function call: `solver(Q,b)` which returns `a`.

# Outputs
- `h` : a vector of linear-phase FIR filter coefficients.
"""
function firls_design(filter_order::Integer, bands_DW::Matrix, D::Matrix, W::Matrix, antisymmetric::Bool; fs::Real = 1, solver::Function = \)
    validate_inputs(filter_order, bands_DW, D, W, fs)
    filter_length, M = get_flength_M(filter_order)
    fir_type = infer_fir_type(isodd(filter_length), antisymmetric)
    a = calc_amplitude_coeff(M, bands_DW, D, W, solver, fir_type)
    h = _to_impulse_response(a, fir_type)
end
@doc """
    firls_design(filter_order::Integer, bands_DW::Matrix, D::Union{Vector,Matrix}, W::Union{Vector,Matrix}, antisymmetric::Bool; fs::Real = 1, solver::Function = \\)

# Arguments
- `filter_order::Integer`   : the order of the FIR filter.
- `bands_DW::Matrix`        : a matrix of size `(N,2)` which contains rows of sequential frequency bands, spanning [0, fs/2].
- `D::Union{Vector,Matrix}` : a matrix of size `(N,2)`, or a vector of size `(N,)`, which amplitude responses for the frequency bands in `bands_DW`. In the case of a matrix, the first and second columns indicate the amplitude response at the lower and upper bound of the frequency bands, interpolated linearly in between. In the case of a vector the elements the amplitude response is piecewise constant.
- `W::Union{Vector,Matrix}` : a matrix of size `(N,2)`, or a vector of size `(N,)`, which contains weighting function values for the frequency bands in `bands_DW`. In the case of a matrix, the first and second columns indicate the weighting function values at the lower and upper bound of the frequency bands, interpolated linearly in between. In the case of a vector the elements the weighting function is piecewise constant.
- `antisymmetric::Bool`     : a Boolean that signifies whether the filter coefficients will be anti-symmetric, as used in type III and IV FIR filters.
- `fs::Real`                : the sampling frequency.
- `solver::Function`        : the function that is called to solve the equation ``Qa = b``, with the function call: `solver(Q,b)` which returns `a`.

# Outputs
- `h` : a vector of linear-phase FIR filter coefficients.
"""
function firls_design(filter_order::Integer, bands_DW::Matrix, D::Union{Vector,Matrix}, W::Union{Vector,Matrix}, antisymmetric::Bool; fs::Real = 1, solver::Function = \)
    firls_design(filter_order, bands_DW, to_matrix_simple(D), to_matrix_simple(W), antisymmetric; fs = fs, solver = solver)
end

@doc """
    firls_design(filter_order::Integer, knotpoints_DW::Vector, D::Vector, W::Vector, antisymmetric::Bool; fs::Real = 1, solver::Function = \\)

# Arguments
- `filter_order::Integer`   : the order of the FIR filter.
- `knotpoints_DW::Vector`   : a vector of size `(N,)` which contains frequency knotpoints, spanning [0, fs/2].
- `D::Vector`               : a vector of size `(N,)` which contains amplitude response values for the frequency knotpoints in `knotpoints_DW`. 
- `W::Vector`               : a vector of size `(N,)` which contains weighting function values for the frequency knotpoints in `knotpoints_DW`. 
- `antisymmetric::Bool`     : a Boolean that signifies whether the filter coefficients will be anti-symmetric, as used in type III and IV FIR filters.
- `fs::Real`                : the sampling frequency.
- `solver::Function`        : the function that is called to solve the equation ``Qa = b``, with the function call: `solver(Q,b)` which returns `a`.

# Outputs
- `h` : a vector of linear-phase FIR filter coefficients.
"""
function firls_design(filter_order::Integer, knotpoints_DW::Vector, D::Vector, W::Vector, antisymmetric::Bool; fs::Real = 1, solver::Function = \)
    firls_design(filter_order, knotpoints_to_matrix(knotpoints_DW), knotpoints_to_matrix(D), knotpoints_to_matrix(W), antisymmetric; fs = fs, solver = solver)
end

@doc """
    firls_design(filter_order::Integer, bands_DW::Matrix, D::Union{Vector,Matrix}, antisymmetric::Bool; fs::Real = 1, solver::Function = \\)

# Arguments
- `filter_order::Integer`   : the order of the FIR filter.
- `bands_DW::Matrix`        : a matrix of size `(N,2)` which contains rows of sequential frequency bands, spanning [0, fs/2].
- `D::Union{Vector,Matrix}` : a matrix of size `(N,2)`, or a vector of size `(N,)`, which amplitude responses for the frequency bands in `bands_DW`. In the case of a matrix, the first and second columns indicate the amplitude response at the lower and upper bound of the frequency bands, interpolated linearly in between. In the case of a vector the elements the amplitude response is piecewise constant.
- `antisymmetric::Bool`     : a Boolean that signifies whether the filter coefficients will be anti-symmetric, as used in type III and IV FIR filters.
- `fs::Real`                : the sampling frequency.
- `solver::Function`        : the function that is called to solve the equation ``Qa = b``, with the function call: `solver(Q,b)` which returns `a`.

# Outputs
- `h` : a vector of linear-phase FIR filter coefficients.
"""
function firls_design(filter_order::Integer, bands_DW::Matrix, D::Union{Vector,Matrix}, antisymmetric::Bool; fs::Real = 1, solver::Function = \)
    D = to_matrix_simple(D)
    validate_inputs(filter_order, bands_DW, D, fs)
    filter_length, M = get_flength_M(filter_order)
    fir_type = infer_fir_type(isodd(filter_length), antisymmetric)
    a = calc_amplitude_coeff(M, bands_DW, D, solver, fir_type)
    h = _to_impulse_response(a, fir_type)
end

@doc """
    firls_design(filter_order::Integer, knotpoints_D::Vector, D::Vector, antisymmetric::Bool; fs::Real = 1, solver::Function = \\)

# Arguments
- `filter_order::Integer`   : the order of the FIR filter.
- `knotpoints_D::Vector`    : a vector of size `(N,)` which contains frequency knotpoints, spanning [0, fs/2].
- `D::Vector`               : a vector of size `(N,)` which contains amplitude response values for the frequency knotpoints in `knotpoints_D::Vector`.
- `antisymmetric::Bool`     : a Boolean that signifies whether the filter coefficients will be anti-symmetric, as used in type III and IV FIR filters.
- `fs::Real`                : the sampling frequency.
- `solver::Function`        : the function that is called to solve the equation ``Qa = b``, with the function call: `solver(Q,b)` which returns `a`.

# Outputs
- `h` : a vector of linear-phase FIR filter coefficients.
"""
function firls_design(filter_order::Integer, knotpoints_D::Vector, D::Vector, antisymmetric::Bool; fs::Real = 1, solver::Function = \)
    firls_design(filter_order, knotpoints_to_matrix(knotpoints_D), knotpoints_to_matrix(D), antisymmetric; fs = fs, solver = solver)
end

function calc_amplitude_coeff(M, bands_DW, D, W, solver, fir_type)
    Q = get_Q(M, bands_DW, W, fir_type)
    b = get_b(M, bands_DW, D, W, fir_type)
    a = solve(Q, b, solver, fir_type)
end
function calc_amplitude_coeff(M, bands_DW, D, solver, fir_type)
    Q = get_Q(M, fir_type)
    b = get_b(M, bands_DW, D, ones(size(D)), fir_type)
    a = solve(Q, b, solver, fir_type)
end

solve(Q, b, solver, fir_type) = solver(Q, b)
function solve(Q, b, solver, fir_type::FIR_III)
    a = zeros(length(b))
    if length(a) > 1
        @views a[2:end] .= solver(Q[2:end,2:end], b[2:end])
    end
    return a
end
function solve(Q::UniformScaling, b, solver, fir_type::FIR_III)
    a = zeros(length(b))
    if length(a) > 1
        @views a[2:end] .= solver(Q, b[2:end])
    end
    return a
end

"""
    get_Q(M, f, W, fir_type)

Constructs the matrix ``Q`` used in the equation ``Qa = b``, based on a set of weights.

# Arguments
- `M::Integer`      : indicator of the amount of elements needed.
- `f::Matrix`       : a matrix of size `(N,2)` which contains rows of sequential frequency bands, spanning [0, fs/2].
- `W::Matrix`       : a matrix of size `(N,2)` which contains rows of weighting coefficients for the frequency bands in `f`. The first and second columns indicate the weighting at the lower and upper bound of the frequency bands, interpolated linearly in between.
- `fir_type::FIR`   : indicates the type of FIR filter.

# Outputs
- `Q::Matrix` : the matrix ``Q`` used in the equation ``Qa = b``.
"""
function get_Q(M, f, W, fir_type)
    q = get_q(M, f, W, fir_type)
    Q1, Q2 = q_to_Q1Q2(q, M, fir_type)
    Q = Q1Q2_to_Q(Q1, Q2, fir_type)
end
"""
    get_Q(M, fir_type)

Constructs the matrix ``Q`` used in the equation ``Qa = b``, when there are no weights. Which results in ``Q`` being the identity matrix.

# Arguments
- `M::Integer`      : indicator of the amount of elements needed.
- `fir_type::FIR`   : indicates the type of FIR filter.

# Outputs
- `Q::Matrix` : the matrix ``Q`` used in the equation ``Qa = b``.
"""
get_Q(M, fir_type) = I
"""
    get_Q(M, fir_type::FIR_I)

Constructs the matrix ``Q`` used in the equation ``Qa = b``, when there are no weights and the FIR filter is of type I.

# Arguments
- `M::Integer`      : indicator of the amount of elements needed.
- `fir_type::FIR_I`   : indicates the type of FIR filter is I.

# Outputs
- `Q::Matrix` : the matrix ``Q`` used in the equation ``Qa = b``.
"""
function get_Q(M, fir_type::FIR_I)
    v_diag = fill(1., M+1)
    v_diag[1] *= 2
    Q = Diagonal(v_diag)
end

function q_to_Q1Q2(q, M, fir_type::Union{FIR_I,FIR_III})
    Q1 = @views to_toeplitz(q[1:M+1])
    Q2 = @views to_hankel(q[1:M+1], q[M+1:end])
    return Q1, Q2
end
function q_to_Q1Q2(q, M, fir_type::Union{FIR_II,FIR_IV})
    Q1 = @views to_toeplitz(q[1:M+1])
    Q2 = @views to_hankel(q[1+1:M+1+1], q[M+1+1:end])
    return Q1, Q2
end

Q1Q2_to_Q(Q1, Q2, fir_type::Union{FIR_I,FIR_II}) = Q1 .+ Q2
Q1Q2_to_Q(Q1, Q2, fir_type::Union{FIR_III,FIR_IV}) = Q1 .- Q2

"""
    get_q(M, f, W, fir_type)

Finds the vector ``q`` which is used to populate the matrix ``Q``.
...
# Arguments
- `M::Integer`      : indicator of the amount of elements needed.
- `f::Matrix`       : a matrix of size `(N,2)` which contains rows of sequential frequency bands, spanning [0, fs/2].
- `W::Matrix`       : a matrix of size `(N,2)` which contains rows of weighting coefficients for the frequency bands in `f`. The first and second columns indicate the weighting at the lower and upper bound of the frequency bands, interpolated linearly in between.
- `fir_type::FIR`   : indicates the type of FIR filter.

# Outputs
- `q_out::Vector`   : a vector of q-values that are used to fill in the Q-matrix.
...
"""
function get_q(M, f, W, fir_type)
    a, b, α, β, γ, k = constants_q(f, W)
    _αn, _βn⁻¹, _γ⁻² = copy(α), copy(β), copy(γ)
    q_out, _qn = allocate_q(M, fir_type), zeros(size(f))
    q_out[1] = qn!(_qn, k, f, a, b)
    for idx = 2:length(q_out)
        _αn .= α; _βn⁻¹ .= β; _γ⁻² .= γ
        n = idx2n_q(idx, fir_type)
        q_out[idx] = qn!(_qn, n, k, _αn, _βn⁻¹, _γ⁻²)
    end
    return q_out
end

function constants_q(f, W)
    fs = 2f[end,2]; k = 2/fs
    a = @views @. (W[:,2] - W[:,1]) / (f[:,2] - f[:,1])
    b = @views @. W[:,1] - f[:,1]*a
    α = @. π*k*f
    β = @. (1/(π*k)) * (a*f + b)
    γ = @. (1/(π*k)^2) * a
    return a, b, α, β, γ, k
end

allocate_q(M, fir_type::Union{FIR_I,FIR_III}) = q_out = zeros(2M+1)
allocate_q(M, fir_type::Union{FIR_II,FIR_IV}) = q_out = zeros(2M+1+1)

idx2n_q(idx, fir_type::Union{FIR_I,FIR_III}) = idx-1.
idx2n_q(idx, fir_type::Union{FIR_II,FIR_IV}) = idx-1.

function qn!(_qn, k, f, a, b)
    # Fallback for when n = 0
    # _qn = (af²/2 + bf)
    #_qn = 1/2 * (af² + 2bf)
    # bn = k * (_qn[:,2] - _qn[:,1]) 
    @. _qn = b*f
    @. _qn *= 2
    @. _qn += a*f^2
    nan2zero!(_qn)
    bn = @views k/2 * (sum(_qn[:,2]) - sum(_qn[:,1]))
end

function qn!(_qn, n, k, _αn, _βn⁻¹, _γn⁻²)
    _αn .*= n
    @. _βn⁻¹ *= 1/n
    @. _γn⁻² *= 1/n^2
    @. _qn = _βn⁻¹ * sin(_αn) + _γn⁻² * cos(_αn)
    nan2zero!(_qn)
    bn = @views k * (sum(_qn[:,2]) - sum(_qn[:,1]))
end

"""
    get_b(M, f, D, W, fir_type)

Finds the vector ``b`` used in the equation ``Qa = b``.
...
# Arguments
- `M::Integer`      : Size of the `b` vector is `M+1` .
- `f::Matrix`       : a matrix of size `(N,2)` which contains rows of sequential frequency bands, spanning the interval [0, fs/2].
- `D::Matrix`       : a matrix of size `(N,2)` which contains rows of desired frequency response values for the frequency bands in `f`. The first and second columns indicate the desired response at the lower and upper bound of the frequency bands, interpolated linearly in between.
- `W::Matrix`       : a matrix of size `(N,2)` which contains rows of weighting coefficients for the frequency bands in `f`. The first and second columns indicate the weighting at the lower and upper bound of the frequency bands, interpolated linearly in between.
- `fir_type::FIR`   : indicates the type of FIR filter.

# Outputs
- `b_out::Vector`   : a vector of size `(M+1,)`, the b-vector used in the equation ``Qa = b``.
...
"""
function get_b(M, f, D, W, fir_type)
    a, b, c, d, α, β, γ, δ, k = constants_b(f, D, W)
    _αn, _βn², _δn = copy(α), copy(β), copy(δ)
    b_out, _bn = zeros(M+1), zeros(size(f))
    b_out[1] = bn_n0!(_bn, idx2n_b(1, fir_type), k, f, a, b, c, d, _αn, _βn², γ, _δn, fir_type)
    for idx = 2:length(b_out)
        _αn .= α; _βn² .= β; _δn .= δ
        n = idx2n_b(idx, fir_type)
        b_out[idx] = bn!(_bn, n, k, _αn, _βn², γ, _δn, fir_type)
    end
    return b_out
end

@doc raw"""
    constants_b(f, D, W)

Calculates data that is reused at every evaluation of [`bn!`](@ref).

...
# Arguments
- `f::Matrix` : a matrix of size `(N,2)` which contains rows of sequential frequency bands, spanning the interval [0, fs/2].
- `D::Matrix` : a matrix of size `(N,2)` which contains rows of desired frequency response values for the frequency bands in `f`. The first and second columns indicate the desired response at the lower and upper bound of the frequency bands, interpolated linearly in between.
- `W::Matrix` : a matrix of size `(N,2)` which contains rows of weighting coefficients for the frequency bands in `f`. The first and second columns indicate the weighting at the lower and upper bound of the frequency bands, interpolated linearly in between.
...
"""
function constants_b(f, D, W)
    fs = 2f[end,2]; k = 2/fs
    a = @views @. (D[:,2] - D[:,1]) / (f[:,2] - f[:,1])
    b = @views @. D[:,1] - f[:,1]*a
    c = @views @. (W[:,2] - W[:,1]) / (f[:,2] - f[:,1])
    d = @views @. W[:,1] - f[:,1]*c
    α = @. k*π*f
    β = @. (k*π)^2 * (a*c*f^2 + (a*d+b*c)*f + b*d)
    γ = @. -2a*c
    δ = @. k*π * (2a*c*f + a*d + b*c) 
    return a, b, c, d, α, β, γ, δ, k
end

@doc raw"""
    bn_n0!(_bn, n, k, f, a, b, c, d, _αn, _βn², γ, _δn, fir_type)

Dispatches to the correct function to calculate the first element of the `b` vector, based on the type of FIR filter.
Needed because for type I FIR filters the value of `n` at the first iteration is `0`.
"""
function bn_n0!(_bn, n, k, f, a, b, c, d, _αn, _βn², γ, _δn, fir_type)
    bn!(_bn, n, k, _αn, _βn², γ, _δn, fir_type)
end
function bn_n0!(_bn, n, k, f, a, b, c, d, _αn, _βn², γ, _δn, fir_type::FIR_I)
    bn!(_bn, k, f, a, b, c, d)
end

@doc raw"""
    idx2n_b(idx, fir_type::Union{FIR_I,FIR_III})

Determine the value of ``n`` based on the index in the b-vector.
For type I and III FIR filters the following holds: ``n = i - 1``.
Based on page 10 and 12 of [this](https://eeweb.engineering.nyu.edu/iselesni/EL713/zoom/linphase.pdf).

...
# Arguments
- `i::Integer` : index in the b-vector.
- `fir_type::Union{FIR_I,FIR_III}`
# Outputs
- `n::Real`
...
"""
idx2n_b(i, fir_type::Union{FIR_I,FIR_III}) = n = i-1.
@doc raw"""
    idx2n_b(idx, fir_type::Union{FIR_II,FIR_IV})

Determine the value of ``n`` based on the index in the b-vector, Where ``n`` is used in 
For type II and IV FIR filters the following holds: ``n = i - 1/2``.
Based on page 11 and 13 of [this](https://eeweb.engineering.nyu.edu/iselesni/EL713/zoom/linphase.pdf).

...
# Arguments
- `i::Integer` : index in the b-vector.
- `fir_type::Union{FIR_II,FIR_IV}`

# Outputs
- `n::Real`
...
"""
idx2n_b(i, fir_type::Union{FIR_II,FIR_IV}) = n = i-0.5

@doc raw"""
    bn!(_bn, n, k, _αn, _βn², γ, _δn, fir_type)

Calculates the elements of the b-vector, which are equal to: 
```math 
b[i] = \frac{2}{f_s} \int_0^{f_s/2} W(f) D(f) cos(\pi \frac{2}{f_s} n f) df, \quad i = 1, 2, \cdots, M+1
```

Both ``D(f)`` and ``W(f)`` are piecewise linear functions, and ``n`` is calculated by [`idx2n_b`](@ref). Using ``k = 2/f_s``, this integral becomes:
```math
b[i] = k \sum_{j=1}^J \int_{F_{j,1}}^{F_{j,2}} g_j(f,n) df = k \sum_{j=1}^J \int_{F_{j,1}}^{F_{j,2}} \big(c_j f+d_j\big) \big(a_j f+b_j\big) \cos(\pi k n f) df, \quad i = 1, 2, \cdots, M+1
```
Where:
* ``F_{j,1}`` and ``F_{j,2}`` are the lower and upper bound of the ``j^{th}`` frequency band.
* ``a_j`` and ``b_j`` are the parameters of the linear function that describes the desired frequency response in the ``j^{th}`` frequency band.
* ``c_j`` and ``d_j`` are the parameters of the linear function that describes the error weighting function in the ``j^{th}`` frequency band.

The antiderivative of ``g(f,n)`` is equal to:
```math 
G_j(f,n) = \frac{1}{\pi^3 k^3 n^3} \bigg(\sin\big(\alpha(f) n\big)\Big(\beta_j(f) n^2 + \gamma_j(f)\Big) + \delta_j(f) n \cos\big(\alpha(f) n\big)\bigg) + constants
```
Where:
*  ``\alpha(f) = \pi kf``, 
* ``\beta(f) = \pi^2 k^2 \big(acf^2 + (ad+bc)f + bd\big)``, 
* ``\gamma(f) = -2ac``, 
* ``\delta(f) = \pi k\big(2acf + ad + bc\big)``.
Note that the subscripts have been dropped here for clarity.

...
# Arguments
- `_bn::Vector`   : a vector of size `(J,)` that is used to store the intermediate values.
- `n::Integer`    : integer denoting the current cosine mode.
- `k::Real`       : equal to `2/fs`.
- `_αn::Vector`   : a vector of size `(J,)` that holds the values of ``\alpha n``.
- `_βn²::Vector`  : a vector of size `(J,)` that holds the values of ``\beta n^2``.
- `γ::Vector`     : a vector of size `(J,)` that holds the values of ``\gamma``.
- `_δn::Vector`   : a vector of size `(J,)` that holds the values of ``\delta n``.
- `fir_type::FIR` : indicates the type of FIR filter.

# Outputs
- `bn`            : ``n^{th}`` element in the ``b`` vector.
...
"""
function bn!(_bn, n, k, _αn, _βn², γ, _δn, fir_type)
    # if antisymmetric is false: _bn = 1/(π³k³n³) * (sin(αn) * (βn² + γ) + δncos(αn))
    # if antisymmetric is true:  _bn = 1/(π³k³n³) * (-cos(αn) * (βn² + γ) + δnsin(αn))
    #                 equal to:  _bn = 1/(π³k³n³) * (sin(αn-π/2) * (βn² + γ) + δncos(αn-π/2))
    # see _update_trig_arg_b
    #  bn = k * (_bn[:,2] - _bn[:,1])
    _update_trig_arg_b!(_αn, n, fir_type)
    @. _βn² *= n^2
    @. _δn *= n
    @. _bn = sin(_αn)*(_βn² + γ) + _δn*cos(_αn)
    nan2zero!(_bn)
    bn = @views k/(π*n*k)^3 * (sum(_bn[:,2]) - sum(_bn[:,1]))
end

@doc raw"""
    bn!(_bn, k, f, a, b, c, d)

Special case for when ``n = 0``, since then the integral is simplified:
```math 
g_j(f,0) = (c_j f+d_j) (a_j f+b_j) \cos(\pi k f 0) = (c_j f+d_j) (a_j f+b_j)
```
And the antiderivative becomes:
```math
G_j(f,0) = a_j c_j \frac{f^3}{3}  + (a_j d_j + b_j c_j)\frac{f^2}{2} + b_j d_j f
```

...
# Arguments
- `_bn::Vector` : a vector of size `(J,)` that is used to store the intermdediate values.
- `k::Real` : equal to ``2/f_s``.
- `f::Matrix` : a matrix of size `(J,2)` which contains rows of sequential frequency bands, spanning [0, fs/2].
- `a::Vector` : a vector of size `(J,)` with the ``a_j`` values in the equation ``a_j f + b_j`` that equates to the linear function that describes the desired frequency response in the ``j^{th}`` frequency band.
- `b::Vector` : a vector of size `(J,)` with the ``b_j`` values in the equation ``a_j f + b_j`` that equates to the linear function that describes the desired frequency response in the ``j^{th}`` frequency band.
- `c::Vector` : a vector of size `(J,)` with the ``c_j`` values in the equation ``c_j f + d_j`` that equates to the linear function that describes the error weighting function in the ``j^{th}`` frequency band.
- `d::Vector` : a vector of size `(J,)` with the ``d_j`` values in the equation ``c_j f + d_j`` that equates to the linear function that describes the error weighting function in the ``j^{th}`` frequency band.
...
"""
function bn!(_bn, k, f, a, b, c, d)
    # _bn = acf³/3 + (ad+bc)f²/2 + bdf
    # bn = k * (_bn[:,2] - _bn[:,1]) 
    @. _bn = b*d*f
    @. _bn *= 2
    @. _bn += (a*d+b*c)*f^2
    @. _bn *= 3/2
    @. _bn += a*c*f^3
    nan2zero!(_bn)
    bn = @views k/3 * (sum(_bn[:,2]) - sum(_bn[:,1]))
end

@doc raw"""
    _update_trig_arg_b!(_αn, n, fir_type::Union{FIR_I,FIR_II})

Updates the argument of the trigonometric functions in [`bn!`](@ref) by multiplying with the current `n`.

...
# Arguments
- `_αn::Vector` : a vector of size `(J,)`
- `n::Real`
- `fir_type{Union{FIR_I,FIR_II}}`
...
"""
_update_trig_arg_b!(_αn, n, fir_type::Union{FIR_I,FIR_II}) = _αn .*= n
@doc raw"""
    _update_trig_arg_b!(_αn, n, fir_type::Union{FIR_III,FIR_IV})
    
Updates the argument of the trigonometric functions in [`bn!`](@ref) by multiplying with ``n`` and subtracting ``\pi``.
The subtraction of ``\pi`` is necessary because when the filter is antisymmetric (type III and IV FIR filters), the filter response is a sum of sines instead of cosines and ``\sin(x) = \cos(x - \pi/2)`` (see page 12 and 13 of [this](https://eeweb.engineering.nyu.edu/iselesni/EL713/zoom/linphase.pdf)).

...
# Arguments
- `_αn::Vector` : a vector of size `(J,)`
- `n::Real`
- `fir_type::Union{FIR_III,FIR_IV}`
...
"""
function _update_trig_arg_b!(_αn, n, fir_type::Union{FIR_III,FIR_IV})
    _αn .*= n
    _αn .-= π/2
end

@doc raw"""
    _to_impulse_response(a, fir_type)

Creates a linear phase FIR filter based on `fir_type` and the coefficients in vector ``a``, which was obtained by solving the linear equation ``Qa = b``.

...
# Arguments
-`a`::Vector      : a vector of size `(M+1,)` with coefficients.   
- `fir_type::FIR` : indicates the type of FIR filter.

# Outputs
-`h`::Vector      : a vector of size (filter_order+1,) with the filter coefficients.
...

For type I FIR filters:
```math
    h = \bigg[ a[M+1] \quad a[M] \quad \cdots \quad a[2] \quad a[1] \quad a[2] \quad \cdots \quad a[M+1] \bigg]^T
```
"""
function _to_impulse_response(a, fir_type::FIR_I)
    filter_length = 2length(a) - 1
    N_a, h = length(a), zeros(eltype(a), filter_length)
    h[1:N_a-1] .= @view(a[end:-1:2])
    h[N_a] = 2a[1]
    h[N_a+1:end] .= @view(a[2:end])
    return h
end
@doc raw"""
For type II FIR filters:
```math
    h = \bigg[ a[M+1] \quad a[M] \quad \cdots \quad a[2] \quad a[1] \quad a[1] \quad a[2] \quad \cdots \quad a[M+1] \bigg]^T
```
"""
function _to_impulse_response(a, fir_type::FIR_II)
    filter_length = 2length(a)
    N_a, h = length(a), zeros(eltype(a), filter_length)
    h[1:N_a] .= @view(a[end:-1:1])
    h[N_a+1:end] .= a
    return h
end
@doc raw"""
For type III FIR filters:
```math
    h = \bigg[ a[M+1] \quad a[M] \quad \cdots \quad a[2] \quad 0 \quad -a[2] \quad \cdots \quad -a[M+1] \bigg]^T
```
"""
function _to_impulse_response(a, fir_type::FIR_III)
    filter_length = 2length(a)-1
    N_a, h = length(a), zeros(eltype(a), filter_length)
    h[1:N_a-1] .= @view(a[end:-1:2])
    h[N_a+1:end] .-= @view(a[2:end])
    return h
end
@doc raw"""
For type IV FIR filters:
```math
    h = \bigg[ a[M+1] \quad a[M] \quad \cdots \quad a[2] \quad 0 \quad 0 \quad -a[2] \quad \cdots \quad -a[M+1] \bigg]^T
```
"""
function _to_impulse_response(a, fir_type::FIR_IV)
    filter_length = 2length(a)
    N_a, h = length(a), zeros(eltype(a), filter_length)
    h[1:N_a] .= @view(a[end:-1:1])
    h[N_a+1:end] .-= a
    return h
end

end


