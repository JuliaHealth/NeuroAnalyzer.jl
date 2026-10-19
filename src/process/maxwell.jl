export maxwell_filter
export maxwell_filter!
export read_fine_calibration
export read_cross_talk

# ===========================================================================
# Signal-space separation (SSS), temporal SSS (tSSS) and movement 
# compensation for Neuromag/MEGIN Vectorview/TRIUX data.
#
# Implementation follows the published method (Taulu & Kajola 2005, Taulu &
# Simola 2006) and the MNE-Python implementation (maxwell_filter) as the
# reference for conventions: real orthonormal spherical harmonics, μ0 field
# scaling, magnetometer scaling (mag_scale), column-normalized pseudoinverse,
# information-based regularization of the internal basis, tSSS overlap
# projector, fine calibration and cross-talk correction.
# ===========================================================================

const _MU0 = 4.0e-7 * π

# ---------------------------------------------------------------------------
# coil integration points ("accurate" definitions, MNE coil_def.dat)
# points: (weight, x, y, z) in the coil frame (ex, ey, ez); normal = ez
# ---------------------------------------------------------------------------
const _VV_PLANAR_ACC = [
    14.9858 1.079e-02 6.713e-03 3.0e-04;
    14.9858 5.891e-03 6.713e-03 3.0e-04;
    14.9858 5.891e-03 -6.713e-03 3.0e-04;
    14.9858 1.079e-02 -6.713e-03 3.0e-04;
    -14.9858 -1.079e-02 6.713e-03 3.0e-04;
    -14.9858 -5.891e-03 6.713e-03 3.0e-04;
    -14.9858 -5.891e-03 -6.713e-03 3.0e-04;
    -14.9858 -1.079e-02 -6.713e-03 3.0e-04
]

function _vv_mag_acc(a::Float64)::Matrix{Float64}
    # 4 × 4 grid at ±a/... spacing as in coil_def.dat (a = 3/8 × coil size)
    c = [-a, -a / 3, a / 3, a]
    return reduce(vcat, [[0.0625 x y 3.0e-4] for x in c for y in c])
end

const _COIL_DEFS = Dict{Int64, Tuple{Symbol, Matrix{Float64}}}(
    2000 => (:mag, [1.0 0.0 0.0 0.0]),
    3012 => (:grad, _VV_PLANAR_ACC),
    3013 => (:grad, _VV_PLANAR_ACC),
    3014 => (:grad, _VV_PLANAR_ACC),
    3015 => (:grad, _VV_PLANAR_ACC),
    3022 => (:mag, _vv_mag_acc(9.675e-03)),
    3023 => (:mag, _vv_mag_acc(9.675e-03)),
    3024 => (:mag, _vv_mag_acc(7.875e-03)),
    3025 => (:mag, _vv_mag_acc(1.050e-02)),
)

# ---------------------------------------------------------------------------
# minimal forward-mode dual numbers (value + gradient in x, y, z)
# ---------------------------------------------------------------------------
struct _D3
    v::Float64
    gx::Float64
    gy::Float64
    gz::Float64
end
_D3(v::Real) = _D3(Float64(v), 0.0, 0.0, 0.0)
Base.:+(a::_D3, b::_D3) = _D3(a.v + b.v, a.gx + b.gx, a.gy + b.gy, a.gz + b.gz)
Base.:-(a::_D3, b::_D3) = _D3(a.v - b.v, a.gx - b.gx, a.gy - b.gy, a.gz - b.gz)
Base.:*(a::_D3, b::_D3) = _D3(
    a.v * b.v,
    a.gx * b.v + b.gx * a.v,
    a.gy * b.v + b.gy * a.v,
    a.gz * b.v + b.gz * a.v,
)
Base.:*(a::_D3, s::Real) = _D3(a.v * s, a.gx * s, a.gy * s, a.gz * s)
Base.:*(s::Real, a::_D3) = a * s
Base.:/(a::_D3, b::_D3) = _D3(
    a.v / b.v,
    (a.gx * b.v - b.gx * a.v) / b.v^2,
    (a.gy * b.v - b.gy * a.v) / b.v^2,
    (a.gz * b.v - b.gz * a.v) / b.v^2,
)
function _d3sqrt(a::_D3)::_D3
    s = sqrt(a.v)
    return _D3(s, a.gx / (2s), a.gy / (2s), a.gz / (2s))
end

# number of SSS components for expansion order L
_sss_ncomp(L::Int64)::Int64 = (L + 1)^2 - 1
# column index of (l, m), l ≥ 1, -l ≤ m ≤ l
_sss_idx(l::Int64, m::Int64)::Int64 = l^2 + l + m

# degree and order of each column
function _sss_degrees_orders(L::Int64)::Tuple{Vector{Int64}, Vector{Int64}}
    d = Int64[]
    o = Int64[]
    for l = 1:L, m = (-l):l
        push!(d, l)
        push!(o, m)
    end
    return d, o
end

"""
    _sss_harmonics(p, L)

Real orthonormal regular solid harmonics `R_lm = r^l Y_lm` and irregular solid harmonics `I_lm = Y_lm / r^(l+1)` (l = 1..L) at point `p` (relative to the expansion origin), as dual numbers carrying their Cartesian gradients.

Uses `r^l P_l^m(cos θ) e^{imφ} = Π_l^m(z, r²) (x + iy)^m` with the recurrence `(l-m) Π_l^m = (2l-1) z Π_{l-1}^m - (l+m-1) r² Π_{l-2}^m`, which is free of singularities on the z-axis.
"""
function _sss_harmonics(p::AbstractVector{Float64}, L::Int64)::Tuple{Vector{_D3}, Vector{_D3}}
    x = _D3(p[1], 1.0, 0.0, 0.0)
    y = _D3(p[2], 0.0, 1.0, 0.0)
    z = _D3(p[3], 0.0, 0.0, 1.0)
    r2 = x * x + y * y + z * z
    one = _D3(1.0)

    # (x + iy)^m = C_m + i S_m
    C = Vector{_D3}(undef, L + 1)
    S = Vector{_D3}(undef, L + 1)
    C[1] = one
    S[1] = _D3(0.0)
    for m = 1:L
        C[m + 1] = x * C[m] - y * S[m]
        S[m + 1] = x * S[m] + y * C[m]
    end

    # Π[l+1, m+1] = Π_l^m
    Π = Matrix{_D3}(undef, L + 1, L + 1)
    for m = 0:L
        dfac = prod(Float64[k for k = 1:2:(2m - 1)]; init = 1.0)   # (2m-1)!!
        Π[m + 1, m + 1] = _D3(dfac)
        m + 1 <= L && (Π[m + 2, m + 1] = z * Π[m + 1, m + 1] * (2m + 1))
        for l = (m + 2):L
            Π[l + 1, m + 1] =
                (z * Π[l, m + 1] * (2l - 1) - r2 * Π[l - 1, m + 1] * (l + m - 1)) *
                (1.0 / (l - m))
        end
    end

    n = _sss_ncomp(L)
    R = Vector{_D3}(undef, n)
    Iv = Vector{_D3}(undef, n)
    r = _d3sqrt(r2)
    # r^(2l+1), starting with l = 0
    rpow = r
    for l = 1:L
        rpow = rpow * r2
        for m = (-l):l
            am = abs(m)
            nrm = sqrt((2l + 1) / (4π) * exp(_lfact(l - am) - _lfact(l + am)))
            m != 0 && (nrm *= sqrt(2.0))
            ang = m >= 0 ? C[am + 1] : S[am + 1]
            rv = Π[l + 1, am + 1] * ang * nrm
            k = _sss_idx(l, m)
            R[k] = rv
            Iv[k] = rv / rpow
        end
    end
    return R, Iv
end

# log factorial
_lfact(n::Int64)::Float64 = sum(log, 1:n; init = 0.0)

"""
    _sss_coil_points(loc, coil_id, T)

Integration points of MEG coils. `loc` is n × 12 (r0, ex, ey, ez; device frame), `T` a 4 × 4 transform applied to positions and normals (device → target frame).

Returns points (np × 3), weighted normals (np × 3) and the coil index of each point.
"""
function _sss_coil_points(
    loc::AbstractMatrix{Float64},
    coil_id::AbstractVector{Int64},
    T::AbstractMatrix{Float64},
)::Tuple{Matrix{Float64}, Matrix{Float64}, Vector{Int64}}
    R = T[1:3, 1:3]
    t = T[1:3, 4]
    pts = Vector{Vector{Float64}}()
    wn = Vector{Vector{Float64}}()
    cidx = Int64[]
    for k in axes(loc, 1)
        haskey(_COIL_DEFS, coil_id[k]) || throw(
            ArgumentError(
                "Coil type $(coil_id[k]) is not supported by Maxwell filtering (Neuromag/MEGIN Vectorview/TRIUX only).",
            ),
        )
        def = _COIL_DEFS[coil_id[k]][2]
        r0 = loc[k, 1:3]
        ex = loc[k, 4:6]
        ey = loc[k, 7:9]
        ez = loc[k, 10:12]
        for row in axes(def, 1)
            p = r0 .+ def[row, 2] .* ex .+ def[row, 3] .* ey .+ def[row, 4] .* ez
            push!(pts, R * p .+ t)
            push!(wn, def[row, 1] .* (R * ez))
            push!(cidx, k)
        end
    end
    return permutedims(reduce(hcat, pts)), permutedims(reduce(hcat, wn)), cidx
end

"""
    _sss_basis(pts, wn, cidx, n_coils, origin, int_order, ext_order)

Internal (n_coils × n_in) and external (n_coils × n_out) SSS bases: field of each real orthonormal multipole, `B = -μ0 ∇V`, projected on the weighted coil normals and summed over integration points. No magnetometer scaling and no column normalization.
"""
function _sss_basis(
    pts::AbstractMatrix{Float64},
    wn::AbstractMatrix{Float64},
    cidx::AbstractVector{Int64},
    n_coils::Int64,
    origin::AbstractVector{Float64},
    int_order::Int64,
    ext_order::Int64,
)::Tuple{Matrix{Float64}, Matrix{Float64}}
    L = max(int_order, ext_order)
    n_in = _sss_ncomp(int_order)
    n_out = _sss_ncomp(ext_order)
    S_in = zeros(n_coils, n_in)
    S_out = zeros(n_coils, n_out)
    for j in axes(pts, 1)
        R, Iv = _sss_harmonics(pts[j, :] .- origin, L)
        nx, ny, nz = wn[j, 1], wn[j, 2], wn[j, 3]
        c = cidx[j]
        for k = 1:n_in
            S_in[c, k] -= _MU0 * (Iv[k].gx * nx + Iv[k].gy * ny + Iv[k].gz * nz)
        end
        for k = 1:n_out
            S_out[c, k] -= _MU0 * (R[k].gx * nx + R[k].gy * ny + R[k].gz * nz)
        end
    end
    return S_in, S_out
end

# pseudoinverse with column normalization (returned in the original scaling)
function _col_norm_pinv(S::AbstractMatrix{Float64})::Tuple{Matrix{Float64}, Vector{Float64}}
    nrm = [norm(view(S, :, k)) for k in axes(S, 2)]
    F = svd(S ./ nrm')
    pS = (F.V ./ F.S') * F.U'
    return pS ./ nrm, F.S
end

"""
    _sss_regularize_in(S, n_in, int_order, keep_out)

Information-based regularization of the internal basis (as in MNE-Python, `regularize="in"`): iteratively remove the internal component with the lowest estimated SNR and keep the set that reaches 98% of the peak total information.

Returns the indices of internal components to keep.
"""
function _sss_regularize_in(
    S::AbstractMatrix{Float64},
    n_in::Int64,
    int_order::Int64,
    keep_out::Vector{Int64},
)::Vector{Int64}
    degrees, orders = _sss_degrees_orders(int_order)
    # expected power of a randomly activated sphere of currents per degree
    ρ_i = 5.91107375632e-07
    r_in = 0.080
    a_l = [ρ_i^2 * l * r_in^(2l + 4) / ((2l + 1)^2 * (l + 1)) for l = 1:int_order]
    a_lm = a_l[degrees]
    noise = (5.0e-13)^2

    nrm = [norm(view(S, :, k)) for k in axes(S, 2)]
    Sn = S ./ nrm'
    keepers = collect(1:n_in)
    remove_order = Int64[]
    info = zeros(n_in)

    for ii = 1:n_in
        cols = vcat(keepers, keep_out)
        F = svd(Sn[:, cols])
        v_in = F.V[1:length(keepers), :] ./ nrm[keepers]
        P = (v_in ./ F.S') * F.U'
        η = noise .* vec(sum(abs2, P; dims = 2))
        η .*= [orders[k] == 0 ? 2.0 : 1.0 for k in keepers] .* 0.0025
        snr = a_lm[keepers] ./ η
        info[ii] = 0.5 * sum(log2.(snr .+ 1.0))
        kmin = argmin(snr)
        push!(remove_order, keepers[kmin])
        deleteat!(keepers, kmin)
        isempty(keepers) && break
        mx = maximum(info)
        ii > 10 && info[ii - 1] < 0.95 * mx && info[ii] < 0.95 * mx && break
    end

    mx = maximum(info)
    lim = findfirst(>=(0.98 * mx), info)
    removed = remove_order[1:(lim - 1)]
    return setdiff(1:n_in, removed)
end

# orthonormal basis of the column space (rank-revealing SVD)
function _orth(A::AbstractMatrix{Float64})::Matrix{Float64}
    F = svd(A)
    isempty(F.S) && return zeros(size(A, 1), 0)
    tol = maximum(size(A)) * F.S[1] * eps(Float64)
    return F.U[:, F.S .> tol]
end

"""
    _tsss_projector(data_int, data_res, corr)

Temporal projector of tSSS (Taulu & Simola 2006, as in MaxFilter and MNE-Python): time courses shared by the internal part and the residual (the part explained by neither the internal nor the external expansion: nearby interference) whose subspace correlation is ≥ `corr`. Inputs are not scaled (T, T/m). Inputs are channels × samples; returns samples × k orthonormal time vectors.
"""
function _tsss_projector(
    data_int::AbstractMatrix{Float64},
    data_res::AbstractMatrix{Float64},
    corr::Float64,
)::Matrix{Float64}
    n_int = norm(data_int)
    n_res = norm(data_res)
    (n_int == 0 || n_res == 0) && return zeros(size(data_int, 2), 0)
    Q_int = _orth(permutedims(data_int ./ n_int))   # samples × k1
    Q_res = _orth(permutedims(data_res ./ n_res))   # samples × k2
    (isempty(Q_int) || isempty(Q_res)) && return zeros(size(data_int, 2), 0)
    F = svd(Q_int' * Q_res)
    keep = F.S .>= corr
    return Q_res * F.V[:, keep]
end

# rotation matrix ↔ unit quaternion (q0 ≥ 0; q = (q1, q2, q3))
function _quat2rot(q::AbstractVector{<:Real})::Matrix{Float64}
    q1, q2, q3 = q
    q0 = sqrt(max(0.0, 1.0 - q1^2 - q2^2 - q3^2))
    return [
        q0^2+q1^2-q2^2-q3^2 2(q1 * q2 - q0 * q3) 2(q1 * q3 + q0 * q2)
        2(q1 * q2 + q0 * q3) q0^2-q1^2+q2^2-q3^2 2(q2 * q3 - q0 * q1)
        2(q1 * q3 - q0 * q2) 2(q2 * q3 + q0 * q1) q0^2-q1^2-q2^2+q3^2
    ]
end

function _rot2quat(R::AbstractMatrix{<:Real})::Vector{Float64}
    tr = R[1, 1] + R[2, 2] + R[3, 3]
    q = if tr > 0
        s = 2sqrt(tr + 1.0)
        [0.25s, (R[3, 2] - R[2, 3]) / s, (R[1, 3] - R[3, 1]) / s, (R[2, 1] - R[1, 2]) / s]
    elseif R[1, 1] > R[2, 2] && R[1, 1] > R[3, 3]
        s = 2sqrt(1.0 + R[1, 1] - R[2, 2] - R[3, 3])
        [(R[3, 2] - R[2, 3]) / s, 0.25s, (R[1, 2] + R[2, 1]) / s, (R[1, 3] + R[3, 1]) / s]
    elseif R[2, 2] > R[3, 3]
        s = 2sqrt(1.0 + R[2, 2] - R[1, 1] - R[3, 3])
        [(R[1, 3] - R[3, 1]) / s, (R[1, 2] + R[2, 1]) / s, 0.25s, (R[2, 3] + R[3, 2]) / s]
    else
        s = 2sqrt(1.0 + R[3, 3] - R[1, 1] - R[2, 2])
        [(R[2, 1] - R[1, 2]) / s, (R[1, 3] + R[3, 1]) / s, (R[2, 3] + R[3, 2]) / s, 0.25s]
    end
    q ./= norm(q)
    q[1] < 0 && (q .*= -1)
    return q[2:4]
end

_hom(R::AbstractMatrix, t::AbstractVector)::Matrix{Float64} = [R t; 0.0 0.0 0.0 1.0]

# channel name normalization for matching calibration / cross-talk files: imported labels
# are "MEG 113", the files use "MEG0113" (both → "MEG113")
_sss_name(s::AbstractString)::String =
    replace(uppercase(replace(s, r"[^A-Za-z0-9]" => "")), r"^MEG0+(?=\d)" => "MEG")

"""
    read_fine_calibration(file_name)

Read a Neuromag/MEGIN fine-calibration file (`sss_cal.dat`).

# Arguments

- `file_name::String`

# Returns

- `Dict{String, NamedTuple}`: normalized channel name (e.g. `"MEG113"` for MEG0113) → `(loc = 12-element vector (r0, ex, ey, ez; device frame, m), imb_cal = calibration factor (magnetometers) or 1/3 imbalance terms (gradiometers))`
"""
function read_fine_calibration(file_name::String)::Dict{String, NamedTuple}
    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))
    cal = Dict{String, NamedTuple}()
    for (ln, line) in enumerate(eachline(file_name))
        s = strip(line)
        (isempty(s) || startswith(s, "#")) && continue
        tok = split(s)
        length(tok) in (14, 16) ||
            throw(ArgumentError("$file_name, line $ln: expected 14 or 16 values, found $(length(tok))."))
        name = String(tok[1])
        if length(name) in (3, 4) && all(isdigit, name)
            name = "MEG" * lpad(name, 4, '0')
        end
        v = parse.(Float64, tok[2:end])
        cal[_sss_name(name)] = (loc = v[1:12], imb_cal = v[13:end])
    end
    return cal
end

"""
    read_cross_talk(file_name)

Read a Neuromag/MEGIN cross-talk compensation file (`ct_sparse.fif`).

# Arguments

- `file_name::String`

# Returns

Named tuple:

- `ch_names::Vector{String}`: normalized channel names (e.g. `"MEG113"` for MEG0113)
- `decoupler::Matrix{Float64}`: cross-talk (decoupler) matrix
"""
function read_cross_talk(file_name::String)::@NamedTuple{
    ch_names::Vector{String},
    decoupler::Matrix{Float64},
}
    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))
    buf, fb, _ = open(file_name, "r") do fid
        _create_fiff_block(fid)
    end
    k_mat = _find_fiff_tag("decoupler_matrix")
    k_names = _find_fiff_tag("proj_item_ch_name_list")
    mat = nothing
    names = nothing
    for i in eachindex(buf)
        if fb[i, 2] == k_mat && isnothing(mat)
            mat = _fiff_decode(fb[i, 3], buf[i])
        elseif fb[i, 2] == k_names && isnothing(names)
            names = _fiff_decode(fb[i, 3], buf[i])
        end
    end
    (isnothing(mat) || isnothing(names)) &&
        throw(ArgumentError("$file_name does not contain a cross-talk decoupler matrix."))
    ch = [_sss_name(String(n)) for n in split(names, ':')]
    size(mat) == (length(ch), length(ch)) ||
        throw(ArgumentError("$file_name: decoupler matrix size does not match channel list."))
    return (ch_names = ch, decoupler = Matrix{Float64}(mat))
end

# rotation matrix mapping unit vector a onto unit vector b (Rodrigues' formula, minimal
# rotation about a × b; antiparallel vectors: half turn about an axis perpendicular to a)
function _vec_rotation(a::AbstractVector{Float64}, b::AbstractVector{Float64})::Matrix{Float64}
    v = cross(a, b)
    s = dot(v, v)
    c = dot(a, b)
    if s < floatmin(Float64)
        c > 0 && return Matrix{Float64}(I(3))
        k = zeros(3)
        k[argmin(abs.(a))] = 1.0
        k = normalize(k .- dot(k, a) .* a)
        return 2.0 .* (k * k') .- I(3)
    end
    vx = [0.0 -v[3] v[2]; v[3] 0.0 -v[1]; -v[2] v[1] 0.0]
    return I(3) .+ vx .+ (vx * vx) .* ((1.0 - c) / s)
end

# sensor orientation (ex, ey, ez; 9 elements) after fine calibration, as in MNE-Python:
# calibrated orientations with unit length; if ex / ey of the calibration file are not
# orthogonal to its ez (files that only rotate ez), the original coil frame is rotated
# onto the calibrated ez
function _fine_cal_orientation(
    o0::AbstractVector{Float64},
    oc::AbstractVector{Float64},
)::Vector{Float64}
    R0 = reshape(Vector{Float64}(o0), 3, 3)
    R0 ./= norm.(eachcol(R0))'
    Rc = reshape(Vector{Float64}(oc), 3, 3)
    Rc ./= norm.(eachcol(Rc))'
    if maximum(abs.(Rc[:, 3]' * Rc[:, 1:2])) > 1.0e-5
        Rc = _vec_rotation(R0[:, 3], Rc[:, 3]) * R0
        Rc ./= norm.(eachcol(Rc))'
    end
    return vec(Rc)
end

# Hann window of n points (as scipy.signal.get_window("hann", n, fftbins = iseven(n)))
function _sss_hann(n::Int64)::Vector{Float64}
    n == 1 && return [1.0]
    m = iseven(n) ? n : n - 1
    return [0.5 - 0.5 * cos(2π * k / m) for k = 0:(n - 1)]
end

# device → head transform of a tSSS buffer (samples b1:b2-1, 0-based): the average of the
# head positions in effect, weighted by their duration (quaternion average: largest
# eigenvector of the weighted sum of outer products, Markley et al. 2007)
function _sss_avg_trans(
    cp::Vector{Int64},
    cp_T::Vector{Matrix{Float64}},
    cp_q::Vector{Vector{Float64}},
    b1::Int64,
    b2::Int64,
)::Matrix{Float64}
    k1 = findlast(<=(b1), cp)
    k2 = findlast(<(b2), cp)
    k1 == k2 && return cp_T[k1]
    w = [min(k == length(cp) ? b2 : cp[k + 1], b2) - max(cp[k], b1) for k = k1:k2]
    w = w ./ sum(w)
    A = zeros(4, 4)
    tr = zeros(3)
    for (wi, k) in zip(w, k1:k2)
        q = cp_q[k]
        q4 = vcat(sqrt(max(0.0, 1.0 - sum(abs2, q[1:3]))), q[1:3])
        A .+= wi .* (q4 * q4')
        tr .+= wi .* q[4:6]
    end
    qa = eigen(Symmetric(A)).vectors[:, end]
    qa .*= sign(qa[findfirst(!=(0), qa)])
    return _hom(_quat2rot(qa[2:4]), tr)
end

"""
    maxwell_filter(obj; <keyword arguments>)

Apply Maxwell filtering (signal-space separation, SSS), optionally temporal SSS (tSSS) and head movement compensation to Neuromag/MEGIN Vectorview/TRIUX MEG data.

The measured field is expanded into internal (brain) and external (environment) multipole components; data are reconstructed from the regularized internal components only. Bad MEG channels are excluded from the expansion and reconstructed from it (they are marked as good afterwards). Non-MEG channels are not modified. SSP projectors present before SSS are removed (they are not valid after SSS).

The processing follows MNE-Python (`mne.preprocessing.maxwell_filter`, defaults `st_fixed=True`, `st_only=False`): cross-talk compensation of the good channels, tSSS on the original data with the average head position of each buffer, then SSS reconstruction (with movement compensation) at `destination`. `st_overlap = false` and `mc_interp = :zero` give MaxFilter-like non-overlapping tSSS buffers and step changes between head positions.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object (continuous MEG data, 1 epoch)
- `int_order::Int64=8`: order of the internal expansion
- `ext_order::Int64=3`: order of the external expansion
- `origin::Union{Symbol, Vector{<:Real}}=:auto`: expansion origin (m) in `coord_frame`; `:auto` uses the head origin fitted to digitization points at import (`coord_frame = :head`) or (0, 0, 0) (`coord_frame = :meg`)
- `coord_frame::Symbol=:head`: coordinate frame of the expansion, `:head` or `:meg` (use `:meg` for empty-room data)
- `regularize::Symbol=:in`: `:in` (information-based removal of internal components) or `:none`
- `mag_scale::Real=100.0`: magnetometer scaling relative to gradiometers (T vs T/m)
- `calibration::Union{Nothing, String}=nothing`: fine-calibration file (`sss_cal.dat`)
- `cross_talk::Union{Nothing, String}=nothing`: cross-talk compensation file (`ct_sparse.fif`)
- `st_duration::Union{Nothing, Real}=nothing`: tSSS buffer length (s); `nothing` = SSS only
- `st_correlation::Real=0.98`: tSSS subspace correlation limit
- `st_overlap::Bool=true`: tSSS buffers overlap by 50% and are combined with a Hann window; `false`: non-overlapping buffers (a short last buffer is merged into the previous one)
- `head_pos::Union{Nothing, AbstractMatrix{<:Real}}=nothing`: head positions for movement compensation, as returned by `chpi_head_positions()` (rows: t (s, from the start of the recording), q1, q2, q3, x, y, z, gof, err, v; device → head); the recording's device → head transform is used before the first position; times in MaxFilter / MNE-Python head position files (`.pos`) are from the start of the acquisition: subtract `obj.header.recording[:first_time]`
- `mc_interp::Symbol=:hann`: interpolation of the reconstruction between consecutive head positions: `:hann` (cos²), `:linear` or `:zero` (step change)
- `destination::Union{Nothing, AbstractMatrix{<:Real}, AbstractVector{<:Real}}=nothing`: head position to transform data to: 4 × 4 device → head transform or 3-element translation (m, without rotation); default: the recording's device → head transform

# Returns

- `NeuroAnalyzer.NEURO`

# References

1. Taulu S, Kajola M. Presentation of electromagnetic multichannel data: the signal space separation method. J Appl Phys. 2005;97(12):124905.
2. Taulu S, Simola J. Spatiotemporal signal space separation method for rejecting nearby interference in MEG measurements. Phys Med Biol. 2006;51(7):1759–1768.
"""
function maxwell_filter(
    obj::NeuroAnalyzer.NEURO;
    int_order::Int64 = 8,
    ext_order::Int64 = 3,
    origin::Union{Symbol, Vector{<:Real}} = :auto,
    coord_frame::Symbol = :head,
    regularize::Symbol = :in,
    mag_scale::Real = 100.0,
    calibration::Union{Nothing, String} = nothing,
    cross_talk::Union{Nothing, String} = nothing,
    st_duration::Union{Nothing, Real} = nothing,
    st_correlation::Real = 0.98,
    st_overlap::Bool = true,
    head_pos::Union{Nothing, AbstractMatrix{<:Real}} = nothing,
    mc_interp::Symbol = :hann,
    destination::Union{Nothing, AbstractMatrix{<:Real}, AbstractVector{<:Real}} = nothing,
)::NeuroAnalyzer.NEURO
    # ------------------------------------------------------------------ #
    # validate                                                           #
    # ------------------------------------------------------------------ #
    _check_datatype(obj, "meg")
    nepochs(obj) == 1 ||
        throw(ArgumentError("Maxwell filtering requires continuous (non-epoched) data."))
    rec = obj.header.recording
    get(rec, :sss, false) && throw(ArgumentError("OBJ has already been Maxwell filtered."))
    1 <= int_order <= 12 || throw(ArgumentError("int_order must be in [1, 12]."))
    1 <= ext_order <= 8 || throw(ArgumentError("ext_order must be in [1, 8]."))
    _check_var(coord_frame, [:head, :meg], "coord_frame")
    _check_var(regularize, [:in, :none], "regularize")
    _check_var(mc_interp, [:hann, :linear, :zero], "mc_interp")
    mag_scale > 0 || throw(ArgumentError("mag_scale must be > 0."))
    0 < st_correlation <= 1 || throw(ArgumentError("st_correlation must be in (0, 1]."))
    !isnothing(st_duration) && st_duration <= 0 &&
        throw(ArgumentError("st_duration must be > 0."))

    for k in (:coil_loc, :coil_id)
        (haskey(rec, k) && !isempty(rec[k])) || throw(
            ArgumentError(
                "OBJ has no MEG sensor geometry ($k); re-import the data with import_fiff().",
            ),
        )
    end

    meg = sort(vcat(rec[:magnetometers], rec[:gradiometers]))
    isempty(meg) && throw(ArgumentError("OBJ contains no MEG sensors."))
    n_meg = length(meg)
    is_mag = [i in rec[:magnetometers] for i in meg]
    coil_id = rec[:coil_id][meg]
    loc = Matrix{Float64}(rec[:coil_loc][meg, :])
    ch_names = _sss_name.(labels(obj)[meg])
    bad = rec[:bad_channel][meg]
    good = findall(.!bad)
    length(good) > _sss_ncomp(int_order) + _sss_ncomp(ext_order) || throw(
        ArgumentError(
            "Too few good MEG channels ($(length(good))) for int_order = $int_order and ext_order = $ext_order.",
        ),
    )

    # device → head transform of the recording
    dev_head = get(rec, :dev_head_t, Matrix{Float64}(undef, 0, 0))
    if coord_frame === :head
        size(dev_head) == (4, 4) || throw(
            ArgumentError(
                "No device → head transform; use coord_frame = :meg (e.g. for empty-room data).",
            ),
        )
    else
        isnothing(head_pos) ||
            throw(ArgumentError("Movement compensation requires coord_frame = :head."))
        isnothing(destination) ||
            throw(ArgumentError("destination requires coord_frame = :head."))
    end

    # expansion origin
    org = if origin === :auto
        if coord_frame === :head
            ho = get(rec, :head_origin, Float64[])
            length(ho) == 3 || throw(
                ArgumentError(
                    "No head origin available (too few digitization points); specify origin explicitly, e.g. [0.0, 0.0, 0.04].",
                ),
            )
            Vector{Float64}(ho)
        else
            zeros(3)
        end
    elseif origin isa Symbol
        throw(ArgumentError("origin must be :auto or a 3-element vector."))
    else
        length(origin) == 3 || throw(ArgumentError("origin must have 3 elements."))
        Vector{Float64}(origin)
    end
    norm(org) > 0.1 &&
        _warn("Origin is $(round(norm(org) * 1000; digits = 1)) mm from the frame center; check its units (m).")

    # destination (reconstruction) transform
    dest = if isnothing(destination)
        coord_frame === :head ? Matrix{Float64}(dev_head) : Matrix{Float64}(I(4))
    elseif destination isa AbstractVector
        length(destination) == 3 || throw(ArgumentError("destination vector must have 3 elements."))
        _hom(Matrix{Float64}(I(3)), Vector{Float64}(destination))
    else
        size(destination) == (4, 4) || throw(ArgumentError("destination matrix must be 4 × 4."))
        Matrix{Float64}(destination)
    end
    T_rec = coord_frame === :head ? Matrix{Float64}(dev_head) : Matrix{Float64}(I(4))

    # ------------------------------------------------------------------ #
    # fine calibration                                                   #
    # ------------------------------------------------------------------ #
    # sensor orientations from the calibration file (positions are not changed),
    # magnetometer calibration factors and gradiometer imbalances (modeled by point
    # magnetometers at the original gradiometer centers and orientations)
    loc0 = copy(loc)
    mag_cal = ones(n_meg)
    imb = [Float64[] for _ = 1:n_meg]
    if !isnothing(calibration)
        cal = read_fine_calibration(calibration)
        n_found = 0
        for k = 1:n_meg
            haskey(cal, ch_names[k]) || continue
            n_found += 1
            loc[k, 4:12] = _fine_cal_orientation(loc0[k, 4:12], cal[ch_names[k]].loc[4:12])
            if is_mag[k]
                mag_cal[k] = cal[ch_names[k]].imb_cal[1]
            else
                imb[k] = cal[ch_names[k]].imb_cal
            end
        end
        n_found < n_meg && _warn(
            "Fine calibration: $(n_meg - n_found) of $n_meg MEG channels not found in $calibration.",
        )
        _info("Fine calibration applied to $n_found channels")
    end

    # ------------------------------------------------------------------ #
    # data (SI units: T, T/m)                                            #
    # ------------------------------------------------------------------ #
    units = rec[:unit][meg]
    d_scale = [_ssp_unit_scale(u) for u in units]
    data = obj.data[meg, :, 1] ./ d_scale

    # cross-talk compensation of the good channels (bad channels are reconstructed)
    if !isnothing(cross_talk)
        ct = read_cross_talk(cross_talk)
        ct_idx = [findfirst(isequal(n), ct.ch_names) for n in ch_names]
        have = findall(!isnothing, ct_idx)
        isempty(have) && throw(ArgumentError("No MEG channels found in $cross_talk."))
        length(have) < n_meg &&
            _warn("Cross-talk: $(n_meg - length(have)) MEG channels not found in $cross_talk.")
        have = intersect(have, good)
        ci = Int64[ct_idx[k] for k in have]
        data[have, :] = ct.decoupler[ci, ci] * data[have, :]
        _info("Cross-talk compensation applied")
    end

    coil_scale = [m ? Float64(mag_scale) : 1.0 for m in is_mag]
    n_in = _sss_ncomp(int_order)
    n_out = _sss_ncomp(ext_order)
    # without magnetometers and fine calibration the homogeneous field (first 3 external
    # components) is invisible to the gradiometers
    out_keep = collect((n_in + 1 + (any(is_mag) || !isnothing(calibration) ? 0 : 3)):(n_in + n_out))

    # ------------------------------------------------------------------ #
    # basis                                                              #
    # ------------------------------------------------------------------ #
    # gradiometer imbalance: point magnetometers (scaled as magnetometers) oriented along
    # ex / ey / ez of the original gradiometer
    function _imbalance!(S::Matrix{Float64}, T::Matrix{Float64})
        for k = 1:n_meg
            isempty(imb[k]) && continue
            dirs = length(imb[k]) == 1 ? (4:6,) : (4:6, 7:9, 10:12)
            for (j, cols) in enumerate(dirs)
                l1 = zeros(1, 12)
                l1[1, 1:3] = loc0[k, 1:3]
                l1[1, 10:12] = loc0[k, cols]
                p, w, c = _sss_coil_points(l1, [2000], T)
                Pi, Po = _sss_basis(p, w, c, 1, org, int_order, ext_order)
                S[k, :] .+= (mag_scale * imb[k][j]) .* vec(hcat(Pi, Po))
            end
        end
        return nothing
    end

    # full decomposition basis (n_meg × (n_in + n_out), coil-scaled) for a device → frame
    # transform
    function _decomp_basis(T::Matrix{Float64})::Matrix{Float64}
        p, w, c = _sss_coil_points(loc, coil_id, T)
        Si, So = _sss_basis(p, w, c, n_meg, org, int_order, ext_order)
        S = hcat(Si, So) .* coil_scale
        _imbalance!(S, T)
        # calibration of magnetometers
        S ./= mag_cal
        return S
    end

    # decomposition at one head position (regularization depends on it): kept internal
    # components, pseudoinverse (applied to unscaled data) and the kept basis at the good
    # sensors (unscaled); cached by transform
    cond_max = Ref(0.0)
    decomps = Dict{Matrix{Float64}, NamedTuple}()
    function _decomp(T::Matrix{Float64})::NamedTuple
        haskey(decomps, T) && return decomps[T]
        S = _decomp_basis(T)[good, :]
        in_keep = regularize === :in ? _sss_regularize_in(S, n_in, int_order, out_keep) :
            collect(1:n_in)
        mom = vcat(in_keep, out_keep)
        pS, sing = _col_norm_pinv(S[:, mom])
        cond_max[] = max(cond_max[], sing[1] / sing[end])
        pS .*= coil_scale[good]'
        d = (in_keep = in_keep, n = length(in_keep), pS = pS, S = S[:, mom] ./ coil_scale[good])
        decomps[T] = d
        return d
    end

    # internal reconstruction basis at the destination (ideal sensors: original geometry,
    # without fine calibration)
    p, w, c = _sss_coil_points(loc0, coil_id, dest)
    S_recon, _ = _sss_basis(p, w, c, n_meg, org, int_order, 1)

    # SSS operator (good channels → all MEG channels) for one head position
    function _op_sss(T::Matrix{Float64})::Matrix{Float64}
        d = _decomp(T)
        return S_recon[:, d.in_keep] * d.pS[1:d.n, :]
    end

    # ------------------------------------------------------------------ #
    # head positions                                                     #
    # ------------------------------------------------------------------ #
    # control points: first samples (0-based) of the head positions and their transforms
    n_s = size(data, 2)
    fs = rec[:sampling_rate]
    cp = [0]
    cp_T = [T_rec]
    cp_q = [vcat(_rot2quat(T_rec[1:3, 1:3]), T_rec[1:3, 4])]
    if !isnothing(head_pos)
        size(head_pos, 2) == 10 ||
            throw(ArgumentError("head_pos must have 10 columns (t, q1-q3, x, y, z, gof, err, v)."))
        t = Float64.(head_pos[:, 1])
        issorted(t; lt = <=) || throw(ArgumentError("head_pos times must be unique and ascending."))
        t[1] < -1.0e-3 && throw(ArgumentError("head_pos times must be ≥ 0 s (from the start of the recording)."))
        maximum(norm.(eachrow(head_pos[:, 5:7]))) > 1.0 &&
            _warn("head_pos contains translations > 1 m; head positions may be invalid.")
        # the recording's transform applies before the first position
        if round(Int64, t[1] * fs) < 1
            empty!(cp)
            empty!(cp_T)
            empty!(cp_q)
        end
        for r in axes(head_pos, 1)
            s = max(0, round(Int64, t[r] * fs))
            q = Vector{Float64}(head_pos[r, 2:7])
            Tn = _hom(_quat2rot(q[1:3]), q[4:6])
            # positions rounded to the same sample: the last one is used
            if !isempty(cp) && s == cp[end]
                cp_T[end] = Tn
                cp_q[end] = q
            else
                push!(cp, s)
                push!(cp_T, Tn)
                push!(cp_q, q)
            end
        end
        cp[1] = 0
        _info("Movement compensation: $(count(<(n_s), cp)) head position(s)")
    end

    # ------------------------------------------------------------------ #
    # tSSS                                                               #
    # ------------------------------------------------------------------ #
    # applied to the original data of the good channels, buffer by buffer, using the
    # average head position of each buffer
    n_proj = Int64[]
    if !isnothing(st_duration)
        n = round(Int64, st_duration * fs)
        n > n_s && throw(ArgumentError("st_duration is longer than the recording."))
        n_ovl = st_overlap ? n ÷ 2 : 0
        step = n - n_ovl
        win = st_overlap ? _sss_hann(n) : ones(n)
        # constant overlap-add normalization
        bins = sum(win[(i * step + 1):((i + 1) * step)] for i = 0:(n ÷ step - 1))
        n % step != 0 && (bins[1:(n % step)] .+= win[(end - n % step + 1):end])
        win ./= median(bins)
        starts = collect(0:step:(n_s - n))
        stops = starts .+ n
        # the remainder is merged into the last buffer
        stops[end] = n_s
        d_good = data[good, :]
        st = zeros(size(d_good))
        for (b, (b1, b2)) in enumerate(zip(starts, stops))
            len = b2 - b1
            w_b = vcat(win, zeros(len - n))
            # the first and the last buffer are not overlapped on one side
            if b == length(starts)
                for shift = step:step:(len - 1)
                    w_b[(shift + 1):len] .+= win[1:(len - shift)]
                end
            end
            if b == 1
                for n_use = (n - step):(-step):1
                    w_b[1:n_use] .+= win[(end - n_use + 1):end]
                end
            end
            D = d_good[:, (b1 + 1):b2]
            d = _decomp(_sss_avg_trans(cp, cp_T, cp_q, b1, b2))
            op_in = d.S[:, 1:d.n] * d.pS[1:d.n, :]
            D_in = op_in * D
            D_res = D .- D_in .- d.S[:, (d.n + 1):end] * (d.pS[(d.n + 1):end, :] * D)
            P = _tsss_projector(D_in, D_res, Float64(st_correlation))
            push!(n_proj, size(P, 2))
            st[:, (b1 + 1):b2] .+= (D .- (D * P) * P') .* w_b'
        end
        data[good, :] = st
        _info(
            "tSSS: $(length(starts)) buffer(s) of $(st_duration) s$(st_overlap ? " (50% overlap)" : ""), $(minimum(n_proj))–$(maximum(n_proj)) temporal component(s) removed per buffer",
        )
    end

    # ------------------------------------------------------------------ #
    # SSS reconstruction                                                 #
    # ------------------------------------------------------------------ #
    out = similar(data)
    d_good = data[good, :]
    k_last = findlast(<(n_s), cp)
    op_r = _op_sss(cp_T[1])
    for k = 1:k_last
        op_l = op_r
        a = cp[k] + 1
        if k == length(cp)
            out[:, a:n_s] = op_l * d_good[:, a:n_s]
            break
        end
        op_r = _op_sss(cp_T[k + 1])
        b = min(cp[k + 1], n_s)
        D = d_good[:, a:b]
        if mc_interp === :zero
            out[:, a:b] = op_l * D
        else
            span = cp[k + 1] - cp[k]
            x = (0:(b - a)) ./ span
            wl = mc_interp === :linear ? 1.0 .- x : cos.(x .* (π / 2)) .^ 2
            out[:, a:b] = (op_l * D) .* wl' .+ (op_r * D) .* (1.0 .- wl')
        end
    end
    n_keep_in = _decomp(cp_T[1]).n
    _info(
        "SSS: using $n_keep_in/$n_in internal and $(length(out_keep))/$n_out external components (int_order = $int_order, ext_order = $ext_order)",
    )
    cond_max[] > 1000 && _warn(
        "SSS basis is ill-conditioned (condition number up to $(round(Int64, cond_max[]))).",
    )

    # ------------------------------------------------------------------ #
    # output object                                                      #
    # ------------------------------------------------------------------ #
    obj_tmp = deepcopy(obj)
    obj_tmp.data[meg, :, 1] = out .* d_scale
    r = obj_tmp.header.recording

    # bad MEG channels were reconstructed
    n_bad = count(bad)
    if n_bad > 0
        r[:bad_channel][meg] .= false
        _info("$n_bad bad MEG channel(s) reconstructed and marked as good.")
    end

    # pre-SSS projectors are no longer valid
    if haskey(r, :ssp_data) && !isempty(r[:ssp_data])
        _info("$(size(r[:ssp_data], 1)) SSP projector(s) removed (not valid after SSS)")
        r[:ssp_data] = Matrix{Float64}(undef, 0, 0)
        r[:ssp_channels] = zeros(Bool, length(r[:ssp_channels]))
        r[:ssp_labels] = String[]
    end

    # calibrated sensor geometry, head position after transformation
    r[:coil_loc][meg, :] = loc
    coord_frame === :head && (r[:dev_head_t] = dest)
    r[:maxshield] = false
    r[:sss] = true
    r[:tsss] = !isnothing(st_duration)
    r[:sss_info] = Dict{Symbol, Any}(
        :int_order => int_order,
        :ext_order => ext_order,
        :origin => org,
        :coord_frame => coord_frame,
        :regularize => regularize,
        :n_in_used => n_keep_in,
        :n_out_used => length(out_keep),
        :rank => n_keep_in,
        :mag_scale => Float64(mag_scale),
        :calibration => isnothing(calibration) ? "" : calibration,
        :cross_talk => isnothing(cross_talk) ? "" : cross_talk,
        :st_duration => isnothing(st_duration) ? 0.0 : Float64(st_duration),
        :st_correlation => Float64(st_correlation),
        :st_overlap => st_overlap,
        :st_components_removed => sum(n_proj; init = 0),
        :movement_compensation => !isnothing(head_pos),
        :mc_interp => mc_interp,
    )

    push!(
        obj_tmp.history,
        "maxwell_filter(obj; int_order=$int_order, ext_order=$ext_order, origin=$org, coord_frame=$coord_frame, regularize=$regularize, mag_scale=$mag_scale, calibration=$calibration, cross_talk=$cross_talk, st_duration=$st_duration, st_correlation=$st_correlation, st_overlap=$st_overlap, head_pos=$(isnothing(head_pos) ? nothing : "[…]"), mc_interp=$mc_interp, destination=$(isnothing(destination) ? nothing : "[…]"))",
    )

    return obj_tmp
end

"""
    maxwell_filter!(obj; <keyword arguments>)

Apply Maxwell filtering (SSS / tSSS / movement compensation) in-place. See `maxwell_filter()` for arguments.

# Returns

- `Nothing`
"""
function maxwell_filter!(obj::NeuroAnalyzer.NEURO; kwargs...)::Nothing
    obj_tmp = maxwell_filter(obj; kwargs...)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    obj.history = obj_tmp.history

    return nothing
end
