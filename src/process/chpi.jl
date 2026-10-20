export chpi_head_positions

# ---------------------------------------------------------------------------
# "normal" accuracy coil definitions (MNE coil_def.dat), used for the many
# dipole fits of continuous HPI; points: (weight, x, y, z) in the coil frame
# ---------------------------------------------------------------------------
const _VV_PLANAR_NORMAL = [
    29.7619 8.4e-03 6.713e-03 3.0e-04;
    29.7619 8.4e-03 -6.713e-03 3.0e-04;
    -29.7619 -8.4e-03 6.713e-03 3.0e-04;
    -29.7619 -8.4e-03 -6.713e-03 3.0e-04
]
_vv_mag_normal(a::Float64)::Matrix{Float64} =
    [0.25 -a -a 3.0e-4; 0.25 -a a 3.0e-4; 0.25 a -a 3.0e-4; 0.25 a a 3.0e-4]

const _COIL_DEFS_NORMAL = Dict{Int64, Matrix{Float64}}(
    3012 => _VV_PLANAR_NORMAL,
    3013 => _VV_PLANAR_NORMAL,
    3014 => _VV_PLANAR_NORMAL,
    3015 => _VV_PLANAR_NORMAL,
    3022 => _vv_mag_normal(6.45e-03),
    3023 => _vv_mag_normal(6.45e-03),
    3024 => _vv_mag_normal(5.25e-03),
    3025 => _vv_mag_normal(7.0e-03),
)

# integration points in the device frame (points, weighted normals, coil index)
function _chpi_coil_points(
    loc::AbstractMatrix{Float64},
    coil_id::AbstractVector{Int64},
)::Tuple{Matrix{Float64}, Matrix{Float64}, Vector{Int64}}
    pts = Vector{Vector{Float64}}()
    wn = Vector{Vector{Float64}}()
    cidx = Int64[]
    for k in axes(loc, 1)
        haskey(_COIL_DEFS_NORMAL, coil_id[k]) || throw(
            ArgumentError(
                "Coil type $(coil_id[k]) is not supported (Neuromag/MEGIN Vectorview/TRIUX only).",
            ),
        )
        def = _COIL_DEFS_NORMAL[coil_id[k]]
        r0, ex, ey, ez = loc[k, 1:3], loc[k, 4:6], loc[k, 7:9], loc[k, 10:12]
        for row in axes(def, 1)
            push!(pts, r0 .+ def[row, 2] .* ex .+ def[row, 3] .* ey .+ def[row, 4] .* ez)
            push!(wn, def[row, 1] .* ez)
            push!(cidx, k)
        end
    end
    return permutedims(reduce(hcat, pts)), permutedims(reduce(hcat, wn)), cidx
end

# lead field (n_coils × 3) of a magnetic dipole at rq (device frame, m)
function _chpi_mag_dipole_fwd(
    rq::AbstractVector{Float64},
    pts::Matrix{Float64},
    wn::Matrix{Float64},
    cidx::Vector{Int64},
    n_coils::Int64,
)::Matrix{Float64}
    G = zeros(n_coils, 3)
    c = _MU0 / (4π)
    for j in axes(pts, 1)
        dx = pts[j, 1] - rq[1]
        dy = pts[j, 2] - rq[2]
        dz = pts[j, 3] - rq[3]
        r2 = dx^2 + dy^2 + dz^2
        r = sqrt(r2)
        r3 = r2 * r
        r5 = r3 * r2
        # B for unit moments along x, y, z, projected on the weighted normal
        dn = dx * wn[j, 1] + dy * wn[j, 2] + dz * wn[j, 3]
        k = cidx[j]
        G[k, 1] += c * (3 * dx * dn / r5 - wn[j, 1] / r3)
        G[k, 2] += c * (3 * dy * dn / r5 - wn[j, 2] / r3)
        G[k, 3] += c * (3 * dz * dn / r5 - wn[j, 3] / r3)
    end
    return G
end

# goodness of fit of the best magnetic dipole at rq for field pattern b
function _chpi_gof(rq::AbstractVector{Float64}, b::Vector{Float64}, fwd::Function)::Float64
    G = fwd(rq)
    m = G \ b
    res = b - G * m
    return 1.0 - sum(abs2, res) / sum(abs2, b)
end

# Nelder-Mead minimization (no external dependencies)
function _nelder_mead(
    f::Function,
    x0::Vector{Float64};
    step::Float64 = 0.005,
    xtol::Float64 = 1.0e-5,
    maxiter::Int64 = 500,
)::Tuple{Vector{Float64}, Float64}
    n = length(x0)
    X = [copy(x0)]
    for i = 1:n
        x = copy(x0)
        x[i] += step
        push!(X, x)
    end
    F = [f(x) for x in X]
    for _ = 1:maxiter
        o = sortperm(F)
        X = X[o]
        F = F[o]
        maximum(norm(x - X[1]) for x in X) < xtol && break
        c = sum(X[1:n]) ./ n
        xr = c .+ (c .- X[end])
        fr = f(xr)
        if fr < F[1]
            xe = c .+ 2 .* (c .- X[end])
            fe = f(xe)
            if fe < fr
                X[end], F[end] = xe, fe
            else
                X[end], F[end] = xr, fr
            end
        elseif fr < F[n]
            X[end], F[end] = xr, fr
        else
            xc = fr < F[end] ? c .+ 0.5 .* (xr .- c) : c .+ 0.5 .* (X[end] .- c)
            fc = f(xc)
            if fc < min(fr, F[end])
                X[end], F[end] = xc, fc
            else
                for i = 2:(n + 1)
                    X[i] = X[1] .+ 0.5 .* (X[i] .- X[1])
                    F[i] = f(X[i])
                end
            end
        end
    end
    i = argmin(F)
    return X[i], F[i]
end

# rigid transform (R, t) with B ≈ R A + t (rows are points), Kabsch algorithm
function _kabsch(A::AbstractMatrix{Float64}, B::AbstractMatrix{Float64})::Tuple{Matrix{Float64}, Vector{Float64}}
    ca = vec(sum(A; dims = 1)) ./ size(A, 1)
    cb = vec(sum(B; dims = 1)) ./ size(B, 1)
    H = (A .- ca')' * (B .- cb')
    F = svd(H)
    d = sign(det(F.V * F.U'))
    R = F.V * Diagonal([1.0, 1.0, d]) * F.U'
    return R, cb - R * ca
end

"""
    chpi_head_positions(obj; <keyword arguments>)

Estimate continuous head positions from continuous head position indicator (cHPI) signals.

For each time window, the amplitude of every HPI coil frequency is estimated by least squares (together with line-frequency harmonics, offset and linear trend), the coil is localized by fitting a magnetic dipole to its field pattern, and the rigid device → head transform is fitted to the digitized coil positions.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object (continuous MEG data with cHPI, imported with `import_fiff()`)
- `t_window::Real=0.2`: window length (s)
- `t_step::Real=0.01`: step between windows (s)
- `gof_limit::Real=0.98`: minimum goodness of fit of a coil to be used
- `dist_limit::Real=0.005`: maximum distance (m) between fitted and digitized coil positions after the transform

# Returns

- `Matrix{Float64}`: one row per accepted window with columns: time (s, window center, from the start of the data), q1, q2, q3 (rotation quaternion), x, y, z (translation, m) of the device → head transform, mean goodness of fit, mean coil position error (m), velocity (m/s). Usable as `head_pos` in `maxwell_filter()`.
"""
function chpi_head_positions(
    obj::NeuroAnalyzer.NEURO;
    t_window::Real = 0.2,
    t_step::Real = 0.01,
    gof_limit::Real = 0.98,
    dist_limit::Real = 0.005,
)::Matrix{Float64}
    # ------------------------------------------------------------------ #
    # validate                                                           #
    # ------------------------------------------------------------------ #
    _check_datatype(obj, "meg")
    nepochs(obj) == 1 || throw(ArgumentError("cHPI fitting requires continuous (non-epoched) data."))
    t_window > 0 || throw(ArgumentError("t_window must be > 0."))
    t_step > 0 || throw(ArgumentError("t_step must be > 0."))
    0 < gof_limit < 1 || throw(ArgumentError("gof_limit must be in (0, 1)."))
    dist_limit > 0 || throw(ArgumentError("dist_limit must be > 0."))

    rec = obj.header.recording
    get(rec, :sss, false) &&
        _warn("OBJ has been Maxwell filtered; cHPI signals may have been altered.")
    hpi = get(rec, :hpi, Dict{Symbol, Any}())
    freqs = Vector{Float64}(get(hpi, :freqs, Float64[]))
    hpi_dig = get(hpi, :dig, Matrix{Float64}(undef, 0, 3))
    n_hpi = length(freqs)
    n_hpi >= 3 || throw(ArgumentError("OBJ contains no cHPI information (at least 3 coil frequencies required)."))
    size(hpi_dig, 1) == n_hpi ||
        throw(ArgumentError("Number of digitized HPI coils ($(size(hpi_dig, 1))) does not match coil frequencies ($n_hpi)."))
    dev_head = get(rec, :dev_head_t, Matrix{Float64}(undef, 0, 0))
    size(dev_head) == (4, 4) || throw(ArgumentError("OBJ has no device → head transform."))
    for k in (:coil_loc, :coil_id)
        (haskey(rec, k) && !isempty(rec[k])) ||
            throw(ArgumentError("OBJ has no MEG sensor geometry ($k); re-import the data with import_fiff()."))
    end

    fs = rec[:sampling_rate]
    any(>=(fs / 2), freqs) && throw(
        ArgumentError("HPI frequencies $(freqs) Hz are not below the Nyquist frequency ($(fs / 2) Hz)."),
    )
    get(rec, :lowpass, nothing) isa Real && 0 < rec[:lowpass] < maximum(freqs) &&
        _warn("Data are low-pass filtered below the highest HPI frequency.")

    # ------------------------------------------------------------------ #
    # sensors and data (SI units)                                        #
    # ------------------------------------------------------------------ #
    meg = sort(vcat(rec[:magnetometers], rec[:gradiometers]))
    meg = meg[.!rec[:bad_channel][meg]]
    n_meg = length(meg)
    n_meg >= 3 * n_hpi + 10 || throw(ArgumentError("Too few good MEG channels."))
    is_mag = [i in rec[:magnetometers] for i in meg]
    coil_scale = [m ? 100.0 : 1.0 for m in is_mag]
    d_scale = [_ssp_unit_scale(u) for u in rec[:unit][meg]]
    data = obj.data[meg, :, 1] ./ d_scale

    pts, wn, cidx = _chpi_coil_points(Matrix{Float64}(rec[:coil_loc][meg, :]), rec[:coil_id][meg])
    fwd = rq -> _chpi_mag_dipole_fwd(rq, pts, wn, cidx, n_meg) .* coil_scale

    # ------------------------------------------------------------------ #
    # design matrix: HPI cos/sin, line harmonics, offset, slope          #
    # ------------------------------------------------------------------ #
    nw = round(Int64, t_window * fs)
    nw >= 10 || throw(ArgumentError("t_window is too short."))
    n_s = size(data, 2)
    nw <= n_s || throw(ArgumentError("t_window is longer than the recording."))
    t = collect(0:(nw - 1)) ./ fs
    lf = get(rec, :line_frequency, 0)
    line_f = lf > 0 ? [h * lf for h = 1:floor(Int64, (fs / 2 - 1) / lf)] : Float64[]
    # do not model line harmonics coinciding with HPI frequencies
    line_f = Base.filter(f -> all(abs.(freqs .- f) .> 1.0 / t_window), line_f)
    cols = Vector{Vector{Float64}}()
    for f in freqs
        push!(cols, cos.(2π * f .* t))
        push!(cols, sin.(2π * f .* t))
    end
    for f in line_f
        push!(cols, cos.(2π * f .* t))
        push!(cols, sin.(2π * f .* t))
    end
    push!(cols, ones(nw))
    push!(cols, t .- mean(t))
    X = reduce(hcat, cols)
    pX = pinv(X)

    # ------------------------------------------------------------------ #
    # windows                                                            #
    # ------------------------------------------------------------------ #
    step = max(1, round(Int64, t_step * fs))
    starts = collect(1:step:(n_s - nw + 1))

    # initial coil positions in the device frame from digitization
    R_dh = dev_head[1:3, 1:3]
    t_dh = dev_head[1:3, 4]
    pos = [R_dh' * (hpi_dig[k, :] .- t_dh) for k = 1:n_hpi]

    rows = Vector{Vector{Float64}}()
    n_skipped = 0
    last_trans = nothing
    last_time = 0.0

    for s in starts
        Y = data[:, s:(s + nw - 1)]
        coef = pX * Y'                         # n_cols × n_meg

        fitted = Matrix{Float64}(undef, n_hpi, 3)
        gofs = zeros(n_hpi)
        for k = 1:n_hpi
            A = permutedims(coef[(2k - 1):(2k), :]) .* coil_scale    # n_meg × 2 (cos, sin)
            F = svd(A)
            b = F.U[:, 1] .* F.S[1]                                    # dominant phase pattern
            x, fv = _nelder_mead(rq -> 1.0 - _chpi_gof(rq, b, fwd), pos[k]; step = 0.002)
            fitted[k, :] = x
            gofs[k] = 1.0 - fv
        end

        use = findall(>=(gof_limit), gofs)
        # rigid transform head → device from good coils; drop the worst coil while errors are too large
        accepted = false
        R = zeros(3, 3)
        tv = zeros(3)
        err = Inf
        while length(use) >= 3
            R, tv = _kabsch(hpi_dig[use, :], fitted[use, :])
            d = [norm(R * hpi_dig[k, :] .+ tv .- fitted[k, :]) for k in use]
            err = mean(d)
            if maximum(d) <= dist_limit
                accepted = true
                break
            end
            length(use) == 3 && break
            deleteat!(use, argmax(d))
        end

        if !accepted
            n_skipped += 1
            continue
        end

        # use accepted fits as starting points for the next window
        for k in use
            pos[k] = fitted[k, :]
        end

        # device → head = inverse of the fitted head → device transform
        R_out = R'
        t_out = -R' * tv
        tc = (s - 1 + nw / 2) / fs
        v = isnothing(last_trans) ? 0.0 : norm(t_out .- last_trans) / (tc - last_time)
        last_trans = t_out
        last_time = tc
        push!(rows, vcat(tc, _rot2quat(R_out), t_out, mean(gofs[use]), err, v))
    end

    n_skipped > 0 && _warn(
        "$n_skipped of $(length(starts)) window(s) skipped (fewer than 3 coils with gof ≥ $gof_limit or coil distance > $(dist_limit * 1000) mm).",
    )
    isempty(rows) && throw(ArgumentError("No head position could be estimated."))
    _info("Head positions estimated for $(length(rows)) window(s)")

    return permutedims(reduce(hcat, rows))
end
