export denoise_cwd
export denoise_cwd!
export denoise_dwd
export denoise_dwd!

"""
    denoise_cwd(s; <keyword arguments>)

Perform denoising by zeroing a frequency band in the continuous wavelet domain.

Computes the CWD, zeros all coefficients in the band `[nf − w, nf + w]` Hz, then reconstructs the signal via the inverse CWD.

# Arguments

- `s::AbstractVector`: signal vector
- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets
- `nf::Real`: noise frequency to suppress in Hz; must be in `[1, fs/2]`
- `w::Int64=5`: half-width of the suppression band in Hz; the zeroed region is `[nf − w, nf + w]`
- `type::Symbol=:pd`: reconstruction method:
    - `:pd`: PenroseDelta (default; generally most accurate)
    - `:nd`: NaiveDelta
    - `:df`: DualFrames

# Returns

- `Vector{Float64}`: denoised signal of the same length as `s`

# Throws

- `ArgumentError`: if `fs < 1`, `nf < 1`, or `nf > fs/2`

# See also

[`denoise_cwd(::AbstractArray)`](@ref), [`denoise_cwd(::NeuroAnalyzer.NEURO)`](@ref)
"""
function denoise_cwd(
    s::AbstractVector;
    fs::Int64,
    wt::T = wavelet(Morlet(2π), β = 2),
    nf::Real,
    w::Int64 = 5,
    type::Symbol = :nd
)::Vector{Float64} where {T <: CWT}

    fs >= 1 || throw(ArgumentError("fs must be ≥ 1."))
    nf >= 1 || throw(ArgumentError("nf must be ≥ 1."))
    nf <= fs / 2 || throw(ArgumentError("nf must be ≤ $(fs / 2)."))

    # perform CWD and zero out noise frequency band
    s_cwd = cwd(s, wt = wt)
    f = cwtfrq(s, fs = fs, wt = wt)
    f_idx1 = vsearch(nf - w, f)
    f_idx2 = vsearch(nf + w, f)
    s_cwd[f_idx1:f_idx2, :] .= 0

    # reconstruct signal
    return vec(icwd(s_cwd, wt = wt, type = type))

end

"""
    denoise_cwd(s; <keyword arguments>)

Perform denoising by zeroing a frequency band in the continuous wavelet domain.

# Arguments

- `s::AbstractArray`: signal array, shape `(channels, samples, epochs)`
- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets
- `nf::Real`: noise frequency to suppress in Hz; must be in `[1, fs/2]`
- `w::Int64=5`: half-width of the suppression band in Hz; the zeroed region is `[nf − w, nf + w]`
- `type::Symbol=:pd`: reconstruction method:
    - `:pd`: PenroseDelta (default; generally most accurate)
    - `:nd`: NaiveDelta
    - `:df`: DualFrames

# Returns

- `Array{Float64, 3}`: denoised signal array

# Throws

- `ArgumentError`: if `s` is not a 3D array or if `fs`, `nf`, or `w` are invalid
"""
function denoise_cwd(
    s::AbstractArray;
    fs::Int64,
    wt::T = wavelet(Morlet(2π), β = 2),
    nf::Real,
    w::Int64 = 5,
    type::Symbol = :nd
)::Array{Float64, 3} where {T <: CWT}

    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # pre-allocate output
    s_new = similar(s, Float64)

    _log_off()
    f = cwtfrq(s[1, :, 1]; fs = fs, wt = wt)
    _log_on()
    f_idx = vsearch(nf, f)
    f_idx1 = vsearch(nf - w, f)
    f_idx2 = vsearch(nf + w, f)
    _info("Noise at: $(f[f_idx]) Hz")
    _info("Noise width: $(f[f_idx1]) to $(f[f_idx2]) Hz")

    _log_off()

    # initialize progress bar
    progbar = Progress(ep_n * ch_n, dt = 1, barlen = 20, color = :white, enabled = progress_bar)

    # calculate over channel and epochs
    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        s_new[ch_idx, :, ep_idx] = denoise_cwd(@view(
            s[ch_idx, :, ep_idx]), fs = fs, wt = wt, nf = nf, w = w, type = type
        )
        # update progress bar
        progress_bar && next!(progbar)
    end

    _log_on()

    return s_new

end

"""
    denoise_cwd(obj; <keyword arguments>)

Perform denoising by zeroing a frequency band in the continuous wavelet domain on selected channels of a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets
- `nf::Real`: noise frequency to suppress in Hz; must be in `[1, fs/2]`
- `w::Int64=5`: half-width of the suppression band in Hz; the zeroed region is `[nf − w, nf + w]`
- `type::Symbol=:pd`: reconstruction method:
    - `:pd`: PenroseDelta (default; generally most accurate)
    - `:nd`: NaiveDelta
    - `:df`: DualFrames

# Returns

- `NeuroAnalyzer.NEURO`: new object with denoised channels

# See also

[`denoise_cwd!`](@ref), [`denoise_cwd(::AbstractArray)`](@ref)
"""
function denoise_cwd(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    wt::T = wavelet(Morlet(2π), β = 2),
    nf::Real,
    w::Int64 = 5,
    type::Symbol = :nd,
)::NeuroAnalyzer.NEURO where {T <: CWT}

    # resolve channel names to integer indices
    ch = get_channel(obj, ch = ch)

    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = denoise_cwd(
        @view(obj.data[ch, :, :]),
        fs = sr(obj),
        wt = wt,
        nf = nf,
        w = w,
        type = type
    )
    push!(obj_new.history, "denoise_cwd(OBJ, ch=$ch, wt=$wt, nf=$nf, w=$w, type=$type)")

    return obj_new

end

"""
    denoise_cwd!(obj; <keyword arguments>)

Perform denoising by zeroing a frequency band in the continuous wavelet domain in-place on selected channels of a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; modified in-place
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=2)`, see ContinuousWavelets.jl documentation for the list of available wavelets
- `nf::Real`: noise frequency to suppress in Hz; must be in `[1, fs/2]`
- `w::Int64=5`: half-width of the suppression band in Hz; the zeroed region is `[nf − w, nf + w]`
- `type::Symbol=:pd`: reconstruction method:
    - `:pd`: PenroseDelta (default; generally most accurate)
    - `:nd`: NaiveDelta
    - `:df`: DualFrames

# Returns

- `Nothing`

# See also
[`denoise_cwd`](@ref)
"""
function denoise_cwd!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    wt::T = wavelet(Morlet(2π), β = 2),
    nf::Real,
    type::Symbol = :nd,
)::Nothing where {T <: CWT}

    obj_new = denoise_cwd(obj, ch = ch, wt = wt, nf = nf, type = type)
    obj.data = obj_new.data
    obj.history = obj_new.history

    return nothing

end

"""
    denoise_dwd(s; <keyword arguments>)

Perform threshold denoising using discrete wavelet decomposition (DWD).

# Arguments

- `s::AbstractVector`: signal vector
- `wt<:DiscreteWavelet=wavelet(WT.haar)`: discrete wavelet, see Wavelets.jl documentation for the list of available wavelets
- `l::Int64=0`: decomposition level; `0` uses the maximum available level (`maxtransformlevels(s)`); must be ≤ `maxtransformlevels(s)`
- `dnt<:DNF=RelErrorShrink(SoftTH())`: denoise type, see WaveletsExt.jl documentation for detailed description of available denoising functions
- `smooth::Symbol=:regular`: smoothing method:
    - `:regular`: threshold all coefficients
    - `:undersmooth`: do not threshold the lowest-frequency subspace

# Returns

- `Vector{Float64}`: denoised signal

# Throws

- `ArgumentError`: if `l > maxtransformlevels(s)` or `smooth` is invalid

# See also

[`denoise_dwd(::AbstractArray)`](@ref), [`denoise_dwd(::NeuroAnalyzer.NEURO)`](@ref)
"""
function denoise_dwd(
    s::AbstractVector;
    wt::T1 = wavelet(WT.haar),
    l::Int64 = 0,
    dnt::T2 = RelErrorShrink(SoftTH()),
    smooth::Symbol = :regular,
)::Vector{Float64} where {T1 <: DiscreteWavelet, T2 <: DNFT}

    _check_var(smooth, [:regular, :undersmooth], "smooth")

    l <= maxtransformlevels(s) || throw(ArgumentError("l must be ≤ $(maxtransformlevels(s))."))
    l == 0 && (l = maxtransformlevels(s))

    if isdyadic(length(s))
        s_new = Wavelets.denoise(s, :sig, wt; L = l, dnt = dnt, smooth = smooth)
    else
        s_new = Wavelets.denoise(s, wt; L = l, dnt = dnt)
    end

    return s_new

end

"""
    denoise_dwd(s; <keyword arguments>)

Perform threshold denoising using discrete wavelet decomposition (DWD).

# Arguments

- `s::AbstractArray`: signal array, shape `(channels, samples, epochs)`
- `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`; see Wavelets.jl documentation for the list of available wavelets
- `l::Int64=0`: decomposition level; `0` uses the maximum available level (`maxtransformlevels(s)`); must be ≤ `maxtransformlevels(s)`
- `dnt<:DNF=RelErrorShrink(SoftTH())`: denoise type, see WaveletsExt.jl documentation for detailed description of available denoising functions
- `smooth::Symbol=:regular`: smoothing method:
    - `:regular`: threshold all coefficients
    - `:undersmooth`: do not threshold the lowest-frequency subspace

# Returns

# Returns

- `Array{Float64, 3}`: denoised signal array

# Throws

- `ArgumentError`: if `s` is not a 3D array or if `fs`, `nf`, or `w` are invalid
"""
function denoise_dwd(
    s::AbstractArray;
    wt::T1 = wavelet(WT.haar),
    l::Int64 = 0,
    dnt::T2 = RelErrorShrink(SoftTH()),
    smooth::Symbol = :regular,
)::Array{Float64, 3} where {T1 <: DiscreteWavelet, T2 <: DNFT}

    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # pre-allocate output
    s_new = similar(s, Float64)

    # calculate over channel and epochs
    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        s_new[ch_idx, :, ep_idx] = denoise_dwd(
            @view(s[ch_idx, :, ep_idx]),
            wt = wt,
            l = l,
            dnt = dnt,
            smooth = smooth
        )
    end

    return s_new

end

"""
    denoise_dwd(obj; <keyword arguments>)

Perform denoising using discrete wavelet decomposition (DWD) on selected channels of a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`; see Wavelets.jl documentation for the list of available wavelets
- `l::Int64=0`: decomposition level; `0` uses the maximum available level (`maxtransformlevels(s)`); must be ≤ `maxtransformlevels(s)`
- `dnt<:DNF=RelErrorShrink(SoftTH())`: denoise type, see WaveletsExt.jl documentation for detailed description of available denoising functions
- `smooth::Symbol=:regular`: smoothing method:
    - `:regular`: threshold all coefficients
    - `:undersmooth`: do not threshold the lowest-frequency subspace

# Returns

- `NeuroAnalyzer.NEURO`: new object with denoised channels

# See also

[`denoise_dwd!`](@ref), [`denoise_dwd(::AbstractArray)`](@ref)
"""
function denoise_dwd(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    wt::T1 = wavelet(WT.haar),
    l::Int64 = 0,
    dnt::T2 = RelErrorShrink(SoftTH()),
    smooth::Symbol = :regular,
)::NeuroAnalyzer.NEURO where {T1 <: DiscreteWavelet, T2 <: DNFT}

    if l == 0
        l = maxtransformlevels(obj.data[1, :, 1])
        _info("Calculating DWD using maximum level: $l")
    end

    ch = get_channel(obj, ch = ch)
    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = @views denoise_dwd(obj.data[ch, :, :], wt = wt, l = l, dnt = dnt, smooth = smooth)
    push!(obj_new.history, "denoise_dwd(OBJ, ch=$ch, wt=$wt, l=$l, dnt=$dnt, smooth=$smooth))")

    return obj_new

end

"""
    denoise_dwd!(obj; <keyword arguments>)

Perform denoising using discrete wavelet decomposition (DWD) in-place on selected channels of a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; modified in-place
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `wt<:DiscreteWavelet`: discrete wavelet, e.g. `wt = wavelet(WT.haar)`; see Wavelets.jl documentation for the list of available wavelets
- `l::Int64=0`: number of levels, default maximum number of levels available or total transformation
- `dnt<:DNF=RelErrorShrink(SoftTH())`: denoise type, see WaveletsExt.jl documentation for detailed description of available denoising functions
- `smooth::Symbol=:regular`: the smoothing method:
    - `:regular`: smoothing thresholds all given coefficients
    - `:undersmooth`: smoothing does not threshold the lowest frequency subspace node of the wavelet transform

# Returns

- `Nothing`

# See also

[`denoise_dwd`](@ref)
"""
function denoise_dwd!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    wt::T1 = wavelet(WT.haar),
    l::Int64 = 0,
    dnt::T2 = RelErrorShrink(SoftTH()),
    smooth::Symbol = :regular,
)::Nothing where {T1 <: DiscreteWavelet, T2 <: DNFT}

    obj_new = denoise_dwd(obj, ch = ch, wt = wt, l = l, dnt = dnt, smooth = smooth)
    obj.data = obj_new.data
    obj.history = obj_new.history

    return nothing

end
