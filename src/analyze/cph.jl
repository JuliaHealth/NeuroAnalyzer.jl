export cph

"""
    cph(s1, s2; <keyword arguments>)

Computes the instantaneous phase of the cross-power spectrum between two signals using multi-taper estimation.

# Arguments

- `s1::AbstractVector`: signal vector
- `s2::AbstractVector`: signal vector
- `fs::Int64`: sampling rate

# Returns

Named tuple containing:

- `ph::Vector{Float64}`: cross-power spectrum phase in radians
- `f::Vector{Float64}`: cross-power spectrum frequencies
"""
function cph(s1::AbstractVector, s2::AbstractVector; fs::Int64)::@NamedTuple{ph::Vector{Float64}, f::Vector{Float64}}

    @assert fs >= 1 "fs must be ≥ 1."
    @assert length(s1) == length(s2) "s1 and s2 must have the same length."

    # stack signals as rows and compute the multi-taper cross-power spectrum
    # mt_cross_power_spectra returns a complex-valued (channels × channels × freq) object
    ps = mt_cross_power_spectra(hcat(s1, s2)', fs = fs)

    # the phase angle of a complex number z is angle(z) = atan(imag(z), real(z))
    # angle() must be applied to the full complex value, not just its imaginary part
    ph = DSP.angle.(ps.power[1, 2, :])

    # get the vector of frequencies
    f = Vector(ps.freq)

    return (ph = ph, f = f)

end

"""
    cph(s; <keyword arguments>)

Computes the instantaneous phase of the cross-power spectrum between all channel pairs using multi-taper estimation.

# Arguments

- `s::AbstractArray`: signal array (channels × samples × epochs)
- `fs::Int64`: sampling rate

# Returns

Named tuple containing:

- `ph::Array{Float64, 4}`: cross-power spectrum phase (in radians), shape `(channels, channels, frequencies, epochs)`
- `f::Vector{Float64}`: cross-power spectrum frequencies
"""
function cph(s::AbstractArray; fs::Int64)::@NamedTuple{ph::Array{Float64, 4}, f::Vector{Float64}}

    # validate that the input is a proper 3-D array (channels × samples × epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # pilot call to determine the frequency vector length — uses the first
    # channel pair
    # f is independent of signal values.
    cph_data = cph(@view(s[1, :, 1]), @view(s[1, :, 1]), fs = fs)
    f = cph_data.f

    # pre-allocate
    ph = zeros(ch_n, ch_n, length(f), ep_n)

    # initialize progress bar
    progbar = Progress(ep_n * ch_n, dt = 1, barlen = 20, color = :white, enabled = progress_bar)

    # calculate over channel and epochs
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx1, ep_idx = idx[1], idx[2]
        for ch_idx2 in 1:ch_idx1
            # @view avoids copying the (samples,) slices per thread.
            ph[ch_idx1, ch_idx2, :, ep_idx], _ = cph(
                @view(s[ch_idx1, :, ep_idx]),
                @view(s[ch_idx2, :, ep_idx]),
                fs = fs,
            )
        end
        progress_bar && next!(progbar)
    end

    # mirror the lower triangle to the upper triangle to produce the full symmetric matrix
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((length(f), ep_n))
        f_idx, ep_idx = idx[1], idx[2]
        for ch_idx1 in 1:(ch_n - 1)
            for ch_idx2 in (ch_idx1 + 1):ch_n
                ph[ch_idx1, ch_idx2, f_idx, ep_idx] = ph[ch_idx2, ch_idx1, f_idx, ep_idx]
            end
        end
    end

    return (ph = ph, f = f)

end

"""
    cph(s1, s2; <keyword arguments>)

Calculate cross-phases between paired channels of two arrays.

# Arguments

- `s1::AbstractArray`: signal array (channels × samples × epochs)
- `s2::AbstractArray`: signal array (channels × samples × epochs)
- `fs::Int64`: sampling rate

# Returns

Named tuple containing:

- `ph::Array{Float64, 3}`: cross-power spectrum phase in radians, shape `(channels, frequencies, epochs)`
- `f::Vector{Float64}`: cross-power spectrum frequencies
"""
function cph(
    s1::AbstractArray,
    s2::AbstractArray;
    fs::Int64,
)::@NamedTuple{ph::Array{Float64, 3}, f::Vector{Float64}}

    # validate that the input is a proper 3-D array (channels × samples × epochs)
    @assert size(s1) == size(s2) "s1 and s2 must have the same size."
    _chk3d(s1)
    _chk3d(s2)

    # number of channels
    ch_n = size(s1, 1)
    # number of epochs
    ep_n = size(s1, 3)

    # pilot call to determine the frequency vector length — uses the first channel pair
    # f is independent of signal values.
    cph_data = cph(@view(s1[1, :, 1]), @view(s2[1, :, 1]), fs = fs)
    f = cph_data.f

    # pre-allocate output
    ph = zeros(ch_n, ch_n, length(f), ep_n)

    # calculate over channel and epochs
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        ph[ch_idx, :, ep_idx] = cph(
            @view(s1[ch_idx, :, ep_idx]),
            @view(s2[ch_idx, :, ep_idx]),
            fs = fs,
        ).ph
    end

    return (ph = ph, f = f)

end

"""
    cph(obj; <keyword arguments>)

Calculate cross-phases between all channel pairs.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)

# Returns

Named tuple containing:

- `ph::Array{Float64, 4}`: cross-power spectrum phase in radians, shape `(channels, channels, frequencies, epochs)`
- `f::Vector{Float64}`: cross-power spectrum frequencies
"""
function cph(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex}
)::@NamedTuple{ph::Array{Float64, 4}, f::Vector{Float64}}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")

    cph_data = cph(@view(obj.data[ch, :, :]), fs = sr(obj))

    return cph_data

end

"""
    cph(obj1, obj2; <keyword arguments>)

Calculate cross-phases between paired channels of two objects.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `ch1::Union{String, Vector{String}, Regex}: channel name(s)
- `ch2::Union{String, Vector{String}, Regex}: channel name(s)
- `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1))`: epoch number(s)
- `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2))`: epoch number(s)

# Returns

Named tuple containing:

- `ph::Array{Float64, 3}`: cross-power spectrum phase in radians, shape `(channels, frequencies, epochs)`
- `f::Vector{Float64}`: cross-power spectrum frequencies
"""
function cph(
    obj1::NeuroAnalyzer.NEURO,
    obj2::NeuroAnalyzer.NEURO;
    ch1::Union{String, Vector{String}, Regex},
    ch2::Union{String, Vector{String}, Regex},
    ep1::Union{Int64, Vector{Int64}, AbstractRange} = _c(nepochs(obj1)),
    ep2::Union{Int64, Vector{Int64}, AbstractRange} = _c(nepochs(obj2)),
)::@NamedTuple{ph::Array{Float64, 3}, f::Vector{Float64}}

    # validate objects
    @assert sr(obj1) == sr(obj2) "OBJ1 and OBJ2 must have the same sampling rate."

    # resolve channel names to integer indices, optionally skipping bad channels
    ch1 = exclude_bads ? get_channel(obj1, ch = ch1, exclude = "bad") : get_channel(obj1, ch = ch1, exclude = "")
    ch2 = exclude_bads ? get_channel(obj2, ch = ch2, exclude = "bad") : get_channel(obj2, ch = ch2, exclude = "")
    @assert length(ch1) == length(ch2) "Lengths of ch1 ($(length(ch1))) and ch2 ($(length(ch2))) must be equal."

    # validate epoch indices and ensure both objects have matching epoch structure
    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    # normalize scalar epoch arguments to vectors so indexing is uniform
    isa(ep1, Int64) && (ep1 = [ep1])
    isa(ep2, Int64) && (ep2 = [ep2])
    @assert length(ep1) == length(ep2) "Lengths of ep1 ($(length(ep1))) and ep2 ($(length(ep2))) must be equal."
    @assert epoch_len(obj1) == epoch_len(obj2) "OBJ1 and OBJ2 must have the same epoch lengths."

    cph_data = cph(
        @view(obj1.data[ch1, :, ep1]),
        @view(obj2.data[ch2, :, ep2]),
        fs = sr(obj1),
    )

    return cph_data

end
