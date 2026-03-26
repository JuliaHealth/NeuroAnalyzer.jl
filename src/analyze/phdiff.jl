export phdiff

"""
    phdiff(s1, s2; <keyword arguments>)

Calculate phase difference between two 1-D signal vectors.

# Arguments

- `s1::AbstractVector`: signal vector
- `s2::AbstractVector`: signal vector
- `pad::Int64=0`: number of zeros to append
- `h::Bool=false`: use Hilbert transform, otherwise use Fourier transform

# Returns

Named tuple:

- `Vector{Float64}`: phase differences in radians
"""
function phdiff(
    s1::AbstractVector,
    s2::AbstractVector;
    pad::Int64 = 0,
    h::Bool = false
)::Vector{Float64}

    h1 = h ? NeuroAnalyzer.htransform(s1) : NeuroAnalyzer.ftransform(s1, pad = pad)
    h2 = h ? NeuroAnalyzer.htransform(s2) : NeuroAnalyzer.ftransform(s2, pad = pad)
    ph1 = h1.ph
    ph2 = h2.ph

    phd = ph1 - ph2

    return phd

end

"""
    phdiff(s; <keyword arguments>)

Calculate phase difference between channels and mean phase of reference `ch` for a 3-D signal array.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `ch::Union{Int64, Vector{Int64}}=_c(size(s, 1))`: index of reference channels, default is all  channels except the analyzed one
- `avg::Symbol=:phase`: method of averaging:
    - `:phase`: phase is calculated for each reference channel separately and then averaged
    - `:signal`: signals are averaged prior to phase calculation
- `pad::Int64=0`: number of zeros to append
- `h::Bool=false`: use Hilbert transform, otherwise use Fourier transform

# Returns

- `Array{Float64, 3}`
"""
function phdiff(
    s::AbstractArray;
    ch::Union{Int64, Vector{Int64}} = _c(size(s, 1)),
    avg::Symbol = :phase,
    pad::Int64 = 0,
    h::Bool = false
)::Array{Float64, 3}

    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # validate
    _check_var(avg, [:phase, :signal], "avg")
    _check_channels(s, ch)

    # number of channels
    ch_n = size(s, 1)
    # epoch length
    ep_len = h ? size(s, 2) : div(size(s, 2) + pad, 2) + 1
    # number of epochs
    ep_n = size(s, 3)

    # pre-allocate output
    phd = zeros(ch_n, ep_len, ep_n)

    # calculate over channel and epochs
    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        if avg === :phase

            ref_channels = setdiff(ch, ch_idx)
            ph_ref = zeros(length(ref_channels), ep_len)
            for ref_idx in eachindex(ref_channels)
                h_data = h ? NeuroAnalyzer.htransform(@view(s[ref_channels[ref_idx], :, ep_idx])) :
                             NeuroAnalyzer.ftransform(@view(s[ref_channels[ref_idx], :, ep_idx]), pad = pad)
                ph_ref[ref_idx, :] = h_data.ph
            end
            ph_ref = vec(mean(ph_ref, dims = 1))

            h_data = h ? NeuroAnalyzer.htransform(@view(s[ch[ch_idx], :, ep_idx])) :
                         NeuroAnalyzer.ftransform(@view(s[ch[ch_idx], :, ep_idx]), pad = pad)
            phd[ch_idx, :, ep_idx] = h_data.ph - ph_ref

        elseif avg === :signal

            ref_channels = setdiff(ch, ch_idx)
            signal_m = vec(mean(@view(s[ref_channels, :, ep_idx]), dims = 1))
            phd[ch_idx, :, ep_idx] = phdiff(@view(s[ch[ch_idx], :, ep_idx]), signal_m)

        end
    end

    return phd

end

"""
    phdiff(obj; <keyword arguments>)

Calculate phase difference between channels and mean phase of reference `ch`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: index of reference channels
- `avg::Symbol=:phase`: method of averaging:
    - `:phase`: phase is calculated for each reference channel separately and then averaged
    - `:signal`: signals are averaged prior to phase calculation
- `pad::Int64=0`: number of zeros to append
- `h::Bool=false`: use Hilbert transform, otherwise use Fourier transform

# Returns

- `Array{Float64, 3}`
"""
function phdiff(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    avg::Symbol = :phase,
    pad::Int64 = 0,
    h::Bool = false
)::Array{Float64, 3}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")

    return phdiff(
        @view(obj.data[ch, :, :]),
        avg = avg,
        pad = pad,
        h = h
    )

end
