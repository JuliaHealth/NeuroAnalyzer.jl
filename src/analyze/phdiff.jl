export phdiff

"""
    phdiff(s1, s2; pad, h)

Calculate phase difference between signals.

# Arguments

- `s1::AbstractVector`
- `s2::AbstractVector`
- `pad::Int64=0`: number of zeros to add
- `h::Bool=false`: use FFT or Hilbert transformation (if h=true)

# Returns

Named tuple containing:
- `phd::Vector{Float64}`: phase differences in radians
"""
function phdiff(s1::AbstractVector, s2::AbstractVector; pad::Int64=0, h::Bool=false)

    if h
        _, _, _, ph1 = hspectrum(s1, pad=pad)
        _, _, _, ph2 = hspectrum(s2, pad=pad)
    else
        _, _, _, ph1 = spectrum(s1, pad=pad)
        _, _, _, ph2 = spectrum(s2, pad=pad)
    end

    phd = ph1 - ph2

    return phd

end

"""
    phdiff(s; ch, pad, h)

Calculate phase difference between channels and mean phase of reference `ch`.

# Arguments

- `s::AbstractArray`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=size(s, 1)`: index of reference channels, default is all  channels except the analyzed one
- `avg::Symbol=:phase`: method of averaging:
    - `:phase`: phase is calculated for each reference channel separately and then averaged
    - `:signal`: signals are averaged prior to phase calculation
- `pad::Int64=0`: pad signals with 0s
- `h::Bool=false`: use FFT or Hilbert transformation

# Returns

- `phd::Array{Float64, 3}`
"""
function phdiff(s::AbstractArray; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=1:size(s, 1), avg::Symbol=:phase, pad::Int64=0, h::Bool=false)

    _check_var(avg, [:phase, :signal], "avg")

    _check_channels(s, ch)

    ch_n = size(s, 1)
    if h
        ep_len = size(s, 2)
    else
        ep_len = div(size(s, 2) + pad, 2) + 1
    end
    ep_n = size(s, 3)

    if h
        phd = similar(s)
    else
        phd = zeros(ch_n, ep_len, ep_n)
    end

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            if avg === :phase
                ref_channels = setdiff(ch, ch_idx)
                ph_ref = zeros(length(ref_channels), ep_len)

                for ref_idx in eachindex(ref_channels)
                    if h
                        _, _, _, ph = @views hspectrum(s[ref_channels[ref_idx], :, ep_idx], pad=pad)
                    else
                        _, _, _, ph = @views spectrum(s[ref_channels[ref_idx], :, ep_idx], pad=pad)
                    end
                    ph_ref[ref_idx, :] = ph
                end

                ph_ref = vec(mean(ph_ref, dims=1))

                if h
                    _, _, _, ph = @views hspectrum(s[ch[ch_idx], :, ep_idx], pad=pad)
                else
                    _, _, _, ph = @views spectrum(s[ch[ch_idx], :, ep_idx], pad=pad)
                end

                phd[ch_idx, :, ep_idx] = ph - ph_ref

            else
                ref_channels = setdiff(ch, ch_idx)

                signal_m = @views vec(mean(s[ref_channels, :, ep_idx], dims=1))

                phd[ch_idx, :, ep_idx] = @views phdiff(s[ch[ch_idx], :, ep_idx], signal_m)
            end
        end
    end

    return phd

end

"""
    phdiff(obj; ch, pad, h)

Calculate phase difference between channels and mean phase of reference `ch`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of reference channels, default is all signal channels except the analyzed one
- `avg::Symbol=:phase`: method of averaging:
    - `:phase`: phase is calculated for each reference channel separately and then averaged
    - `:signal`: signals are averaged prior to phase calculation
- `pad::Int64=0`: pad signals with 0s
- `h::Bool=false`: use FFT or Hilbert transformation

# Returns

- `phd::Array{Float64, 3}`
"""
function phdiff(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), avg::Symbol=:phase, pad::Int64=0, h::Bool=false)

    _check_channels(obj, ch)
    length(ch) == 1 && (ch = [ch])

    phd = @views phdiff(obj.data[ch, :, :], avg=avg, pad=pad, h=h)

    return phd

end
