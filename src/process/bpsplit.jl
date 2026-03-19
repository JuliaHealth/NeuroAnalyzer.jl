export bpsplit

"""
    bpsplit(obj; <keyword arguments>)

Split a signal into standard EEG frequency bands using FIR band-pass filters.

One filter is designed per band and applied to every selected channel and epoch.

The following bands are extracted: `:delta`, `:theta`, `:alpha`, `:alpha_lower`, `:alpha_higher`, `:beta`, `:beta_lower`, `:beta_higher`, `:gamma`, `:gamma_1`, `:gamma_2`, `:gamma_lower`, `:gamma_higher`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `order::Int64=91`: number of taps for the FIR band-pass filter; must be odd for band-pass filters
- `w::Union{Nothing, AbstractVector, <:Real}=nothing`: window for the FIR filter; `nothing` applies the default Hamming window of length `order`

# Returns

Named tuple:

- `s::Array{Float64, 4}`: band-split signal of shape `(n_bands, n_channels, epoch_len, n_epochs)`
- `bn::Vector{Symbol}`: band name symbols (in the order used in `s`)
- `bf::Vector{Tuple{Real, Real}}`: frequency limits `(f_low, f_high)` in Hz for each band, in the same order as `bn`

# Throws

- `ArgumentError`: if `order` is even (band-pass FIR requires odd tap count), or if any band frequency exceeds the Nyquist frequency for `obj`

# See also

[`filter_create`](@ref), [`filter_apply`](@ref), [`band_frq`](@ref)
"""
function bpsplit(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    order::Int64 = 91,
    w::Union{Nothing, AbstractVector, <:Real} = nothing
)::@NamedTuple{
    s::Array{Float64, 4},
    bn::Vector{Symbol},
    bf::Vector{Tuple{Real, Real}}
}

    bn = [
        :delta,
        :theta,
        :alpha,
        :alpha_lower,
        :alpha_higher,
        :beta,
        :beta_lower,
        :beta_higher,
        :gamma,
        :gamma_1,
        :gamma_2,
        :gamma_lower,
        :gamma_higher,
    ]

    # resolve channel names to integer indices
    ch = get_channel(obj, ch=ch)

    # number of channels
    ch_n = length(ch)
    # number of epochs
    ep_n = nepochs(obj)
    # sampling rate
    fs = sr(obj)
    # epoch lengths
    el = epoch_len(obj)

    s  = zeros(length(const_bn), ch_n, el, ep_n)

    # pre-allocate output
    bf = Vector{Tuple{Real, Real}}(undef, length(const_bn))

    # design one filter per band, then apply it to all channels and epochs.
    # the outer band loop is sequential (each band uses a different filter);
    # the inner channel loop is parallelized
    @inbounds for band_idx in eachindex(const_bn)
        band_f = band_frq(obj; band=const_bn[band_idx])
        bf[band_idx] = band_f
        flt = filter_create(
            fs=fs,
            fprototype=:fir,
            ftype=:bp,
            cutoff=band_f,
            order=order,
            w=w,
        )

        # calculate over channel and epochs
        @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
            ch_idx, ep_idx = idx[1], idx[2]
            s[band_idx, ch_idx, :, ep_idx] = filter_apply(
                @view(obj.data[ch[ch_idx], :, ep_idx]),
                flt=flt
            )
        end

    end

    return (s=s, bn=const_bn, bf=bf)

end
