export eros

"""
    eros(obj; <keyword arguments>)

Calculate ERO (Event-Related Oscillations) spectrogram. If `obj` is ERP, `ero()` returns two epochs: ERP spectrogram (`ero_s[:, :, 1]`) and averaged spectrograms of all ERP epochs (`ero_s[:, :, 2]`). Otherwise, `ero()` returns averaged spectrograms of all `obj` epochs (`ero_s[:, :, 1]`)

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Int64`: channel to analyze
- `method::Symbol=:stft`: method of calculating spectrogram:
    - `:stft`: short-time Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian and Hilbert transform
    - `:cwt`: continuous wavelet transformation
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length, default is 4 seconds
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window
- `pad::Int64=0`: number of zeros to add
- `frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2)`: frequency limits
- `frq_n::Int64=_tlength(frq_lim)`: number of frequencies
- `norm::Bool=true`: normalize powers to dB
- `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
- `gw::Real=5`: Gaussian width in Hz
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number o cycles is used per frequency: `ncyc = logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc = linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`
- `wt<:CWT=wavelet(Morlet(2π), β=32, Q=128)`: continuous wavelet, e.g. `wt = wavelet(Morlet(2π), β=32, Q=128)`, see ContinuousWavelets.jl documentation for the list of available wavelets

# Returns

Named tuple containing:
- `ero_s::Array{Float64, 3}`: spectrogram(s)
- `ero_f::Vector{Float64}`: frequencies
- `ero_t::Vector{Float64}`: time
"""
function eros(obj::NeuroAnalyzer.NEURO; ch::Int64, pad::Int64=0, frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2), frq_n::Int64=_tlength(frq_lim), method::Symbol=:stft, nt::Int64=7, wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, norm::Bool=true, frq::Symbol=:log, gw::Real=5, ncyc::Union{Int64, Tuple{Int64, Int64}}=32, wt::T=wavelet(Morlet(2π), β=32, Q=128)) where {T <: CWT}

    _check_channels(obj, ch)
    _check_var(method, [:stft, :mt, :mw, :gh, :cwt], "method")

    ero_s, ero_f, ero_t = NeuroAnalyzer.spectrogram(obj, ch=ch, method=method, nt=nt, pad=pad, frq_lim=frq_lim, frq_n=frq_n, norm=norm, frq=frq, gw=gw, ncyc=ncyc, wt=wt, wlen=wlen, woverlap=woverlap, w=w)

    ero_s = ero_s[:, :, 1, :]

    if method in [:mt, :stft]
        f1_idx = vsearch(frq_lim[1], ero_f)
        f2_idx = vsearch(frq_lim[2], ero_f)
        ero_f = ero_f[f1_idx:f2_idx]
        ero_s = ero_s[f1_idx:f2_idx, :, :]
    end

    if datatype(obj) == "erp"
        ero_s = cat(ero_s[:, :, 1], mean(ero_s, dims=3), dims=3)
    else
        ero_s = mean(ero_s, dims=3)
    end

    return (ero_s=ero_s, ero_f=ero_f, ero_t=ero_t)

end
