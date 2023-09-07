export erop

"""
    erop(obj; <keyword arguments>)

Calculate ERO (Event-Related Oscillations) power-spectrum. If `obj` is ERP, `ero()` returns two epochs: ERP power-spectrum (`ero_s[:, :, 1]`) and averaged power-spectra of all ERP epochs (`ero_s[:, :, 2]`). Otherwise, `ero()` returns averaged power-spectra of all `obj` epochs (`ero_s[:, :, 1]`)

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Int64`: channel to analyze
- `method::Symbol=:welch`: method of calculating power-spectrum:
    - `:welch`: Welch periodogram
    - `:stft`: short time Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:mw`: Morlet wavelet convolution
- `nt::Int64=8`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window for Welch and STFT
- `pad::Int64=0`: number of zeros to add
- `frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2)`: frequency limits
- `frq_n::Int64=_tlength(frq_lim)`: number of frequencies
- `norm::Bool=true`: normalize powers to dB
- `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet, for tuple a variable number o cycles is used per frequency: `ncyc = logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc = linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

# Returns

Named tuple containing:
- `ero_p::Array{Float64, 3}`: powers
- `ero_f::Vector{Float64}`: frequencies
"""
function erop(obj::NeuroAnalyzer.NEURO; ch::Int64, nt::Int64=8, wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, pad::Int64=0, frq_lim::Tuple{Real, Real}=(0, sr(obj) / 2), frq_n::Int64=_tlength(frq_lim), method::Symbol=:welch, norm::Bool=true, frq::Symbol=:log, ncyc::Union{Int64, Tuple{Int64, Int64}}=6)

    _check_channels(obj, ch)
    _check_var(method, [:welch, :stft, :mt, :mw], "method")

    frq_lim = tuple_order(frq_lim)
    @assert !(frq_lim[1] < 0 || frq_lim[2] < 0 || frq_lim[1] > sr(obj) / 2 || frq_lim[2] > sr(obj) / 2) "frq_lim must be in [0, $(sr(obj) / 2)]."
    frq_lim[1] == 0 && (frq_lim = (0.1, frq_lim[2]))

    if method === :welch
        ero_p, ero_f = psd(obj, ch=ch, norm=norm, nt=nt, wlen=wlen, woverlap=woverlap, w=w)
    elseif method === :stft
        ero_p, ero_f = psd(obj, ch=ch, norm=norm, st=true, nt=nt, wlen=wlen, woverlap=woverlap, w=w)
    elseif method === :mt
        ero_p, ero_f = psd(obj, ch=ch, norm=norm, mt=true, nt=nt, wlen=wlen, woverlap=woverlap, w=w)
    else
        ero_p, ero_f = psd_mw(obj, ch=ch, pad=pad, norm=norm, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc)
    end

    ero_p = ero_p[1, :, :]

    if method in [:welch, :stft, :mt]
        f1_idx = vsearch(frq_lim[1], ero_f)
        f2_idx = vsearch(frq_lim[2], ero_f)
        ero_f = ero_f[f1_idx:f2_idx]
        ero_p = ero_p[f1_idx:f2_idx, :]
    end

    if obj.header.recording[:data_type] == "erp"
        ero_p = cat(ero_p[:, 1], mean(ero_p, dims=2), dims=2)
    else
        ero_p = mean(ero_p, dims=2)
    end

    return (ero_p=ero_p, ero_f=ero_f)

end
