export erop

"""
    erop(obj; <keyword arguments>)

Calculate ERO (Event-Related Oscillations) power-spectrum. If `obj` is ERP, `ero()` returns two epochs: ERP power-spectrum (`ero_s[:, :, 1]`) and averaged power-spectra of all ERP epochs (`ero_s[:, :, 2]`). Otherwise, `ero()` returns averaged power-spectra of all `obj` epochs (`ero_s[:, :, 1]`)

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Int64`: channel to analyze
- `method::Symbol=:welch`: method of calculating power-spectrum:
    - `:welch`: Welch's periodogram
    - `:stft`: short time Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:fft`: Fast Fourier transform
    - `:mw`: Morlet wavelet convolution
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window
- `db::Bool=true`: normalize powers to dB
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(sr(obj) / 2)`

# Returns

Named tuple containing:
- `ero_p::Array{Float64, 3}`: powers
- `ero_f::Vector{Float64}`: frequencies
"""
function erop(obj::NeuroAnalyzer.NEURO; ch::Int64, nt::Int64=7, wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, method::Symbol=:welch, db::Bool=true, ncyc::Union{Int64, Tuple{Int64, Int64}}=32)

    _check_channels(obj, ch)
    _check_var(method, [:welch, :stft, :fft, :mt, :mw], "method")

    ero_p, ero_f = psd(obj, ch=ch, db=db, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc)

    ero_p = ero_p[1, :, :]

    if datatype(obj) == "erp"
        ero_p = cat(ero_p[:, 1], mean(ero_p, dims=2), dims=2)[:, :]
    else
        ero_p = mean(ero_p, dims=2)[:, :]
    end

    return (ero_p=ero_p, ero_f=ero_f)

end
