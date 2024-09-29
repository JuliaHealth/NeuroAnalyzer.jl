export erop

"""
    erop(obj; <keyword arguments>)

Calculate ERO (Event-Related Oscillations) power-spectrum. If `obj` is ERP, `erop()` returns two epochs: ERP power-spectrum (`s[:, 1]`) and averaged power-spectra of all ERP epochs (`s[:, 2]`). Otherwise, `erop()` returns averaged power-spectra of all `obj` epochs (`s[:, 1]`)

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::String`: channel to analyze
- `method::Symbol=:welch`: method of calculating power-spectrum:
    - `:welch`: Welch's periodogram
    - `:stft`: short time Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:fft`: Fast Fourier transform
    - `:mw`: Morlet wavelet convolution
    - `:gh`: Gaussian and Hilbert transform
    - `:cwt`: continuous wavelet transformation
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window
- `db::Bool=true`: normalize powers to dB
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(sr(obj) / 2)`
- `gw::Real=5`: Gaussian width in Hz
- `wt::T where {T <: CWT}=wavelet(Morlet(2π), β=32, Q=128)`: continuous wavelet, see ContinuousWavelets.jl documentation for the list of available wavelets

# Returns

Named tuple containing:
- `p::Matrix{Float64}`: powers
- `f::Vector{Float64}`: frequencies
"""
function erop(obj::NeuroAnalyzer.NEURO; ch::String, nt::Int64=7, wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, method::Symbol=:welch, db::Bool=true, ncyc::Union{Int64, Tuple{Int64, Int64}}=32, gw::Real=5, wt::T=wavelet(Morlet(2π), β=32, Q=128))::NamedTuple{p::Matrix{Float64}, f::Vector{Float64}} where {T <: CWT}

    p, f = psd(obj, ch=ch, db=db, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc, gw=gw, wt=wt)

    p = p[1, :, :]

    if datatype(obj) == "erp"
        p = cat(p[:, 1], mean(p, dims=2), dims=2)[:, :]
    else
        p = mean(p, dims=2)[:, :]
    end

    return (p=p, f=f)

end
