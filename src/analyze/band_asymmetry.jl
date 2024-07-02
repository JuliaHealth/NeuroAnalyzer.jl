export band_asymmetry

"""
    band_asymmetry(obj; ch1, ch2, frq_lim, method, nt, wlen, woverlap, ncyc)

Calculate band asymmetry: ln(channel 1 band power) - ln(channel 2 band power).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch1::Union{Int64, Vector{Int64}, <:AbstractRange}`: index of channels, e.g. left frontal channels
- `ch2::Union{Int64, Vector{Int64}, <:AbstractRange}`: index of channels, e.g. right frontal channels
- `frq_lim::Tuple{Real, Real}`: lower and upper frequency bounds
- `method::Symbol=:welch`: method used to calculate PSD:
    - `:welch`: Welch's periodogram
    - `:fft`: fast Fourier transform
    - `:mt`: multi-tapered periodogram
    - `:stft`: short time Fourier transform
- `nt::Int64=7`: number of Slepian tapers
- `wlen::Int64=sr(obj)`: window length (in samples), default is 1 second
- `woverlap::Int64=round(Int64, wlen * 0.97)`: window overlap (in samples)
- `w::Bool=true`: if true, apply Hanning window
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=32`: number of cycles for Morlet wavelet, for tuple a variable number of cycles is used per frequency: `ncyc=linspace(ncyc[1], ncyc[2], frq_n)`, where `frq_n` is the length of `0:(sr(obj) / 2)`

# Returns

- `ba::Float64`: band asymmetry
- `ba_norm::Float64`: normalized band asymmetry
"""
function band_asymmetry(obj::NeuroAnalyzer.NEURO; ch1::Union{Int64, Vector{Int64}, <:AbstractRange}, ch2::Union{Int64, Vector{Int64}, <:AbstractRange}, frq_lim::Tuple{Real, Real}, method::Symbol=:welch, nt::Int64=7, wlen::Int64=sr(obj), woverlap::Int64=round(Int64, wlen * 0.97), w::Bool=true, ncyc::Union{Int64, Tuple{Int64, Int64}}=32)

    _check_channels(obj, ch1)
    _check_channels(obj, ch2)
    length(ch1) == 1 && (ch1 = [ch1])
    length(ch2) == 1 && (ch2 = [ch2])

    bp1 = @views band_power(obj.data[ch1, :, :], fs=sr(obj), frq_lim=frq_lim, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc)
    bp2 = @views band_power(obj.data[ch2, :, :], fs=sr(obj), frq_lim=frq_lim, method=method, nt=nt, wlen=wlen, woverlap=woverlap, w=w, ncyc=ncyc)

    ba = log(mean(bp1)) - log(mean(bp2))
    ba_norm = (mean(bp1) - mean(bp2)) / (mean(bp1) + mean(bp2))

    return (ba=ba, ba_norm=ba_norm)

end
