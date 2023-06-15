export erop

"""
    erop(obj; <keyword arguments>)

Calculate ERO (Event-Related Oscillations) power-spectrum. If `obj` is ERP, `ero()` returns two epochs: ERP power-spectrum (`ero_s[:, :, 1]`) and averaged power-spectra of all ERP epochs (`ero_s[:, :, 2]`). Otherwise, `ero()` returns averaged power-spectra of all `obj` epochs (`ero_s[:, :, 1]`)

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Int64`: channel to analyze
- `method::Symbol=:standard`: method of calculating power-spectrum:
    - `:standard`: standard
    - `:mt`: multi-tapered periodogram
    - `:mw`: Morlet wavelet convolution
- `nt::Int64=8`: number of Slepian tapers
- `pad::Int64=0`: number of zeros to add
- `frq_lim::Tuple{Real, Real}=(0, sr(obj) ÷ 2)`: frequency limits
- `frq_n::Int64=_tlength(frq_lim)`: number of frequencies
- `norm::Bool=true`: normalize powers to dB
- `frq::Symbol=:log`: linear (`:lin`) or logarithmic (`:log`) frequencies
- `ncyc::Union{Int64, Tuple{Int64, Int64}}=6`: number of cycles for Morlet wavelet, for tuple a variable number o cycles is used per frequency: `ncyc = logspace(log10(ncyc[1]), log10(ncyc[2]), frq_n)` for `frq = :log` or `ncyc = linspace(ncyc[1], ncyc[2], frq_n)` for `frq = :lin`

# Returns

Named tuple containing:
- `ero_p::Array{Float64, 3}`: powers
- `ero_f::Vector{Float64}`: frequencies
"""
function erop(obj::NeuroAnalyzer.NEURO; ch::Int64, nt::Int64=8, pad::Int64=0, frq_lim::Tuple{Real, Real}=(0, sr(obj) ÷ 2), frq_n::Int64=_tlength(frq_lim), method::Symbol=:standard, norm::Bool=true, frq::Symbol=:log, ncyc::Union{Int64, Tuple{Int64, Int64}}=6)

    _check_channels(obj, ch)
    _check_var(method, [:standard, :mt, :mw], "method")

    frq_lim = tuple_order(frq_lim)
    frq_lim[1] < 0 && throw(ArgumentError("Lower frequency bound must be ≥ 0."))
    frq_lim[2] > sr(obj) ÷ 2 && throw(ArgumentError("Upper frequency bound must be ≤ $(sr(obj) ÷ 2)."))
    frq_n < 2 && throw(ArgumentError("frq_n frequency bound must be ≥ 2."))
    frq_lim[1] == 0 && (frq_lim = (0.1, frq_lim[2]))

    if method === :standard
        ero_p, ero_f = psd(obj, ch=ch, norm=norm, mt=false, nt=nt)
    elseif method === :mt
        ero_p, ero_f = psd(obj, ch=ch, norm=norm, mt=true, nt=nt)
    else
        ero_p, ero_f = psd_mw(obj, ch=ch, pad=pad, norm=norm, frq_lim=frq_lim, frq_n=frq_n, frq=frq, ncyc=ncyc)
    end

    ero_p = ero_p[1, :, :]

    if method in [:standard, :mt]
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
