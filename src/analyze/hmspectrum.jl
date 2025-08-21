export hmspectrum

"""
Hilbert marginal spectrum

"""
function hmspectrum(obj; ch)
    p, _, t = NeuroAnalyzer.spectrogram(obj, ch=ch, method=:hht, db=false)
    p = dropdims(sum(p, dims=1), dims=3)
    # Plots.plot(t, p[1, :, 1])
    return p, t
end
