using Pkg
Pkg.add("Wavelets")
Pkg.add("ContinuousWavelets")
Pkg.add("Artifacts")

using Wavelets
using ContinuousWavelets
using Artifacts

#=
using BenchmarkTools

BenchmarkTools.DEFAULT_PARAMETERS.samples = 1
b = @benchmarkable import_edf(joinpath(testfiles_path, "eeg-test-edf.edf")) evals=20
run(b)
=#

@info "Reporting system data"
println("# System info")
println()
if Sys.isunix() || Sys.isapple()
    print("OS: $(read(`uname -a`, String))")
elseif Sys.iswindows()
    print("OS: Windows $(Sys.windows_version())")
end
println("CPU: $(Sys.cpu_info()[1].model) $(length(Sys.cpu_info()) ÷ 2) cores ($(round(Sys.cpu_info()[1].speed / 1024, digits=2)) GHz)")
println("RAM: $(round((Int64(Sys.free_memory())) / 1024^3, digits=1)) GB free / $(round((Int64(Sys.total_memory())) / 1024^3, digits=1)) GB")
println()
@info "Loading NeuroAnalyzer.jl"
println("# Loading NeuroAnalyzer.jl")
println()
# @time_imports using NeuroAnalyzer
@time using NeuroAnalyzer
println()
println("# NeuroAnalyzer info")
println()
na_info()

global testfiles_path = joinpath(artifact"NeuroAnalyzer_test-files", "neuroanalyzer-test-files")

@info "Benchmarking: IO"
println()
println("# IO")
println()

print(rpad("Import EDF", 36))
eeg = import_edf(joinpath(testfiles_path, "eeg-test-edf.edf"))
@time import_edf(joinpath(testfiles_path, "eeg-test-edf.edf"))
print(rpad("Import EDF+", 36))
import_edf(joinpath(testfiles_path, "eeg-test-edfplus.edf"))
@time import_edf(joinpath(testfiles_path, "eeg-test-edfplus.edf"))
print(rpad("Import BDF+", 36))
import_bdf(joinpath(testfiles_path, "eeg-test-bdf.bdf"))
@time import_bdf(joinpath(testfiles_path, "eeg-test-bdf.bdf"))
print(rpad("Import GDF 1.25", 36))
import_gdf(joinpath(testfiles_path, "eeg-test-gdf_1.25.gdf"))
@time import_gdf(joinpath(testfiles_path, "eeg-test-gdf_1.25.gdf"))
print(rpad("Import GDF 2.20", 36))
import_gdf(joinpath(testfiles_path, "eeg-test-gdf_2.20.gdf"))
@time import_gdf(joinpath(testfiles_path, "eeg-test-gdf_2.20.gdf"))
print(rpad("Import Digitrack", 36))
import_digitrack(joinpath(testfiles_path, "eeg-test-digitrack.txt"))
@time import_digitrack(joinpath(testfiles_path, "eeg-test-digitrack.txt"))
print(rpad("Import BrainVision", 36))
import_bv(joinpath(testfiles_path, "eeg-test-bv.vhdr"))
@time import_bv(joinpath(testfiles_path, "eeg-test-bv.vhdr"))
print(rpad("Import CSV ch×t", 36))
import_csv(joinpath(testfiles_path, "eeg-test_chxt.csv.gz"))
@time import_csv(joinpath(testfiles_path, "eeg-test_chxt.csv.gz"))
print(rpad("Import CSV t×ch", 36))
import_csv(joinpath(testfiles_path, "eeg-test_txch.csv.gz"))
@time import_csv(joinpath(testfiles_path, "eeg-test_txch.csv.gz"))
print(rpad("Import EEGLAB dataset", 36))
import_set(joinpath(testfiles_path, "eeg-test-eeglab.set"))
@time import_set(joinpath(testfiles_path, "eeg-test-eeglab.set"))
print(rpad("Import SNIRF", 36))
import_snirf(joinpath(testfiles_path, "fnirs-test-snirf.snirf"))
@time import_snirf(joinpath(testfiles_path, "fnirs-test-snirf.snirf"))
print(rpad("Import NIRS", 36))
import_nirs(joinpath(testfiles_path, "fnirs-test-nirs.nirs"))
@time import_nirs(joinpath(testfiles_path, "fnirs-test-nirs.nirs"))
print(rpad("Import NIRX", 36))
import_nirx(joinpath(testfiles_path, "nirx/NIRS-2020-08-18_001.hdr"))
@time import_nirx(joinpath(testfiles_path, "nirx/NIRS-2020-08-18_001.hdr"))
print(rpad("Import FIFF", 36))
import_fiff(joinpath(testfiles_path, "meg-test-fiff.fif"))
@time import_fiff(joinpath(testfiles_path, "meg-test-fiff.fif"))
print(rpad("Save HDF5", 36))
tmp = tempname() * ".hdf"
NeuroAnalyzer.save(eeg, file_name=tmp)
tmp = tempname() * ".hdf"
@time NeuroAnalyzer.save(eeg, file_name=tmp)
print(rpad("Load HDF5", 36))
NeuroAnalyzer.load(tmp)
@time NeuroAnalyzer.load(tmp)

eeg = import_edf(joinpath(testfiles_path, "eeg-test-edf.edf"))
NeuroAnalyzer.filter!(eeg, ch=labels(eeg), fprototype=:fir, ftype=:lp, cutoff=40, order=8)
NeuroAnalyzer.filter!(eeg, ch=labels(eeg), fprototype=:fir, ftype=:hp, cutoff=1, order=8)
e10 = epoch(eeg, ep_len=10)
keep_epoch!(e10, ep=1:20)
load_locs!(e10, file_name=joinpath(NeuroAnalyzer.PATH, "locs", "standard-10-20-cap19-elmiko.ced"))
n = import_nirs(joinpath(testfiles_path, "fnirs-test-nirs.nirs"))

@info "Benchmarking: EDIT"
println()
println("# EDIT")
println()

e10_tmp = deepcopy(e10)
print(rpad("Make epochs", 36))
epoch(e10_tmp, ep_len=10)
@time epoch(e10_tmp, ep_len=10)
e10_tmp = deepcopy(e10)
print(rpad("Delete channel", 36))
delete_channel(e10_tmp, ch="F3")
@time delete_channel(e10_tmp, ch="F3")
e10_tmp = deepcopy(e10)
print(rpad("Keep channel", 36))
keep_channel(e10_tmp, ch="F3")
@time keep_channel(e10_tmp, ch="F3")
e10_tmp = deepcopy(e10)
print(rpad("Delete epoch", 36))
delete_epoch(e10_tmp, ep=1)
@time delete_epoch(e10_tmp, ep=1)
e10_tmp = deepcopy(e10)
print(rpad("Keep epoch", 36))
keep_epoch(e10_tmp, ep=1)
@time keep_epoch(e10_tmp, ep=1)
print(rpad("Virtual channel", 36))
vch(e10, f="fp1 + fp2")
@time vch(e10, f="fp1 + fp2")
print(rpad("Detect bad channels and epochs", 36))
detect_bad(e10, ch=get_channel(e10, type="eeg"))
@time detect_bad(e10, ch=get_channel(e10, type="eeg"))

@info "Benchmarking: PROCESS"
println()
println("# PROCESS")
println()

print(rpad("apply()", 36))
apply(e10, f="mean(obj, dims=1)", ch=["F3", "F4"])
@time apply(e10, f="mean(obj, dims=1)", ch=["F3", "F4"])
print(rpad("Add signal", 36))
x = rand(epoch_len(e10))
add_signal(e10, s=x, ch="all")
@time add_signal(e10, s=x, ch="all")
print(rpad("Average", 36))
average(e10, e10)
@time average(e10, e10)
print(rpad("CBP", 36))
cbp(e10, ch="all", frq=4)
@time cbp(e10, ch="all", frq=4)
print(rpad("FFT denoising", 36))
denoise_fft(e10, ch="all")
@time denoise_fft(e10, ch="all")
print(rpad("Wiener denoising", 36))
denoise_wien(e10, ch="all")
@time denoise_wien(e10, ch="all")
print(rpad("Wavelet denoising: DWT", 36))
denoise_dwt(e10, ch="all", wt=wavelet(WT.haar))
@time denoise_dwt(e10, ch="all", wt=wavelet(WT.haar))
print(rpad("Wavelet denoising: CWT", 36))
denoise_cwt(e10, ch="all", nf=50)
@time denoise_cwt(e10, ch="all", nf=50)
print(rpad("Derivative", 36))
derivative(e10, ch="all")
@time derivative(e10, ch="all")
print(rpad("Detrend: ls", 36))
detrend(e10, ch="all", type=:ls)
@time detrend(e10, ch="all", type=:ls)
print(rpad("Detrend: linear", 36))
detrend(e10, ch="all", type=:linear)
@time detrend(e10, ch="all", type=:linear)
print(rpad("Detrend: constant", 36))
detrend(e10, ch="all", type=:constant)
@time detrend(e10, ch="all", type=:constant)
print(rpad("Detrend: poly", 36))
detrend(e10, ch="all", type=:poly)
@time detrend(e10, ch="all", type=:poly)
print(rpad("Detrend: loess", 36))
detrend(e10, ch="all", type=:loess)
@time detrend(e10, ch="all", type=:loess)
print(rpad("DWT: sdwt", 36))
dw_trans(e10, ch="all", wt=wavelet(WT.haar), type=:sdwt)
@time dw_trans(e10, ch="all", wt=wavelet(WT.haar), type=:sdwt)
print(rpad("DWT: acdwt", 36))
dw_trans(e10, ch="all", wt=wavelet(WT.haar), type=:acdwt)
@time dw_trans(e10, ch="all", wt=wavelet(WT.haar), type=:acdwt)
print(rpad("CWT", 36))
cw_trans(e10, ch="all", wt=wavelet(Morlet(π), β=2))
@time cw_trans(e10, ch="all", wt=wavelet(Morlet(π), β=2))
print(rpad("Split using DWT", 36))
dw_split(e10, ch="all", wt=wavelet(WT.db2), type=:sdwt, n=5)
@time dw_split(e10, ch="all", wt=wavelet(WT.db2), type=:sdwt, n=5)
print(rpad("Split using BP", 36))
bpsplit(e10, ch="all")
@time bpsplit(e10, ch="all")
print(rpad("ERP", 36))
average_epochs(e10)
@time average_epochs(e10)
print(rpad("Frequency convolution", 36))
fconv(e10, ch="all", kernel=generate_morlet(256, 1, 32, complex=true))
@time fconv(e10, ch="all", kernel=generate_morlet(256, 1, 32, complex=true))
print(rpad("Filter: moving average", 36))
NeuroAnalyzer.filter_mavg(e10, ch="all", k=2)
@time NeuroAnalyzer.filter_mavg(e10, ch="all", k=2)
print(rpad("Filter: moving median", 36))
NeuroAnalyzer.filter_mmed(e10, ch="all", k=2)
@time NeuroAnalyzer.filter_mmed(e10, ch="all", k=2)
print(rpad("Filter: poly", 36))
NeuroAnalyzer.filter_poly(e10, ch="all")
@time NeuroAnalyzer.filter_poly(e10, ch="all")
print(rpad("Filter: SG", 36))
NeuroAnalyzer.filter_sg(e10, ch="all")
@time NeuroAnalyzer.filter_sg(e10, ch="all")
print(rpad("Filter: FIR LP", 36))
NeuroAnalyzer.filter(e10, ch="all", fprototype=:fir, ftype=:lp, cutoff=45.0, order=8)
@time NeuroAnalyzer.filter(e10, ch="all", fprototype=:fir, ftype=:lp, cutoff=45.0, order=8)
print(rpad("Filter: FIR HP", 36))
NeuroAnalyzer.filter(e10, ch="all", fprototype=:fir, ftype=:hp, cutoff=0.1, order=8)
@time NeuroAnalyzer.filter(e10, ch="all", fprototype=:fir, ftype=:hp, cutoff=0.1, order=8)
print(rpad("Filter: IIR LP", 36))
NeuroAnalyzer.filter(e10, ch="all", fprototype=:chebyshev1, ftype=:lp, cutoff=45.0, order=8, rs=45)
@time NeuroAnalyzer.filter(e10, ch="all", fprototype=:butterworth, ftype=:lp, cutoff=45.0, order=8, rs=45)
print(rpad("Filter: IIR HP", 36))
NeuroAnalyzer.filter(e10, ch="all", fprototype=:butterworth, ftype=:hp, cutoff=0.1, order=8)
@time NeuroAnalyzer.filter(e10, ch="all", fprototype=:butterworth, ftype=:hp, cutoff=0.1, order=8)
print(rpad("Filter: IIR notch", 36))
NeuroAnalyzer.filter(e10, ch="all", fprototype=:iirnotch, cutoff=50, bw=2)
@time NeuroAnalyzer.filter(e10, ch="all", fprototype=:iirnotch, cutoff=50, bw=2)
print(rpad("Filter: Remez", 36))
NeuroAnalyzer.filter(e10, ch="all", fprototype=:remez, ftype=:lp, order=128, cutoff=20, bw=0.5)
@time NeuroAnalyzer.filter(e10, ch="all", fprototype=:remez, ftype=:lp, order=128, cutoff=20, bw=0.5)
print(rpad("Filter: Gaussian", 36))
NeuroAnalyzer.filter_g(e10, ch="all", f=20)
@time NeuroAnalyzer.filter_g(e10, ch="all", f=20)
print(rpad("Interpolate: LR", 36))
lrinterpolate_channel(e10, ch="Fp1", ep=1)
@time lrinterpolate_channel(e10, ch="Fp1", ep=1)
print(rpad("Interpolate: PL", 36))
plinterpolate_channel(e10, ch="Fp1", ep=1)
@time plinterpolate_channel(e10, ch="Fp1", ep=1)
print(rpad("Normalize", 36))
NeuroAnalyzer.normalize(e10, ch="all", method=:zscore)
@time NeuroAnalyzer.normalize(e10, ch="all", method=:zscore)
print(rpad("Remove DC", 36))
remove_dc(e10, ch="all")
@time remove_dc(e10, ch="all")
print(rpad("Channel referencing", 36))
reference_ce(e10, ch="F3")
@time reference_ce(e10, ch="F3")
print(rpad("A referencing", 36))
reference_a(e10)
@time reference_a(e10)
print(rpad("M referencing", 36))
a1 = labels(e10)[findfirst(isequal("a1"), lowercase.(labels(e10)))]
a2 = labels(e10)[findfirst(isequal("a2"), lowercase.(labels(e10)))]
edit_channel!(e10, ch="A1", field=:label, value="M1")
edit_channel!(e10, ch="A2", field=:label, value="M2")
reference_m(e10)
@time reference_m(e10)
print(rpad("AVG referencing", 36))
reference_avg(e10)
@time reference_avg(e10)
print(rpad("PLAP referencing", 36))
reference_plap(e10)
@time reference_plap(e10)
print(rpad("Custom referencing", 36))
reference_custom(e10)
@time reference_custom(e10)
print(rpad("CSD", 36))
csd(e10)
@time csd(e10)
print(rpad("Standardize", 36))
standardize(e10, ch="all")
@time standardize(e10, ch="all")
print(rpad("Taper", 36))
taper(e10, ch="all", t=generate_window(:hann, epoch_len(e10)))
@time taper(e10, ch="all", t=generate_window(:hann, epoch_len(e10)))
print(rpad("Time convolution", 36))
tconv(e10, ch="all", kernel=[0.2, 0.1, 0.2])
@time tconv(e10, ch="all", kernel=[0.2, 0.1, 0.2])
print(rpad("Wavelet band-pass filtering", 36))
wbp(e10, ch="all", frq=4)
@time wbp(e10, ch="all", frq=4)
print(rpad("Remove pops", 36))
remove_pops(eeg, ch="all")
@time remove_pops(eeg, ch="all")
print(rpad("Remove power line noise: IIR", 36))
eeg = import_edf(joinpath(testfiles_path, "eeg-test-edf.edf"))
e10_tmp = epoch(eeg, ep_len=10)
keep_epoch!(e10_tmp, ep=1)
remove_powerline(e10_tmp, ch="all", pl_frq=50, method=:iir);
@time remove_powerline(e10_tmp, ch="all", pl_frq=50, method=:iir);
print(rpad("Normalize power", 36))
normpower(e10, ch="all");
@time normpower(e10, ch="all");

print(rpad("NIRS: int2od", 36))
n_tmp = intensity2od(n)
@time intensity2od(n)
print(rpad("NIRS: od2conc", 36))
od2conc(n_tmp)
@time od2conc(n_tmp)

@info "Benchmarking: ANALYZE"
println()
println("# ANALYZE")
println()

e10_2 = deepcopy(e10)
e10_2.data .*= 0.75

print(rpad("Amplitude difference", 36))
ampdiff(e10, ch="all")
@time ampdiff(e10, ch="all")
print(rpad("Auto-covariance", 36))
acov(e10, ch="all")
@time acov(e10, ch="all")
print(rpad("Auto-correlation", 36))
acor(e10, ch="all")
@time acor(e10, ch="all")
print(rpad("Band power: welch", 36))
band_power(e10, ch="all", frq_lim=(10, 20))
@time band_power(e10, ch="all", frq_lim=(10, 20))
print(rpad("Band power: fft", 36))
band_power(e10, ch="all", frq_lim=(10, 20), method=:fft)
@time band_power(e10, ch="all", frq_lim=(10, 20), method=:fft)
print(rpad("Band power: stft", 36))
band_power(e10, ch="all", frq_lim=(10, 20), method=:stft)
@time band_power(e10, ch="all", frq_lim=(10, 20), method=:stft)
print(rpad("Band power: mt", 36))
band_power(e10, ch="all", frq_lim=(10, 20), method=:mt)
@time band_power(e10, ch="all", frq_lim=(10, 20), method=:mt)
print(rpad("Band power: mw", 36))
band_power(e10, ch="all", frq_lim=(10, 20), method=:mw)
@time band_power(e10, ch="all", frq_lim=(10, 20), method=:mw)
print(rpad("Band mpower: welch", 36))
band_mpower(e10, ch="all", frq_lim=(10, 20))
@time band_mpower(e10, ch="all", frq_lim=(10, 20))
print(rpad("Band mpower: fft", 36))
band_mpower(e10, ch="all", frq_lim=(10, 20), method=:fft)
@time band_mpower(e10, ch="all", frq_lim=(10, 20), method=:fft)
print(rpad("Band mpower: stft", 36))
band_mpower(e10, ch="all", frq_lim=(10, 20), method=:stft)
@time band_mpower(e10, ch="all", frq_lim=(10, 20), method=:stft)
print(rpad("Band mpower: mt", 36))
band_mpower(e10, ch="all", frq_lim=(10, 20), method=:mt)
@time band_mpower(e10, ch="all", frq_lim=(10, 20), method=:mt)
print(rpad("Band mpower: mw", 36))
band_mpower(e10, ch="all", frq_lim=(10, 20), method=:mw)
@time band_mpower(e10, ch="all", frq_lim=(10, 20), method=:mw)
print(rpad("Correlation matrix", 36))
corm(e10, ch="all")
@time corm(e10, ch="all")
print(rpad("Covariance matrix", 36))
covm(e10, ch="all")
@time covm(e10, ch="all")
print(rpad("Cross-phases 1", 36))
cph(e10, ch="all")
@time cph(e10, ch="all")
print(rpad("Cross-phases 2", 36))
cph(e10, e10, ch1="F3", ch2="F4", ep1=1, ep2=1)
@time cph(e10, e10, ch1="F3", ch2="F4", ep1=1, ep2=1)
print(rpad("DISS", 36))
diss(e10, ch="all")
@time diss(e10, ch="all")
print(rpad("Entropy", 36))
entropy(e10, ch="all")
@time entropy(e10, ch="all")
print(rpad("Negentropy", 36))
negentropy(e10, ch="all")
@time negentropy(e10, ch="all")
print(rpad("Temporal envelope", 36))
tenv(e10, ch="all")
@time tenv(e10, ch="all")
print(rpad("Temporal envelope: mean", 36))
tenv_mean(e10, ch="all", dims=1)
@time tenv_mean(e10, ch="all", dims=1)
print(rpad("Temporal envelope: median", 36))
tenv_median(e10, ch="all", dims=1)
@time tenv_median(e10, ch="all", dims=1)
print(rpad("Power envelope", 36))
penv(e10, ch="all")
@time penv(e10, ch="all")
print(rpad("Power envelope: mean", 36))
penv_mean(e10, ch="all", dims=1)
@time penv_mean(e10, ch="all", dims=1)
print(rpad("Power envelope: median", 36))
penv_median(e10, ch="all", dims=1)
@time penv_median(e10, ch="all", dims=1)
print(rpad("Spectral envelope", 36))
senv(e10, ch="all")
@time senv(e10, ch="all")
print(rpad("Spectral envelope: mean", 36))
senv_mean(e10, ch="all", dims=1)
@time senv_mean(e10, ch="all", dims=1)
print(rpad("Spectral envelope: median", 36))
senv_median(e10, ch="all", dims=1)
@time senv_median(e10, ch="all", dims=1)
print(rpad("H-amplitude envelope", 36))
henv(e10, ch="all")
@time henv(e10, ch="all")
print(rpad("H-amplitude envelope: mean", 36))
henv_mean(e10, ch="all", dims=1)
@time henv_mean(e10, ch="all", dims=1)
print(rpad("H-amplitude envelope: median", 36))
henv_median(e10, ch="all", dims=1)
@time henv_median(e10, ch="all", dims=1)
print(rpad("Envelope correlations", 36))
env1, _ = tenv(e10, ch="F3")
env2, _ = tenv(e10, ch="F4")
env_cor(env1, env2)
@time env_cor(env1, env2)
print(rpad("ERP peaks", 36))
e10_erp = average_epochs(e10)
erp_peaks(e10_erp)
@time erp_peaks(e10_erp)
print(rpad("Coherence: MT", 36))
coherence(e10, e10, ch1="F3", ch2="F4", ep1=1, ep2=1, method=:mt)
@time coherence(e10, e10, ch1="F3", ch2="F4", ep1=1, ep2=1, method=:mt)
print(rpad("Coherence: FFT", 36))
coherence(e10, e10, ch1="F3", ch2="F4", ep1=1, ep2=1, method=:fft)
@time coherence(e10, e10, ch1="F3", ch2="F4", ep1=1, ep2=1, method=:fft)
print(rpad("Instant frequency", 36))
frqinst(e10)
@time frqinst(e10)
print(rpad("GED", 36))
ged(e10, e10, ch1=["F3", "F4"], ch2=["F3", "F4"], ep1=1:10, ep2=11:20)
@time ged(e10, e10, ch1=["F3", "F4"], ch2=["F3", "F4"], ep1=1:10, ep2=11:20)
print(rpad("Generate ICA", 36))
ica_decompose(eeg, n=15, iter=10)
@time ica_decompose(eeg, n=15, iter=10)
ic, ic_mw = ica_decompose(eeg, n=15, iter=10)
print(rpad("Reconstruct ICA", 36))
ica_reconstruct(eeg, ic, ic_mw, ic_idx=1)
@time ica_reconstruct(eeg, ic, ic_mw, ic_idx=1)
print(rpad("Remove ICA", 36))
ica_remove(eeg, ic, ic_mw, ic_idx=1)
@time ica_remove(eeg, ic, ic_mw, ic_idx=1)
print(rpad("ISPC 1", 36))
ispc(e10)
@time ispc(e10)
print(rpad("ISPC 2", 36))
ispc(e10, e10, ch1="F3", ch2="F4", ep1=1, ep2=1)
@time ispc(e10, e10, ch1="F3", ch2="F4", ep1=1, ep2=1)
print(rpad("ITPC", 36))
itpc(e10, ch="F3", t=256)
@time itpc(e10, ch="F3", t=256)
print(rpad("ITPC spectrogram", 36))
itpc_spec(e10, ch="F3", frq_lim=(10, 20), frq_n=11)
@time itpc_spec(e10, ch="F3", frq_lim=(10, 20), frq_n=11)
print(rpad("Signal difference: absdiff", 36))
mdiff(e10, e10, ch1="F3", ch2="F4", ep1=1, ep2=1, method=:absdiff)
@time mdiff(e10, e10, ch1="F3", ch2="F4", ep1=1, ep2=1, method=:absdiff)
print(rpad("Signal difference: diff2int", 36))
mdiff(e10, e10, ch1="F3", ch2="F4", ep1=1, ep2=1, method=:diff2int)
@time mdiff(e10, e10, ch1="F3", ch2="F4", ep1=1, ep2=1, method=:diff2int)
print(rpad("Mutual information 1", 36))
mutual_information(e10, ch="all")
@time mutual_information(e10, ch="all")
print(rpad("Mutual information 2", 36))
mutual_information(e10, e10, ch1="F3", ch2="F4", ep1=1, ep2=1)
@time mutual_information(e10, e10, ch1="F3", ch2="F4", ep1=1, ep2=1)
print(rpad("MSCI95: normal", 36))
msci95(e10, ch="all")
@time msci95(e10, ch="all")
print(rpad("MSCI95: boot", 36))
msci95(e10, ch="all", method=:boot)
@time msci95(e10, ch="all")
print(rpad("MSCI95: obj1, obj2", 36))
msci95(e10, e10, ch1="F3", ch2="F4", ep1=1, ep2=1)
@time msci95(e10, e10, ch1="F3", ch2="F4", ep1=1, ep2=1)
print(rpad("Generate 4 PCs", 36))
pca_decompose(e10, ch="all", n=4)
@time pca_decompose(e10, ch="all", n=4)
print(rpad("Reconstruct using 4 PCs", 36))
pc, pcv, pcm, pcmodel = pca_decompose(e10, ch="all", n=4)
pca_reconstruct(e10, ch="all", pc, pcmodel);
@time pca_reconstruct(e10, pc, ch="all", pcmodel);
print(rpad("Phase difference: phase", 36))
phdiff(e10, ch="all", avg=:phase)
@time phdiff(e10, ch="all", avg=:phase)
print(rpad("Phase difference: signal", 36))
phdiff(e10, ch="all", avg=:signal)
@time phdiff(e10, ch="all", avg=:signal)
print(rpad("PLI 1", 36))
pli(e10, ch="all")
@time pli(e10, ch="all")
print(rpad("PLI 2", 36))
pli(e10, e10, ch1="F3", ch2="F4", ep1=1, ep2=1)
@time pli(e10, e10, ch1="F3", ch2="F4", ep1=1, ep2=1)
print(rpad("PSD: welch", 36))
psd(e10, ch="all")
@time psd(e10, ch="all")
print(rpad("PSD: fft", 36))
psd(e10, ch="all", method=:fft)
@time psd(e10, ch="all", method=:fft)
print(rpad("PSD: stft", 36))
psd(e10, ch="all", method=:stft)
@time psd(e10, ch="all", method=:stft)
print(rpad("PSD: mt", 36))
psd(e10, ch="all", method=:mt)
@time psd(e10, ch="all", method=:mt)
print(rpad("PSD: mw", 36))
psd(e10, ch="all", method=:mw)
@time psd(e10, ch="all", method=:mw)
print(rpad("PSD: gh", 36))
psd(e10, ch="all", method=:gh)
@time psd(e10, ch="all", method=:gh)
print(rpad("PSD: cwt", 36))
psd(e10, ch="all", method=:cwt)
@time psd(e10, ch="all", method=:cwt)
print(rpad("PSD: relative", 36))
psd_rel(e10, ch="all", frq_lim=(10, 20))
@time psd_rel(e10, ch="all", frq_lim=(10, 20))
print(rpad("PSD slope", 36))
psd_slope(e10, ch="all")
@time psd_slope(e10, ch="all")
print(rpad("AMP", 36))
amp(e10, ch="all")
@time amp(e10, ch="all")
print(rpad("SNR: mean", 36))
snr(e10, ch="all", type=:mean)
@time snr(e10, ch="all", type=:mean)
print(rpad("SNR: RMS", 36))
snr(e10, ch="all", type=:rms)
@time snr(e10, ch="all", type=:rms)
print(rpad("Spectrogram: stft", 36))
NeuroAnalyzer.spectrogram(e10, ch="all")
@time NeuroAnalyzer.spectrogram(e10, ch="all")
print(rpad("Spectrogram: mt", 36))
NeuroAnalyzer.spectrogram(e10, ch="all", method=:mt)
@time NeuroAnalyzer.spectrogram(e10, ch="all", method=:mt)
print(rpad("Spectrogram: mw", 36))
NeuroAnalyzer.spectrogram(e10, ch="all", method=:mw)
@time NeuroAnalyzer.spectrogram(e10, ch="all", method=:mw)
print(rpad("Spectrogram: gh", 36))
NeuroAnalyzer.spectrogram(e10, ch="all", method=:gh)
@time NeuroAnalyzer.spectrogram(e10, ch="all", method=:gh)
print(rpad("Spectrogram: cwt", 36))
NeuroAnalyzer.spectrogram(e10, ch="all", method=:cwt)
@time NeuroAnalyzer.spectrogram(e10, ch="all", method=:cwt)
print(rpad("Spectrum: FFT", 36))
spectrum(e10, ch="all")
@time spectrum(e10, ch="all")
print(rpad("Spectrum: Hilbert", 36))
spectrum(e10, ch="all", h=true)
@time spectrum(e10, ch="all", h=true)
print(rpad("Stationarity: mean", 36))
stationarity(e10, ch="all", method=:mean)
@time stationarity(e10, ch="all", method=:mean)
print(rpad("Stationarity: var", 36))
stationarity(e10, ch="all", method=:var)
@time stationarity(e10, ch="all", method=:var)
print(rpad("Stationarity: cov", 36))
stationarity(e10, ch="all", method=:cov)
@time stationarity(e10, ch="all", method=:cov)
print(rpad("Stationarity: hilbert", 36))
stationarity(e10, ch="all", method=:hilbert)
@time stationarity(e10, ch="all", method=:hilbert)
print(rpad("Stationarity: adf", 36))
stationarity(e10, ch="all", method=:adf)
@time stationarity(e10, ch="all", method=:adf)
print(rpad("Epoch stats", 36))
epoch_stats(e10, ch="all")
@time epoch_stats(e10, ch="all")
print(rpad("Channel stats", 36))
channel_stats(e10, ch="all")
@time channel_stats(e10, ch="all")
print(rpad("Coherence: MT", 36))
coherence(e10, e10, ch1="F3", ch2="F4", ep1=1, ep2=2, method=:mt)
@time coherence(e10, e10, ch1="F3", ch2="F4", ep1=1, ep2=2, method=:mt)
print(rpad("Coherence: FFT", 36))
coherence(e10, e10, ch1="F3", ch2="F4", ep1=1, ep2=2, method=:fft)
@time coherence(e10, e10, ch1="F3", ch2="F4", ep1=1, ep2=2, method=:fft)
print(rpad("TKEO", 36))
tkeo(e10, ch="all")
@time tkeo(e10, ch="all")
print(rpad("Total power: welch", 36))
total_power(e10, ch="all")
@time total_power(e10, ch="all")
print(rpad("Total power: fft", 36))
total_power(e10, ch="all", method=:fft)
@time total_power(e10, ch="all", method=:fft)
print(rpad("Total power: stft", 36))
total_power(e10, ch="all", ch="all", method=:stft)
@time total_power(e10, ch="all", ch="all", method=:stft)
print(rpad("Total power: mt", 36))
total_power(e10, ch="all", ch="all", ch="all", method=:mt)
@time total_power(e10, ch="all", ch="all", method=:mt)
print(rpad("Total power: mw", 36))
total_power(e10, ch="all", ch="all", method=:mw)
@time total_power(e10, ch="all", ch="all", method=:mw)
print(rpad("F-test 1", 36))
vartest(e10, ch="all")
@time vartest(e10, ch="all")
print(rpad("F-test 2", 36))
vartest(e10, e10, ch1="F3", ch2="F4", ep1=1, ep2=1)
@time vartest(e10, e10, ch1="F3", ch2="F4", ep1=1, ep2=1)
print(rpad("Cross-covariance", 36))
x = xcov(e10, e10, ch1="F3", ch2="F4", ep1=1, ep2=1)
@time xcov(e10, e10, ch1="F3", ch2="F4", ep1=1, ep2=1)
print(rpad("Cross-correlation", 36))
x = xcor(e10, e10, ch1="F3", ch2="F4", ep1=1, ep2=1)
@time xcor(e10, e10, ch1="F3", ch2="F4", ep1=1, ep2=1)
print(rpad("Cross-power spectral density", 36))
pxy, f = cpsd(e10, e10, ch1="F3", ch2="F4", ep1=1, ep2=1)
@time cpsd(e10, e10, ch1="F3", ch2="F4", ep1=1, ep2=1)

@info "Benchmarking completed"
