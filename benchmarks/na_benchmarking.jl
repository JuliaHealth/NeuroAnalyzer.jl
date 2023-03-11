using Wavelets
using ContinuousWavelets

#=
using BenchmarkTools

BenchmarkTools.DEFAULT_PARAMETERS.samples = 1
b = @benchmarkable import_edf("test/eeg-test-edf.edf") evals=20
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

@info "Benchmarking: io.jl"
println()
println("# IO")
println()

print(rpad("Import EDF", 36))
eeg = import_edf("test/eeg-test-edf.edf");
@time import_edf("test/eeg-test-edf.edf");
print(rpad("Import EDF+", 36))
import_edf("test/eeg-test-edfplus.edf");
@time import_edf("test/eeg-test-edfplus.edf");
print(rpad("Import BDF+", 36))
import_bdf("test/eeg-test-bdf.bdf");
@time import_bdf("test/eeg-test-bdf.bdf");
print(rpad("Import Digitrack", 36))
import_digitrack("test/eeg-test-digitrack.txt");
@time import_digitrack("test/eeg-test-digitrack.txt");
print(rpad("Import BrainVision", 36))
import_bv("test/eeg-test-bv.vhdr");
@time import_bv("test/eeg-test-bv.vhdr");
print(rpad("Import CSV ch×t", 36))
import_csv("test/eeg-test_chxt.csv.gz");
@time import_csv("test/eeg-test_chxt.csv.gz");
print(rpad("Import CSV t×ch", 36))
import_csv("test/eeg-test_txch.csv.gz");
@time import_csv("test/eeg-test_txch.csv.gz");
print(rpad("Import EEGLAB dataset", 36))
import_set("test/eeg-test.set");
@time import_set("test/eeg-test.set");
print(rpad("Save HDF5", 36))
tmp = tempname()
save(eeg, file_name=tmp);
tmp = tempname()
@time save(eeg, file_name=tmp);
print(rpad("Load HDF5", 36))
load(tmp);
@time load(tmp);

eeg = import_edf("test/eeg-test-edf.edf");
delete_channel!(eeg, channel=20:24);
e10 = epoch(eeg, ep_len=2560);
load_locs!(e10, file_name="locs/standard-10-20-cap19-elmiko.ced")

@info "Benchmarking: edit.jl"
println()
println("# EDIT")
println()

e10_tmp = deepcopy(e10);
print(rpad("Delete channel", 36))
delete_channel(e10_tmp, channel=1);
@time delete_channel(e10_tmp, channel=1);
e10_tmp = deepcopy(e10);
print(rpad("Keep channel", 36))
keep_channel(e10_tmp, channel=1);
@time keep_channel(e10_tmp, channel=1);
e10_tmp = deepcopy(e10);
print(rpad("Delete epoch", 36))
delete_epoch(e10_tmp, epoch=1);
@time delete_epoch(e10_tmp, epoch=1);
e10_tmp = deepcopy(e10);
print(rpad("Keep epoch", 36))
keep_epoch(e10_tmp, epoch=1);
@time keep_epoch(e10_tmp, epoch=1);
print(rpad("Virtual channel", 36))
vch(e10, f="fp1 + fp2");
@time vch(e10, f="fp1 + fp2");

@info "Benchmarking: process.jl"
println()
println("# PROCESS")
println()

print(rpad("A referencing", 36))
edf_am = import_edf("test/eeg-test-edf.edf")
delete_channel!(edf_am, channel=[22, 23, 24])
e10_am = epoch(edf_am, ep_len=2560)
reference_a(e10_am)
@time reference_a(e10_am);
print(rpad("M referencing", 36))
edit_channel!(e10_am, channel=17, field=:labels, value="M1")
edit_channel!(e10_am, channel=18, field=:labels, value="M2")
reference_m(e10_am)
@time reference_m(e10_am);
print(rpad("CAR referencing", 36))
reference_car(e10);
@time reference_car(e10);
print(rpad("Channel referencing", 36))
reference_ch(e10, channel=1);
@time reference_ch(e10, channel=1);
print(rpad("Filter: FIR LP", 36))
NeuroAnalyzer.filter(e10, fprototype=:fir, ftype=:lp, cutoff=45.0, order=8);
@time NeuroAnalyzer.filter(e10, fprototype=:fir, ftype=:lp, cutoff=45.0, order=8);
print(rpad("Filter: FIR HP", 36))
NeuroAnalyzer.filter(e10, fprototype=:fir, ftype=:hp, cutoff=0.1, order=8);
@time NeuroAnalyzer.filter(e10, fprototype=:fir, ftype=:hp, cutoff=0.1, order=8);
print(rpad("Filter: IIR LP", 36))
NeuroAnalyzer.filter(e10, fprototype=:butterworth, ftype=:lp, cutoff=45.0, order=8);
@time NeuroAnalyzer.filter(e10, fprototype=:butterworth, ftype=:lp, cutoff=45.0, order=8);
print(rpad("Filter: IIR HP", 36))
NeuroAnalyzer.filter(e10, fprototype=:butterworth, ftype=:hp, cutoff=0.1, order=8);
@time NeuroAnalyzer.filter(e10, fprototype=:butterworth, ftype=:hp, cutoff=0.1, order=8);
print(rpad("Filter: IIR notch", 36))
NeuroAnalyzer.filter(e10, fprototype=:iirnotch, cutoff=50, bw=2);
@time NeuroAnalyzer.filter(e10, fprototype=:iirnotch, cutoff=50, bw=2);
print(rpad("Filter: Remez", 36))
NeuroAnalyzer.filter(e10, fprototype=:remez, ftype=:lp, order=128, cutoff=20, bw=0.5);
@time NeuroAnalyzer.filter(e10, fprototype=:remez, ftype=:lp, order=128, cutoff=20, bw=0.5);
print(rpad("Filter: polynomial (8)", 36))
NeuroAnalyzer.filter(e10, fprototype=:poly, order=8, window=256);
@time NeuroAnalyzer.filter(e10, fprototype=:poly, order=8, window=256);
print(rpad("Filter: moving average", 36))
NeuroAnalyzer.filter(e10, fprototype=:mavg);
@time NeuroAnalyzer.filter(e10, fprototype=:mavg);
print(rpad("Interpolate: PL", 36))
plinterpolate_channel(e10, channel=1, epoch=1);
@time plinterpolate_channel(e10, channel=1, epoch=1);
print(rpad("Interpolate: LR", 36))
lrinterpolate_channel(e10, channel=1, epoch=1);
@time lrinterpolate_channel(e10, channel=1, epoch=1);
print(rpad("Surface Laplacian", 36))
slaplacian(e10);
@time slaplacian(e10);

@info "Benchmarking: analyze.jl"
println()
println("# ANALYZE")
println()

print(rpad("Total power", 36))
total_power(e10);
@time total_power(e10);
print(rpad("Total power: mt", 36))
total_power(e10, mt=true);
@time total_power(e10, mt=true);
print(rpad("Band power", 36))
band_power(e10, f=(10, 20));
@time band_power(e10, f=(10, 20));
print(rpad("Band power: mt", 36))
band_power(e10, f=(10, 20), mt=true);
@time band_power(e10, f=(10, 20), mt=true);
print(rpad("Covariance matrix", 36))
covm(e10);
@time covm(e10);
print(rpad("Correlation matrix", 36))
corm(e10);
@time corm(e10);
print(rpad("Auto-covariance", 36))
acov(e10);
@time acov(e10);
print(rpad("Cross-covariance 1", 36))
xcov(e10, lag=10, demean=true);
@time xcov(e10, lag=10, demean=true);
print(rpad("Cross-covariance 2", 36))
xcov(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=2, lag=10, demean=true);
@time xcov(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=2, lag=10, demean=true);
print(rpad("PSD", 36))
psd(e10);
@time psd(e10);
print(rpad("PSD: mt", 36))
psd(e10, mt=true);
@time psd(e10, mt=true);
print(rpad("PSD: mw ", 36))
psd_mw(e10);
@time psd_mw(e10);
print(rpad("Stationarity: mean", 36))
stationarity(e10, method=:mean);
@time stationarity(e10, method=:mean);
print(rpad("Stationarity: var", 36))
stationarity(e10, method=:var);
@time stationarity(e10, method=:var);
print(rpad("Stationarity: cov", 36))
stationarity(e10, method=:cov);
@time stationarity(e10, method=:cov);
print(rpad("Stationarity: hilbert", 36))
stationarity(e10, method=:hilbert);
@time stationarity(e10, method=:hilbert);
print(rpad("Stationarity: adf", 36))
stationarity(e10, method=:adf);
@time stationarity(e10, method=:adf);
print(rpad("Mutual information 1", 36))
mutual_information(e10);
@time mutual_information(e10);
print(rpad("Mutual information 2", 36))
mutual_information(e10, e10);
@time mutual_information(e10, e10);
print(rpad("Entropy", 36))
entropy(e10);
@time entropy(e10);
print(rpad("Negentropy", 36))
negentropy(e10);
@time negentropy(e10);
print(rpad("Time coherence", 36))
tcoherence(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=2);
@time tcoherence(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=2);
print(rpad("Signal difference: absdiff", 36))
difference(e10, e10; method=:absdiff);
@time difference(e10, e10; method=:absdiff);
print(rpad("Signal difference: diff2int", 36))
difference(e10, e10; method=:diff2int);
@time difference(e10, e10; method=:diff2int);
print(rpad("Epoch stats", 36))
epoch_stats(e10);
@time epoch_stats(e10);
print(rpad("Spectrogram: standard", 36))
NeuroAnalyzer.spectrogram(e10);
@time NeuroAnalyzer.spectrogram(e10);
print(rpad("Spectrogram: mt", 36))
NeuroAnalyzer.spectrogram(e10, method=:mt);
@time NeuroAnalyzer.spectrogram(e10, method=:mt);
print(rpad("Spectrogram: stft", 36))
NeuroAnalyzer.spectrogram(e10, method=:stft);
@time NeuroAnalyzer.spectrogram(e10, method=:stft);
print(rpad("Spectrogram: mw", 36))
NeuroAnalyzer.spectrogram(e10, method=:mw);
@time NeuroAnalyzer.spectrogram(e10, method=:mw);
print(rpad("Spectrogram: cwt", 36))
NeuroAnalyzer.spectrogram(e10, method=:cwt);
@time NeuroAnalyzer.spectrogram(e10, method=:cwt);
print(rpad("Spectrum: FFT", 36))
spectrum(e10);
@time spectrum(e10);
print(rpad("Spectrum: Hilbert", 36))
spectrum(e10, h=true);
@time spectrum(e10, h=true);
print(rpad("Channel stats", 36))
channel_stats(e10);
@time channel_stats(e10);
print(rpad("SNR", 36))
snr(e10);
@time snr(e10);
print(rpad("Standardize", 36))
standardize(e10);
@time standardize(e10);
print(rpad("Frequency convolution", 36))
fconv(e10, kernel=generate_morlet(256, 1, 32, complex=true));
@time fconv(e10, kernel=generate_morlet(256, 1, 32, complex=true));
print(rpad("Time convolution", 36))
tconv(e10, kernel=generate_morlet(256, 1, 32, complex=true));
@time tconv(e10, kernel=generate_morlet(256, 1, 32, complex=true));
print(rpad("DFT", 36))
dft(e10);
@time dft(e10);
print(rpad("MSCI95: normal", 36))
msci95(e10);
@time msci95(e10);
print(rpad("MSCI95: boot", 36))
msci95(e10, method=:boot);
@time msci95(e10);
print(rpad("MSCI85: 2 objects", 36))
msci95(e10, e10);
@time mean(e10, e10);
print(rpad("Subtract channels", 36))
chdiff(e10, e10); 
@time chdiff(e10, e10); 
print(rpad("Temporal envelope", 36))
tenv(e10);
@time tenv(e10);
print(rpad("Temporal envelope: mean", 36))
tenv_mean(e10, dims=1);
@time tenv_mean(e10, dims=1);
print(rpad("Temporal envelope: median", 36))
tenv_median(e10, dims=1);
@time tenv_median(e10, dims=1);
print(rpad("Power envelope", 36))
penv(e10);
@time penv(e10);
print(rpad("Power envelope: mean", 36))
penv_mean(e10, dims=1);
@time penv_mean(e10, dims=1);
print(rpad("Power envelope: median", 36))
penv_median(e10, dims=1);
@time penv_median(e10, dims=1);
print(rpad("Spectral envelope", 36))
senv(e10);
@time senv(e10);
print(rpad("Spectral envelope: mean", 36))
senv_mean(e10, dims=1);
@time senv_mean(e10, dims=1);
print(rpad("Spectral envelope: median", 36))
senv_median(e10, dims=1);
@time senv_median(e10, dims=1);
print(rpad("Hilbert amplitude envelope", 36))
henv(e10);
@time henv(e10);
print(rpad("Hilbert amplitude envelope: mean", 36))
henv_mean(e10, dims=1);
@time henv_mean(e10, dims=1);
print(rpad("Hilbert amplitude envelope: median", 36))
henv_median(e10, dims=1);
@time henv_median(e10, dims=1);
print(rpad("ISPC 1", 36))
ispc(e10);
@time ispc(e10);
print(rpad("ISPC 2", 36))
ispc(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=1);
@time ispc(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=1);
print(rpad("ITPC", 36))
itpc(e10, channel=1, t=256);
@time itpc(e10, channel=1, t=256);
print(rpad("ITPC spectrogram", 36))
itpc_s(e10, channel=1, frq_lim=(10, 20), frq_n=11);
@time itpc_s(e10, channel=1, frq_lim=(10, 20), frq_n=11);
print(rpad("PLI 1", 36))
pli(e10);
@time pli(e10);
print(rpad("PLI 2", 36))
pli(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=1);
@time pli(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=1);
print(rpad("Envelope correlation: amp", 36))
env_cor(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=1, type=:amp);
@time env_cor(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=1, type=:amp);
print(rpad("Envelope correlation: pow", 36))
env_cor(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=1, type=:pow);
@time env_cor(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=1, type=:pow);
print(rpad("Envelope correlation: spec", 36))
env_cor(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=1, type=:spec);
@time env_cor(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=1, type=:spec);
print(rpad("Envelope correlation: hamp", 36))
env_cor(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=1, type=:hamp);
@time env_cor(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=1, type=:hamp);
print(rpad("GED", 36))
ged(e10, e10);
@time ged(e10, e10);
print(rpad("Instant frequency", 36))
frqinst(e10);
@time frqinst(e10);
print(rpad("TKEO", 36))
tkeo(e10);
@time tkeo(e10);
print(rpad("Frequency coherence", 36))
fcoherence(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=1);
@time fcoherence(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=1);
print(rpad("F-test 1", 36))
vartest(e10);
@time vartest(e10);
print(rpad("F-test 2", 36))
vartest(e10, e10);
@time vartest(e10, e10);
print(rpad("Band power", 36))
band_mpower(e10; f=(10, 20));
@time band_mpower(e10; f=(10, 20));
print(rpad("Relative PSD", 36))
psd_rel(e10; f=(10, 20));
@time psd_rel(e10; f=(10, 20));
print(rpad("Frequency band split", 36))
fbsplit(e10);
@time fbsplit(e10);
print(rpad("Channel difference", 36))
chdiff(e10, e10, channel1=1, channel2=2);
@time chdiff(e10, e10, channel1=1, channel2=2);
print(rpad("Cross power spectrum 1", 36))
cps(e10);
@time cps(e10);
print(rpad("Cross power spectrum 2", 36))
cps(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=1);
@time cps(e10, e10, channel1=1, channel2=2, epoch1=1, epoch2=1);
print(rpad("Amplitude difference", 36))
ampdiff(e10);
@time ampdiff(e10);
print(rpad("DWT", 36))
dw_trans(e10, wt=wavelet(WT.haar), type=:sdwt);
@time dw_trans(e10, wt=wavelet(WT.haar), type=:sdwt);
print(rpad("CWT", 36))
cw_trans(e10, wt=wavelet(Morlet(π), β=2));
@time cw_trans(e10, wt=wavelet(Morlet(π), β=2));
print(rpad("PSD slope", 36))
psdslope(e10);
@time psdslope(e10);
print(rpad("apply()", 36))
apply(e10, f="mean(eeg, dims=1)", channel=1:4);
@time apply(e10, f="mean(eeg, dims=1)", channel=1:4);
print(rpad("Normalize", 36))
normalize(e10, method=:zscore);
@time normalize(e10, method=:zscore);
print(rpad("Remove DC", 36))
remove_dc(e10);
@time remove_dc(e10);
print(rpad("Taper", 36))
taper(e10, t=generate_window(:hann, epoch_len(e10)));
@time taper(e10, t=generate_window(:hann, epoch_len(e10)));
print(rpad("Derivative", 36))
derivative(e10);
@time derivative(e10);
print(rpad("Detrend", 36))
detrend(e10, type=:constant);
@time detrend(e10, type=:constant);
print(rpad("Wiener denoising", 36))
denoise_wien(e10);
@time denoise_wien(e10);
print(rpad("FFT denoising", 36))
denoise_fft(e10);
@time denoise_fft(e10);
print(rpad("Generate PCA", 36))
pca(e10, n=4);
@time pca(e10, n=4);
print(rpad("Generate ICA", 36))
ica(e10, n=15, tol=1.0);
@time ica(e10, n=15, tol=1.0);
i, i_mw = ica(e10, n=15, tol=1.0)
e10_ica = add_component(e10, c=:ica, v=i);
add_component!(e10_ica, c=:ica_mw, v=i_mw);
print(rpad("Remove ICA", 36))
ica_reconstruct(e10_ica, ic=1);
@time ica_reconstruct(e10_ica, ic=1);
print(rpad("Split using DWT", 36))
bands_dwt(e10, channel=1, wt=wavelet(WT.db2), type=:sdwt, n=5);
@time bands_dwt(e10, channel=1, wt=wavelet(WT.db2), type=:sdwt, n=5);