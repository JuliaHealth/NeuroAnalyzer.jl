using NeuroJ
using Test

edf = eeg_import_edf("eeg-test-edf.edf")

edf1 = eeg_delete_channel(edf, 1)
@test edf1.eeg_header[:channels_no] == 18

edf1 = eeg_keep_channel(edf, 1)
@test edf1.eeg_header[:channels_no] == 1

edf1 = eeg_derivative(edf)
@test size(edf1.eeg_signals) == (19, 354816, 1)

tbp = eeg_total_power(edf)
@test size(tbp) == (19, 1)

abp = eeg_band_power(edf, f1=2, f2=4)
@test size(abp) == (19, 1)

edf1 = eeg_detrend(edf)
@test size(edf1.eeg_signals) == (19, 354816, 1)

edf1 = eeg_reference_channel(edf, 1)
@test size(edf1.eeg_signals) == (19, 354816, 1)

edf1 = eeg_reference_car(edf)
@test size(edf1.eeg_signals) == (19, 354816, 1)

edf1 = eeg_extract_channel(edf, "Cz")
@test size(edf1) == (354816, )

edf1 = eeg_extract_channel(edf, 18)
@test size(edf1) == (354816, )

@test eeg_get_channel(edf, 1) == "Fp1"
@test eeg_get_channel(edf, "Fp1") == 1

edf1 = eeg_rename_channel(edf, "Cz", "CZ")
@test edf1.eeg_header[:labels][18] == "CZ"
edf1 = eeg_rename_channel(edf, 1, "FP1")
@test edf1.eeg_header[:labels][1] == "FP1"

edf1 = eeg_taper(edf, edf.eeg_signals[1, :, 1])
@test size(edf1.eeg_signals) == (19, 354816, 1)

edf1 = eeg_demean(edf)
@test size(edf1.eeg_signals) == (19, 354816, 1)

edf1 = eeg_normalize_zscore(edf)
@test size(edf1.eeg_signals) == (19, 354816, 1)

edf1 = eeg_normalize_minmax(edf)
@test size(edf1.eeg_signals) == (19, 354816, 1)

cov_m = eeg_cov(edf)
@test size(cov_m) == (19, 19)

cor_m = eeg_cor(edf)
@test size(cor_m) == (19, 19)

edf1 = eeg_upsample(edf, new_sr=512)
@test size(edf1.eeg_signals) == (19, 709631, 1)

@test eeg_history(edf) == String[]

@test eeg_labels(edf)[1] == "Fp1"

@test eeg_samplingrate(edf) == 256

edf1 = eeg_epochs(edf, epochs_len=10, average=true)
@test size(edf1.eeg_signals) == (19, 10, 1)

edf1 = eeg_extract_epoch(edf, 1)
@test size(edf1.eeg_signals) == (19, 354816, 1)

edf1 = eeg_tconv(edf, kernel=generate_hanning(256))
@test size(edf1.eeg_signals) == (19, 354815, 1)

edf1 = eeg_filter(edf, fprototype=:butterworth, ftype=:lp, cutoff=2, order=8)
@test size(edf1.eeg_signals) == (19, 354816, 1)
edf1 = eeg_filter(edf, fprototype=:mavg, d=10)
@test size(edf1.eeg_signals) == (19, 354816, 1)
edf1 = eeg_filter(edf, fprototype=:mmed, d=10)
@test size(edf1.eeg_signals) == (19, 354816, 1)
edf1 = eeg_filter(edf, fprototype=:mavg, window=generate_gaussian(eeg_samplingrate(edf), 32, 0.01))
@test size(edf1.eeg_signals) == (19, 354816, 1)

edf1 = eeg_downsample(edf, new_sr=128)
@test size(edf1.eeg_signals) == (19, 177408, 1)

acov_m, _ = eeg_autocov(edf)
@test size(acov_m) == (19, 3)

ccov_m, _ = eeg_crosscov(edf)
@test size(ccov_m) == (361, 3)

p, f = eeg_psd(edf1)
@test size(p, 1) == 19

p = eeg_stationarity(edf, method=:mean)
@test size(p) == (19, 10, 1)
p = eeg_stationarity(edf, method=:var)
@test size(p) == (19, 10, 1)
p = eeg_stationarity(edf, method=:hilbert)
@test size(p) == (19, 354815, 1)
p = eeg_stationarity(edf, window=10000, method=:euclid)
@test size(p) == (37, 1)

e = eeg_trim(edf, trim_len=(10 * eeg_samplingrate(edf)), offset=(20 * eeg_samplingrate(edf)), from=:start)
@test size(e.eeg_signals) == (19, 352256, 1)

m = eeg_mi(edf)
@test size(m) == (19, 19)
m = eeg_mi(edf, edf)
@test size(m) == (19, 19)

e = eeg_entropy(edf)
@test size(e) == (19, 1)

a = eeg_band(:alpha)
@test a == (8, 13)

m = eeg_coherence(edf, edf)
@test size(m) == (19, 354816)

hz, _ = eeg_freqs(edf)
@test typeof(hz) == Vector{Float64}

edf1 = eeg_epochs(edf, epochs_len=10*eeg_samplingrate(edf), average=true)
p, v = eeg_pca(edf1, edf1, n=2)
@test size(p) == (2, 2560, 19, 1)
@test size(v) == (2, 19, 1)

true