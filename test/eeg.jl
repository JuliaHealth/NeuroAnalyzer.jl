using NeuroJ
using Test

edf = eeg_import_edf("eeg-test-edf.edf")

edf1 = eeg_delete_channel(edf, 1)
@test edf1.eeg_header[:channels_no] == 18

edf1 = eeg_keep_channel(edf, 1)
@test edf1.eeg_header[:channels_no] == 1

edf1 = eeg_derivative(edf)
@test round(sum(edf1.eeg_signals), digits=2) == -1911.99

tbp = eeg_total_power(edf)
@test round(sum(tbp), digits=2) == 375052.7

abp = eeg_band_power(edf, f1=2, f2=4)
@test round(sum(abp), digits=2) == 1542.35

edf1 = eeg_detrend(edf)
@test round(sum(edf1.eeg_signals), digits=2) == -0.0

edf1 = eeg_reference_channel(edf, 1)
@test round(sum(edf1.eeg_signals), digits=2) == -6.2119013128e8

edf1 = eeg_reference_car(edf)
@test round(sum(edf1.eeg_signals), digits=2) == 0.0

edf1 = eeg_extract_channel(edf, "Cz")
@test round(sum(edf1), digits=2) == -4.839052834e7

edf1 = eeg_extract_channel(edf, 18)
@test round(sum(edf1), digits=2) == -4.839052834e7

@test eeg_get_channel(edf, 1) == "Fp1"
@test eeg_get_channel(edf, "Fp1") == 1

edf1 = eeg_rename_channel(edf, "Cz", "CZ")
@test edf1.eeg_header[:labels][18] == "CZ"
edf1 = eeg_rename_channel(edf, 1, "FP1")
@test edf1.eeg_header[:labels][1] == "FP1"

edf1 = eeg_taper(edf, edf.eeg_signals[1, :, 1])
@test round(sum(edf1.eeg_signals), digits=2) == 6.687805610437e10

edf1 = eeg_demean(edf)
@test round(sum(edf1.eeg_signals), digits=2) == 0.0

edf1 = eeg_normalize_zscore(edf)
@test round(sum(edf1.eeg_signals), digits=2) == 0.0

edf1 = eeg_normalize_minmax(edf)
@test round(sum(edf1.eeg_signals), digits=2) == 3.89948242e6

cov_m = eeg_cov(edf)
@test round(sum(cov_m), digits=2) == 156620.75

cor_m = eeg_cor(edf)
@test round(sum(cor_m), digits=2) == 141.04

edf1 = eeg_upsample(edf, new_sr=512)
@test round(sum(edf1.eeg_signals), digits=2) == -9.5454636605e8

@test eeg_history(edf) == String[]

@test eeg_labels(edf)[1] == "Fp1"

@test eeg_samplingrate(edf) == 256

edf1 = eeg_epochs(edf, epochs_len=10, average=true)
@test round(sum(edf1.eeg_signals), digits=2) == -13451.5

edf1 = eeg_extract_epoch(edf, 1)
@test round(sum(edf1.eeg_signals), digits=2) == -4.7727273831e8

edf1 = eeg_tconv(edf, kernel=generate_hanning(256))
@test round(sum(edf1.eeg_signals), digits=2) == -6.13360628002e10

edf1 = eeg_filter(edf, fprototype=:butterworth, ftype=:lp, cutoff=2, order=8)
@test round(sum(edf1.eeg_signals), digits=2) == -4.7726293721e8
edf1 = eeg_filter(edf, fprototype=:mavg, d=10)
@test round(sum(edf1.eeg_signals), digits=2) == -4.7727283843e8
edf1 = eeg_filter(edf, fprototype=:mmed, d=10)
@test size(edf1.eeg_signals) == (19, 354816, 1)
edf1 = eeg_filter(edf, fprototype=:mavg, window=generate_gaussian(eeg_samplingrate(edf), 32, 0.01))
@test round(sum(edf1.eeg_signals), digits=2) == -2.165654661766e10

edf1 = eeg_downsample(edf, new_sr=128)
@test round(sum(edf1.eeg_signals), digits=2) == -2.3863254176e8

acov_m, _ = eeg_autocov(edf)
@test round(sum(acov_m), digits=2) == 7.6170287087299e11

ccov_m, _ = eeg_crosscov(edf)
@test round(sum(ccov_m), digits=2) == 4.34119627649116e12

p, f = eeg_psd(edf1)
@test size(p, 1) == 19

p = eeg_stationarity(edf, method=:mean)
@test round(sum(p), digits=2) == -13451.5
p = eeg_stationarity(edf, method=:var)
@test round(sum(p), digits=2) == 5.54177262e6
p = eeg_stationarity(edf, method=:hilbert)
@test round(sum(p), digits=2) == 3.53311741e6
p = eeg_stationarity(edf, window=10000, method=:euclid)
@test round(sum(p), digits=2) == 1.0040717529475e11

e = eeg_trim(edf, trim_len=(10 * eeg_samplingrate(edf)), offset=(20 * eeg_samplingrate(edf)), from=:start)
@test size(e.eeg_signals) == (19, 352256, 1)

m = eeg_mi(edf)
@test round(sum(m), digits=2) == 363.64
m = eeg_mi(edf, edf)
@test round(sum(m), digits=2) == 363.64

e = eeg_entropy(edf)
@test round(sum(e), digits=2) == 116.06

a = eeg_band(:alpha)
@test a == (8, 13)

true