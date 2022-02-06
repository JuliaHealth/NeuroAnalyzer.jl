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
@test round(sum(edf1.eeg_signals), digits=2) == -4.7741781538e8

edf1 = eeg_filter(edf, fprototype=:butterworth, ftype=:lp, cutoff=2, order=8)
@test round(sum(edf1.eeg_signals), digits=2) == -4.7726293721e8

edf1 = eeg_downsample(edf, new_sr=128)
@test round(sum(edf1.eeg_signals), digits=2) == -2.3863254176e8

acov_m, _ = eeg_autocov(edf)
@test round(sum(acov_m), digits=2) == 7.6170287087299e11

ccov_m, _ = eeg_crosscov(edf)
@test round(sum(ccov_m), digits=2) == 4.34119627649116e12

p, f = eeg_psd(edf1)
@test size(p, 1) == 19

true