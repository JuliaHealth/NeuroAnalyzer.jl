using NeuroJ
using Test

edf = eeg_import_edf("eeg-test-edf.edf")

edf1 = eeg_delete_channel(edf, channel=1)
@test edf1.eeg_header[:channel_n] == 18

edf1 = eeg_keep_channel(edf, channel=1)
@test edf1.eeg_header[:channel_n] == 1

edf1 = eeg_derivative(edf)
@test size(edf1.eeg_signals) == (19, 354816, 1)

tbp = eeg_total_power(edf)
@test size(tbp) == (19, 1)

abp = eeg_band_power(edf, f1=2, f2=4)
@test size(abp) == (19, 1)

edf1 = eeg_detrend(edf)
@test size(edf1.eeg_signals) == (19, 354816, 1)

edf1 = eeg_reference_channel(edf, channel=1)
@test size(edf1.eeg_signals) == (19, 354816, 1)

edf1 = eeg_reference_car(edf)
@test size(edf1.eeg_signals) == (19, 354816, 1)

edf1 = eeg_extract_channel(edf, channel="Cz")
@test size(edf1) == (354816, )

edf1 = eeg_extract_channel(edf, channel=18)
@test size(edf1) == (354816, )

@test eeg_get_channel(edf, channel=1) == "Fp1"
@test eeg_get_channel(edf, channel="Fp1") == 1

edf1 = eeg_rename_channel(edf, channel="Cz", new_name="CZ")
@test edf1.eeg_header[:labels][18] == "CZ"
edf1 = eeg_rename_channel(edf, channel=1, new_name="FP1")
@test edf1.eeg_header[:labels][1] == "FP1"

edf1 = eeg_taper(edf, taper=edf.eeg_signals[1, :, 1])
@test size(edf1.eeg_signals) == (19, 354816, 1)

edf1 = eeg_demean(edf)
@test size(edf1.eeg_signals) == (19, 354816, 1)

edf1 = eeg_normalize_zscore(edf)
@test size(edf1.eeg_signals) == (19, 354816, 1)

edf1 = eeg_normalize_minmax(edf)
@test size(edf1.eeg_signals) == (19, 354816, 1)

cov_m = eeg_cov(edf)
@test size(cov_m) == (19, 19, 1)

cor_m = eeg_cor(edf)
@test size(cor_m) == (19, 19, 1)

edf1 = eeg_upsample(edf, new_sr=512)
@test size(edf1.eeg_signals) == (19, 709631, 1)

@test eeg_history(edf) == String[]

@test eeg_labels(edf)[1] == "Fp1"

@test eeg_sr(edf) == 256

edf1 = eeg_epochs(edf, epoch_len=10, average=true)
@test size(edf1.eeg_signals) == (19, 10, 1)

edf1 = eeg_extract_epoch(edf, epoch=1)
@test size(edf1.eeg_signals) == (19, 354816, 1)

edf1 = eeg_tconv(edf, kernel=generate_hanning(256))
@test size(edf1.eeg_signals) == (19, 354816, 1)

edf1 = eeg_filter(edf, fprototype=:butterworth, ftype=:lp, cutoff=2, order=8)
@test size(edf1.eeg_signals) == (19, 354816, 1)
edf1 = eeg_filter(edf, fprototype=:mavg, d=10)
@test size(edf1.eeg_signals) == (19, 354816, 1)
edf1 = eeg_filter(edf, fprototype=:mmed, d=10)
@test size(edf1.eeg_signals) == (19, 354816, 1)
edf1 = eeg_filter(edf, fprototype=:mavg, window=generate_gaussian(eeg_sr(edf), 32, 0.01))
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

e = eeg_trim(edf, trim_len=(10 * eeg_sr(edf)), offset=(20 * eeg_sr(edf)), from=:start)
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

e = eeg_fconv(edf, kernel=[1, 2, 3, 4])
@test size(edf.eeg_signals) == (19, 354816, 1)

p, v = eeg_pca(edf, n=2)
@test size(p) == (2, 354816, 1)
@test size(v) == (2, 1)

e = eeg_edit(edf, field=:patient, value="unknown")
@test e.eeg_header[:patient] == "unknown"

e = eeg_epochs(edf, epoch_n=10)
e9 = eeg_delete_epoch(e, epoch=10)
@test size(e9.eeg_signals) == (19, 35481, 9)
e1 = eeg_keep_epoch(e, epoch=1)
@test size(e1.eeg_signals) == (19, 35481, 1)

e = eeg_pick(edf, pick=:left)
@test length(e) == 8

e = eeg_epochs(edf, epoch_len=20*256)
m, s, v = eeg_epochs_stats(e)
@test size(v) == (69, )

e = eeg_epochs(edf, epoch_len=20, average=true)
i = eeg_ica(e, n=5)
@test size(i) == (5, 20, 1)

true