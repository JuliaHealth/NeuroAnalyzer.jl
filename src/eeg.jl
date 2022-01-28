"""
    signal_derivative(signal)

Returns the derivative of the `signal` vector with length same as the signal.

# Arguments

- `signals::Vector{Float64}` - the signal vector to analyze
"""
signal_derivative(signal::Vector{Float64}) = vcat(diff(signal), diff(signal)[end])

"""
    signal_derivative(signal)

Returns the derivative of each the `signal` matrix channels with length same as the signal.

# Arguments

- `signal::Matrix{Float64}` - the signal matrix to analyze (rows: channels, columns: time).
"""
function signal_derivative(signal::Matrix{Float64})
    channels_no = size(signal, 1)
    signal_der = zeros(size(signal))

    for idx in 1:channels_no
        signal_der[idx, :] = signal_derivative(signal[idx, :])
    end

    return signal_der
end

"""
    signal_total_power(signal)

Calculates total power for the `signal` vector.

# Arguments

- `signals::Vector{Float64}` - the signal vector to analyze
"""
function signal_total_power(signal::Vector{Float64}, fs)
    psd = welch_pgram(signal, 4*fs, fs=fs)
    # dx: frequency resolution
    dx = psd.freq[2] - psd.freq[1]
    stp = simpson(psd.power, dx=dx)

    return stp
end

"""
    signal_total_power(signal)

Calculates total power for each the `signal` matrix channels.

# Arguments

- `signal::Matrix{Float64}` - the signal matrix to analyze (rows: channels, columns: time).
"""
function signal_total_power(signal::Matrix{Float64}, fs)
    channels_no = size(signal, 1)
    stp = zeros(size(signal, 1))

    for idx in 1:channels_no
        stp[idx] = signal_total_power(signal[idx, :])
    end

    return stp
end

"""
    signal_band_power(signal, f1, f2)

Calculates absolute band power between frequencies `f1` and `f2` for the `signal` vector.

# Arguments

- `signals::Vector{Float64}` - the signal vector to analyze
- `fs::Float64` - Sampling rate of the signal
- `f1::Float64` - Lower frequency bound
- `f2::Float64` - Upper frequency bound
"""
function signal_band_power(signal::Vector{Float64}, fs::Float64, f1::Float64, f2::Float64)
    psd = welch_pgram(signal, 4*fs, fs=fs)
    frq_idx = [vsearch(Vector(psd.freq), f1), vsearch(Vector(psd.freq), f2)]
    # dx: frequency resolution
    dx = psd.freq[2] - psd.freq[1]
    sbp = simpson(psd.power[frq_idx[1]:frq_idx[2]], psd.freq[frq_idx[1]:frq_idx[2]], dx=dx)

    return sbp
end

"""
    signal_band_power(signal, f1, f2)

Calculates absolute band power between frequencies `f1` and `f2` for each the `signal` matrix channels.

# Arguments

- `signals::Matrix{Float64}` - the signal matrix to analyze
- `fs::Float64` - Sampling rate of the signal
- `f1::Float64` - Lower frequency bound
- `f2::Float64` - Upper frequency bound
"""
function signal_band_power(signal::Matrix{Float64}, fs, f1, f2)
    channels_no = size(signal, 1)
    sbp = zeros(size(signal, 1))

    for idx in 1:channels_no
        sbp[idx] = signal_band_power(signal[idx, :], fs, f1, f2)
    end

    return sbp
end

"""
    make_spectrum(signal, fs)

Returns FFT and DFT sample frequencies for a DFT for the `signal` vector.

# Arguments

- `signals::Vector{Float64}` - the signal vector to analyze
- `fs::Float64` - Sampling rate of the signal
"""
function signal_make_spectrum(signal::Vector{Float64}, fs)
    signal_fft = fft(signal)
    # number of samples
    n = length(signal)
    # time between samples
    d = 1 / fs
    signal_sf = fftfreq(n, d)

    return signal_fft, signal_sf
end

"""
    make_spectrum(signal, fs)

Returns FFT and DFT sample frequencies for a DFT for each the `signal` matrix channels.

# Arguments

- `signals::Matrix{Float64}` - the signal matrix to analyze
- `fs::Float64` - Sampling rate of the signal
"""
function signal_make_spectrum(signal::Matrix{Float64}, fs)
    channels_no = size(signal, 1)
    signal_fft = zeros(ComplexF64, size(signal))
    signal_sf = zeros(size(signal))

    for idx in 1:channels_no
        signal_fft[idx, :], signal_sf[idx, :] = signal_make_spectrum(signal[idx, :], fs)
    end

    return signal_fft, signal_sf
end

"""
    signal_detrend(signal, type=:linear)

Removes linear trend from the `signal` vector.

# Arguments

- `signals::Vector{Float64}` - the signal vector to analyze
- `type::Symbol[:linear, :constant]`, optional
    - `linear` - the result of a linear least-squares fit to `signal` is subtracted from `signal`
    - `constant` - the mean of `signal` is subtracted
"""
function signal_detrend(signal::Vector{Float64}; trend=:linear)
    trend in [:linear, :constant] || throw(ArgumentError("""Trend type must be ":linear" or ":constant"."""))

    if trend == :constant
        signal_det = demean(signal)
    else
        A = ones(length(signal))
        coef = A \ signal
        signal_det = @. signal - dot(A, coef)
    end

    return signal_det
end

"""
    signal_detrend(signal, type=:linear)

Removes linear trend for each the `signal` matrix channels.

# Arguments

- `signals::Matrix{Float64}` the signal matrix to analyze
- `type::Symbol[:linear, :constant]`, optional
    - `linear` - the result of a linear least-squares fit to `signal` is subtracted from `signal`
    - `constant` - the mean of `signal` is subtracted
"""
function signal_detrend(signal::Matrix{Float64}; trend=:linear)
    trend in [:linear, :constant] || throw(ArgumentError("""Trend type must be ":linear" or ":constant"."""))
    channels_no = size(signal, 1)
    signal_det = zeros(size(signal))

    for idx in 1:channels_no
        signal_det[idx, :] = signal_detrend(signal[idx, :], trend=trend)
    end

    return signal_det
end

"""
    signal_ci95(signal, n=3; method=:normal)

Calculates mean, std and 95% confidence interval for each the `signal` matrix channels.

# Arguments

- `signal::Matrix{Float64}` - the signal matrix to analyze
- `n::Int` - number of bootstraps.
- `method::Symbol[:normal, :boot]` - use normal method or `n`-times boostrapping.
"""
function signal_ci95(signal::Matrix{Float64}; n=3, method=:normal)
    method in [:normal, :boot] || throw(ArgumentError("""Method must be ":normal" or ":boot"."""))

    if method === :normal
        signal_mean = mean(signal, dims=1)'
        signal_sd = std(signal, dims=1)' / sqrt(size(signal, 1))
        upper_bound = signal_mean + 1.96 * signal_sd
        lower_bound = signal_mean - 1.96 * signal_sd
    else
        signal_tmp1 = zeros(size(signal, 1) * n, size(signal, 2))
        Threads.@threads for idx1 in 1:size(signal, 1) * n
            signal_tmp2 = zeros(size(signal))
            sample_idx = rand(1:size(signal, 1), size(signal, 1))
            for idx2 in 1:size(signal, 1)
                signal_tmp2[idx2, :] = signal[sample_idx[idx2], :]'
            end
            signal_tmp1[idx1, :] = mean(signal_tmp2, dims=1)
        end
        signal_mean = mean(signal_tmp1, dims=1)'
        signal_sd = std(signal_tmp1, dims=1)' / sqrt(size(signal_tmp1, 1))
        signal_sorted = sort(signal_tmp1, dims=1)
        lower_bound = signal_sorted[round(Int, 0.025 * size(signal_tmp1, 1)), :]
        upper_bound = signal_sorted[round(Int, 0.975 * size(signal_tmp1, 1)), :]
    end

    return Vector(signals_mean[:, 1]), Vector(signals_sd[:, 1]), Vector(upper_bound[:, 1]), Vector(lower_bound[:, 1])
end

"""
    signal_mean(signal1, signal2)

Calculates mean and 95% confidence interval for 2 signals.

# Arguments

- `signal1::Matrix{Float64}` - the signal 1 matrix to analyze
- `signal2:Matrix{Float64}` - the signal 2 matrix to analyze
"""
function signal_mean(signal1::Matrix{Float64}, signal2::Matrix{Float64})
    signal1_mean = mean(signal1, dims=1)'
    signal2_mean = mean(signal2, dims=1)'
    signals_mean = signal1_mean - signal2_mean
    signal1_sd = std(signal1, dims=1) / sqrt(size(signal1, 1))
    signal2_sd = std(signal2, dims=1) / sqrt(size(signal2, 1))
    signals_mean_sd = sqrt.(signal1_sd.^2 .+ signal2_sd.^2)'

    return Vector(signals_mean[:, 1]), Vector(signals_mean_sd[:, 1]), Vector((signals_mean + 1.96 * signals_mean_sd)[:, 1]), Vector((signals_mean - 1.96 * signals_mean_sd)[:, 1])
end

"""
    signal_difference(signal1::Matrix, signal2::Matrix, n=3; method=:absdiff)

Calculates mean difference and 95% confidence interval for 2 signals.

# Arguments

- `signal1::Matrix` - the signal 1 matrix to analyze
- `signal2:Matrix` - the signal 2 matrix to analyze
- `n::Int` - number of bootstraps.
- `method::Symbol[:absdiff, :diff2int]`
    - `:absdiff` - maximum difference
    - `:diff2int` - integrated area of the squared difference
"""
function signal_difference(signal1::Matrix, signal2::Matrix; n=3, method=:absdiff)
    method in [:absdiff, :diff2int] || throw(ArgumentError("""Method must be ":absdiff" or ":diff2int"."""))
    signal1_mean = mean(signal1, dims=1)'
    signal2_mean = mean(signal2, dims=1)'

    if method === :absdiff
        # statistic: maximum difference
        signals_diff = signal1_mean - signal2_mean
        signals_statistic_single = findmax(abs.(signals_diff))[1]
    else
        # statistic: integrated area of the squared difference
        signals_diff_squared = (signal1_mean - signal2_mean).^2
        signals_statistic_single = simpson(signals_diff_squared)
    end

    signals = [signal1; signal2]
    signals_statistic = zeros(size(signal1, 1) * n)

    Threads.@threads for idx1 in 1:(size(signal1, 1) * n)
        signals_tmp1 = zeros(size(signal1, 1), size(signal1, 2))
        sample_idx = rand(1:size(signals, 1), size(signals, 1))
        # sample_idx = sample_idx[1:1000]
        for idx2 in 1:size(signal1, 1)
            signals_tmp1[idx2, :] = signals[sample_idx[idx2], :]'
        end
        signal1_mean = mean(signals_tmp1, dims=1)
        signals_tmp1 = zeros(size(signal1, 1), size(signal1, 2))
        sample_idx = rand(1:size(signals, 1), size(signals, 1))
        # sample_idx = sample_idx[1:1000]
        for idx2 in 1:size(signal1, 1)
            signals_tmp1[idx2, :] = signals[sample_idx[idx2], :]'
        end
        signal2_mean = mean(signals_tmp1, dims=1)
        if method === :absdiff
            # statistic: maximum difference
            signals_diff = signal1_mean - signal2_mean
            signals_statistic[idx1] = findmax(abs.(signals_diff))[1]
        else
            # statistic: integrated area of the squared difference
            signals_diff_squared = (signal1_mean - signal2_mean).^2
            signals_statistic[idx1] = simpson(signals_diff_squared)
        end
    end

    p = length(signals_statistic[signals_statistic .> signals_statistic_single]) / size(signal1, 1) * n

    return signals_statistic, signals_statistic_single, p
end

"""
   signal_autocov(signal, lag=1, demean=false, normalize=false)

Calculates autocovariance of the `signal` vector.

# Arguments

- `signal::Vector{Float64}` - the signal vector to analyze
- `lag::Int` - lags range is `-lag:lag`
- `demean::Bool[true, false]` - demean signal
- `normalize::Bool[true, false]` - normalize signal
"""
function signal_autocov(signal::Vector{Float64}; lag=1, demean=false, normalize=false)
    signal_lags = collect(-lag:lag)

    if demean == true
        signal_demeaned = signal .- mean(signal)
    else
        signal_demeaned = signal
    end

    signal_ac = zeros(length(signal_lags))

    for idx in 1:length(signal_lags)
        if signal_lags[idx] == 0
            # no lag
            signal_lagged = signal_demeaned
            signals_mul = signal_demeaned .* signal_lagged
        elseif signal_lags[idx] > 0
            # positive lag
            signal_lagged = signal_demeaned[1:(end - signal_lags[idx])]
            signals_mul = signal_demeaned[(1 + signal_lags[idx]):end] .* signal_lagged
        elseif signal_lags[idx] < 0
            # negative lag
            signal_lagged = signal_demeaned[(1 + abs(signal_lags[idx])):end]
            signals_mul = signal_demeaned[1:(end - abs(signal_lags[idx]))] .* signal_lagged
        end
        signals_sum = sum(signals_mul)
        if normalize == true
            signal_ac[idx] = signals_sum / length(signal)
        else
            signal_ac[idx] = signals_sum
        end
    end

    return signal_ac, signal_lags
end

"""
   signal_autocov(signal, lag=1, demean=false, normalize=false)

Calculates autocovariance of each the `signal` matrix channels.

# Arguments

- `signal::Matrix{Float64}` - the signal vector to analyze
- `lag::Int` - lags range is `-lag:lag`
- `demean::Bool[true, false]` - demean signals
- `normalize::Bool[true, false]` - normalize signals
"""
function signal_autocov(signal::Matrix{Float64}; lag=1, demean=false, normalize=false)
    signal_lags = collect(-lag:lag)
    channels_no = size(signal, 1)
    signal_ac = zeros(channels_no, length(signal_lags))

    for idx in 1:channels_no
        signal_ac[idx, :], _ = signal_autocov(signal[idx, :], lag=lag, demean=demean, normalize=normalize)
    end

    return signal_ac, signal_lags
end

"""
   signal_crosscov(signal1, signal2, lag=1, demean=false, normalize=false)

Calculates cross-covariance between `signal1` and `signal2` vectors.

# Arguments

- `signal1::Vector{Float64}` - the signal 1 vector to analyze
- `signal2::Vector{Float64}` - the signal 2 vector to analyze
- `lag::Int` - lags range is `-lag:lag`
- `demean::Bool[true, false]` - demean signals
- `normalize::Bool[true, false]` - normalize signals
"""
function signal_crosscov(signal1::Vector{Float64}, signal2::Vector{Float64}; lag=1, demean=false, normalize=false)
    lags = collect(-lag:lag)

    if demean == true
        signal_demeaned1 = signal1 .- mean(signal1)
        signal_demeaned2 = signal2 .- mean(signal2)
    else
        signal_demeaned1 = signal1
        signal_demeaned2 = signal2
    end

    ac = zeros(length(lags))

    for idx in 1:length(lags)
        if lags[idx] == 0
            # no lag
            signal_lagged = signal_demeaned2
            signals_mul = signal_demeaned1 .* signal_lagged
        elseif lags[idx] > 0
            # positive lag
            signal_lagged = signal_demeaned2[1:(end - lags[idx])]
            signals_mul = signal_demeaned1[(1 + lags[idx]):end] .* signal_lagged
        elseif lags[idx] < 0
            # negative lag
            signal_lagged = signal_demeaned2[(1 + abs(lags[idx])):end]
            signals_mul = signal_demeaned1[1:(end - abs(lags[idx]))] .* signal_lagged
        end
        signals_sum = sum(signals_mul)
        if normalize == true
            ac[idx] = signals_sum / length(signal1)
        else
            ac[idx] = signals_sum
        end
    end

    return ac, lags
end

"""
   signal_crosscov(signal1, signal2, lag=1, demean=false, normalize=false)

Calculates cross-covariance between same channels in `signal1` and `signal2` matrices.

# Arguments

- `signal1::Matrix{Float64}` - the signal 1 matrix to analyze
- `signal2::Matrix{Float64}` - the signal 2 matrix to analyze
- `lag::Int` - lags range is `-lag:lag`
- `demean::Bool[true, false]` - demean signals
- `normalize::Bool[true, false]` - normalize signals
"""
function signal_crosscov(signal1::Matrix{Float64}, signal2::Matrix{Float64}; lag=1, demean=false, normalize=false)
    signal_lags = collect(-lag:lag)
    channels_no = size(signal, 1)
    signal_ac = zeros(channels_no, length(signal_lags))

    for idx in 1:channels_no
        signal_ac[idx, :], _ = signal_crosscov(signal1[idx, :], signal2[idx, :], lag=lag, demean=demean, normalize=normalize)
    end

    return signal_ac, signal_lags
end

"""
   signal_crosscov(signal1, lag=1, demean=false, normalize=false)

Calculates cross-covariance for all channels in `signal` matrix.

# Arguments

- `signal::Matrix{Float64}` - the signal matrix to analyze
- `lag::Int` - lags range is `-lag:lag`
- `demean::Bool[true, false]` - demean signals
- `normalize::Bool[true, false]` - normalize signals
"""
function signal_crosscov(signal::Matrix{Float64}; lag=1, demean=false, normalize=false)
    signal_lags = collect(-lag:lag)
    channels_no = size(signal, 1)
    signal_ac = zeros(channels_no^2, length(signal_lags))

    for idx in 1:channels_no^2
        for idx1 in 1:channels_no, idx2 in 1:channels_no
            signal_ac[idx, :], _ = signal_crosscov(signal[idx1, :], signal[idx2, :], lag=lag, demean=demean, normalize=normalize)
        end
    end

    return signal_ac, signal_lags
end

"""
    signal_spectrum(signal, pad=0, remove_dc=false, detrend=false, taper=nothing)
"""
function signal_spectrum(signal::Vector{Float64}, pad::Int=0, remove_dc=false, detrend=false, derivative=false, taper=nothing)
    pad < 0 && throw(ArgumentError("""Value of "pad" cannot be negative."""))
    remove_dc == true && signal = demean(signal)
    detrend == true && (signal = signal_detrend(signal))
    derivative == true && (signal = signal_derivative(signal))
    taper != nothing && (signal = signal .* taper)

    if pad == 0
        signal_fft = fft(signal)
    else
        signal_fft = fft0(signal, pad)
    end

    # normalize
    signal_fft ./= length(signal)

    # amplitudes
    signal_amplitudes = @. 2 * abs(signal_fft)

    # power
    signal_powers = signal_amplitudes.^2

    # phases
    signal_phases = atan.(imag(signal_fft), real(signal_fft))

    return signal_fft, signal_amplitudes, signal_powers, signal_phases
end

"""
    eeg_load(in_file, read_annotations=true, header_only=false, clean_labels=true)

Loads EDF/EDFPlus/BDF/BDFPlus file and returns eeg_signal::eeg object.

 8 ascii : version of this data format (0)
80 ascii : local patient identification (mind item 3 of the additional EDF+ specs)
80 ascii : local recording identification (mind item 4 of the additional EDF+ specs)
 8 ascii : startdate of recording (dd.mm.yy) (mind item 2 of the additional EDF+ specs)
 8 ascii : starttime of recording (hh.mm.ss)
 8 ascii : number of bytes in header record
44 ascii : reserved
 8 ascii : number of data records (-1 if unknown, obey item 10 of the additional EDF+ specs)
 8 ascii : duration of a data record, in seconds
 4 ascii : number of signals (ns) in data record

ns * 16 ascii : ns * label (e.g. EEG Fpz-Cz or Body temp) (mind item 9 of the additional EDF+ specs)
ns * 80 ascii : ns * transducer type (e.g. AgAgCl electrode)
ns * 8 ascii : ns * physical dimension (e.g. uV or degreeC)
ns * 8 ascii : ns * physical minimum (e.g. -500 or 34)
ns * 8 ascii : ns * physical maximum (e.g. 500 or 40)
ns * 8 ascii : ns * digital minimum (e.g. -2048)
ns * 8 ascii : ns * digital maximum (e.g. 2047)
ns * 80 ascii : ns * prefiltering (e.g. HP:0.1Hz LP:75Hz)
ns * 8 ascii : ns * nr of samples in each data record
ns * 32 ascii : ns * reserved

DATA RECORD
nr of samples[1] * integer : first signal in the data record
nr of samples[2] * integer : second signal
..
..
nr of samples[ns] * integer : last signal 

samplingrate <- n.samples / data.record.duration
gain = (physicalmaximum - physicalminimum) / (digitalmaximum - digitalminimum)
value = (value - digitalminimum ) * gain + physicalminimum

Source: Kemp B, Värri A, Rosa AC, Nielsen KD, Gade J. A simple format for exchange of digitized polygraphic recordings. Electroencephalography and Clinical Neurophysiology. 1992 May;82(5):391–3. 

"""
function eeg_load(in_file, read_annotations=true, header_only=false, clean_labels=false)
    fid = open(in_file)

    header = zeros(UInt8, 256)
    readbytes!(fid, header, 256)
    header = String(Char.(header))

    version = parse(Int, rstrip(header[1:8]))
    version == 0 && (eeg_filetype = "EDF")
    patient = rstrip(header[9:88])
    recording = rstrip(header[89:168])
    recording_date = header[169:176]
    recording_time = header[177:184]
    data_offset = parse(Int, rstrip(header[185:192]))
    reserved  = header[193:236]
    data_records = parse(Int, rstrip(header[237:244]))
    data_records_duration  = parse(Float64, rstrip(header[245:252]))
    channels_no  = parse(Int, rstrip(header[253:256]))

    labels = Vector{String}(undef, channels_no)
    transducers = Vector{String}(undef, channels_no)
    physical_dimension = Vector{String}(undef, channels_no)
    physical_minimum = Vector{Float64}(undef, channels_no)
    physical_maximum = Vector{Float64}(undef, channels_no)
    digital_minimum = Vector{Float64}(undef, channels_no)
    digital_maximum = Vector{Float64}(undef, channels_no)
    prefiltering = Vector{String}(undef, channels_no)
    samples_per_datarecord = Vector{Int64}(undef, channels_no)

    header = zeros(UInt8, channels_no * 16)
    readbytes!(fid, header, channels_no * 16)
    header = String(Char.(header))
    for idx in 1:channels_no
        labels[idx] = rstrip(header[1 + ((idx - 1) * 16):(idx * 16)])
    end
    header = zeros(UInt8, channels_no * 80)
    readbytes!(fid, header, channels_no * 80)
    header = String(Char.(header))
    for idx in 1:channels_no
        transducers[idx] = rstrip(header[1 + ((idx - 1) * 80):(idx * 80)])
    end
    header = zeros(UInt8, channels_no * 8)
    readbytes!(fid, header, channels_no * 8)
    header = String(Char.(header))
    for idx in 1:channels_no
        physical_dimension[idx] = rstrip(header[1 + ((idx - 1) * 8):(idx * 8)])
    end
    header = zeros(UInt8, channels_no * 8)
    readbytes!(fid, header, channels_no * 8)
    header = String(Char.(header))
    for idx in 1:channels_no
        physical_minimum[idx] = parse(Float64, rstrip(header[1 + ((idx - 1) * 8):(idx * 8)]))
    end
    header = zeros(UInt8, channels_no * 8)
    readbytes!(fid, header, channels_no * 8)
    header = String(Char.(header))
    for idx in 1:channels_no
        physical_maximum[idx] = parse(Float64, rstrip(header[1 + ((idx - 1) * 8):(idx * 8)]))
    end
    header = zeros(UInt8, channels_no * 8)
    readbytes!(fid, header, channels_no * 8)
    header = String(Char.(header))
    for idx in 1:channels_no
        digital_minimum[idx] = parse(Float64, rstrip(header[1 + ((idx - 1) * 8):(idx * 8)]))
    end
    header = zeros(UInt8, channels_no * 8)
    readbytes!(fid, header, channels_no * 8)
    header = String(Char.(header))
    for idx in 1:channels_no
        digital_maximum[idx] = parse(Float64, rstrip(header[1 + ((idx - 1) * 8):(idx * 8)]))
    end
    header = zeros(UInt8, channels_no * 80)
    readbytes!(fid, header, channels_no * 80)
    header = String(Char.(header))
    for idx in 1:channels_no
        prefiltering[idx] = rstrip(header[1 + ((idx - 1) * 80):(idx * 80)])
    end
    header = zeros(UInt8, channels_no * 8)
    readbytes!(fid, header, channels_no * 8)
    header = String(Char.(header))
    for idx in 1:channels_no
        samples_per_datarecord[idx] = parse(Int, rstrip(header[1 + ((idx - 1) * 8):(idx * 8)]))
    end

    sampling_rate = Vector{Float64}(undef, channels_no)
    gain = Vector{Float64}(undef, channels_no)
    for idx in 1:channels_no
        sampling_rate[idx] = samples_per_datarecord[idx] / data_records_duration
        gain[idx] = (physical_maximum[idx] - physical_minimum[idx]) / (digital_maximum[idx] - digital_minimum[idx])
    end

    if clean_labels == true
        labels = replace.(labels, "EEG " => "")
        labels = replace.(labels, "ECG " => "")
    end

    close(fid)

    fid = open(in_file)
    header = zeros(UInt8, data_offset)
    readbytes!(fid, header, data_offset)
    eeg_signals = zeros(channels_no, samples_per_datarecord[1] * data_records)

    for idx1 in 1:data_records
        for idx2 in 1:channels_no
            signal = zeros(UInt8, samples_per_datarecord[idx2] * 2);
            readbytes!(fid, signal, samples_per_datarecord[idx2] * 2);
            signal = map(ltoh, reinterpret(Int16, signal));
            eeg_signals[idx2, ((idx1 - 1) * samples_per_datarecord[idx2] + 1):(idx1 * samples_per_datarecord[idx2])] = @. (signal - digital_minimum[idx2]) * gain[idx2] + physical_minimum[idx2];
        end
    end

    close(fid)

    eeg_file_header = Dict(:version => version, :eeg_filetype => eeg_filetype, :patient => patient, :recording => recording, :recording_date => recording_date, :recording_time => recording_time, :data_records => data_records, :data_records_duration => data_records_duration, :channels_no => channels_no)
    eeg_signal_header = Dict(:labels => labels, :transducers => transducers, :physical_dimension => physical_dimension, :physical_minimum => physical_minimum, :physical_maximum => physical_maximum, :digital_minimum => digital_minimum, :digital_maximum => digital_maximum, :prefiltering => prefiltering, :samples_per_datarecord => samples_per_datarecord, :sampling_rate => sampling_rate, :gain => gain)
    eeg = EEG(eeg_file_header, eeg_signal_header, eeg_signals)
    return eeg
end

"""
    signals_epoch_avg(signals, n)

Divides `signals` matrix into `n`-samples long averaged epochs and average them.
"""
function signals_epoch_avg(signals::Matrix, n)
    channels = size(signals, 1)
    m = size(signals, 2)
    no_epochs = m ÷ n
    epochs = zeros(channels, n)
    for idx1 in 1:channels
        epoch_tmp = zeros(no_epochs, n)
        for idx2 in 1:no_epochs
            epoch_tmp[idx2, :] = signals[idx1, (((idx2 - 1) * n) + 1):(idx2 * n)]
        end
        epochs[idx1, :] = mean(epoch_tmp, dims=1)
    end
    return epochs
end

"""
    signal_epoch_avg(signals, n)

Divides `signal` vector into `n`-samples long averaged epochs and average them.
"""
function signal_epoch_avg(signal::Vector{Float64}, n)
    m = length(signal)
    no_epochs = m ÷ n
    epochs = zeros(n)
    epoch_tmp = zeros(no_epochs, n)
    for idx in 1:no_epochs
        epoch_tmp[idx, :] = signal[(((idx - 1) * n) + 1):(idx * n)]
    end
    epochs = Vector(mean(epoch_tmp, dims=1)[1, :])
    return epochs
end

"""
    signal_filter_butter(signal, type, cutoff, fs, poles=8)

Filters `signal` vector using filter `type`=[:lp, :hp, :bp, :bs], `cutoff` in Hz, `fs` sampling rate and `poles`-pole Butterworth filter.
"""
function signal_filter_butter(signal::Vector{Float64}, type, cutoff, fs, poles=8)
    if type == :lp
        responsetype = Lowpass(cutoff; fs=fs)
        prototype = Butterworth(poles)
    elseif type == :hp
        responsetype = Highpass(cutoff; fs=fs)
        prototype = Butterworth(poles)
    elseif type == :bp
        responsetype = Bandpass(cutoff[1], cutoff[2]; fs=fs)
        prototype = Butterworth(poles)
    elseif type == :bs
        responsetype = Bandstop(cutoff[1], cutoff[2]; fs=fs)
        prototype = Butterworth(poles)
    end
    filter = digitalfilter(responsetype, prototype)
    signal_filtered = filt(filter, signal)
    return signal_filtered
end

"""
    signals_filter_butter(signal, type, cutoff, fs, poles=8)

Filters `signals` matrix using filter `type` [:lp, :hp, :bp, :bs], `cutoff` in Hz, `fs` sampling rate and `poles`-pole Butterworth filter.
"""
function signals_filter_butter(signals::Matrix, type, cutoff, fs, poles=8)
    no_channels = size(signals, 1)
    signals_filtered = zeros(size(signals))
    for idx in 1:no_channels
        signals_filtered[idx, :] = signal_filter(signals[idx, :], type, cutoff, fs, poles)
    end
    return signals_filtered
end

"""
    signal_plot(time, signal)

Plots `signal` vector.
"""
function signal_plot(time, signal)
    amplitude_min, _ = findmin(signal)
    amplitude_min = ceil(Int64, amplitude_min)
    amplitude_max, _ = findmax(signal)
    amplitude_max = floor(Int64, amplitude_max)
    plot(time, signal, ylim=(amplitude_min, amplitude_max))
end

"""
    signal_plot(t, signals; labels=[], rescale=false, xlabel="Time [s]", ylabel="Amplitude [μV]", yamp=100)

Plots `signal` vector.
"""
function signal_plot(t, signal; labels=[], rescale=false, xlabel="Time [s]", ylabel="Amplitude [μV]", yamp=100)
    if rescale == true
        # rescale
        signal = (signal .- mean(signal)) ./ std(signal)
        yamp = ceil(findmax(signal)[1])
    end
    p = plot(t, signal, xlabel="Time [s]", ylabel="Amplitude [μV]", legend=false, t=:line, c=:black, ylims=(-yamp, yamp))
    return p
end

"""
    signals_plot(t, signals; labels=[], rescale=false, xlabel="Time [s]", ylabel="Channels")

Plots `signals` matrix.
"""
function signals_plot(t, signals; labels=[], rescale=true, xlabel="Time [s]", ylabel="Channels")
    no_channels = size(signals, 1)

    # reverse so 1st channel is on top
    signals = reverse(signals, dims = 1)

    if rescale == true
        # rescale and shift so all channels are visible
        variances = var(signals; dims=2)
        mean_variance = mean(variances)
        for idx in 1:no_channels
            signals[idx, :] = (signals[idx, :] .- mean(signals[idx, :])) ./ mean_variance .+ (idx - 1)
        end
    end

    p = plot(xlabel=xlabel, ylabel=ylabel, ylim=(-0.5, no_channels-0.5))
    for idx in 1:no_channels
        # Rescale and shift so all channels are visible
        p = plot!(t, signals[idx, :], legend=false, t=:line, c=:black)
    end
    p = plot!(p, yticks = (no_channels-1:-1:0, labels))
    return p
end

"""
    eeg_drop_channel(eeg, channels)

Removes `channels` from the `eeg` set.
"""
function eeg_drop_channel(eeg, channels)
    if typeof(channels) == UnitRange{Int64}
        channels = collect(channels)
    end
    channels = sort!(channels, rev=true)

    eeg_file_header = eeg.eeg_file_header
    eeg_signal_header = eeg.eeg_signal_header
    eeg_signals = eeg.eeg_signals

    channels_no = eeg_file_header[:channels_no]

    # update headers
    eeg_file_header[:channels_no] = channels_no - length(channels)
    for idx1 in 1:length(channels)
        for idx2 in 1:channels_no
            if idx2 == channels[idx1]
                deleteat!(eeg_signal_header[:labels], idx2)
                deleteat!(eeg_signal_header[:transducers], idx2)
                deleteat!(eeg_signal_header[:physical_dimension], idx2)
                deleteat!(eeg_signal_header[:physical_minimum], idx2)
                deleteat!(eeg_signal_header[:physical_maximum], idx2)
                deleteat!(eeg_signal_header[:digital_minimum], idx2)
                deleteat!(eeg_signal_header[:digital_maximum], idx2)
                deleteat!(eeg_signal_header[:prefiltering], idx2)
                deleteat!(eeg_signal_header[:samples_per_datarecord], idx2)
                deleteat!(eeg_signal_header[:sampling_rate], idx2)
                deleteat!(eeg_signal_header[:gain], idx2)
            end
        end 
    end

    # remove channels
    eeg_signals = eeg_signals[setdiff(1:end, (channels)), :]

    # create new dataset    
    eeg_new = EEG(eeg_file_header, eeg_signal_header, eeg_signals)
    return eeg_new
end