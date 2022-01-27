"""
    signal_derivative(signal)

Returns the derivative of the `signal` with length same as the signal
"""
signal_derivative(signal::Vector) = vcat(diff(signal), diff(signal)[end])

"""
    band_power(psd, f1, f2)

Calculates absolute band power between frequencies `f1` and `f2` for the `psd` power spectrum.
"""
function signal_band_power(psd, f1, f2)
    frq_idx = [vsearch(psd.freq, f1), vsearch(psd.freq, f2)]
    # dx: frequency resolution
    dx = psd.freq[2] - psd.freq[1]
    result = simpson(psd.power[frq_idx[1]:frq_idx[2]], start=frq_idx[1], stop=frq_idx[1], dx=dx)
    return result
end

"""
    make_spectrum(signal, fs)

Returns FFT and DFT sample frequencies for a DFT of the `x` signal 
"""
function signal_make_spectrum(signal::Vector, fs)
    hs = fft(signal)
    n = length(signal)               # number of samples
    d = 1 / fs                    # time between samples
    fs = fftfreq(n, d)
    return hs, fs
end

"""
    signal_detrend(signal, type=:linear)

Removes linear trend of `signal`.

# Arguments
- `x::Vector` - the signal to detrend.
- `type::Symbol[:linear, :constant]`, optional.
    linear: the result of a linear least-squares fit to `y` is subtracted from `y`
    constant: the mean of `y` is subtracted.
"""
function signal_detrend(signal::Vector; trend=:linear)
    trend in [:linear, :constant] || throw(ArgumentError("""Trend type must be ":linear" or ":constant"."""))
    if trend == :constant
        result = signal .- mean(signal)
    else
        A = ones(length(signal))
        coef = A \ signal
        result = @. signal - dot(A, coef)
    end
    return result
end

"""
    signals_ci95(signals::Matrix, n=3; method=:normal)

Calculates mean, std and 95% confidence interval for the matrix of signals .

# Arguments

- `signals::Matrix` - the signal matrix to analyze (rows: trials, columns: time).
- `n::Int` - number of bootstraps.
- `method::Symbol[:normal, :boot]` - use normal method or `n`-times boostrapping.
"""
function signals_ci95(signals::Matrix, n=3; method=:normal)
    method in [:normal, :boot] || throw(ArgumentError("""Method must be ":normal" or ":boot"."""))
    if method === :normal
        signals_mean = mean(signals, dims=1)'
        signals_sd = std(signals, dims=1)' / sqrt(size(signals, 1))
        upper_bound = signals_mean + 1.96 * signals_sd
        lower_bound = signals_mean - 1.96 * signals_sd
    else
        signals_tmp1 = zeros(size(signals, 1) * n, size(signals, 2))
        Threads.@threads for idx1 in 1:size(signals, 1) * n
            signals_tmp2 = zeros(size(signals))
            sample_idx = rand(1:size(signals, 1), size(signals, 1))
            for idx2 in 1:size(signals, 1)
                signals_tmp2[idx2, :] = signals[sample_idx[idx2], :]'
            end
            signals_tmp1[idx1, :] = mean(signals_tmp2, dims=1)
        end
        signals_mean = mean(signals_tmp1, dims=1)'
        signals_sd = std(signals_tmp1, dims=1)' / sqrt(size(signals_tmp1, 1))
        signal_sorted = sort(signals_tmp1, dims=1)
        lower_bound = signal_sorted[round(Int, 0.025 * size(signals_tmp1, 1)), :]
        upper_bound = signal_sorted[round(Int, 0.975 * size(signals_tmp1, 1)), :]
    end
    return signals_mean, signals_sd, upper_bound, lower_bound
end

"""
    signals_mean(signals1::Matrix, signals2::Matrix)

Calculates mean difference and 95% confidence interval for 2 signals.

# Arguments

- `signal1::Matrix` - the signal matrix to analyze (rows: trials, columns: time).
- `signal2:Matrix` - the signal matrix to analyze (rows: trials, columns: time).
"""
function signals_mean(signals1::Matrix, signals2::Matrix)
    signals1_mean = mean(signals1, dims=1)'
    signals2_mean = mean(signals2, dims=1)'
    signals_mean = signals1_mean - signals2_mean
    signals1_sd = std(signals1, dims=1) / sqrt(size(signals1, 1))
    signals2_sd = std(signals2, dims=1) / sqrt(size(signals2, 1))
    signals_mean_sd = sqrt.(signals1_sd.^2 .+ signals2_sd.^2)'

    return signals_mean, signals_mean_sd, signals_mean + 1.96 * signals_mean_sd, signals_mean - 1.96 * signals_mean_sd
end

"""
    signals_difference(signals1::Matrix, signals2::Matrix, n=3; method=:absdiff)

Calculates mean difference and 95% confidence interval for 2 signals.

# Arguments

- `signal1::Matrix` - the signal matrix to analyze (rows: trials, columns: time).
- `signal2:Matrix` - the signal matrix to analyze (rows: trials, columns: time).
- `n::Int` - number of bootstraps.
- `method::Symbol[:absdiff, :diff2int]`
"""
function signals_difference(signals1::Matrix, signals2::Matrix, n::Int=3; method=:absdiff)
    method in [:absdiff, :diff2int] || throw(ArgumentError("""Method must be ":absdiff" or ":diff2int"."""))
    signals1_mean = mean(signals1, dims=1)'
    signals2_mean = mean(signals2, dims=1)'

    if method === :absdiff
        # statistic: maximum difference
        signals_diff = signals1_mean - signals2_mean
        signals_statistic_single = findmax(abs.(signals_diff))[1]
    else
        # statistic: integrated area of the squared difference
        signals_diff_squared = (signals1_mean - signals2_mean).^2
        signals_statistic_single = simpson(signals_diff_squared)
    end

    signals = [signals1; signals2]
    signals_statistic = zeros(size(signals1, 1) * n)

    Threads.@threads for idx1 in 1:(size(signals1, 1) * n)
        signals_tmp1 = zeros(size(signals1, 1), size(signals1, 2))
        sample_idx = rand(1:size(signals, 1), size(signals, 1))
        sample_idx = sample_idx[1:1000]
        for idx2 in 1:size(signals1, 1)
            signals_tmp1[idx2, :] = signals[sample_idx[idx2], :]'
        end
        signals1_mean = mean(signals_tmp1, dims=1)
        signals_tmp1 = zeros(size(signals1, 1), size(signals1, 2))
        sample_idx = rand(1:size(signals, 1), size(signals, 1))
        sample_idx = sample_idx[1:1000]
        for idx2 in 1:size(signals1, 1)
            signals_tmp1[idx2, :] = signals[sample_idx[idx2], :]'
        end
        signals2_mean = mean(signals_tmp1, dims=1)
        if method === :absdiff
            # statistic: maximum difference
            signals_diff = signals1_mean - signals2_mean
            signals_statistic[idx1] = findmax(abs.(signals_diff))[1]
        else
            # statistic: integrated area of the squared difference
            signals_diff_squared = (signals1_mean - signals2_mean).^2
            signals_statistic[idx1] = simpson(signals_diff_squared)
        end
    end

    p = length(signals_statistic[signals_statistic .> signals_statistic_single]) / size(signals1, 1) * n

    return signals_statistic, signals_statistic_single, p
end

"""
   signal_autocov(signal, lag=1, demean=false, normalize=false)

Calculates autocovariance of the `signal` vector for lags = -lag:lag.
"""
function signal_autocov(signal::Vector, lag=1, demean=false, normalize=false)
    lags = collect(-lag:lag)

    if demean == true
        signal_demeaned = signal .- mean(signal)
    else
        signal_demeaned = signal
    end

    ac = zeros(length(lags))

    for idx in 1:length(lags)
        if lags[idx] == 0
            # no lag
            signal_lagged = signal_demeaned
            signals_mul = signal_demeaned .* signal_lagged
        elseif lags[idx] > 0
            # positive lag
            signal_lagged = signal_demeaned[1:(end - lags[idx])]
            signals_mul = signal_demeaned[(1 + lags[idx]):end] .* signal_lagged
        elseif lags[idx] < 0
            # negative lag
            signal_lagged = signal_demeaned[(1 + abs(lags[idx])):end]
            signals_mul = signal_demeaned[1:(end - abs(lags[idx]))] .* signal_lagged
        end
        signals_sum = sum(signals_mul)
        if normalize == true
            ac[idx] = signals_sum / length(signal)
        else
            ac[idx] = signals_sum
        end
    end

    return ac, lags
end

"""
   signals_crosscov(signal1, signal2, lag=1, demean=false, normalize=false)

Calculates cross-covariance between `signal1` and `signal2` vectors for lags = -lag:lag.
"""
function signals_crosscov(signal1::Vector, signal2::Vector, lag=1, demean=false, normalize=false)
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
    signal_spectrum(signal, pad=0, remove_dc=false, detrend=false, taper=nothing)
"""
function signal_spectrum(signal::Vector, pad::Int=0, remove_dc=false, detrend=false, derivative=false, taper=nothing)
    pad < 0 && throw(ArgumentError("""Value of "pad" cannot be negative."""))
    remove_dc == true && (signal = signal .- mean(signal))
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
function eeg_load(in_file, read_annotations=true, header_only=false, clean_labels=true)
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
    eeg_signals = zeros(channels_no, samples_per_datarecord[idx] * data_records)

    for idx1 in 1:data_records
        for idx2 in 1:channels_no
            signal = zeros(UInt8, samples_per_datarecord[idx2] * 2);
            readbytes!(fid, signal, samples_per_datarecord[idx2] * 2);
            signal = map(ltoh, reinterpret(Int16, signal));
            eeg_signals[idx2, ((idx1 - 1) * samples_per_datarecord[idx2] + 1):(idx1 * samples_per_datarecord[idx2])] = @. (signal - digital_minimum[idx2]) * gain[idx2] + physical_minimum[idx2];
        end
    end

    close(fid)

    eeg_file_header = [version, eeg_filetype, patient, recording, recording_date, recording_time, data_records,     data_records_duration, channels_no]
    eeg_signal_header = DataFrame(:labels => labels, :transducers => transducers, :physical_dimension => physical_dimension, :physical_minimum => physical_minimum, :physical_maximum => physical_maximum, :digital_minimum => digital_minimum, :digital_maximum => digital_maximum, :prefiltering => prefiltering, :samples_per_datarecord => samples_per_datarecord)
    eeg = EEG(eeg_file_header, eeg_signal_header, eeg_signals)
    return eeg
end