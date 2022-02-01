function clean_eeg(eeg_channel, lp = 100, hp = 0.1, ac = 50, fs = 256)
    % CLEAN_EEG(eeg_channel, low_pass, high_pass, ac_frequency)
    % Filters single EEG channel (eeg_channel) at AC frequency (ac_frequency [50 Hz]), low pass at low_pass and high pass at high_pass frequencies.

    # filter BS 50 / 60 Hz
    if ac != 0
        responsetype = Bandstop(ac-10, ac+10; fs);
        designmethod = Butterworth(8);
        eeg_filtered = filt(digitalfilter(responsetype, designmethod), eeg_channel);
    else
        eeg_filtered = eeg_channel
    end

    # filter LP 20 Hz
    if lp != 0
        responsetype = Lowpass(lp; fs);
        designmethod = Butterworth(8);
        eeg_filtered = filt(digitalfilter(responsetype, designmethod), eeg_filtered);
    end

    # filter HP 0.1 Hz
    if hp != 0
        responsetype = Highpass(hp; fs);
        designmethod = Butterworth(8);
        eeg_filtered = filt(digitalfilter(responsetype, designmethod), eeg_filtered);
    end

    return eeg_filtered

end
