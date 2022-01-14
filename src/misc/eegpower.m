function [total_power, delta_band, theta_band, alpha_band, alpha1_band, alpha2_band, beta_band, beta1_band, beta2_band, gamma_band, gamma1_band, gamma2_band] = eegpower(eeg_data, fs)

    nyquist_freq = fs/2;

    total_power = bandpower(eeg_data, fs, [0 nyquist_freq]);
    delta_band = bandpower(eeg_data, fs, [0 4]);
    theta_band = bandpower(eeg_data, fs, [4 8]);
    alpha_band = bandpower(eeg_data, fs, [8 13]);
    alpha1_band = bandpower(eeg_data, fs, [8 10]);
    alpha2_band = bandpower(eeg_data, fs, [10 13]);
    beta_band = bandpower(eeg_data, fs, [13 31]);
    beta1_band = bandpower(eeg_data, fs, [13 18]);
    beta2_band = bandpower(eeg_data, fs, [18 31]);
    gamma_band = bandpower(eeg_data, fs, [31 nyquist_freq]);
    gamma1_band = bandpower(eeg_data, fs, [31 40.99]);
    gamma2_band = bandpower(eeg_data, fs, [41 nyquist_freq]);
end