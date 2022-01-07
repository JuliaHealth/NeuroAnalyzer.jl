import os
from autoreject import AutoReject
import mne
import numpy as np

eeg_file = os.path.expanduser('~/Documents')+"/Data/test-data/eeg/eeg-test.edf"
eeg_file = os.path.expanduser('~/Documents')+"/Data/tdcs-avh/eeg/29-06.edf"

print("Loading EDF file..")
eeg_data = mne.io.read_raw_edf(eeg_file, preload=True)

# montage
eeg_data.rename_channels({'EEG Fp1': 'Fp1', 'EEG Fp2': 'Fp2', 'EEG F3': 'F3', 'EEG F4': 'F4', 'EEG C3': 'C3', 'EEG C4': 'C4', 'EEG P3': 'P3', 'EEG P4': 'P4', 'EEG O1': 'O1', 'EEG O2': 'O2', 'EEG F7': 'F7', 'EEG F8': 'F8', 'EEG T3': 'T3', 'EEG T4': 'T4', 'EEG T5': 'T5', 'EEG T6': 'T6', 'EEG Fz': 'Fz', 'EEG Cz': 'Cz', 'EEG Pz': 'Pz'})
eeg_data.rename_channels({'EEG Fp1': 'Fp1', 'EEG Fp2': 'Fp2', 'EEG F3': 'F3', 'EEG F4': 'F4', 'EEG C3': 'C3', 'EEG C4': 'C4', 'EEG P3': 'P3', 'EEG P4': 'P4', 'EEG O1': 'O1', 'EEG O2': 'O2', 'EEG F7': 'F7', 'EEG F8': 'F8', 'EEG T3': 'T3', 'EEG T4': 'T4', 'EEG T5': 'T5', 'EEG T6': 'T6', 'EEG Fz': 'Fz', 'EEG Cz': 'Cz', 'EEG Pz': 'Pz', 'ECG S1': 'ECG'})

eeg_data.set_channel_types({'EEG A1': 'misc', 'EEG A2': 'misc', 'ECG': 'ecg'})

ecg_data = eeg_data.pick_channels(['ECG'])

montage = mne.channels.make_standard_montage("standard_1020")
eeg_data.set_montage(montage)

# average referencing
eeg_data = eeg_data.copy().set_eeg_reference(ref_channels='average')

# REST referencing
eeg_data.del_proj()  # remove our average reference projector first
sphere = mne.make_sphere_model('auto', 'auto', eeg_data.info)
src = mne.setup_volume_source_space(sphere=sphere, exclude=30., pos=15.)
forward = mne.make_forward_solution(eeg_data.info, trans=None, src=src, bem=sphere)
eeg_data = eeg_data.copy().set_eeg_reference('REST', forward=forward)

# power line NF
eeg_picks = mne.pick_types(eeg_data.info, eeg=True)
power_line_freqs = (50, 100)
eeg_data = eeg_data.copy().notch_filter(freqs=power_line_freqs, picks=eeg_picks)

# 1 Hz HP
eeg_data = eeg_data.copy().filter(l_freq=1, h_freq=40)

# remove eye blinks and eye movements: ICA
ica = mne.preprocessing.ICA(n_components=12, random_state=1, max_iter=800)
ica.fit(eeg_data)
# check bad components
ica.plot_components(outlines='skirt')
# store bad components
bad_idx = [2, 6, 7, 8, 10, 11]
bad_idx = [0, 4, 8, 10]
ica.plot_properties(eeg_data, picks=bad_idx)
# automatic algorithm: identify eye artifacts based on channel
bad_idx, scores = ica.find_bads_eog(eeg_data, ['Fp1', 'Fp2'], threshold=1.5)
print(bad_idx)
eeg_data = ica.apply(eeg_data.copy(), exclude=bad_idx)

eeg_data.plot_psd(fmax=50, average=True)
eeg_data.plot()
eeg_data.plot_psd_topo(fmax=128)

# 10s epochs
epochs = mne.make_fixed_length_epochs(eeg_data, duration=5, preload=True)
ar = AutoReject()
epochs_clean = ar.fit_transform(epochs)
epochs_clean.plot(block=True)

from autoreject import Ransac
rsc = Ransac()
epochs_clean = rsc.fit_transform(epochs)

# automatic algorithm: identify eye artifacts based on channel
bad_idx, scores = ica.find_bads_eog(epochs, ['Fp1', 'Fp2'], threshold=1.5)
print(bad_idx)
epochs = ica.apply(epochs.copy(), exclude=bad_idx)

# ECG artifacts
ecg_picks = mne.pick_types(ecg_data.info, eeg=True)
power_line_freqs = (50, 100)
ecg_data = ecg_data.copy().notch_filter(freqs=power_line_freqs, picks=ecg_picks)

# 1 Hz HP
# convert μV to mV
ecg_data.data = ecg_data.get_data()*1000
ecg_data = eeg_data.copy().filter(l_freq=1, h_freq=40)
ecg_data.plot()
import matplotlib.pyplot as plt
plt.plot(ecg_data.times, ecg_data.get_data().flatten())
plt.show()
average_ecg = mne.preprocessing.create_ecg_epochs(eeg_data).average()
print('We found %i ECG events' % average_ecg.nave)
average_ecg.plot_joint()