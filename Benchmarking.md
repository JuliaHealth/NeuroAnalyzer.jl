# NeuroAnalyzer.jl benchmarking

```julia
@time_imports using NeuroAnalyzer
using BenchmarkTools

# IO
@benchmark edf = eeg_import_edf("test/eeg-test-edfplus.edf")
@benchmark bdf = eeg_import_edf("test/eeg-test-bdf.bdf")
@benchmark edf = eeg_import_edf("test/eeg-test-edf.edf")
@benchmark eeg_delete_channel!(edf, channel=[17, 18, 22, 23, 24])

# EDIT
@benchmark edf_epoched = eeg_epochs(edf, epoch_len=10*eeg_sr(edf))

# PROCESS
@benchmark e10 = eeg_reference_car(edf_epoched)
@benchmark eeg_filter!(e10, fprototype=:iirnotch, cutoff=50, bw=2)
@benchmark eeg_filter!(e10, fprototype=:butterworth, ftype=:lp, cutoff=45.0, order=8)
@benchmark eeg_filter!(e10, fprototype=:butterworth, ftype=:hp, cutoff=0.1, order=8)

# ANALYZE
tbp = eeg_total_power(e10)
ac = eeg_acov(e10, norm=false)
cc = eeg_xcov(e10, lag=10, demean=true)
mconv = eeg_tconv(e10, kernel=generate_morlet(256, 1, 32, complex=true))

# PLOT
```

## Results

### Workstation

Specs: AMD Threadripper 3960X, RAM: 128 GB, NVIDIA Quadro P2200

### Laptop

Specs: ThinkPad T14: AMD Ryzen 5 Pro, 12 cores, RAM 20 GB