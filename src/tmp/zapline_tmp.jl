
"""
    signal_zapline(signal::Array{Float64, 3}; <keyword arguments>)

Automatic removal of power-line artifacts.

# Arguments

- `signal::Array{Float64, 3}`
- `fs::Int64`: sampling rate
- `fline::Int64=50`: line frequency [Hz]
- `nremove::Int64=1`: number of components to remove
- `nfft::Int64=1024`: size of FFT
- `keep::Int64=0`: number of components to keep in DSS, default is all
- `iter::Int64=1`: number of iterations for smoothing filter

# Returns

- `s_clean::AbstractArray`

# Source

de Cheveigné A. ZapLine: A simple and effective method to remove power line artifacts. NeuroImage. 2020 Feb;207:116356.
"""
function signal_zapline(signal::Array{Float64, 3}; fs::Int64, fline::Int64=50, nremove::Int64=1, nfft::Int64=1024, keep::Int64=0, iter::Int64=1)

    channel_n, signal_len, epoch_n = size(signal)
    signal_len < nfft && (nfft = floor(Int64, length(signal), 2))

    s_clean = deepcopy(signal)
    # smooth by convolution with square window
    square_window = rect(floor(Int64, 1 / (fline / fs))) ./ (1 / (fline / fs))
    if iter > 1
        for idx in 1:(iter - 1)
            square_window = conv(square_window, square_window)
        end
    end
    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx2 in 1:channel_n
            s = @view s_clean[idx2, :, epoch]
            s_clean[idx2, :, epoch] = filtfilt(square_window, s)
        end
    end

    # reduce dimensionality to avoid overfitting
    keep == 0 && (keep = channel_n)
    s_pca = signal - s_clean
    @inbounds @simd for epoch in 1:epoch_n
        M = MultivariateStats.fit(PCA, s_clean[:, :, epoch], maxoutdim=keep)
        Threads.@threads for idx in 1:keep
            MultivariateStats.predict(M, signal_pca)[idx, :]
        end
    end
    # xxxx=nt_pca(x-xx,[],p.nkeep); % reduce dimensionality to avoid overfitting

    for idx in 1:n
        pc_var[idx, epoch] = v[idx]
        # pc[idx, :, epoch] = (eig_vec[:, idx] .* s)[idx, :]
        pc[idx, :, epoch] = MultivariateStats.predict(M, s)[idx, :]
    end

    return signal_clean::AbstractArray
end