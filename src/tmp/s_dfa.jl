"""
    s_dfa(signal; fs)

Perform detrended fluctuation analysis of the `signal`.

# Arguments

- `signal::Array{<:Real, 3}`
- `fs::Int64`: sampling rate

# Returns

- `signal_cs::Vector{Float64}`
"""
function s_dfa(signal::Array{<:Real, 3}; fs::Int64)
    
    # INCOMPLETE, DOES NOT WORK, DO NOT USE
    # also, eeg_dfa() is missing

    signal_dm = s_demean(signal)
    signal_cs = s_cums(signal_dm)

    scale_order = collect(1:20)
    scale_duration = log.(scale_order)
    scatter(scale_order, scale_duration)

    @inbounds @simd for idx in length(scale_duration):-1:1
        if scale_duration[idx] > length(signal) ÷ 2
            popat!(scale_duration, idx)
            popat!(scale_order, idx)
        end
    end
    @inbounds @simd for idx in length(scale_duration):-1:2
        if scale_duration[idx] == scale_duration[idx - 1]
            popat!(scale_duration, idx)
            popat!(scale_order, idx)
        end
    end

    epoch_rms = zeros(length(scale_duration))
    @inbounds @simd for order_idx in 1:length(scale_duration)
        epochs_n = Int(length(signal_cs) ÷ scale_duration[order_idx])
        epochs = zeros(epochs_n, scale_duration[order_idx])
        for epoch_idx in 1:epochs_n
            epochs[epoch_idx, :] = signal_cs[(((epoch_idx - 1) * scale_duration[order_idx]) + 1):(epoch_idx * scale_duration[order_idx])]
        end
        epochs_dt = s_detrend(epochs, type=:ls)
        epochs_rms = zeros(epochs_n)
        for idx in 1:size(epochs, 1)
            epochs_rms[idx] = s_rms(epochs_dt[idx, :])
        end
        epoch_rms[order_idx] = mean(epochs_rms)
    end
    scale_duration = log.(scale_duration)
    epoch_rms = log.(epoch_rms)
    scatter(scale_order, epoch_rms)
    atilde = pinv(log.(scale_order)) * log.(epoch_rms)
    Plots.plot!(atilde .* log.(scale_order))

    return signal_cs
end