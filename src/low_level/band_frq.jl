export band_frq

"""
    band(obj, band)

Return frequency limits for a `band`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `band::Symbol`: band range name:
    - `:list`
    - `:total`
    - `:delta`
    - `:theta`
    - `:alpha`
    - `:beta`
    - `:beta_high`
    - `:gamma`
    - `:gamma_1`
    - `:gamma_2`
    - `:gamma_lower`
    - `:gamma_higher`.
# Returns

- `band_frequency::Tuple{Real, Real}`
"""
function band_frq(obj::NeuroAnalyzer.NEURO; band::Symbol)

    bands = [:list, :total, :delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher]
    _check_var(band, [:list, :total, :delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher], "band")
    if band === :list
        print("Available band names: ")
        [print(":$(bands[x]), ") for x in 2:(length(bands) - 1)]
        println(":$(bands[end])")
        return nothing
    end

    band === :total && (band_frq = (0.1, round(sr(obj) / 2, digits=1)))
    band === :delta && (band_frq = (0.1, 4.0))
    band === :theta && (band_frq = (4.0, 8.0))
    band === :alpha && (band_frq = (8.0, 13.0))
    band === :beta && (band_frq = (14.0, 30.0))
    band === :beta_high && (band_frq = (25.0, 30.0))
    band === :gamma && (band_frq = (30.0, 150.0))
    band === :gamma_1 && (band_frq = (31.0, 40.0))
    band === :gamma_2 && (band_frq = (41.0, 50.0))
    band === :gamma_lower && (band_frq = (30.0, 80.0))
    band === :gamma_higher && (band_frq = (80.0, 150.0))
    
    if band_frq[1] > sr(obj) / 2
        _info("Nyquist frequency based on sampling rate ($(sr(obj) / 2)) is lower than $band range: $band_frq, band frequency truncated to: ($(sr(obj) / 2 - 0.2), $(sr(obj) / 2 - 0.1))")
        band_frq = (sr(obj) / 2 - 0.2, sr(obj) / 2 - 0.1)
    end
    if band_frq[2] > sr(obj) / 2
        _info("Nyquist frequency based on sampling rate ($(sr(obj) / 2)) is lower than $band range: $band_frq, band frequency truncated to: ($(band_frq[1]), $(sr(obj) / 2 - 0.1))")
        band_frq = (band_frq[1], sr(obj) / 2 - 0.1)
    end

    return band_frq
end

"""
    band_frq(fs, band)

Return frequency limits of a `band`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `band::Symbol`: band range name:
    - `:list`
    - `:total`
    - `:delta`
    - `:theta`
    - `:alpha`
    - `:beta`
    - `:beta_high`
    - `:gamma`
    - `:gamma_1`
    - `:gamma_2`
    - `:gamma_lower`
    - `:gamma_higher`.

# Returns

- `band_frq::Tuple{Real, Real}`
"""
function band_frq(fs::Int64; band::Symbol)

    bands = [:list, :total, :delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher]
    _check_var(band, [:list, :total, :delta, :theta, :alpha, :beta, :beta_high, :gamma, :gamma_1, :gamma_2, :gamma_lower, :gamma_higher], "band")
    if band === :list
        print("Available band names: ")
        [print(":$(bands[x]), ") for x in 2:(length(bands) - 1)]
        println(":$(bands[end])")
        return nothing
    end

    band === :total && (band_frq = (0.1, round(fs / 2, digits=1)))
    band === :delta && (band_frq = (0.1, 4.0))
    band === :theta && (band_frq = (4.0, 8.0))
    band === :alpha && (band_frq = (8.0, 13.0))
    band === :beta && (band_frq = (14.0, 30.0))
    band === :beta_high && (band_frq = (25.0, 30.0))
    band === :gamma && (band_frq = (30.0, 150.0))
    band === :gamma_1 && (band_frq = (31.0, 40.0))
    band === :gamma_2 && (band_frq = (41.0, 50.0))
    band === :gamma_lower && (band_frq = (30.0, 80.0))
    band === :gamma_higher && (band_frq = (80.0, 150.0))
    
    if band_frq[1] > fs / 2
        _info("Nyquist frequency based on EEG sampling rate ($(fs / 2)) is lower than $band range: $band_frq, band frequency truncated to: ($(fs / 2 - 0.2), $(fs / 2 - 0.1))")
        band_frq = (fs / 2 - 0.2, fs / 2 - 0.1)
    end
    if band_frq[2] > fs / 2
        _info("Nyquist frequency based on EEG sampling rate ($(fs / 2)) is lower than $band range: $band_frq, band frequency truncated to: ($(band_frq[1]), $(fs / 2 - 0.1))")
        band_frq = (band_frq[1], fs / 2 - 0.1)
    end

    return band_frq
end

