# NeuroAnalyzer test plugin - also a template for new plugins.
#
# Plugin layout: <plugins_path>/<plugin_name>/src/*.jl
# Every .jl file in src/ is included into the NeuroAnalyzer module by plugins_reload()
# (called at start-up and after plugins_install() / plugins_update()), therefore:
#
# - public and internal NeuroAnalyzer functions (e.g. _info(), _zeros(), get_channel()) can be called directly
# - names must not clash with NeuroAnalyzer functions (or other plugins): use a common prefix
# - functions are called as NeuroAnalyzer.<name>() (export only takes effect for `using NeuroAnalyzer`
#   executed after the plugin was loaded)
# - follow the NeuroAnalyzer conventions: a private array method `_name(s)` and a public method `name(obj; ch)`
#   for NEURO objects, NEURO.data is Array{Float64, 3} (channels × samples × epochs)

export test_plugin
export test_plugin_zcr

"""
    _test_plugin_zcr(s; fs)

Zero-crossing rate of a signal.

# Arguments

- `s::AbstractVector`: signal vector
- `fs::Int64`: sampling rate in Hz

# Returns

- `Float64`: number of zero crossings per second
"""
function _test_plugin_zcr(s::AbstractVector; fs::Int64)::Float64
    # validate
    fs >= 1 || throw(ArgumentError("fs must be ≥ 1."))
    length(s) >= 2 || throw(ArgumentError("Signal must contain at least 2 samples."))

    return _zeros(s) / (length(s) / fs)
end

"""
    _test_plugin_zcr(s; fs)

Zero-crossing rate of a signal array.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `fs::Int64`: sampling rate in Hz

# Returns

- `Matrix{Float64}`: number of zero crossings per second, shape (channels, epochs)
"""
function _test_plugin_zcr(s::AbstractArray; fs::Int64)::Matrix{Float64}
    _chk3d(s)

    zcr = zeros(size(s, 1), size(s, 3))
    @inbounds Threads.@threads :static for idx in CartesianIndices(zcr)
        zcr[idx] = _test_plugin_zcr(@view(s[idx[1], :, idx[2]]); fs = fs)
    end

    return zcr
end

"""
    test_plugin_zcr(obj; <keyword arguments>)

Zero-crossing rate (an example plugin function). For a sine of frequency `f` the rate is `2f`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)

# Returns

- `Matrix{Float64}`: number of zero crossings per second, shape (channels, epochs)
"""
function test_plugin_zcr(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::Matrix{Float64}
    ch = exclude_bads ? get_channel(obj; ch = ch, exclude = "bad") : get_channel(obj; ch = ch, exclude = "")
    isempty(ch) && throw(ArgumentError("No channels selected."))

    return _test_plugin_zcr(@view(obj.data[ch, :, :]); fs = sr(obj))
end

"""
    test_plugin()

Check that plugins are loaded and work: `test_plugin_zcr()` is applied to a synthetic 10 Hz sine (expected zero-crossing rate: 20 per second).

# Arguments

Nothing

# Returns

- `Bool`: `true` if the self-test passed
"""
function test_plugin()::Bool
    _info("NeuroAnalyzer test plugin")

    # synthetic object: one EEG channel, 10 Hz sine (phase-shifted so that no sample is exactly 0 and
    # all 20 crossings of each 1 s epoch fall between samples), 256 Hz, 2 epochs × 1 s
    fs = 256
    t = (0:(fs - 1)) ./ fs
    s = sin.(2π * 10 .* t .- 0.1)
    obj = create_object(; data_type = "eeg")
    obj.header.recording[:sampling_rate] = fs
    add_channel!(obj; data = reshape(repeat(s, 2), 1, fs, 2), label = "sine", type = "eeg", unit = "μV")
    create_time!(obj; fs = fs)

    zcr = test_plugin_zcr(obj; ch = "sine")
    passed = size(zcr) == (1, 2) && all(zcr .== 20)
    passed ? _info("Test passed") : _warn("Test failed: zero-crossing rate $(vec(zcr)), expected 20")

    return passed
end
