# ---------------------------------------------------------------------------
# low-pass and high-pass frequencies of the channels
# ---------------------------------------------------------------------------

"""
    _prefilter_freqs(s, fs)

Low-pass and high-pass frequency (Hz) from the description of the filters of a channel stored by a file format, e.g. `"HP:0.1Hz LP:70Hz N:50Hz"` (EDF), `"HP: DC; LP: 500 Hz"`, `"LP: 172.2 Hz; HP: 0.1 Hz"`, `"HP:0,18Hz LP:104,0Hz"` (decimal comma). A high-pass given in seconds (time constant τ) is converted to 1 / (2πτ) Hz. Missing, non-positive or non-numeric values (e.g. `DC`, `n/a`, `0`): no filter, i.e. the low-pass frequency is the Nyquist frequency and the high-pass frequency is 0 (as in MNE-Python); the low-pass frequency is at most the Nyquist frequency.

# Returns

- `Tuple{Float64, Float64}`: low-pass and high-pass frequency in Hz
"""
function _prefilter_freqs(s::AbstractString, fs::Real)::Tuple{Float64, Float64}
    nyq = fs / 2
    # decimal comma (e.g. "0,18Hz")
    t = replace(string(s), r"(\d),(\d)" => s"\1.\2")
    function val(re)
        m = match(re, t)
        isnothing(m) && return nothing
        v = tryparse(Float64, m[1])
        (isnothing(v) || !isfinite(v) || v <= 0) && return nothing
        return (v, lowercase(something(m[2], "")))
    end
    num = raw"([0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)\s*(hz|s)?"
    lp = val(Regex("(?:\\bLP|low[\\s_-]*pass)\\s*[:=]?\\s*" * num, "i"))
    hp = val(Regex("(?:\\bHP|high[\\s_-]*pass)\\s*[:=]?\\s*" * num, "i"))
    f_lp = isnothing(lp) ? nyq : min(lp[1], nyq)
    f_hp = if isnothing(hp)
        0.0
    else
        # time constant in seconds → cut-off frequency
        hp[2] == "s" ? 1 / (2π * hp[1]) : hp[1]
    end
    return f_lp, min(f_hp, nyq)
end

"""
    _prefilter_notch(s, fs)

Notch filter frequencies (Hz) from the description of the filters of a channel stored by a file format, e.g. `"HP:0.1Hz LP:70Hz N:50Hz"` (EDF) or `"Notch: 50 Hz"`; frequencies that are not positive or not below the Nyquist frequency are ignored.

# Returns

- `Vector{Float64}`: unique notch frequencies, sorted
"""
function _prefilter_notch(s::AbstractString, fs::Real)::Vector{Float64}
    t = replace(string(s), r"(\d),(\d)" => s"\1.\2")
    re = r"(?:\bN|notch)\s*[:=]?\s*([0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)"i
    f = Float64[]
    for m in eachmatch(re, t)
        v = tryparse(Float64, m[1])
        (isnothing(v) || !isfinite(v) || v <= 0 || v >= fs / 2) && continue
        push!(f, v)
    end
    return sort(unique(f))
end

"""
    _rec_filters!(rec)

Make sure that the recording header `rec` has the per-channel low-pass, high-pass and notch filter frequencies (`:lowpass`, `:highpass`: one value per channel; `:notch`: a vector of frequencies per channel): missing or wrongly sized fields are created from the `:prefiltering` descriptions of objects saved before v0.26.12 (which are then removed) or set to no filter (low-pass: Nyquist frequency, high-pass: 0, notch: none).

# Returns

- `Nothing`
"""
function _rec_filters!(rec::Dict)::Nothing
    n = length(rec[:label])
    fs = rec[:sampling_rate]
    pf = get(rec, :prefiltering, nothing)
    ok(k) = haskey(rec, k) && rec[k] isa AbstractVector{<:Real} && length(rec[k]) == n
    if !(ok(:lowpass) && ok(:highpass))
        f = if pf isa AbstractVector && length(pf) == n
            [_prefilter_freqs(string(x), fs) for x in pf]
        else
            fill((fs / 2, 0.0), n)
        end
        rec[:lowpass] = Float64[x[1] for x in f]
        rec[:highpass] = Float64[x[2] for x in f]
    end
    if !(haskey(rec, :notch) && rec[:notch] isa AbstractVector && length(rec[:notch]) == n)
        rec[:notch] = if pf isa AbstractVector && length(pf) == n
            Vector{Float64}[_prefilter_notch(string(x), fs) for x in pf]
        else
            [Float64[] for _ = 1:n]
        end
    end
    delete!(rec, :prefiltering)
    return nothing
end

# ---------------------------------------------------------------------------
# internal helper: build the fields shared by every recording type
# ---------------------------------------------------------------------------
function _common_recording_fields(;
    data_type::String,
    file_name::String,
    file_size_mb::Real,
    file_type::String,
    recording::String,
    recording_date::String,
    recording_time::String,
    recording_notes::String,
    channel_type::Vector{String},
    channel_order::Vector{Int64},
    clabels::Vector{String},
    units::Vector{String},
    sampling_rate::Int64,
    bad_channels::Vector{Bool},
    prefiltering::Vector{String} = String[],
)::Dict
    # low-pass and high-pass frequency of each channel (Hz), from the filter descriptions
    # of the file (`prefiltering`); no filter: Nyquist frequency and 0
    f = length(prefiltering) == length(clabels) ?
        [_prefilter_freqs(p, sampling_rate) for p in prefiltering] :
        fill((sampling_rate / 2, 0.0), length(clabels))
    notch = length(prefiltering) == length(clabels) ?
        Vector{Float64}[_prefilter_notch(p, sampling_rate) for p in prefiltering] :
        [Float64[] for _ in clabels]
    return Dict(
        :lowpass => Float64[x[1] for x in f],
        :highpass => Float64[x[2] for x in f],
        :notch => notch,
        :data_type => data_type,
        :file_name => file_name,
        :file_size_mb => file_size_mb,
        :file_type => file_type,
        :recording => recording,
        :recording_date => recording_date,
        :recording_time => recording_time,
        :recording_notes => recording_notes,
        :channel_type => channel_type,
        :channel_order => channel_order,
        :label => clabels,
        :unit => units,
        :sampling_rate => sampling_rate,
        :bad_channel => bad_channels,
    )
end

# ---------------------------------------------------------------------------
# subject
# ---------------------------------------------------------------------------

function _create_subject(;
    id::String,
    first_name::String,
    middle_name::String,
    last_name::String,
    handedness::String,
    head_circumference::Real,
    weight::Real,
    height::Real,
)::Dict
    return Dict(
        :id => id,
        :first_name => first_name,
        :middle_name => middle_name,
        :last_name => last_name,
        :handedness => handedness,
        :head_circumference => head_circumference,
        :weight => weight,
        :height => height,
    )
end

# ---------------------------------------------------------------------------
# EEG / sEEG / ECoG  (identical schema; data_type distinguishes them)
# ---------------------------------------------------------------------------

function _create_recording_eeg_like(;
    data_type::String,
    file_name::String,
    file_size_mb::Real,
    file_type::String,
    recording::String,
    recording_date::String,
    recording_time::String,
    recording_notes::String,
    channel_type::Vector{String},
    channel_order::Vector{Int64},
    reference::String,
    clabels::Vector{String},
    transducers::Vector{String},
    units::Vector{String},
    prefiltering::Vector{String},
    line_frequency::Real,
    sampling_rate::Int64,
    gain::Vector{Float64},
    bad_channels::Vector{Bool},
)::Dict
    d = _common_recording_fields(;
        data_type, file_name, file_size_mb, file_type, recording,
        recording_date, recording_time, recording_notes,
        channel_type, channel_order, clabels, units, sampling_rate, bad_channels,
        prefiltering,
    )
    merge!(
        d,
        Dict(
            :reference => reference,
            :transducers => transducers,
            :line_frequency => line_frequency,
            :gain => gain,
        ),
    )
    return d
end

# public thin wrappers - callers use the modality-specific name; all three delegate to the shared implementation above.
_create_recording_eeg(; kwargs...) = _create_recording_eeg_like(; kwargs...)
_create_recording_seeg(; kwargs...) = _create_recording_eeg_like(; kwargs...)
_create_recording_ecog(; kwargs...) = _create_recording_eeg_like(; kwargs...)

# ---------------------------------------------------------------------------
# MEG
# ---------------------------------------------------------------------------

function _create_recording_meg(;
    data_type::String,
    file_name::String,
    file_size_mb::Real,
    file_type::String,
    recording::String,
    recording_date::String,
    recording_time::String,
    recording_notes::String,
    channel_type::Vector{String},
    channel_order::Vector{Int64},
    reference::String,
    clabels::Vector{String},
    units::Vector{String},
    prefiltering::Vector{String},
    line_frequency::Real,
    sampling_rate::Int64,
    magnetometers::Vector{Int64},
    gradiometers::Vector{Int64},
    coil_type::Vector{String},
    bad_channels::Vector{Bool},
    ssp_labels::Vector{String},
    ssp_channels::Vector{Bool},
    ssp_data::Matrix{Float64},
    maxshield::Bool = false,
    coil_loc::Matrix{Float64} = Matrix{Float64}(undef, 0, 12),
    dev_head_t::Matrix{Float64} = Matrix{Float64}(undef, 0, 0),
    coil_id::Vector{Int64} = Int64[],
    head_origin::Vector{Float64} = Float64[],
    hpi::Dict{Symbol, Any} = Dict{Symbol, Any}(),
    first_time::Real = 0.0,
)::Dict
    # sensor geometry (FIFF/MNE `loc` convention), one row per channel:
    # r0 (1:3, m), ex (4:6), ey (7:9), ez (10:12; coil normal)
    # MEG channels: device frame; EEG channels: head frame (4:6 = reference electrode)
    (isempty(coil_loc) || size(coil_loc) == (length(clabels), 12)) || throw(
        ArgumentError("coil_loc must be an n_channels × 12 matrix."),
    )
    (isempty(dev_head_t) || size(dev_head_t) == (4, 4)) ||
        throw(ArgumentError("dev_head_t must be a 4 × 4 matrix."))
    # FIFF coil type code per channel (selects coil integration points)
    (isempty(coil_id) || length(coil_id) == length(clabels)) ||
        throw(ArgumentError("coil_id must have one entry per channel."))
    # head origin (sphere fitted to digitization, head frame, m)
    (isempty(head_origin) || length(head_origin) == 3) ||
        throw(ArgumentError("head_origin must have 3 elements."))
    # hpi: :freqs (Hz), :dig (n_coils × 3, head frame, m), :coil_no
    # first_time: time (s) of the first sample from the start of the acquisition (FIFF
    # first_samp / sampling rate); MaxFilter / MNE-Python head positions are relative to it
    # normalize time separator (MEG files sometimes use '.' instead of ':')
    recording_time = replace(recording_time, '.' => ':')
    d = _common_recording_fields(;
        data_type, file_name, file_size_mb, file_type, recording,
        recording_date, recording_time, recording_notes,
        channel_type, channel_order, clabels, units, sampling_rate, bad_channels,
        prefiltering,
    )
    merge!(
        d,
        Dict(
            :reference => reference,
            :line_frequency => line_frequency,
            :magnetometers => magnetometers,
            :gradiometers => gradiometers,
            :coil_type => coil_type,
            :ssp_labels => ssp_labels,
            :ssp_channels => ssp_channels,
            :ssp_data => ssp_data,
            :maxshield => maxshield,
            :coil_loc => coil_loc,
            :dev_head_t => dev_head_t,
            :coil_id => coil_id,
            :head_origin => head_origin,
            :hpi => hpi,
            :first_time => Float64(first_time),
            :sss => false,
            :tsss => false,
        ),
    )
    return d
end

# ---------------------------------------------------------------------------
# NIRS
# ---------------------------------------------------------------------------

function _create_recording_nirs(;
    data_type::String,
    file_name::String,
    file_size_mb::Real,
    file_type::String,
    recording::String,
    recording_date::String,
    recording_time::String,
    recording_notes::String,
    wavelengths::Vector{Float64},
    wavelength_index::Vector{Int64},
    optode_pairs::Matrix{Int64},
    channel_type::Vector{String},
    channel_order::Vector{Int64},
    clabels::Vector{String},
    units::Vector{String},
    src_labels::Vector{String},
    det_labels::Vector{String},
    opt_labels::Vector{String},
    sampling_rate::Int64,
    bad_channels::Vector{Bool},
    length_unit::String = "",
)::Dict
    d = _common_recording_fields(;
        data_type, file_name, file_size_mb, file_type, recording,
        recording_date, recording_time, recording_notes,
        channel_type, channel_order, clabels, units, sampling_rate, bad_channels,
    )
    merge!(
        d,
        Dict(
            :wavelengths => wavelengths,
            # unit of the optode positions (e.g. "mm", "cm"); "" if unknown
            :length_unit => length_unit,
            :wavelength_index => wavelength_index,
            :optode_pairs => optode_pairs,
            :src_labels => src_labels,
            :det_labels => det_labels,
            :optode_labels => opt_labels,
        ),
    )
    return d
end

# ---------------------------------------------------------------------------
# sensors / EDA / TPT  (identical schema; data_type distinguishes them)
# ---------------------------------------------------------------------------

function _create_recording_prefiltered(;
    data_type::String,
    file_name::String,
    file_size_mb::Real,
    file_type::String,
    recording::String,
    recording_date::String,
    recording_time::String,
    recording_notes::String,
    channel_type::Vector{String},
    channel_order::Vector{Int64},
    clabels::Vector{String},
    units::Vector{String},
    prefiltering::Vector{String},
    sampling_rate::Int64,
    bad_channels::Vector{Bool},
)::Dict
    d = _common_recording_fields(;
        data_type, file_name, file_size_mb, file_type, recording,
        recording_date, recording_time, recording_notes,
        channel_type, channel_order, clabels, units, sampling_rate, bad_channels,
        prefiltering,
    )
    return d
end

_create_recording_sensors(; kwargs...) = _create_recording_prefiltered(; kwargs...)
_create_recording_eda(; kwargs...) = _create_recording_prefiltered(; kwargs...)
_create_recording_tpt(; kwargs...) = _create_recording_prefiltered(; kwargs...)

# ---------------------------------------------------------------------------
# MEP
# ---------------------------------------------------------------------------

function _create_recording_mep(;
    data_type::String,
    file_name::String,
    file_size_mb::Real,
    file_type::String,
    recording::String,
    recording_date::String,
    recording_time::String,
    recording_notes::String,
    channel_type::Vector{String},
    channel_order::Vector{Int64},
    clabels::Vector{String},
    units::Vector{String},
    sampling_rate::Int64,
    stimulation_intensity::Vector{Int64},
    coil_type::Vector{String},
    stimulation_sample::Vector{Int64},
    markers_pos::Vector{Int64},
    markers_neg::Vector{Int64},
    bad_channels::Vector{Bool},
)::Dict
    d = _common_recording_fields(;
        data_type, file_name, file_size_mb, file_type, recording,
        recording_date, recording_time, recording_notes,
        channel_type, channel_order, clabels, units, sampling_rate, bad_channels,
    )
    merge!(
        d,
        Dict(
            :stimulation_intensity => stimulation_intensity,
            :coil_type => coil_type,
            :stimulation_sample => stimulation_sample,
            :markers_pos => markers_pos,
            :markers_neg => markers_neg,
        ),
    )
    return d
end

# ---------------------------------------------------------------------------
# header
# ---------------------------------------------------------------------------

function _create_experiment(; name::String, notes::String, design::String)::Dict
    return Dict(:name => name, :notes => notes, :design => design)
end

function _create_header(;
    subject::Dict,
    recording::Dict,
    experiment::Dict,
)::NeuroAnalyzer.HEADER
    return NeuroAnalyzer.HEADER(subject, recording, experiment)
end
