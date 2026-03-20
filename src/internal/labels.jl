# ------------------------------------------------------------------ #
# internal helper: apply a sequence of replacements in a single pass,
# avoiding repeated intermediate allocations.
# ------------------------------------------------------------------ #

"""
    _clean_labels(clabels)

Normalize raw channel labels read from an EDF/BDF file.

Removes format prefixes (EEG, EDF, BDF), normalizes duplicate-label patterns (EOG EOG → EOG, ECG EKG → ECG), strips stray dots, and normalizes MEG labels to the form "MEG <number>" with no leading zeros.

# Arguments

- `clabels::Vector{String}`: raw channel labels from file header

# Returns

- `Vector{String}`: cleaned labels (same length, same order)
"""
function _clean_labels(clabels::Vector{String})::Vector{String}
    # string is immutable in Julia — shallow copy is sufficient
    l = copy(clabels)

    # remove EEG prefix (case-insensitive, with any trailing whitespace)
    # four separate replace calls collapsed into one regex pass
    l = replace.(l, r"(?i)EEG\s*" => "")

    # normalize duplicate-label patterns (case-insensitive)
    l = replace.(l, r"(?i)EOG\s+EOG" => "EOG")
    l = replace.(l, r"(?i)ECG\s+EKG" => "ECG")

    # remove EDF/BDF format prefixes with any trailing whitespace
    l = replace.(l, r"(?i)EDF\s*" => "")
    l = replace.(l, r"(?i)BDF\s*" => "")

    # strip dots - longest pattern first to avoid leaving isolated dots
    l = replace.(l, ".." => "", "." => "")

    # normalise MEG labels: ensure exactly one space after "MEG" and strip ALL leading zeros
    l = replace.(l, r"MEG\s*0*" => "MEG ")

    # collapse any double spaces introduced by prefix removal
    l = replace.(l, "  " => " ")

    return strip.(l)

end

"""
    _clean_meg_labels(clabels)

Normalize channel labels from an MEG file.

Ensures a single space between each modality prefix (MEG, EEG, EOG, EMG) and its channel number, and strips ALL leading zeros from the number.

# Arguments

- `clabels::Vector{String}`: raw channel labels from MEG file header

# Returns

- `Vector{String}`: cleaned labels (same length, same order)
"""
function _clean_meg_labels(clabels::Vector{String})::Vector{String}
    # string is immutable in Julia — shallow copy is sufficient
    l = copy(clabels)

    # ensure one space after each modality prefix and strip ALL leading zeros
    for prefix in ("MEG", "EEG", "EOG", "EMG")
        l = replace.(l, Regex("$prefix\\s*0*") => "$prefix ")
    end

    # collapse any double spaces introduced above
    l = replace.(l, "  " => " ")

    return strip.(l)

end

"""
    _clean_eeg_labels(clabels)

Normalize channel labels from an EEG file.

Removes the "EEG " prefix, normalises duplicate-label patterns, and strips ALL leading zeros from EOG/EMG channel numbers.

# Arguments

- `clabels::Vector{String}`: raw channel labels from EEG file header

# Returns

- `Vector{String}`: cleaned labels (same length, same order)
"""
function _clean_eeg_labels(clabels::Vector{String})::Vector{String}
    # string is immutable in Julia — shallow copy is sufficient
    l = copy(clabels)

    # remove "EEG " prefix (with trailing space) — a plain label like "EEG Fp1" becomes "Fp1"
    l = replace.(l, "EEG " => "")

    # normalise duplicate-label patterns
    l = replace.(l, "EOG EOG" => "EOG")
    l = replace.(l, "ECG EKG" => "ECG")

    # strip ALL leading zeros from EOG/EMG channel numbers
    for prefix in ("EEG", "EMG", "EOG")
        l = replace.(l, Regex("$prefix\\s*0+") => "$prefix ")
    end

    # collapse any double spaces
    l = replace.(l, "  " => " ")

    return strip.(l)

end
