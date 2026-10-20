# standard 10-20 / 10-10 electrode names (case-insensitive): prefix + number or "z"
const _LABEL_1010 = r"^(fp|af|ft|fc|tp|cp|po|f|c|t|p|o|i|n|a|m)(\d{1,2}|z)$"i

# numbered channels of MEG systems (≥ 3 digits: MEG0113, EEG 001, EOG061, STI 014); "ECG1" is not one
const _LABEL_NUMBERED = r"^(MEG|EEG|EOG|ECG|EKG|EMG|STI|MISC|IAS|SYS|CHPI)\s*(\d{3,})$"i

"""
    _clean_label(label)

Normalize a single channel label (see `_clean_labels()`).
"""
function _clean_label(label::AbstractString)::String
    raw = strip(String(label))
    l = replace(raw, r"\s+" => " ")

    # numbered channels: one space between the prefix and the number, leading zeros kept
    # ("MEG0113" → "MEG 0113", "EOG061" → "EOG 061", "EEG 001" unchanged)
    m = match(_LABEL_NUMBERED, l)
    m !== nothing && return m[1] * " " * m[2]

    # format prefixes before a channel name ("EEG Fp1", "EDF C3", "EEG-Fp1"), not before a number
    l = replace(l, r"^(EEG|EDF|BDF)[\s_-]+(?=[^\d\s])"i => "")

    # duplicated type prefixes
    l = replace(l, r"^EOG\s+EOG"i => "EOG")
    l = replace(l, r"^ECG\s+EKG"i => "ECG")

    # trailing dots used as padding (PhysioNet: "Fc5.", "T9..")
    l = replace(l, r"\.+$" => "")

    # reference suffixes (TUH: "FP1-REF", "C3-LE")
    l = replace(l, r"-(REF|LE|AR|AVG)$"i => "")

    # standard capitalization of 10-20 / 10-10 names: FP1 → Fp1, CZ → Cz, afz → AFz, Fc5 → FC5
    m = match(_LABEL_1010, l)
    if m !== nothing
        prefix = lowercase(m[1]) == "fp" ? "Fp" : uppercase(m[1])
        l = prefix * lowercase(m[2])
    end

    return isempty(l) ? raw : l
end

"""
    _clean_labels(clabels)

Normalize channel labels read from a file header.

- numbered channels of MEG systems (MEG, EEG, EOG, ECG, EMG, STI, MISC, … with a ≥ 3-digit number): exactly one space between the prefix and the number; leading zeros are kept (`"MEG0113"` → `"MEG 0113"`, as in MNE-Python and the MEGIN naming)
- format prefixes before a channel name are removed (`"EEG Fp1"`, `"EDF C3"` → `"Fp1"`, `"C3"`)
- duplicated type prefixes are merged (`"EOG EOG"` → `"EOG"`, `"ECG EKG"` → `"ECG"`)
- trailing padding dots are removed (`"Fc5."` → `"FC5"`)
- reference suffixes are removed (`"FP1-REF"`, `"C3-LE"` → `"Fp1"`, `"C3"`)
- standard 10-20 / 10-10 names get the standard capitalization (`"FP1"` → `"Fp1"`, `"CZ"` → `"Cz"`)

Labels that would become identical after cleaning keep their original form.

# Arguments

- `clabels::AbstractVector{<:AbstractString}`: raw channel labels

# Returns

- `Vector{String}`: cleaned labels (same length, same order)
"""
function _clean_labels(clabels::AbstractVector{<:AbstractString})::Vector{String}
    raw = [String(strip(l)) for l in clabels]
    l = _clean_label.(clabels)
    # cleaning must not create duplicates: revert the colliding labels
    for lbl in unique(l)
        idx = findall(isequal(lbl), l)
        length(idx) > 1 && length(unique(raw[idx])) > 1 && (l[idx] = raw[idx])
    end
    return l
end

# label compared without separators, case and leading zeros of the number ("MEG 0113" ≡ "MEG113" ≡ "meg0113")
_label_key(s::AbstractString)::String =
    replace(uppercase(replace(s, r"[^A-Za-z0-9]" => "")), r"^([A-Z]+)0+(?=\d)" => s"\1")

"""
    _match_labels(labels, clabels)

Rename `labels` (e.g. of a bundled layout) to the matching recording labels `clabels`; labels without a match are kept.
"""
function _match_labels(labels::AbstractVector{<:AbstractString}, clabels::AbstractVector{<:AbstractString})::Vector{String}
    key = Dict(_label_key(c) => String(c) for c in clabels)
    return [get(key, _label_key(l), String(l)) for l in labels]
end
