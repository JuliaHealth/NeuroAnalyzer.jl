"""
    _len(obj::NeuroAnalyzer.NEURO, len::Int64, def_l::Int64)::Int64

Determine appropriate segment length for analysis, returning either the requested length or a default length based on the object's epoch duration.  

# Arguments
- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `len::Int64`: requested length in samples (0 means use default length)
- `def_l::Int64`: default length in seconds to use when `len=0`

# Returns
- `Int64`: Appropriate segment length in samples:
    - if `len > 0`: Returns `len` (requested length)
    - if `len == 0`:
        - if epoch duration > `def_l` seconds: returns `def_l * sr(obj)`
        - otherwise: returns the full epoch length
"""
function _len(obj::NeuroAnalyzer.NEURO, len::Int64, def_l::Int64)::Int64
    len >= 0 || throw(ArgumentError("Length must be non-negative."))
    def_l > 0 || throw(ArgumentError("Default length must be positive."))
    # return requested length if specified (len > 0)
    len > 0 && return len
    # calculate default length in samples
    default_samples = def_l * sr(obj)
    # return the appropriate length based on epoch duration
    if epoch_len(obj) > default_samples
        return default_samples
    else
        return epoch_len(obj)
    end
end
