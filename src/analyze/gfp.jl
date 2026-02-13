export gfp

"""
    gfp(obj)

Calculate Global Field Power (GFP).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names

# Returns

- `gfp::Vector{Float64}`

# Notes

Global field power is a measure of agreement of the signals picked up by all sensors across the entire scalp: if all sensors have the same value at a given time point, the GFP will be zero at that time point; if the signals differ, the GFP will be non-zero at that time point. GFP peaks may reflect “interesting” brain activity, warranting further investigation. Mathematically, the GFP is the population standard deviation across all sensors, calculated separately for every time point.
"""
function gfp(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::Vector{Float64}

    _check_datatype(obj, ["erp", "erf"])

    ch = exclude_bads ? get_channel(obj, ch=ch, exclude="bad") : get_channel(obj, ch=ch, exclude="")
    @assert length(ch) > 1 "More than 1 channel must be selected."
    s = @views obj.data[ch, :, 1]

    return std(s, dims=1)[:]

end