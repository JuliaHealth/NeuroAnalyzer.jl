export tes_model

"""
    tes_model(; <keyword arguments>)

Create model of TES stimulation.

The electric potential at every electrode of the standard 10-10 (47 channels) cap is the superposition of the anode and cathode contributions, `k × I / (r + 1)` (`k = 1`, `r`: distance from the electrode in the plane of the cap). It is plotted as a topographical map.

# Arguments

- `anode::String`: anode location (label of the 10-10 system, e.g. `"F3"`)
- `cathode::String`: cathode location (e.g. `"Fp2"`)
- `anode_curr::Real=2.0`: anode current [mA]; must be > 0
- `cathode_curr::Real=-2.0`: cathode current [mA]; must be < 0

# Returns

- `GLMakie.Figure`

# Notes

This is a very initial version, simplified model - just superficial spread of the electric field

# To do

Model spread of electric field at the cortical surface - reduce charge for skull resistance
"""
function tes_model(;
    anode::String,
    cathode::String,
    anode_curr::Real = 2.0,
    cathode_curr::Real = -2.0,
)::GLMakie.Figure
    _wip()

    # validate
    anode_curr > 0 || throw(ArgumentError("anode_curr must be > 0."))
    cathode_curr < 0 || throw(ArgumentError("cathode_curr must be < 0."))

    locs = import_locs(joinpath(NeuroAnalyzer.PATH, "locs", "standard-10-10-cap47.ced"))
    clabels = String.(locs[!, :label])
    ch_n = length(clabels)

    anode_ch = findfirst(==(lowercase(anode)), lowercase.(clabels))
    cathode_ch = findfirst(==(lowercase(cathode)), lowercase.(clabels))
    isnothing(anode_ch) && throw(ArgumentError("Anode location $anode not found."))
    isnothing(cathode_ch) && throw(ArgumentError("Cathode location $cathode not found."))
    anode_ch != cathode_ch || throw(ArgumentError("anode and cathode must be different."))

    loc_x = locs[!, :loc_x]
    loc_y = locs[!, :loc_y]

    # kₑ = 8.988 * 10^9 # for air / vacuum
    kₑ = 1

    # distances from the anode and the cathode (+1: finite at the electrode itself)
    r_a = [hypot(loc_x[idx] - loc_x[anode_ch], loc_y[idx] - loc_y[anode_ch]) for idx = 1:ch_n] .+ 1
    r_c = [hypot(loc_x[idx] - loc_x[cathode_ch], loc_y[idx] - loc_y[cathode_ch]) for idx = 1:ch_n] .+ 1

    E = (kₑ * anode_curr) ./ r_a .+ (kₑ * cathode_curr) ./ r_c
    E[anode_ch] = anode_curr
    E[cathode_ch] = cathode_curr

    # EEG-type object holding the potential as a constant 2-sample signal
    obj = create_object(; data_type = "eeg")
    create_data!(obj; data = reshape(repeat(E, 1, 2), ch_n, 2, 1), fs = 1, type = "eeg")
    obj.header.recording[:label] = clabels
    obj.locs = locs

    return plot_topo(
        obj;
        ch = "eeg",
        tpos = 0,
        title = "TES model: $anode (+) / $cathode (−)",
        nmethod = :none,
    )
end
