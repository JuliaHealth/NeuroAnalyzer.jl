export tdcs_dose
export tacs_dose
export tpcs_dose
export tes_protocol

"""
    tdcs_dose(; <keyword arguments>)

Calculate charge, current density, and charge density for tDCS stimulation.

# Arguments

- `current::Real`: stimulation current in mA; must be > 0
- `pad_area::Real`: electrode pad area in cm²; must be > 0
- `duration::Int64`: stimulation duration in seconds; must be > 0

# Returns

Named tuple:

- `charge::Float64`: total delivered charge in C
- `current_density::Float64`: current density in A/m²
- `charge_density::Float64`: delivered charge density in kC/m²

# Throws

- `ArgumentError`: if `current`, `pad_area`, or `duration` are ≤ 0

# References

Chhatbar PY, George MS, Kautz SA, Feng W. Quantitative reassessment of safety limits of tDCS for two animal studies. Brain Stimulation. 2017;10(5):1011–2.

# See also

[`tacs_dose`](@ref), [`tpcs_dose`](@ref), [`tes_protocol`](@ref)
"""
function tdcs_dose(;
    current::Real,
    pad_area::Real,
    duration::Int64
)::@NamedTuple{
    charge::Float64,
    current_density::Float64,
    charge_density::Float64
}

    !(current  > 0) && throw(ArgumentError("current must be > 0."))
    !(pad_area > 0) && throw(ArgumentError("pad_area must be > 0."))
    !(duration > 0) && throw(ArgumentError("duration must be > 0."))

    # Unit conversions:
    #   current:   mA → A    (÷ 1_000)
    #   pad_area:  cm² → m²  (÷ 10_000)   was: ÷ 1_000 (wrong)
    #   charge:    C  → kC   (÷ 1_000)
    i_A = current  / 1_000
    a_m2 = pad_area / 10_000

    # [A × s = C]
    charge = i_A * duration
    # [A/m²]
    current_density = i_A / a_m2
    # [kC/m²]
    charge_density  = (charge / 1_000) / a_m2

    return (charge=charge, current_density=current_density, charge_density=charge_density)
end

"""
    tacs_dose(; <keyword arguments>)

Calculate charge, current density, and charge density for tACS stimulation.

The effective current is computed by integrating one cycle of the sinusoidal waveform (including the DC offset) over a 1 ms grid using Simpson's rule, then scaling by the number of cycles.

# Arguments

- `current::Real`: peak-to-peak stimulation current in mA; must be > 0
- `pad_area::Real`: electrode pad area in cm²; must be > 0
- `duration::Real`: stimulation duration in seconds; must be > 0
- `offset::Real`: DC current offset in μA
- `frequency::Real`: sinusoidal frequency in Hz; must be > 0
- `phase::Real`: phase shift in degrees

# Returns

Named tuple:

- `charge::Float64`: total delivered charge in C
- `current_density::Float64`: current density in A/m²
- `charge_density::Float64`: delivered charge density in kC/m²

# Throws

- `ArgumentError`: if `current`, `pad_area`, `duration`, or `frequency` are ≤ 0

# See also

[`tdcs_dose`](@ref), [`tpcs_dose`](@ref), [`tes_protocol`](@ref)
"""
function tacs_dose(;
    current::Real,
    pad_area::Real,
    duration::Int64,
    offset::Real,
    frequency::Real,
    phase::Real
)::@NamedTuple{
    charge::Float64,
    current_density::Float64,
    charge_density::Float64
}

    !(current > 0) && throw(ArgumentError("current must be > 0."))
    !(pad_area > 0) && throw(ArgumentError("pad_area must be > 0."))
    !(duration > 0) && throw(ArgumentError("duration must be > 0."))
    !(frequency > 0) && throw(ArgumentError("frequency must be > 0."))

    # integrate the rectified sinusoid over one cycle to get the effective current
    t = collect(0:0.001:1)
    i_cycle = abs.(generate_sine(frequency, t, current, phase) .+ offset)
    # effective mA over one cycle
    eff_current = simpson(i_cycle, t)

    cycles = frequency * duration
    _info("Number of cycles: $cycles")
    _info("Effective current: $eff_current mA")

    # unit conversions: mA → A (÷ 1_000), cm² → m² (÷ 10_000)
    i_A  = eff_current / 1_000
    a_m2 = pad_area / 10_000

    charge = i_A * duration
    current_density = i_A / a_m2
    charge_density = (charge / 1_000) / a_m2

    return (charge=charge, current_density=current_density, charge_density=charge_density)
end

"""
    tpcs_dose(; <keyword arguments>)

Calculate charge, current density, and charge density for tPCS stimulation.

The effective current is computed as `current × (pw / isi)`, i.e. scaled by the duty cycle.

# Arguments

# Arguments
- `current::Real`: peak-to-peak stimulation current in mA; must be > 0
- `pad_area::Real`: electrode pad area in cm²; must be > 0
- `duration::Real`: stimulation duration in seconds; must be > 0
- `pw::Real`: pulse width in ms; must be > 0 and < `isi`
- `isi::Real`: inter-stimulus interval in ms; must be > `pw`

# Returns

Named tuple:

- `charge::Float64`: total delivered charge in C
- `current_density::Float64`: current density in A/m²
- `charge_density::Float64`: delivered charge density in kC/m²

# Throws

- `ArgumentError`: if `current`, `pad_area`, `duration`, or `pw` are ≤ 0, or if `isi ≤ pw`

# See also

[`tdcs_dose`](@ref), [`tacs_dose`](@ref), [`tes_protocol`](@ref)
"""
function tpcs_dose(;
    current::Real,
    pad_area::Real,
    duration::Real,
    pw::Real,
    isi::Real
)::@NamedTuple{
    charge::Float64,
    current_density::Float64,
    charge_density::Float64}

    !(current > 0) && throw(ArgumentError("current must be > 0."))
    !(pad_area > 0) && throw(ArgumentError("pad_area must be > 0."))
    !(duration > 0) && throw(ArgumentError("duration must be > 0."))
    !(pw > 0) && throw(ArgumentError("pw must be > 0."))
    !(isi > pw) && throw(ArgumentError("isi must be > pw."))

    # convert pulse timings from ms → s
    pw_s  = pw  / 1_000
    isi_s = isi / 1_000

    cycles = duration / isi_s
    # duty-cycle–weighted effective current [mA]
    eff_current = current * (pw_s / isi_s)
    _info("Number of cycles: $cycles")
    _info("Effective current: $eff_current mA")

    # unit conversions: mA → A (÷ 1_000), cm² → m² (÷ 10_000)
    i_A = eff_current / 1_000
    a_m2 = pad_area / 10_000

    charge = i_A * duration
    current_density = i_A / a_m2
    charge_density = (charge / 1_000) / a_m2

    return (charge=charge, current_density=current_density, charge_density=charge_density)
end

"""
    tes_protocol(; <keyword arguments>)

Create a TES (tDCS/tACS/tRNS/tPCS) stimulation protocol dictionary.

# Arguments

- `type::Symbol`: stimulation type; one of `:tDCS`, `:tACS`, `:tRNS`, `:tPCS`
- `hd::Bool`: if `true`, use high-density electrodes
- `current::Real`: stimulation current in mA; must be > 0
- `frequency::Real=0`: stimulation frequency in Hz; must be > 0 for `:tACS` and `:tRNS`
- `anode_size::Tuple{Int64, Int64}`: anode dimensions `(width, height)` in mm; both values must be > 0
- `cathode_size::Tuple{Int64, Int64}`: cathode dimensions `(width, height)` in mm; both values must be > 0
- `anode_loc::Symbol`: anode location (10-20 Positioning System label)
- `cathode_loc::Symbol`: cathode location (10-20 Positioning System label)
- `duration::Real`: stimulation duration in seconds; must be > 0
- `ramp_in::Real`: ramp-in duration in seconds; must be ≥ 0
- `ramp_out::Real`: ramp-out duration in seconds; must be ≥ 0
- `sham::Bool`: if `true`, the protocol includes sham stimulation

# Returns

- `Dict`: protocol dictionary with all stimulation parameters

# Throws

- `ArgumentError`: if any argument fails its validation check

# See also

[`tdcs_dose`](@ref), [`tacs_dose`](@ref), [`tpcs_dose`](@ref)
"""
function tes_protocol(;
    type::Symbol,
    hd::Bool,
    current::Real,
    frequency::Real = 0,
    anode_size::Tuple{Int64, Int64},
    cathode_size::Tuple{Int64, Int64},
    anode_loc::Symbol,
    cathode_loc::Symbol,
    duration::Real,
    ramp_in::Real,
    ramp_out::Real,
    sham::Bool,
)::Dict

    _check_var(type, [:tDCS, :tACS, :tRNS, :tPCS], "type")
    !(current > 0) && throw(ArgumentError("current must be > 0 mA."))
    if type === :tACS || type === :tRNS
        !(frequency > 0) && throw(ArgumentError("frequency must be > 0 Hz."))
    end
    !(anode_size[1] > 0) && throw(ArgumentError("anode_size width must be > 0 mm."))
    !(anode_size[2] > 0) && throw(ArgumentError("anode_size height must be > 0 mm."))
    !(cathode_size[1] > 0) && throw(ArgumentError("cathode_size width must be > 0 mm."))
    !(cathode_size[2] > 0) && throw(ArgumentError("cathode_size height must be > 0 mm."))
    !(duration > 0) && throw(ArgumentError("duration must be > 0 s."))
    !(ramp_in  >= 0) && throw(ArgumentError("ramp_in must be ≥ 0 s."))
    !(ramp_out >= 0) && throw(ArgumentError("ramp_out must be ≥ 0 s."))

    protocol = Dict(
        :type => type,
        :hd => hd,
        :current => current,
        :frequency => frequency,
        :anode_size => anode_size
        :cathode_size => cathode_size,
        :anode_loc => anode_loc,
        :cathode_loc => cathode_loc,
        :duration => duration,
        :ramp_in => ramp_in,
        :ramp_out => ramp_out,
        :sham => sham,
    )

    return protocol

end
