export tdcs_dose
export tacs_dose
export tpcs_dose
export tes_protocol

"""
    tdcs_dose(; <keyword arguments>)

Calculate `charge`, `current_density` and `charge_ density` for tDCS stimulation.

# Arguments

- `current::Real`: stimulation current [mA]
- `pad_area::Real`: electrode pad area [cm²]
- `duration::Int64`: stimulation duration [s]

# Returns

Named tuple containing:
- `charge::Float64`: charge [C]
- `current_density::Float64`: current density [A/m²]
- `charge_density::Float64`: delivered charge density [kC/m²]

# Source

Chhatbar PY, George MS, Kautz SA, Feng W. Quantitative reassessment of safety limits of tDCS for two animal studies. Brain Stimulation. 2017;10(5):1011–2.
"""
function tdcs_dose(; current::Real, pad_area::Real, duration::Int64)

    charge = (current / 1_000) * duration
    current_density = (current / 1_000) / (pad_area / 1_000)
    charge_density = (charge / 1_000) / (pad_area / 1_000)

    return (charge=charge, current_density=current_density, charge_density=charge_density)

end

"""
    tacs_dose(; <keyword arguments>)

Calculate `charge`, `current_density` and `charge_ density` for tACS stimulation.

# Arguments

- `current::Real`: stimulation current (peak to peak) [mA]
- `pad_area::Real`: electrode pad area [cm²]
- `duration::Real`: stimulation duration [s]
- `offset::Real`: current offset [μA]
- `frequency::Real`: sinus frequency [Hz]
- `phase::Real`: phase shift [degree]

# Returns

Named tuple containing:
- `charge::Float64`: charge [C]
- `current_density::Float64`: current density [A/m²]
- `charge_density::Float64`: delivered charge density [kC/m²]
"""
function tacs_dose(; current::Real, pad_area::Real, duration::Int64, offset::Real, frequency::Real, phase::Real)

    # calculate sine current along one cycle
    t = collect(0:0.001:1)
    current = abs.(generate_sine(frequency, t, current, phase) .+ offset)
    current = simpson(current, t)
    
    cycles = frequency * duration
    _info("Number of cycles: $cycles")
    _info("Effective current: $current")

    charge = (current / 1_000) * duration
    current_density = (current / 1_000) / (pad_area / 1_000)
    charge_density = (charge / 1_000) / (pad_area / 1_000)

    return (charge=charge, current_density=current_density, charge_density=charge_density)

end

"""
    tpcs_dose(; <keyword arguments>)

Calculate `charge`, `current_density` and `charge_ density` for tPCS stimulation.

# Arguments

- `current::Real`: stimulation current (peak to peak) [mA]
- `pad_area::Real`: electrode pad area [cm²]
- `duration::Real`: stimulation duration [s]
- `pw::Real`: pulse width [ms]
- `isi::Real`: interstimulus interval [ms] (pulse width + interval)

# Returns

Named tuple containing:
- `charge::Float64`: charge [C]
- `current_density::Float64`: current density [A/m²]
- `charge_density::Float64`: delivered charge density [kC/m²]
"""
function tpcs_dose(; current::Real, pad_area::Real, duration::Real, pw::Real, isi::Real)

    @assert isi > pw "isi must be > pw."

    # convert to seconds
    pw = pw / 1000
    isi = isi / 1000

    cycles = duration / isi

    # calculate pulse wave current
    current = current * (pw / isi)
    _info("Number of cycles: $cycles")
    _info("Effective current: $current")

    charge = (current / 1_000) * duration
    current_density = (current / 1_000) / (pad_area / 1_000)
    charge_density = (charge / 1_000) / (pad_area / 1_000)

    return (charge=charge, current_density=current_density, charge_density=charge_density)

end

"""
    tes_protocol(; <keyword arguments>)

Create TES (tDCS/tACS/tRNS/tPCS) protocol.

# Arguments

- `type::Symbol`: stimulation type (`:tDCS`, `:tACS`, `:tRNS`, `:tPCS`)
- `hd::Bool`: high-density electrodes
- `current::Real`: stimulation current [mA]
- `frequency::Real=0`: stimulation frequency [mA]
- `anode_size::Tuple{Int64, Int64}`: anode dimensions [mm]
- `cathode_size::Tuple{Int64, Int64}`: cathode dimensions [mm]
- `anode_loc::Symbol`: anode location (according to 10-20 Positioning System)
- `cathode_loc::Symbol`: cathode location (according to 10-20 Positioning System)
- `duration::Real`: stimulation duration [s]
- `ramp_in::Real`: stimulation duration [s]
- `ramp_out::Real`: stimulation duration [s]
- `sham::Bool`: protocol includes sham stimulations

# Returns

- `protocol::Dict`
"""
function tes_protocol(; type::Symbol, hd::Bool, current::Real, frequency::Real=0, anode_size::Tuple{Int64, Int64}, cathode_size::Tuple{Int64, Int64}, anode_loc::Symbol, cathode_loc::Symbol, duration::Real, ramp_in::Real, ramp_out::Real, sham::Bool)

    _check_var(type, [:tDCS, :tACS, :tRNS, :tPCS], "type")
    @assert current > 0 "current must be > 0 mA."
    if type === :tACS || type === :tRNS
        @assert frequency > 0 "frequency must be > 0 mA."
    end
    (anode_size[1] <= 0 || anode_size[2] <= 0) && @error "anode dimensions > 0 mm."
    (cathode_size[1] <= 0 || cathode_size[2] <= 0) && @error "anode dimensions > 0 mm."
    @assert duration > 0 "duration must be > 0 s."
    @assert ramp_in >= 0 "ramp_in must be ≥ 0 s."
    @assert ramp_out >= 0 "ramp_out must be ≥ 0 s."

    protocol = Dict(:type=>type, :hd=>hd, :current=>current, :frequency=>frequency, :cathode_size=>cathode_size, :cathode_size=>cathode_size, :anode_loc=>anode_loc, :cathode_loc=>cathode_loc, :duration=>duration, :ramp_in=>ramp_in, :ramp_out=>ramp_out, :sham=>sham, )

    return protocol

end