using Test

@info "Test: tdcs_dose()"
@test tdcs_dose(current = 2.0, pad_area = 35, duration = 1200)[1] == 2.4

@info "Test: ect_charge()"
@test ect_charge(pw = 0.5, pint = 10, pf = 10, duration = 10) == 0.5

@info "Test: tes_protocol()"
p = tes_protocol(;
    type = :tDCS,
    hd = false,
    current = 2.0,
    anode_size = (50, 70),
    cathode_size = (50, 70),
    anode_loc = :F3,
    cathode_loc = :F4,
    duration = 1200,
    ramp_in = 20,
    ramp_out = 20,
    sham = false,
)
@test p isa Dict{Symbol, Any}

@info "Test: tacs_dose()"
@test round(
    tacs_dose(
        current = 2.0,
        pad_area = 35,
        duration = 1200,
        frequency = 1,
        phase = 0,
        offset = 0,
    )[1],
    digits = 2,
) == 0.76

@info "Test: tpcs_dose()"
@test tpcs_dose(current = 2.0, pad_area = 35, duration = 1200, pw = 10, isi = 100)[1] ==
      0.24


@info "Test: stimulation dose values"
# tDCS 2 mA, 35 cm², 1200 s: charge 2.4 C, current density 2 mA / 35 cm² = 0.571 A/m², charge density 0.686 kC/m²
d = tdcs_dose(current = 2.0, pad_area = 35, duration = 1200)
@test d.charge ≈ 2.4
@test d.current_density ≈ 0.002 / 0.0035
@test d.charge_density ≈ 2.4e-3 / 0.0035
# tACS 2 mA peak-to-peak: effective current = mean |1 mA × sin| = 2/π mA
d = tacs_dose(current = 2.0, pad_area = 35, duration = 1200, frequency = 1, phase = 0, offset = 0)
@test d.charge ≈ 2 / π * 1e-3 * 1200 rtol = 1e-5
@test d.current_density ≈ 2 / π * 1e-3 / 0.0035 rtol = 1e-5
# tPCS: duty cycle pw / isi
d = tpcs_dose(current = 2.0, pad_area = 35, duration = 1200, pw = 10, isi = 100)
@test d.charge ≈ 0.2e-3 * 1200
@test d.current_density ≈ 0.2e-3 / 0.0035
# ECT: charge = pw × current × frequency × duration (mC)
@test ect_charge(pw = 1.0, pint = 800, pf = 50, duration = 5) == 200.0
p = tes_protocol(; type = :tDCS, hd = false, current = 2.0, anode_size = (50, 70), cathode_size = (50, 70),
                 anode_loc = :F3, cathode_loc = :F4, duration = 1200, ramp_in = 20, ramp_out = 20, sham = false)
@test p[:type] == :tDCS && p[:current] == 2.0 && p[:anode_loc] == :F3 && p[:duration] == 1200

true
