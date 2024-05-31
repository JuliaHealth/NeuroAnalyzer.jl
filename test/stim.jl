using Test

@info "Test 1/4: tdcs_dose()"
@test tdcs_dose(current=2.0, pad_area=35, duration=1200)[1] == 2.4

@info "Test 2/4: ect_charge()"
@test ect_charge(pw=0.5, pint=10, pf=10, duration=10) == 500.0

@info "Test 3/4: tes_protocol()"
p = tes_protocol(type=:tDCS, hd=false, current=2.0, anode_size=(50, 70), cathode_size=(50, 70), anode_loc=:F3, cathode_loc=:F4, duration=1200, ramp_in=20, ramp_out=20, sham=false)
@test p isa Dict{Symbol, Any}

@info "Test 4/4: tacs_dose()"
@test tacs_dose(current=2.0, pad_area=35, duration=1200, frequency=1, phase=0, offset=0)[1] == 0.7639437268477123

true