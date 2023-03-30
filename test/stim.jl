using NeuroAnalyzer
using Test

@info "test 1/3: tes_dose()"
@test tes_dose(current=2.0, pad_area=35, duration=1200)[1] == 2.4

@info "test 2/3: ect_charge()"
@test ect_charge(pw=0.5, pint=10, pf=10, duration=10) == 500.0

@info "test 3/3: tes_protocol()"
p = tes_protocol(type=:tDCS, hd=false, current=2.0, anode_size=(50, 70), cathode_size=(50, 70), anode_loc=:F3, cathode_loc=:F4, duration=1200, ramp_in=20, ramp_out=20, sham=false)
@test p isa Dict{Symbol, Any}

true