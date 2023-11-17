export locs_generate
export locs_generate!

"""
    locs_generate(locs)

Generate spherical coordinates according to 10/5 system.

# Arguments

- `locs::DataFrame`

# Returns

- `locs_new::DataFrame`
"""
function locs_generate(locs::DataFrame)

    locs_new = deepcopy(locs)

    lab = locs[!, :labels]

    x = zeros(length(lab))
    y = zeros(length(lab))
    z = zeros(length(lab))

    x[lab .== "Cz"] .= sph2cart(1.0, 0, 90)[1]
    y[lab .== "Cz"] .= sph2cart(1.0, 0, 90)[2]
    z[lab .== "Cz"] .= sph2cart(1.0, 0, 90)[3]

    x[lab .== "Fz"] .= sph2cart(1.0, 90, 45)[1]
    y[lab .== "Fz"] .= sph2cart(1.0, 90, 45)[2]
    z[lab .== "Fz"] .= sph2cart(1.0, 90, 45)[3]

    x[lab .== "Pz"] .= sph2cart(1.0, 270, 45)[1]
    y[lab .== "Pz"] .= sph2cart(1.0, 270, 45)[2]
    z[lab .== "Pz"] .= sph2cart(1.0, 270, 45)[3]

    x[lab .== "T4" .|| lab .== "T8"] .= sph2cart(1.0, 0, 0)[1]
    y[lab .== "T4" .|| lab .== "T8"] .= sph2cart(1.0, 0, 0)[2]
    z[lab .== "T4" .|| lab .== "T8"] .= sph2cart(1.0, 0, 0)[3]
    
    x[lab .== "F8"] .= sph2cart(1.0, 36, 0)[1]
    y[lab .== "F8"] .= sph2cart(1.0, 36, 0)[2]
    z[lab .== "F8"] .= sph2cart(1.0, 36, 0)[3]

    x[lab .== "Fp2"] .= sph2cart(1.0, 72, 0)[1]
    y[lab .== "Fp2"] .= sph2cart(1.0, 72, 0)[2]
    z[lab .== "Fp2"] .= sph2cart(1.0, 72, 0)[3]

    x[lab .== "Fpz"] .= sph2cart(1.0, 90, 0)[1]
    y[lab .== "Fpz"] .= sph2cart(1.0, 90, 0)[2]
    z[lab .== "Fpz"] .= sph2cart(1.0, 90, 0)[3]

    x[lab .== "Fp1"] .= sph2cart(1.0, 108, 0)[1]
    y[lab .== "Fp1"] .= sph2cart(1.0, 108, 0)[2]
    z[lab .== "Fp1"] .= sph2cart(1.0, 108, 0)[3]

    x[lab .== "F7"] .= sph2cart(1.0, 144, 0)[1]
    y[lab .== "F7"] .= sph2cart(1.0, 144, 0)[2]
    z[lab .== "F7"] .= sph2cart(1.0, 144, 0)[3]

    x[lab .== "T3" .|| lab .== "T7"] .= sph2cart(1.0, 180, 0)[1]
    y[lab .== "T3" .|| lab .== "T7"] .= sph2cart(1.0, 180, 0)[2]
    z[lab .== "T3" .|| lab .== "T7"] .= sph2cart(1.0, 180, 0)[3]

    x[lab .== "T5"] .= sph2cart(1.0, 216, 0)[1]
    y[lab .== "T5"] .= sph2cart(1.0, 216, 0)[2]
    z[lab .== "T5"] .= sph2cart(1.0, 216, 0)[3]

    x[lab .== "O1"] .= sph2cart(1.0, 252, 0)[1]
    y[lab .== "O1"] .= sph2cart(1.0, 252, 0)[2]
    z[lab .== "O1"] .= sph2cart(1.0, 252, 0)[3]

    x[lab .== "Oz"] .= sph2cart(1.0, 270, 0)[1]
    y[lab .== "Oz"] .= sph2cart(1.0, 270, 0)[2]
    z[lab .== "Oz"] .= sph2cart(1.0, 270, 0)[3]

    x[lab .== "O2"] .= sph2cart(1.0, 288, 0)[1]
    y[lab .== "O2"] .= sph2cart(1.0, 288, 0)[2]
    z[lab .== "O2"] .= sph2cart(1.0, 288, 0)[3]

    x[lab .== "T6"] .= sph2cart(1.0, 324, 0)[1]
    y[lab .== "T6"] .= sph2cart(1.0, 324, 0)[2]
    z[lab .== "T6"] .= sph2cart(1.0, 324, 0)[3]

    x[lab .== "C3"] .= sph2cart(1.0, 180, 45)[1]
    y[lab .== "C3"] .= sph2cart(1.0, 180, 45)[2]
    z[lab .== "C3"] .= sph2cart(1.0, 180, 45)[3]

    x[lab .== "C4"] .= sph2cart(1.0, 0, 45)[1]
    y[lab .== "C4"] .= sph2cart(1.0, 0, 45)[2]
    z[lab .== "C4"] .= sph2cart(1.0, 0, 45)[3]

    x[lab .== "F4"] .= sph2cart(1.0, 45, 45)[1]
    y[lab .== "F4"] .= sph2cart(1.0, 45, 45)[2]
    z[lab .== "F4"] .= sph2cart(1.0, 45, 45)[3]

    x[lab .== "F3"] .= sph2cart(1.0, 135, 45)[1]
    y[lab .== "F3"] .= sph2cart(1.0, 135, 45)[2]
    z[lab .== "F3"] .= sph2cart(1.0, 135, 45)[3]

    x[lab .== "P3"] .= sph2cart(1.0, 225, 45)[1]
    y[lab .== "P3"] .= sph2cart(1.0, 225, 45)[2]
    z[lab .== "P3"] .= sph2cart(1.0, 225, 45)[3]

    x[lab .== "P4"] .= sph2cart(1.0, 315, 45)[1]
    y[lab .== "P4"] .= sph2cart(1.0, 315, 45)[2]
    z[lab .== "P4"] .= sph2cart(1.0, 315, 45)[3]

    x[lab .== "A1"] .= -0.92
    y[lab .== "A1"] .= -0.23
    z[lab .== "A1"] .= -0.55
    
    x[lab .== "A2"] .= 0.92
    y[lab .== "A2"] .= -0.23
    z[lab .== "A2"] .= -0.55

    x[lab .== "M1"] .= -0.94
    y[lab .== "M1"] .= -0.10
    z[lab .== "M1"] .= -0.30

    x[lab .== "M2"] .= 0.94
    y[lab .== "M2"] .= -0.10
    z[lab .== "M2"] .= -0.30

    x[lab .== "EMG1"] .= -0.70
    y[lab .== "EMG1"] .= 0.70
    z[lab .== "EMG1"] .= -1.10

    x[lab .== "EMG2"] .= 0.70
    y[lab .== "EMG2"] .= 0.70
    z[lab .== "EMG2"] .= -1.10

    x[lab .== "EOG1"] .= -0.87
    y[lab .== "EOG1"] .= 0.51
    z[lab .== "EOG1"] .= -0.37

    x[lab .== "EOG2"] .= 0.87
    y[lab .== "EOG2"] .= 0.51
    z[lab .== "EOG2"] .= -0.37

    x[lab .== "VEOG1"] .= -0.87
    y[lab .== "VEOG1"] .= 0.51
    z[lab .== "VEOG1"] .= -0.37

    x[lab .== "VEOG2"] .= 0.87
    y[lab .== "VEOG2"] .= 0.51
    z[lab .== "VEOG2"] .= -0.37

    x[lab .== "HEOG1"] .= -0.64
    y[lab .== "HEOG1"] .= 0.77
    z[lab .== "HEOG1"] .= -0.04

    x[lab .== "HEOG2"] .= 0.64
    y[lab .== "HEOG2"] .= 0.77
    z[lab .== "HEOG2"] .= -0.04

    x[lab .== "VEOG"] .= 0.87
    y[lab .== "VEOG"] .= 0.51
    z[lab .== "VEOG"] .= -0.37

    x[lab .== "HEOG"] .= 0.64
    y[lab .== "HEOG"] .= 0.77
    z[lab .== "HEOG"] .= -0.04

    locs_new[:, :loc_x] = x
    locs_new[:, :loc_y] = y
    locs_new[:, :loc_z] = z

    locs_cart2sph!(locs_new) 
    locs_sph2pol!(locs_new) 

    return locs_new

end

"""
    locs_generate!(locs)

Generate spherical coordinates according to 10/5 system.

# Arguments

- `locs::DataFrame`
"""
function locs_generate!(locs::DataFrame)

    locs_tmp = locs_generate(locs)

    locs[:, :loc_radius] = locs_tmp[:, :loc_radius]
    locs[:, :loc_theta] = locs_tmp[:, :loc_theta]
    locs[:, :loc_x] = locs_tmp[:, :loc_x]
    locs[:, :loc_y] = locs_tmp[:, :loc_y]
    locs[:, :loc_z] = locs_tmp[:, :loc_z]
    locs[:, :loc_radius_sph] = locs_tmp[:, :loc_radius_sph]
    locs[:, :loc_theta_sph] = locs_tmp[:, :loc_theta_sph]
    locs[:, :loc_phi_sph] = locs_tmp[:, :loc_phi_sph]

    return nothing

end