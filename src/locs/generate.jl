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

    lab = lowercase.(locs[!, :labels])

    x = zeros(length(lab))
    y = zeros(length(lab))
    z = zeros(length(lab))

    r = zeros(length(lab))
    t = zeros(length(lab))

    x[lab .== "afz"] .= sph2cart(1.0, 90, 22.5)[1]
    y[lab .== "afz"] .= sph2cart(1.0, 90, 22.5)[2]
    z[lab .== "afz"] .= sph2cart(1.0, 90, 22.5)[3]

    x[lab .== "fz"] .= sph2cart(1.0, 90, 45)[1]
    y[lab .== "fz"] .= sph2cart(1.0, 90, 45)[2]
    z[lab .== "fz"] .= sph2cart(1.0, 90, 45)[3]

    x[lab .== "fcz"] .= sph2cart(1.0, 90, 67.5)[1]
    y[lab .== "fcz"] .= sph2cart(1.0, 90, 67.5)[2]
    z[lab .== "fcz"] .= sph2cart(1.0, 90, 67.5)[3]

    x[lab .== "cz"] .= sph2cart(1.0, 0, 90)[1]
    y[lab .== "cz"] .= sph2cart(1.0, 0, 90)[2]
    z[lab .== "cz"] .= sph2cart(1.0, 0, 90)[3]

    x[lab .== "cpz"] .= sph2cart(1.0, 270, 67.5)[1]
    y[lab .== "cpz"] .= sph2cart(1.0, 270, 67.5)[2]
    z[lab .== "cpz"] .= sph2cart(1.0, 270, 67.5)[3]

    x[lab .== "pz"] .= sph2cart(1.0, 270, 45)[1]
    y[lab .== "pz"] .= sph2cart(1.0, 270, 45)[2]
    z[lab .== "pz"] .= sph2cart(1.0, 270, 45)[3]

    x[lab .== "poz"] .= sph2cart(1.0, 270, 22.5)[1]
    y[lab .== "poz"] .= sph2cart(1.0, 270, 22.5)[2]
    z[lab .== "poz"] .= sph2cart(1.0, 270, 22.5)[3]

    x[lab .== "t4" .|| lab .== "t8"] .= sph2cart(1.0, 0, 0)[1]
    y[lab .== "t4" .|| lab .== "t8"] .= sph2cart(1.0, 0, 0)[2]
    z[lab .== "t4" .|| lab .== "t8"] .= sph2cart(1.0, 0, 0)[3]
    
    x[lab .== "fc8"] .= sph2cart(1.0, 18, 0)[1]
    y[lab .== "fc8"] .= sph2cart(1.0, 18, 0)[2]
    z[lab .== "fc8"] .= sph2cart(1.0, 18, 0)[3]
    
    x[lab .== "f8"] .= sph2cart(1.0, 36, 0)[1]
    y[lab .== "f8"] .= sph2cart(1.0, 36, 0)[2]
    z[lab .== "f8"] .= sph2cart(1.0, 36, 0)[3]

    x[lab .== "af8"] .= sph2cart(1.0, 54, 0)[1]
    y[lab .== "af8"] .= sph2cart(1.0, 54, 0)[2]
    z[lab .== "af8"] .= sph2cart(1.0, 54, 0)[3]

    x[lab .== "fp2"] .= sph2cart(1.0, 72, 0)[1]
    y[lab .== "fp2"] .= sph2cart(1.0, 72, 0)[2]
    z[lab .== "fp2"] .= sph2cart(1.0, 72, 0)[3]

    x[lab .== "fpz"] .= sph2cart(1.0, 90, 0)[1]
    y[lab .== "fpz"] .= sph2cart(1.0, 90, 0)[2]
    z[lab .== "fpz"] .= sph2cart(1.0, 90, 0)[3]

    x[lab .== "fp1"] .= sph2cart(1.0, 108, 0)[1]
    y[lab .== "fp1"] .= sph2cart(1.0, 108, 0)[2]
    z[lab .== "fp1"] .= sph2cart(1.0, 108, 0)[3]

    x[lab .== "af7"] .= sph2cart(1.0, 126, 0)[1]
    y[lab .== "af7"] .= sph2cart(1.0, 126, 0)[2]
    z[lab .== "af7"] .= sph2cart(1.0, 126, 0)[3]

    x[lab .== "f7"] .= sph2cart(1.0, 144, 0)[1]
    y[lab .== "f7"] .= sph2cart(1.0, 144, 0)[2]
    z[lab .== "f7"] .= sph2cart(1.0, 144, 0)[3]

    x[lab .== "ft7"] .= sph2cart(1.0, 162, 0)[1]
    y[lab .== "ft7"] .= sph2cart(1.0, 162, 0)[2]
    z[lab .== "ft7"] .= sph2cart(1.0, 162, 0)[3]

    x[lab .== "t3" .|| lab .== "t7"] .= sph2cart(1.0, 180, 0)[1]
    y[lab .== "t3" .|| lab .== "t7"] .= sph2cart(1.0, 180, 0)[2]
    z[lab .== "t3" .|| lab .== "t7"] .= sph2cart(1.0, 180, 0)[3]

    y[lab .== "tp7"] .= sph2cart(1.0, 198, 0)[2]
    x[lab .== "tp7"] .= sph2cart(1.0, 198, 0)[1]
    z[lab .== "tp7"] .= sph2cart(1.0, 198, 0)[3]

    y[lab .== "t5" .|| lab .== "p7"] .= sph2cart(1.0, 216, 0)[2]
    x[lab .== "t5" .|| lab .== "p7"] .= sph2cart(1.0, 216, 0)[1]
    z[lab .== "t5" .|| lab .== "p7"] .= sph2cart(1.0, 216, 0)[3]

    y[lab .== "po7"] .= sph2cart(1.0, 234, 0)[2]
    x[lab .== "po7"] .= sph2cart(1.0, 234, 0)[1]
    z[lab .== "po7"] .= sph2cart(1.0, 234, 0)[3]

    x[lab .== "o1"] .= sph2cart(1.0, 252, 0)[1]
    y[lab .== "o1"] .= sph2cart(1.0, 252, 0)[2]
    z[lab .== "o1"] .= sph2cart(1.0, 252, 0)[3]

    x[lab .== "oz"] .= sph2cart(1.0, 270, 0)[1]
    y[lab .== "oz"] .= sph2cart(1.0, 270, 0)[2]
    z[lab .== "oz"] .= sph2cart(1.0, 270, 0)[3]

    x[lab .== "o2"] .= sph2cart(1.0, 288, 0)[1]
    y[lab .== "o2"] .= sph2cart(1.0, 288, 0)[2]
    z[lab .== "o2"] .= sph2cart(1.0, 288, 0)[3]

    x[lab .== "po8"] .= sph2cart(1.0, 306, 0)[1]
    y[lab .== "po8"] .= sph2cart(1.0, 306, 0)[2]
    z[lab .== "po8"] .= sph2cart(1.0, 306, 0)[3]

    x[lab .== "t6" .|| lab .== "p8"] .= sph2cart(1.0, 324, 0)[1]
    y[lab .== "t6" .|| lab .== "p8"] .= sph2cart(1.0, 324, 0)[2]
    z[lab .== "t6" .|| lab .== "p8"] .= sph2cart(1.0, 324, 0)[3]

    x[lab .== "tp8"] .= sph2cart(1.0, 342, 0)[1]
    y[lab .== "tp8"] .= sph2cart(1.0, 342, 0)[2]
    z[lab .== "tp8"] .= sph2cart(1.0, 342, 0)[3]

    x[lab .== "c5"] .= sph2cart(1.0, 180, 22.5)[1]
    y[lab .== "c5"] .= sph2cart(1.0, 180, 22.5)[2]
    z[lab .== "c5"] .= sph2cart(1.0, 180, 22.5)[3]

    x[lab .== "c3"] .= sph2cart(1.0, 180, 45)[1]
    y[lab .== "c3"] .= sph2cart(1.0, 180, 45)[2]
    z[lab .== "c3"] .= sph2cart(1.0, 180, 45)[3]

    x[lab .== "c1"] .= sph2cart(1.0, 180, 67.5)[1]
    y[lab .== "c1"] .= sph2cart(1.0, 180, 67.5)[2]
    z[lab .== "c1"] .= sph2cart(1.0, 180, 67.5)[3]

    x[lab .== "c6"] .= sph2cart(1.0, 0, 22.5)[1]
    y[lab .== "c6"] .= sph2cart(1.0, 0, 22.5)[2]
    z[lab .== "c6"] .= sph2cart(1.0, 0, 22.5)[3]

    x[lab .== "c4"] .= sph2cart(1.0, 0, 45)[1]
    y[lab .== "c4"] .= sph2cart(1.0, 0, 45)[2]
    z[lab .== "c4"] .= sph2cart(1.0, 0, 45)[3]

    x[lab .== "c2"] .= sph2cart(1.0, 0, 67.5)[1]
    y[lab .== "c2"] .= sph2cart(1.0, 0, 67.5)[2]
    z[lab .== "c2"] .= sph2cart(1.0, 0, 67.5)[3]

    x[lab .== "f4"] .= sph2cart(1.0, 45, 45)[1]
    y[lab .== "f4"] .= sph2cart(1.0, 45, 45)[2]
    z[lab .== "f4"] .= sph2cart(1.0, 45, 45)[3]

    x[lab .== "f3"] .= sph2cart(1.0, 135, 45)[1]
    y[lab .== "f3"] .= sph2cart(1.0, 135, 45)[2]
    z[lab .== "f3"] .= sph2cart(1.0, 135, 45)[3]

    x[lab .== "p3"] .= sph2cart(1.0, 225, 45)[1]
    y[lab .== "p3"] .= sph2cart(1.0, 225, 45)[2]
    z[lab .== "p3"] .= sph2cart(1.0, 225, 45)[3]

    x[lab .== "p4"] .= sph2cart(1.0, 315, 45)[1]
    y[lab .== "p4"] .= sph2cart(1.0, 315, 45)[2]
    z[lab .== "p4"] .= sph2cart(1.0, 315, 45)[3]

    x[lab .== "a1"] .= -0.92
    y[lab .== "a1"] .= -0.23
    z[lab .== "a1"] .= -0.55
    
    x[lab .== "a2"] .= 0.92
    y[lab .== "a2"] .= -0.23
    z[lab .== "a2"] .= -0.55

    x[lab .== "m1"] .= -0.94
    y[lab .== "m1"] .= -0.10
    z[lab .== "m1"] .= -0.30

    x[lab .== "m2"] .= 0.94
    y[lab .== "m2"] .= -0.10
    z[lab .== "m2"] .= -0.30

    x[lab .== "emg1"] .= -0.70
    y[lab .== "emg1"] .= 0.70
    z[lab .== "emg1"] .= -1.10

    x[lab .== "emg2"] .= 0.70
    y[lab .== "emg2"] .= 0.70
    z[lab .== "emg2"] .= -1.10

    x[lab .== "eog1"] .= -0.87
    y[lab .== "eog1"] .= 0.51
    z[lab .== "eog1"] .= -0.37

    x[lab .== "eog2"] .= 0.87
    y[lab .== "eog2"] .= 0.51
    z[lab .== "eog2"] .= -0.37

    x[lab .== "veog1"] .= -0.87
    y[lab .== "veog1"] .= 0.51
    z[lab .== "veog1"] .= -0.37

    x[lab .== "veog2"] .= 0.87
    y[lab .== "veog2"] .= 0.51
    z[lab .== "veog2"] .= -0.37

    x[lab .== "heog1"] .= -0.64
    y[lab .== "heog1"] .= 0.77
    z[lab .== "heog1"] .= -0.04

    x[lab .== "heog2"] .= 0.64
    y[lab .== "heog2"] .= 0.77
    z[lab .== "heog2"] .= -0.04

    x[lab .== "veog"] .= 0.87
    y[lab .== "veog"] .= 0.51
    z[lab .== "veog"] .= -0.37

    x[lab .== "heog"] .= 0.64
    y[lab .== "heog"] .= 0.77
    z[lab .== "heog"] .= -0.04

    locs_new[:, :loc_x] = x
    locs_new[:, :loc_y] = y
    locs_new[:, :loc_z] = z

    locs_cart2sph!(locs_new) 
   
    # TO DO: generate polar coordinates
    locs_sph2pol!(locs_new) 

    r[lab .== "afz"] .= 0.75
    t[lab .== "afz"] .= 90

    r[lab .== "fz"] .= 0.5
    t[lab .== "fz"] .= 90

    r[lab .== "fcz"] .= 0.25
    t[lab .== "fcz"] .= 90

    r[lab .== "cz"] .= 0
    t[lab .== "cz"] .= 0

    r[lab .== "cpz"] .= 0.25
    t[lab .== "cpz"] .= 270

    r[lab .== "pz"] .= 0.5
    t[lab .== "pz"] .= 270

    r[lab .== "poz"] .= 0.75
    t[lab .== "poz"] .= 270

    r[lab .== "c2"] .= 0.25
    t[lab .== "c2"] .= 0

    r[lab .== "c4"] .= 0.5
    t[lab .== "c4"] .= 0

    r[lab .== "c6"] .= 0.75
    t[lab .== "c6"] .= 0

    r[lab .== "t4" .|| lab .== "t8"] .= 1
    t[lab .== "t4" .|| lab .== "t8"] .= 0

    r[lab .== "c1"] .= 0.25
    t[lab .== "c1"] .= 180

    r[lab .== "c3"] .= 0.5
    t[lab .== "c3"] .= 180

    r[lab .== "c5"] .= 0.75
    t[lab .== "c5"] .= 180

    r[lab .== "fc8"] .= 1
    t[lab .== "fc8"] .= 18
    
    r[lab .== "f8"] .= 1
    t[lab .== "f8"] .= 36

    r[lab .== "af8"] .= 1
    t[lab .== "af8"] .= 54

    r[lab .== "fp2"] .= 1
    t[lab .== "fp2"] .= 72

    r[lab .== "fpz"] .= 1
    t[lab .== "fpz"] .= 90

    r[lab .== "fp1"] .= 1
    t[lab .== "fp1"] .= 108

    r[lab .== "af7"] .= 1
    t[lab .== "af7"] .= 126

    r[lab .== "f7"] .= 1
    t[lab .== "f7"] .= 144

    r[lab .== "ft7"] .= 1
    t[lab .== "ft7"] .= 162

    r[lab .== "t3" .|| lab .== "t7"] .= 1
    t[lab .== "t3" .|| lab .== "t7"] .= 180

    r[lab .== "tp7"] .= 1
    t[lab .== "tp7"] .= 198

    r[lab .== "t5" .|| lab .== "p7"] .= 1
    t[lab .== "t5" .|| lab .== "p7"] .= 216

    r[lab .== "po7"] .= 1
    t[lab .== "po7"] .= 234

    r[lab .== "o1"] .= 1
    t[lab .== "o1"] .= 252

    r[lab .== "oz"] .= 1
    t[lab .== "oz"] .= 270

    r[lab .== "o2"] .= 1
    t[lab .== "o2"] .= 288

    r[lab .== "po8"] .= 1
    t[lab .== "po8"] .= 306

    r[lab .== "t6" .|| lab .== "p8"] .= 1
    t[lab .== "t6" .|| lab .== "p8"] .= 324

    r[lab .== "tp8"] .= 1
    t[lab .== "tp8"] .= 342

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