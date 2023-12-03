export locs_generate
export locs_generate!

"""
    locs_generate(locs)

Generate spherical coordinates according to 10/10 system.

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

    e_labels = ["cz", "c2", "c4", "c6", "t4", "t8", "t10", "c1", "c3", "c5", "t3", "t7", "t9", "fcz", "fc2", "fc4", "fc6", "fc8", "ft8", "fc10", "ft10", "fc1", "fc3", "fc5", "fc7", "ft7", "fc9", "ft9", "fz", "f2", "f4", "f6", "f8", "f10", "f1", "f3", "f5", "f7", "f9", "afz", "af2", "af4", "af6", "af1", "af3", "af7", "fpz", "fp2", "fp1", "cpz", "cp2", "cp4", "cp6", "cp8", "tp8", "tp10", "cp1", "cp3", "cp5", "cp7", "tp7", "tp9", "pz", "p2", "p4", "p6", "p8", "p10", "t6", "p1", "p3", "p5", "p7", "p9", "t5", "poz", "po2", "po4", "po6", "po8", "po1", "po3", "po5", "po7", "oz", "o2", "o1", "a1", "a2", "m1", "m2", "emg1", "emg2", "eog1", "eog2", "veog1", "veog2", "heog1", "heog2", "veog", "heog"]

    x[lab .== "cz"] .= sph2cart(1.0, 0, 90)[1]
    y[lab .== "cz"] .= sph2cart(1.0, 0, 90)[2]
    z[lab .== "cz"] .= sph2cart(1.0, 0, 90)[3]

    x[lab .== "c2"] .= sph2cart(1.0, 0, 67.5)[1]
    y[lab .== "c2"] .= sph2cart(1.0, 0, 67.5)[2]
    z[lab .== "c2"] .= sph2cart(1.0, 0, 67.5)[3]

    x[lab .== "c4"] .= sph2cart(1.0, 0, 45)[1]
    y[lab .== "c4"] .= sph2cart(1.0, 0, 45)[2]
    z[lab .== "c4"] .= sph2cart(1.0, 0, 45)[3]

    x[lab .== "c6"] .= sph2cart(1.0, 0, 22.5)[1]
    y[lab .== "c6"] .= sph2cart(1.0, 0, 22.5)[2]
    z[lab .== "c6"] .= sph2cart(1.0, 0, 22.5)[3]

    x[lab .== "t4" .|| lab .== "t8"] .= sph2cart(1.0, 0, 0)[1]
    y[lab .== "t4" .|| lab .== "t8"] .= sph2cart(1.0, 0, 0)[2]
    z[lab .== "t4" .|| lab .== "t8"] .= sph2cart(1.0, 0, 0)[3]

    x[lab .== "t10"] .= sph2cart(1.0, -22.5, 0)[1]
    y[lab .== "t10"] .= sph2cart(1.0, -22.5, 0)[2]
    z[lab .== "t10"] .= sph2cart(1.0, -22.5, 0)[3]

    x[lab .== "c1"] .= sph2cart(1.0, 0, 112.5)[1]
    y[lab .== "c1"] .= sph2cart(1.0, 0, 112.5)[2]
    z[lab .== "c1"] .= sph2cart(1.0, 0, 112.5)[3]

    x[lab .== "c3"] .= sph2cart(1.0, 0, 135)[1]
    y[lab .== "c3"] .= sph2cart(1.0, 0, 135)[2]
    z[lab .== "c3"] .= sph2cart(1.0, 0, 135)[3]

    x[lab .== "c5"] .= sph2cart(1.0, 0, 157.5)[1]
    y[lab .== "c5"] .= sph2cart(1.0, 0, 157.5)[2]
    z[lab .== "c5"] .= sph2cart(1.0, 0, 157.5)[3]

    x[lab .== "t3" .|| lab .== "t7"] .= sph2cart(1.0, 0, 180)[1]
    y[lab .== "t3" .|| lab .== "t7"] .= sph2cart(1.0, 0, 180)[2]
    z[lab .== "t3" .|| lab .== "t7"] .= sph2cart(1.0, 0, 180)[3]

    x[lab .== "t9"] .= sph2cart(1.0, -22.5, 180)[1]
    y[lab .== "t9"] .= sph2cart(1.0, -22.5, 180)[2]
    z[lab .== "t9"] .= sph2cart(1.0, -22.5, 180)[3]

    x[lab .== "fcz"] .= sph2cart(1.0, 90, 67.5)[1]
    y[lab .== "fcz"] .= sph2cart(1.0, 90, 67.5)[2]
    z[lab .== "fcz"] .= sph2cart(1.0, 90, 67.5)[3]

    x[lab .== "fc2"] .= sph2cart(cosd(22.5), 0, 67.5)[1]
    y[lab .== "fc2"] .= cosd(22.5) + sph2cart(cosd(22.5), 0, 67.5)[2]
    z[lab .== "fc2"] .= sph2cart(cosd(22.5), 0, 67.5)[3]
    
    x[lab .== "fc4"] .= sph2cart(cosd(22.5), 0, 45)[1]
    y[lab .== "fc4"] .= cosd(22.5) + sph2cart(cosd(22.5), 0, 45)[2]
    z[lab .== "fc4"] .= sph2cart(cosd(22.5), 0, 45)[3]
    
    x[lab .== "fc6"] .= sph2cart(cosd(22.5), 0, 22.5)[1]
    y[lab .== "fc6"] .= cosd(22.5) + sph2cart(cosd(22.5), 0, 22.5)[2]
    z[lab .== "fc6"] .= sph2cart(cosd(22.5), 0, 22.5)[3]

    x[lab .== "fc8" .|| lab .== "ft8"] .= sph2cart(cosd(22.5), 0, 0)[1]
    y[lab .== "fc8" .|| lab .== "ft8"] .= cosd(22.5) + sph2cart(cosd(22.5), 0, 0)[2]
    z[lab .== "fc8" .|| lab .== "ft8"] .= sph2cart(cosd(22.5), 0, 0)[3]

    x[lab .== "fc10" .|| lab .== "ft10"] .= sph2cart(cosd(22.5), -22.5, 0)[1]
    y[lab .== "fc10" .|| lab .== "ft10"] .= cosd(22.5) + sph2cart(cosd(22.5), -22.5, 0)[2]
    z[lab .== "fc10" .|| lab .== "ft10"] .= sph2cart(cosd(22.5), -22.5, 0)[3]

    x[lab .== "fc1"] .= sph2cart(cosd(22.5), 0, 112.5)[1]
    y[lab .== "fc1"] .= cosd(22.5) + sph2cart(cosd(22.5), 0, 112.5)[2]
    z[lab .== "fc1"] .= sph2cart(cosd(22.5), 0, 112.5)[3]
    
    x[lab .== "fc3"] .= sph2cart(cosd(22.5), 0, 135)[1]
    y[lab .== "fc3"] .= cosd(22.5) + sph2cart(cosd(22.5), 0, 135)[2]
    z[lab .== "fc3"] .= sph2cart(cosd(22.5), 0, 135)[3]
    
    x[lab .== "fc5"] .= sph2cart(cosd(22.5), 0, 157.5)[1]
    y[lab .== "fc5"] .= cosd(22.5) + sph2cart(cosd(22.5), 0, 157.5)[2]
    z[lab .== "fc5"] .= sph2cart(cosd(22.5), 0, 157.5)[3]

    x[lab .== "fc7" .|| lab .== "ft6"] .= sph2cart(cosd(22.5), 0, 180)[1]
    y[lab .== "fc7" .|| lab .== "ft6"] .= cosd(22.5) + sph2cart(cosd(22.5), 0, 180)[2]
    z[lab .== "fc7" .|| lab .== "ft6"] .= sph2cart(cosd(22.5), 0, 180)[3]

    x[lab .== "fc9" .|| lab .== "ft9"] .= sph2cart(cosd(22.5), -22.5, 180)[1]
    y[lab .== "fc9" .|| lab .== "ft9"] .= cosd(22.5) + sph2cart(cosd(22.5), -22.5, 180)[2]
    z[lab .== "fc9" .|| lab .== "ft9"] .= sph2cart(cosd(22.5), -22.5, 180)[3]

    x[lab .== "fz"] .= sph2cart(1.0, 90, 45)[1]
    y[lab .== "fz"] .= sph2cart(1.0, 90, 45)[2]
    z[lab .== "fz"] .= sph2cart(1.0, 90, 45)[3]
    
    x[lab .== "f2"] .= sph2cart(cosd(45), 0, 67.5)[1]
    y[lab .== "f2"] .= cosd(45) + sph2cart(cosd(45), 0, 67.5)[2]
    z[lab .== "f2"] .= sph2cart(cosd(45), 0, 67.5)[3]
    
    x[lab .== "f4"] .= sph2cart(cosd(45), 0, 45)[1]
    y[lab .== "f4"] .= cosd(45) + sph2cart(cosd(45), 0, 45)[2]
    z[lab .== "f4"] .= sph2cart(cosd(45), 0, 45)[3]
    
    x[lab .== "f6"] .= sph2cart(cosd(45), 0, 22.5)[1]
    y[lab .== "f6"] .= cosd(45) + sph2cart(cosd(45), 0, 22.5)[2]
    z[lab .== "f6"] .= sph2cart(cosd(45), 0, 22.5)[3]
    
    x[lab .== "f8"] .= sph2cart(cosd(45), 0, 0)[1]
    y[lab .== "f8"] .= cosd(45) + sph2cart(cosd(45), 0, 0)[2]
    z[lab .== "f8"] .= sph2cart(cosd(45), 0, 0)[3]

    x[lab .== "f10"] .= sph2cart(cosd(45), -22.5, 0)[1]
    y[lab .== "f10"] .= cosd(45) + sph2cart(cosd(45), -22.5, 0)[2]
    z[lab .== "f10"] .= sph2cart(cosd(45), -22.5, 0)[3]

    x[lab .== "f1"] .= sph2cart(cosd(45), 0, 112.5)[1]
    y[lab .== "f1"] .= cosd(45) + sph2cart(cosd(45), 0, 112.5)[2]
    z[lab .== "f1"] .= sph2cart(cosd(45), 0, 112.5)[3]
    
    x[lab .== "f3"] .= sph2cart(cosd(45), 0, 135)[1]
    y[lab .== "f3"] .= cosd(45) + sph2cart(cosd(45), 0, 135)[2]
    z[lab .== "f3"] .= sph2cart(cosd(45), 0, 135)[3]
    
    x[lab .== "f5"] .= sph2cart(cosd(45), 0, 157.5)[1]
    y[lab .== "f5"] .= cosd(45) + sph2cart(cosd(45), 0, 157.5)[2]
    z[lab .== "f5"] .= sph2cart(cosd(45), 0, 157.5)[3]
    
    x[lab .== "f7"] .= sph2cart(cosd(45), 0, 180)[1]
    y[lab .== "f7"] .= cosd(45) + sph2cart(cosd(45), 0, 180)[2]
    z[lab .== "f7"] .= sph2cart(cosd(45), 0, 180)[3]
    
    x[lab .== "f9"] .= sph2cart(cosd(45), -22.5, 180)[1]
    y[lab .== "f9"] .= cosd(45) + sph2cart(cosd(45), -22.5, 180)[2]
    z[lab .== "f9"] .= sph2cart(cosd(45), -22.5, 180)[3]
    
    x[lab .== "afz"] .= sph2cart(1.0, 90, 67.5)[1]
    y[lab .== "afz"] .= sph2cart(1.0, 90, 67.5)[2]
    z[lab .== "afz"] .= sph2cart(1.0, 90, 67.5)[3]
        
    x[lab .== "af2"] .= sph2cart(cosd(67.5), 0, 67.5)[1]
    y[lab .== "af2"] .= cosd(67.5) + sph2cart(cosd(67.5), 0, 67.5)[2]
    z[lab .== "af2"] .= sph2cart(cosd(67.5), 0, 67.5)[3]
    
    x[lab .== "af4"] .= sph2cart(cosd(67.5), 0, 45)[1]
    y[lab .== "af4"] .= cosd(67.5) + sph2cart(cosd(67.5), 0, 45)[2]
    z[lab .== "af4"] .= sph2cart(cosd(67.5), 0, 45)[3]
    
    x[lab .== "af6"] .= sph2cart(cosd(67.5), 0, 22.5)[1]
    y[lab .== "af6"] .= cosd(67.5) + sph2cart(cosd(67.5), 0, 22.5)[2]
    z[lab .== "af6"] .= sph2cart(cosd(67.5), 0, 22.5)[3]
    
    x[lab .== "af1"] .= sph2cart(cosd(67.5), 0, 112.5)[1]
    y[lab .== "af1"] .= cosd(67.5) + sph2cart(cosd(67.5), 0, 112.5)[2]
    z[lab .== "af1"] .= sph2cart(cosd(67.5), 0, 112.5)[3]
    
    x[lab .== "af3"] .= sph2cart(cosd(67.5), 0, 135)[1]
    y[lab .== "af3"] .= cosd(67.5) + sph2cart(cosd(67.5), 0, 135)[2]
    z[lab .== "af3"] .= sph2cart(cosd(67.5), 0, 135)[3]
    
    x[lab .== "af7"] .= sph2cart(cosd(67.5), 0, 157.5)[1]
    y[lab .== "af7"] .= cosd(67.5) + sph2cart(cosd(67.5), 0, 157.5)[2]
    z[lab .== "af7"] .= sph2cart(cosd(67.5), 0, 157.5)[3]

    x[lab .== "fpz"] .= sph2cart(1.0, 90, 0)[1]
    y[lab .== "fpz"] .= sph2cart(1.0, 90, 0)[2]
    z[lab .== "fpz"] .= sph2cart(1.0, 90, 0)[3]
        
    x[lab .== "fp2"] .= sph2cart(1.0, 67.5, 0)[1]
    y[lab .== "fp2"] .= sph2cart(1.0, 67.5, 0)[2]
    z[lab .== "fp2"] .= sph2cart(1.0, 67.5, 0)[3]
        
    x[lab .== "fp1"] .= sph2cart(1.0, 112.5, 0)[1]
    y[lab .== "fp1"] .= sph2cart(1.0, 112.5, 0)[2]
    z[lab .== "fp1"] .= sph2cart(1.0, 112.5, 0)[3]

    x[lab .== "cpz"] .= sph2cart(1.0, 270, 67.5)[1]
    y[lab .== "cpz"] .= sph2cart(1.0, 270, 67.5)[2]
    z[lab .== "cpz"] .= sph2cart(1.0, 270, 67.5)[3]

    x[lab .== "cp2"] .= sph2cart(cosd(22.5), 0, 67.5)[1]
    y[lab .== "cp2"] .= -cosd(67.5) + sph2cart(cosd(22.5), 0, 67.5)[2]
    z[lab .== "cp2"] .= sph2cart(cosd(22.5), 0, 67.5)[3]
    
    x[lab .== "cp4"] .= sph2cart(cosd(22.5), 0, 45)[1]
    y[lab .== "cp4"] .= -cosd(67.5) + sph2cart(cosd(22.5), 0, 45)[2]
    z[lab .== "cp4"] .= sph2cart(cosd(22.5), 0, 45)[3]
    
    x[lab .== "cp6"] .= sph2cart(cosd(22.5), 0, 22.5)[1]
    y[lab .== "cp6"] .= -cosd(67.5) + sph2cart(cosd(22.5), 0, 22.5)[2]
    z[lab .== "cp6"] .= sph2cart(cosd(22.5), 0, 22.5)[3]
    
    x[lab .== "cp8" .|| lab .== "tp8"] .= sph2cart(cosd(22.5), 0, 0)[1]
    y[lab .== "cp8" .|| lab .== "tp8"] .= -cosd(67.5) + sph2cart(cosd(22.5), 0, 0)[2]
    z[lab .== "cp8" .|| lab .== "tp8"] .= sph2cart(cosd(22.5), 0, 0)[3]

    x[lab .== "tp10"] .= sph2cart(cosd(22.5), -22.5, 0)[1]
    y[lab .== "tp10"] .= -cosd(67.5) + sph2cart(cosd(22.5), -22.5, 0)[2]
    z[lab .== "tp10"] .= sph2cart(cosd(22.5), -22.5, 0)[3]

    x[lab .== "tp9"] .= sph2cart(cosd(22.5), -22.5, 0)[1]
    y[lab .== "tp9"] .= -cosd(67.5) + sph2cart(cosd(22.5), -22.5, 0)[2]
    z[lab .== "tp9"] .= sph2cart(cosd(22.5), -22.5, 0)[3]

    x[lab .== "cp1"] .= sph2cart(cosd(22.5), 0, 112.5)[1]
    y[lab .== "cp1"] .= -cosd(67.5) + sph2cart(cosd(22.5), 0, 112.5)[2]
    z[lab .== "cp1"] .= sph2cart(cosd(22.5), 0, 112.5)[3]
    
    x[lab .== "cp3"] .= sph2cart(cosd(22.5), 0, 135)[1]
    y[lab .== "cp3"] .= -cosd(67.5) + sph2cart(cosd(22.5), 0, 135)[2]
    z[lab .== "cp3"] .= sph2cart(cosd(22.5), 0, 135)[3]
    
    x[lab .== "cp5"] .= sph2cart(cosd(22.5), 0, 157.5)[1]
    y[lab .== "cp5"] .= -cosd(67.5) + sph2cart(cosd(22.5), 0, 157.5)[2]
    z[lab .== "cp5"] .= sph2cart(cosd(22.5), 0, 157.5)[3]
    
    x[lab .== "cp7" .|| lab .== "tp7"] .= sph2cart(cosd(22.5), 0, 180)[1]
    y[lab .== "cp7" .|| lab .== "tp7"] .= -cosd(67.5) + sph2cart(cosd(22.5), 0, 180)[2]
    z[lab .== "cp7" .|| lab .== "tp7"] .= sph2cart(cosd(22.5), 0, 180)[3]
 
    x[lab .== "pz"] .= sph2cart(1.0, 270, 45)[1]
    y[lab .== "pz"] .= sph2cart(1.0, 270, 45)[2]
    z[lab .== "pz"] .= sph2cart(1.0, 270, 45)[3]

    x[lab .== "p2"] .= sph2cart(cosd(45), 0, 67.5)[1]
    y[lab .== "p2"] .= -cosd(45) + sph2cart(cosd(45), 0, 67.5)[2]
    z[lab .== "p2"] .= sph2cart(cosd(45), 0, 67.5)[3]
    
    x[lab .== "p4"] .= sph2cart(cosd(45), 0, 45)[1]
    y[lab .== "p4"] .= -cosd(45) + sph2cart(cosd(45), 0, 45)[2]
    z[lab .== "p4"] .= sph2cart(cosd(45), 0, 45)[3]
    
    x[lab .== "p6"] .= sph2cart(cosd(45), 0, 22.5)[1]
    y[lab .== "p6"] .= -cosd(45) + sph2cart(cosd(45), 0, 22.5)[2]
    z[lab .== "p6"] .= sph2cart(cosd(45), 0, 22.5)[3]
    
    x[lab .== "p8" .|| lab .== "t6"] .= sph2cart(cosd(45), 0, 0)[1]
    y[lab .== "p8" .|| lab .== "t6"] .= -cosd(45) + sph2cart(cosd(45), 0, 0)[2]
    z[lab .== "p8" .|| lab .== "t6"] .= sph2cart(cosd(45), 0, 0)[3]

    x[lab .== "p10"] .= sph2cart(cosd(45), -22.5, 0)[1]
    y[lab .== "p10"] .= -cosd(45) + sph2cart(cosd(45), -22.5, 0)[2]
    z[lab .== "p10"] .= sph2cart(cosd(45), -22.5, 0)[3]

    x[lab .== "p1"] .= sph2cart(cosd(45), 0, 112.5)[1]
    y[lab .== "p1"] .= -cosd(45) + sph2cart(cosd(45), 0, 112.5)[2]
    z[lab .== "p1"] .= sph2cart(cosd(45), 0, 112.5)[3]
    
    x[lab .== "p3"] .= sph2cart(cosd(45), 0, 135)[1]
    y[lab .== "p3"] .= -cosd(45) + sph2cart(cosd(45), 0, 135)[2]
    z[lab .== "p3"] .= sph2cart(cosd(45), 0, 135)[3]
    
    x[lab .== "p5"] .= sph2cart(cosd(45), 0, 157.5)[1]
    y[lab .== "p5"] .= -cosd(45) + sph2cart(cosd(45), 0, 157.5)[2]
    z[lab .== "p5"] .= sph2cart(cosd(45), 0, 157.5)[3]
    
    x[lab .== "p7" .|| lab .== "t5"] .= sph2cart(cosd(45), 0, 180)[1]
    y[lab .== "p7" .|| lab .== "t5"] .= -cosd(45) + sph2cart(cosd(45), 0, 180)[2]
    z[lab .== "p7" .|| lab .== "t5"] .= sph2cart(cosd(45), 0, 180)[3]

    x[lab .== "p9"] .= sph2cart(cosd(45), -22.5, 180)[1]
    y[lab .== "p9"] .= -cosd(45) + sph2cart(cosd(45), -22.5, 180)[2]
    z[lab .== "p9"] .= sph2cart(cosd(45), -22.5, 180)[3]
 
    x[lab .== "poz"] .= sph2cart(1.0, 270, 22.5)[1]
    y[lab .== "poz"] .= sph2cart(1.0, 270, 22.5)[2]
    z[lab .== "poz"] .= sph2cart(1.0, 270, 22.5)[3]

    x[lab .== "po2"] .= sph2cart(cosd(67.5), 0, 67.5)[1]
    y[lab .== "po2"] .= -cosd(67.5) + sph2cart(cosd(67.5), 0, 67.5)[2]
    z[lab .== "po2"] .= sph2cart(cosd(67.5), 0, 67.5)[3]
    
    x[lab .== "po4"] .= sph2cart(cosd(67.5), 0, 45)[1]
    y[lab .== "po4"] .= -cosd(67.5) + sph2cart(cosd(67.5), 0, 45)[2]
    z[lab .== "po4"] .= sph2cart(cosd(67.5), 0, 45)[3]
    
    x[lab .== "po6"] .= sph2cart(cosd(67.5), 0, 22.5)[1]
    y[lab .== "po6"] .= -cosd(67.5) + sph2cart(cosd(67.5), 0, 22.5)[2]
    z[lab .== "po6"] .= sph2cart(cosd(67.5), 0, 22.5)[3]
    
    x[lab .== "po8"] .= sph2cart(cosd(67.5), 0, 0)[1]
    y[lab .== "po8"] .= -cosd(67.5) + sph2cart(cosd(67.5), 0, 0)[2]
    z[lab .== "po8"] .= sph2cart(cosd(67.5), 0, 0)[3]
    
    x[lab .== "po1"] .= sph2cart(cosd(67.5), 0, 112.5)[1]
    y[lab .== "po1"] .= -cosd(67.5) + sph2cart(cosd(67.5), 0, 112.5)[2]
    z[lab .== "po1"] .= sph2cart(cosd(67.5), 0, 112.5)[3]
    
    x[lab .== "po3"] .= sph2cart(cosd(67.5), 0, 135)[1]
    y[lab .== "po3"] .= -cosd(67.5) + sph2cart(cosd(67.5), 0, 135)[2]
    z[lab .== "po3"] .= sph2cart(cosd(67.5), 0, 135)[3]
    
    x[lab .== "po5"] .= sph2cart(cosd(67.5), 0, 157.5)[1]
    y[lab .== "po5"] .= -cosd(67.5) + sph2cart(cosd(67.5), 0, 157.5)[2]
    z[lab .== "po5"] .= sph2cart(cosd(67.5), 0, 157.5)[3]

    x[lab .== "po7"] .= sph2cart(cosd(67.5), 0, 180)[1]
    y[lab .== "po7"] .= -cosd(67.5) + sph2cart(cosd(67.5), 0, 180)[2]
    z[lab .== "po7"] .= sph2cart(cosd(67.5), 0, 180)[3]

    x[lab .== "oz"] .= sph2cart(1.0, 270, 0)[1]
    y[lab .== "oz"] .= sph2cart(1.0, 270, 0)[2]
    z[lab .== "oz"] .= sph2cart(1.0, 270, 0)[3]

    x[lab .== "o2"] .= sph2cart(1.0, 292.5, 0)[1]
    y[lab .== "o2"] .= sph2cart(1.0, 292.5, 0)[2]
    z[lab .== "o2"] .= sph2cart(1.0, 292.5, 0)[3]
        
    x[lab .== "o1"] .= sph2cart(1.0, 247.5, 0)[1]
    y[lab .== "o1"] .= sph2cart(1.0, 247.5, 0)[2]
    z[lab .== "o1"] .= sph2cart(1.0, 247.5, 0)[3]

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

    x = round.(x, digits=3)
    y = round.(y, digits=3)
    z = round.(z, digits=3)

    locs_new[:, :loc_x] = x
    locs_new[:, :loc_y] = y
    locs_new[:, :loc_z] = z

    locs_cart2sph!(locs_new) 
    locs_sph2pol!(locs_new) 

    f_labels = lowercase.(locs[!, :labels])
    no_match = setdiff(f_labels, e_labels)
    length(no_match) > 0 && _warn("Location$(_pl(no_match)): $(uppercase.(no_match)) could not be generated.")

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

"""
    locs_generate(obj)

Generate spherical coordinates according to 10/5 system.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function locs_generate(obj::NeuroAnalyzer.NEURO)

    obj_new = deepcopy(obj)
    locs = locs_generate(obj.locs)
    obj_new.locs = locs

    # add entry to :history field
    push!(obj_new.history, "locs_generate(OBJ)")

    return obj_new

end

"""
    locs_generate!(obj)

Generate spherical coordinates according to 10/5 system.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
"""
function locs_generate!(obj::NeuroAnalyzer.NEURO)

    obj_new = locs_generate(obj)
    obj.history = obj_new.history
    obj.locs = obj_new.locs

    return nothing

end
