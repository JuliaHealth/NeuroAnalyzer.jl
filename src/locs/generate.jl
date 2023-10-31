export locs_generate

function locs_generate()

    lab = ["Fpz", "Fz", "Cz", "Pz", "Oz", "Fp1", "Fp2", "F7", "F3", "F4", "F8", "T3", "T7", "C3", "C4", "T4", "T5", "T8", "P3", "P4", "T6", "O1", "O2"]

    x = zeros(length(lab))
    y = zeros(length(lab))
    z = zeros(length(lab))

    x[lab .== "Cz"] .= 0.0
    y[lab .== "Cz"] .= 0.0
    z[lab .== "Cz"] .= 1.0

    x[lab .== "Fpz"] .= 0.0
    y[lab .== "Fpz"] .= 1.0
    z[lab .== "Fpz"] .= 0.0

    x[lab .== "Oz"] .= 0.0
    y[lab .== "Oz"] .= -1.0
    z[lab .== "Oz"] .= 0.0

    x[lab .== "T3" .|| lab .== "T7"] .= -1.0
    y[lab .== "T3" .|| lab .== "T7"] .= 0.0
    z[lab .== "T3" .|| lab .== "T7"] .= 0.0

    x[lab .== "T4" .|| lab .== "T8"] .= -1.0
    y[lab .== "T4" .|| lab .== "T8"] .= 0.0
    z[lab .== "T4" .|| lab .== "T8"] .= 0.0

end