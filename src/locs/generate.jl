export locs_generate

"""
    locs_generate()

Generate locs according to 10/5 system.
"""
function locs_generate()

    _wip()

    lab = ["Fpz", "Fz", "Cz", "Pz", "Oz", "Fp1", "Fp2", "F7", "F3", "F4", "F8", "T3", "T7", "C3", "C4", "T4", "T5", "T8", "P3", "P4", "T6", "O1", "O2", "A1", "A2", "M1", "M2", "EMG1", "EMG2", "EOG1", "EOG2"]

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

end