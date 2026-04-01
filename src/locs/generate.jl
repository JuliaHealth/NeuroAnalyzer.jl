export locs_generate
export locs_generate!

"""
    locs_generate(locs)

Generate Cartesian and spherical electrode coordinates according to the 10/10 international system.

Labels are matched case-insensitively. Bipolar labels (e.g. `Fp1-Fz`) are stripped to their first component before matching. Unrecognized labels are reported in a warning but left with their original (zero) coordinates.

# Arguments

- `locs::DataFrame`: channel location data

# Returns

- `DataFrame`: modified channel location data
"""
function locs_generate(locs::DataFrame)::DataFrame
    # create new dataset
    locs_new = deepcopy(locs)

    # work with lowercase labels for case-insensitive matching
    lab = lowercase.(locs[!, :label])

    # strip bipolar reference suffixes: "Fp1-Fz" → "fp1"
    m = match.(r"(.+)\-(.+)", lab)
    for idx in eachindex(lab)
        m[idx] !== nothing && (lab[idx] = m[idx].captures[1])
    end

    # strip concatenated reference labels: "Fp1Fz" → "fp1"
    m = match.(r"([a-z]+[0-9]+[0-9]?)([a-z]+[0-9]+)", lab)
    for idx in eachindex(lab)
        m[idx] !== nothing && (lab[idx] = m[idx].captures[1])
    end

    n = length(lab)
    x = zeros(n)
    y = zeros(n)
    z = zeros(n)

    # complete list of electrode labels handled by this function
    # used at the end to warn about unrecognized labels
    e_labels = [
        "cz", "c2", "c4", "c6", "t4", "t8", "t10",
        "c1", "c3", "c5", "t3", "t7", "t9",
        "fcz", "fc2", "fc4", "fc6", "fc8", "ft8", "fc10", "ft10",
        "fc1", "fc3", "fc5", "fc7", "ft7", "fc9", "ft9",
        "fz", "f2", "f4", "f6", "f8", "f10",
        "f1", "f3", "f5", "f7", "f9",
        "afz", "af2", "af4", "af6", "af1", "af3", "af7",
        "fpz", "fp2", "fp1",
        "cpz", "cp2", "cp4", "cp6", "cp8", "tp8", "tp10",
        "tp9", "cp1", "cp3", "cp5", "cp7", "tp7",
        "pz", "p2", "p4", "p6", "p8", "p10", "t6",
        "p1", "p3", "p5", "p7", "p9", "t5",
        "poz", "po2", "po4", "po6", "po8",
        "po1", "po3", "po5", "po7",
        "oz", "o2", "o1",
        "a1", "a2", "m1", "m2",
        "emg1", "emg2", "eog1", "eog2",
        "veog1", "veog2", "heog1", "heog2",
        "veog", "heog", "reog", "leog",
    ]

    # ------------------------------------------------------------------ #
    # Electrode coordinate assignment                                    #
    #                                                                    #
    # Convention: sph2cart(radius, azimuth_deg, elevation_deg)           #
    #   - radius 1.0 for mid-sagittal and equatorial electrodes          #
    #   - radius cosd(θ) for electrodes at latitude θ above the equator  #
    #   - y-offset ±cosd(θ) shifts the ring to the correct AP position   #
    #   - phi=0 → right hemisphere; phi=180 → left hemisphere            #
    # ------------------------------------------------------------------ #

    # central row (C strip): equatorial plane, z = sin(elevation)
    # Cz sits at the apex of the unit sphere (elevation = 90°)
    let (xi, yi, zi) = sph2cart(1.0, 0, 90)
        x[lab .== "cz"] .= xi
        y[lab .== "cz"] .= yi
        z[lab .== "cz"] .= zi
    end
    let (xi, yi, zi) = sph2cart(1.0, 0, 67.5)
        x[lab .== "c2"] .= xi
        y[lab .== "c2"] .= yi
        z[lab .== "c2"] .= zi
    end
    let (xi, yi, zi) = sph2cart(1.0, 0, 45)
        x[lab .== "c4"] .= xi
        y[lab .== "c4"] .= yi
        z[lab .== "c4"] .= zi
    end
    let (xi, yi, zi) = sph2cart(1.0, 0, 22.5)
        x[lab .== "c6"] .= xi
        y[lab .== "c6"] .= yi
        z[lab .== "c6"] .= zi
    end
    let (xi, yi, zi) = sph2cart(1.0, 0, 0)
        mask = lab .== "t4" .|| lab .== "t8"
        x[mask] .= xi
        y[mask] .= yi
        z[mask] .= zi
    end
    let (xi, yi, zi) = sph2cart(1.0, -22.5, 0)
        x[lab .== "t10"] .= xi
        y[lab .== "t10"] .= yi
        z[lab .== "t10"] .= zi
    end
    let (xi, yi, zi) = sph2cart(1.0, 0, 112.5)
        x[lab .== "c1"] .= xi
        y[lab .== "c1"] .= yi
        z[lab .== "c1"] .= zi
    end
    let (xi, yi, zi) = sph2cart(1.0, 0, 135)
        x[lab .== "c3"] .= xi
        y[lab .== "c3"] .= yi
        z[lab .== "c3"] .= zi
    end
    let (xi, yi, zi) = sph2cart(1.0, 0, 157.5)
        x[lab .== "c5"] .= xi
        y[lab .== "c5"] .= yi
        z[lab .== "c5"] .= zi
    end
    let (xi, yi, zi) = sph2cart(1.0, 0, 180)
        mask = lab .== "t3" .|| lab .== "t7"
        x[mask] .= xi
        y[mask] .= yi
        z[mask] .= zi
    end
    let (xi, yi, zi) = sph2cart(1.0, -22.5, 180)
        x[lab .== "t9"] .= xi
        y[lab .== "t9"] .= yi
        z[lab .== "t9"] .= zi
    end

    # FC row: latitude 22.5° anterior; radius = cosd(22.5), y-offset = +cosd(22.5)
    let r = cosd(22.5), off = cosd(22.5)
        let (xi, yi, zi) = sph2cart(r, 90, 67.5)
            x[lab .== "fcz"] .= xi
            y[lab .== "fcz"] .= yi
            z[lab .== "fcz"] .= zi
        end
        for (lbl, phi) in [("fc2", 67.5), ("fc4", 45), ("fc6", 22.5)]
            let (xi, yi, zi) = sph2cart(r, 0, phi)
                x[lab .== lbl] .= xi
                y[lab .== lbl] .= off .+ yi
                z[lab .== lbl] .= zi
            end
        end
        let (xi, yi, zi) = sph2cart(r, 0, 0)
            mask = lab .== "fc8" .|| lab .== "ft8"
            x[mask] .= xi
            y[mask] .= off .+ yi
            z[mask] .= zi
        end
        let (xi, yi, zi) = sph2cart(r, -22.5, 0)
            mask = lab .== "fc10" .|| lab .== "ft10"
            x[mask] .= xi
            y[mask] .= off .+ yi
            z[mask] .= zi
        end
        for (lbl, phi) in [("fc1", 112.5), ("fc3", 135), ("fc5", 157.5)]
            let (xi, yi, zi) = sph2cart(r, 0, phi)
                x[lab .== lbl] .= xi
                y[lab .== lbl] .= off .+ yi
                z[lab .== lbl] .= zi
            end
        end
        let (xi, yi, zi) = sph2cart(r, 0, 180)
            mask = lab .== "fc7" .|| lab .== "ft7"
            x[mask] .= xi
            y[mask] .= off .+ yi
            z[mask] .= zi
        end
        let (xi, yi, zi) = sph2cart(r, -22.5, 180)
            mask = lab .== "fc9" .|| lab .== "ft9"
            x[mask] .= xi
            y[mask] .= off .+ yi
            z[mask] .= zi
        end
    end

    # F row: latitude 45° anterior; radius = cosd(45), y-offset = +cosd(45)
    let r = cosd(45), off = cosd(45)
        let (xi, yi, zi) = sph2cart(1.0, 90, 45)
            x[lab .== "fz"] .= xi
            y[lab .== "fz"] .= yi
            z[lab .== "fz"] .= zi
        end
        for (lbl, phi) in [("f2", 67.5), ("f4", 45), ("f6", 22.5), ("f8", 0)]
            let (xi, yi, zi) = sph2cart(r, 0, phi)
                x[lab .== lbl] .= xi
                y[lab .== lbl] .= off .+ yi
                z[lab .== lbl] .= zi
            end
        end
        let (xi, yi, zi) = sph2cart(r, -22.5, 0)
            x[lab .== "f10"] .= xi
            y[lab .== "f10"] .= off .+ yi
            z[lab .== "f10"] .= zi
        end
        for (lbl, phi) in [("f1", 112.5), ("f3", 135), ("f5", 157.5), ("f7", 180)]
            let (xi, yi, zi) = sph2cart(r, 0, phi)
                x[lab .== lbl] .= xi
                y[lab .== lbl] .= off .+ yi
                z[lab .== lbl] .= zi
            end
        end
        let (xi, yi, zi) = sph2cart(r, -22.5, 180)
            x[lab .== "f9"] .= xi
            y[lab .== "f9"] .= off .+ yi
            z[lab .== "f9"] .= zi
        end
    end

    # AF row: latitude 67.5° anterior; radius = cosd(67.5), y-offset = +cosd(67.5)
    let r = cosd(67.5), off = cosd(67.5)
        let (xi, yi, zi) = sph2cart(1.0, 90, 67.5)
            x[lab .== "afz"] .= xi
            y[lab .== "afz"] .= yi
            z[lab .== "afz"] .= zi
        end
        for (lbl, phi) in [
            ("af2", 67.5), ("af4", 45), ("af6", 22.5),
            ("af1", 112.5), ("af3", 135), ("af7", 157.5),
        ]
            let (xi, yi, zi) = sph2cart(r, 0, phi)
                x[lab .== lbl] .= xi
                y[lab .== lbl] .= off .+ yi
                z[lab .== lbl] .= zi
            end
        end
    end

    # Fp row: frontopolar, on the unit sphere
    let (xi, yi, zi) = sph2cart(1.0, 90, 0)
        x[lab .== "fpz"] .= xi
        y[lab .== "fpz"] .= yi
        z[lab .== "fpz"] .= zi
    end
    let (xi, yi, zi) = sph2cart(1.0, 67.5, 0)
        x[lab .== "fp2"] .= xi
        y[lab .== "fp2"] .= yi
        z[lab .== "fp2"] .= zi
    end
    let (xi, yi, zi) = sph2cart(1.0, 112.5, 0)
        x[lab .== "fp1"] .= xi
        y[lab .== "fp1"] .= yi
        z[lab .== "fp1"] .= zi
    end

    # CP row: latitude 22.5° posterior; y-offset = −cosd(67.5)
    let r = cosd(22.5), off = -cosd(67.5)
        let (xi, yi, zi) = sph2cart(1.0, 270, 67.5)
            x[lab .== "cpz"] .= xi
            y[lab .== "cpz"] .= yi
            z[lab .== "cpz"] .= zi
        end
        for (lbl, phi) in [("cp2", 67.5), ("cp4", 45), ("cp6", 22.5)]
            let (xi, yi, zi) = sph2cart(r, 0, phi)
                x[lab .== lbl] .= xi
                y[lab .== lbl] .= off .+ yi
                z[lab .== lbl] .= zi
            end
        end
        let (xi, yi, zi) = sph2cart(r, 0, 0)
            mask = lab .== "cp8" .|| lab .== "tp8"
            x[mask] .= xi
            y[mask] .= off .+ yi
            z[mask] .= zi
        end
        let (xi, yi, zi) = sph2cart(r, -22.5, 0)
            x[lab .== "tp10"] .= xi
            y[lab .== "tp10"] .= off .+ yi
            z[lab .== "tp10"] .= zi
        end
        let (xi, yi, zi) = sph2cart(r, -22.5, 180)
            x[lab .== "tp9"] .= xi
            y[lab .== "tp9"] .= off .+ yi
            z[lab .== "tp9"] .= zi
        end
        for (lbl, phi) in [("cp1", 112.5), ("cp3", 135), ("cp5", 157.5)]
            let (xi, yi, zi) = sph2cart(r, 0, phi)
                x[lab .== lbl] .= xi
                y[lab .== lbl] .= off .+ yi
                z[lab .== lbl] .= zi
            end
        end
        let (xi, yi, zi) = sph2cart(r, 0, 180)
            mask = lab .== "cp7" .|| lab .== "tp7"
            x[mask] .= xi
            y[mask] .= off .+ yi
            z[mask] .= zi
        end
    end

    # P row: latitude 45° posterior; y-offset = −cosd(45)
    let r = cosd(45), off = -cosd(45)
        let (xi, yi, zi) = sph2cart(1.0, 270, 45)
            x[lab .== "pz"] .= xi
            y[lab .== "pz"] .= yi
            z[lab .== "pz"] .= zi
        end
        for (lbl, phi) in [("p2", 67.5), ("p4", 45), ("p6", 22.5)]
            let (xi, yi, zi) = sph2cart(r, 0, phi)
                x[lab .== lbl] .= xi
                y[lab .== lbl] .= off .+ yi
                z[lab .== lbl] .= zi
            end
        end
        let (xi, yi, zi) = sph2cart(r, 0, 0)
            mask = lab .== "p8" .|| lab .== "t6"
            x[mask] .= xi
            y[mask] .= off .+ yi
            z[mask] .= zi
        end
        let (xi, yi, zi) = sph2cart(r, -22.5, 0)
            x[lab .== "p10"] .= xi
            y[lab .== "p10"] .= off .+ yi
            z[lab .== "p10"] .= zi
        end
        for (lbl, phi) in [("p1", 112.5), ("p3", 135), ("p5", 157.5)]
            let (xi, yi, zi) = sph2cart(r, 0, phi)
                x[lab .== lbl] .= xi
                y[lab .== lbl] .= off .+ yi
                z[lab .== lbl] .= zi
            end
        end
        let (xi, yi, zi) = sph2cart(r, 0, 180)
            mask = lab .== "p7" .|| lab .== "t5"
            x[mask] .= xi
            y[mask] .= off .+ yi
            z[mask] .= zi
        end
        let (xi, yi, zi) = sph2cart(r, -22.5, 180)
            x[lab .== "p9"] .= xi
            y[lab .== "p9"] .= off .+ yi
            z[lab .== "p9"] .= zi
        end
    end

    # PO row: latitude 67.5° posterior; y-offset = −cosd(67.5)
    let r = cosd(67.5), off = -cosd(67.5)
        let (xi, yi, zi) = sph2cart(1.0, 270, 22.5)
            x[lab .== "poz"] .= xi
            y[lab .== "poz"] .= yi
            z[lab .== "poz"] .= zi
        end
        for (lbl, phi) in [
            ("po2", 67.5), ("po4", 45), ("po6", 22.5), ("po8", 0),
            ("po1", 112.5), ("po3", 135), ("po5", 157.5), ("po7", 180),
        ]
            let (xi, yi, zi) = sph2cart(r, 0, phi)
                x[lab .== lbl] .= xi
                y[lab .== lbl] .= off .+ yi
                z[lab .== lbl] .= zi
            end
        end
    end

    # O row: occipital, on the unit sphere
    for (lbl, az, phi) in [("oz", 270, 0), ("o2", 292.5, 0), ("o1", 247.5, 0)]
        let (xi, yi, zi) = sph2cart(1.0, az, phi)
            x[lab .== lbl] .= xi
            y[lab .== lbl] .= yi
            z[lab .== lbl] .= zi
        end
    end

    # ------------------------------------------------------------------ #
    # non-scalp electrodes: hardcoded Cartesian coordinates              #
    # these are approximations based on typical head-surface placement   #
    # ------------------------------------------------------------------ #

    # ear/mastoid references
    x[lab .== "a1"] .= -0.92
    y[lab .== "a1"] .= -0.23
    z[lab .== "a1"] .= -0.55
    x[lab .== "a2"] .= 0.92
    y[lab .== "a2"] .= -0.23
    z[lab .== "a2"] .= -0.55
    x[lab .== "m1"] .= -0.94
    y[lab .== "m1"] .= -0.1
    z[lab .== "m1"] .= -0.3
    x[lab .== "m2"] .= 0.94
    y[lab .== "m2"] .= -0.1
    z[lab .== "m2"] .= -0.3

    # EMG electrodes (chin/neck placement)
    x[lab .== "emg1"] .= -0.7
    y[lab .== "emg1"] .= 0.7
    z[lab .== "emg1"] .= -1.1
    x[lab .== "emg2"] .= 0.7
    y[lab .== "emg2"] .= 0.7
    z[lab .== "emg2"] .= -1.1

    # EOG electrodes: vertical (VEOG) and horizontal (HEOG) eye movements
    for lbl in ("eog1", "veog1")
        x[lab .== lbl] .= -0.87
        y[lab .== lbl] .= 0.51
        z[lab .== lbl] .= -0.37
    end
    for lbl in ("eog2", "veog2", "veog")
        x[lab .== lbl] .= 0.87
        y[lab .== lbl] .= 0.51
        z[lab .== lbl] .= -0.37
    end
    for lbl in ("heog1", "leog")
        x[lab .== lbl] .= -0.64
        y[lab .== lbl] .= 0.77
        z[lab .== lbl] .= -0.04
    end
    for lbl in ("heog2", "reog", "heog")
        x[lab .== lbl] .= 0.64
        y[lab .== lbl] .= 0.77
        z[lab .== lbl] .= -0.04
    end

    # round to 3 decimal places to avoid floating-point noise
    x = round.(x; digits = 3)
    y = round.(y; digits = 3)
    z = round.(z; digits = 3)

    # write the Cartesian coordinates back into the copied DataFrame
    locs_new[:, :loc_x] = x
    locs_new[:, :loc_y] = y
    locs_new[:, :loc_z] = z

    # derive spherical and polar coordinates from the updated Cartesian values
    locs_cart2sph!(locs_new)
    locs_sph2pol!(locs_new)

    # warn about any labels that had no match in the electrode lookup table
    no_match = setdiff(lab, e_labels)
    if !isempty(no_match)
        _warn(
            "Location$(_pl(length(no_match))): $(uppercase.(no_match)) could not be generated.",
        )
    end

    return locs_new
end

"""
    locs_generate!(locs)

Generate spherical coordinates according to the 10/10 system, modifying `locs` in-place.

# Arguments

- `locs::DataFrame`: channel location data

# Returns

- `Nothing`
"""
function locs_generate!(locs::DataFrame)::Nothing
    locs_tmp = locs_generate(locs)

    # copy all location columns back into the original DataFrame.
    for col in (
        :loc_radius, :loc_theta, :loc_x, :loc_y, :loc_z,
        :loc_radius_sph, :loc_theta_sph, :loc_phi_sph,
    )
        locs[:, col] = locs_tmp[:, col]
    end

    return nothing
end

"""
    locs_generate(obj)

Generate spherical coordinates according to the 10/10 system for all channels in a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function locs_generate(obj::NeuroAnalyzer.NEURO)::NeuroAnalyzer.NEURO

    # create new dataset
    obj_new = deepcopy(obj)

    locs = locs_generate(obj.locs)
    obj_new.locs = locs
    push!(obj_new.history, "locs_generate(OBJ)")

    return obj_new
end

"""
    locs_generate!(obj)

Generate spherical coordinates according to the 10/10 system, modifying `locs` in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object

# Returns

- `Nothing`
"""
function locs_generate!(obj::NeuroAnalyzer.NEURO)::Nothing
    obj_new = locs_generate(obj)
    obj.history = obj_new.history
    obj.locs = obj_new.locs

    return nothing
end
