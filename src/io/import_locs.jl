export import_locs
export import_locs_ced
export import_locs_locs
export import_locs_csd
export import_locs_elc
export import_locs_tsv
export import_locs_sfp
export import_locs_geo
export import_locs_mat

"""
    import_locs(file_name)

Load channel locations. Supported formats:
- CED
- ELC
- LOCS
- TSV
- SFP
- CSD
- GEO
- MAT

This is a meta-function that triggers appropriate `import_locs_*()` function. File format is detected based on file extension.

# Arguments

- `file_name::String`: name of the file to load

# Returns

- `locs::DataFrame`
"""
function import_locs(file_name::String)

    isfile(file_name) || throw(ArgumentError("File $file_name cannot be loaded."))

    _info("Send standard location for your channels to adam.wysokinski@neuroanalyzer.org")
    _info("Nose direction is set at '+Y'.")

    if splitext(file_name)[2] == ".ced"
        locs = import_locs_ced(file_name)
    elseif splitext(file_name)[2] == ".elc"
        locs = import_locs_elc(file_name)
    elseif splitext(file_name)[2] == ".locs"
        locs = import_locs_locs(file_name)
    elseif splitext(file_name)[2] == ".tsv"
        locs = import_locs_tsv(file_name)
    elseif splitext(file_name)[2] == ".sfp"
        locs = import_locs_sfp(file_name)
    elseif splitext(file_name)[2] == ".csd"
        locs = import_locs_csd(file_name)
    elseif splitext(file_name)[2] == ".geo"
        locs = import_locs_geo(file_name)
    elseif splitext(file_name)[2] == ".mat"
        locs = import_locs_mat(file_name)
    else
        throw(ArgumentError("Unknown file format."))
    end

    return locs

end

"""
    import_locs_ced(file_name)

Load channel locations from CED file.

# Arguments

- `file_name::String`

# Returns

- `locs::DataFrame`
"""
function import_locs_ced(file_name::String)

    isfile(file_name) || throw(ArgumentError("$file_name not found."))
    splitext(file_name)[2] == ".ced" || throw(ArgumentError("Not a CED file."))

    locs = CSV.read(file_name, delim="\t", DataFrame)

    colnames = lowercase.(names(locs))
    DataFrames.rename!(locs, Symbol.(colnames))

    clabels = lstrip.(locs[!, "labels"])

    x = zeros(length(clabels))
    y = zeros(length(clabels))
    z = zeros(length(clabels))
    radius = zeros(length(clabels))
    theta = zeros(length(clabels))
    radius_sph = zeros(length(clabels))
    theta_sph = zeros(length(clabels))
    phi_sph = zeros(length(clabels))

    "x" in colnames && (x = Float64.(locs[!, "x"]))
    "y" in colnames && (y = Float64.(locs[!, "y"]))
    "z" in colnames && (z = Float64.(locs[!, "z"]))
    "theta" in colnames && (theta = Float64.(locs[!, "theta"]))
    "radius" in colnames && (radius = Float64.(locs[!, "radius"]))
    "sph_radius" in colnames && (radius_sph = Float64.(locs[!, "sph_radius"]))
    "sph_theta" in colnames && (theta_sph = Float64.(locs[!, "sph_theta"]))
    "sph_phi" in colnames && (phi_sph = Float64.(locs[!, "sph_phi"]))

    locs = DataFrame(:channel=>1:length(clabels), :labels=>clabels, :loc_theta=>theta, :loc_radius=>radius, :loc_x=>x, :loc_y=>y, :loc_z=>z, :loc_radius_sph=>radius_sph, :loc_theta_sph=>theta_sph, :loc_phi_sph=>phi_sph)

    locs = _round_locs(locs)

    locs_swapxy!(locs)
    locs_flipx!(locs, planar=true, spherical=false)

    return locs

end

"""
    import_locs_locs(file_name)

Load channel locations from LOCS file.

# Arguments

- `file_name::String`

# Returns

- `locs::DataFrame`
"""
function import_locs_locs(file_name::String)

    isfile(file_name) || throw(ArgumentError("$file_name not found."))
    splitext(file_name)[2] == ".locs" || throw(ArgumentError("Not a LOCS file."))

    locs = CSV.read(file_name, header=false, delim="\t", DataFrame)

    DataFrames.rename!(locs, [:number, :theta, :radius, :labels])
    clabels = lstrip.(locs[!, "labels"])

    x = zeros(length(clabels))
    y = zeros(length(clabels))
    z = zeros(length(clabels))
    radius = zeros(length(clabels))
    theta = zeros(length(clabels))
    radius_sph = zeros(length(clabels))
    theta_sph = zeros(length(clabels))
    phi_sph = zeros(length(clabels))

    theta = Float64.(locs[!, "theta"])
    radius = Float64.(locs[!, "radius"])

    locs = DataFrame(:channel=>1:length(clabels), :labels=>clabels, :loc_theta=>theta, :loc_radius=>radius, :loc_x=>x, :loc_y=>y, :loc_z=>z, :loc_radius_sph=>radius_sph, :loc_theta_sph=>theta_sph, :loc_phi_sph=>phi_sph)

    locs = _round_locs(locs)

    locs_swapxy!(locs)
    locs_flipx!(locs, planar=true, spherical=false)

    locs[!, :loc_phi_sph] = zeros(nrow(locs))

    return locs

end

"""
    import_locs_elc(file_name)

Load channel locations from ELC file.

# Arguments

- `file_name::String`

# Returns

- `locs::DataFrame`
"""
function import_locs_elc(file_name::String)

    isfile(file_name) || throw(ArgumentError("$file_name not found."))
    splitext(file_name)[2] == ".elc" || throw(ArgumentError("Not a ELC file."))

    f = open(file_name, "r")
    elc_file = readlines(f)
    close(f)

    locs_n = 0
    locs_l = 0
    for idx in eachindex(elc_file)
        if occursin("NumberPositions", elc_file[idx]) == true
            locs_n = parse(Int64, replace(elc_file[idx], "NumberPositions=" => ""))
            locs_l = idx + 2
        end
    end
    clabels = repeat([""], locs_n)

    x = zeros(length(clabels))
    y = zeros(length(clabels))
    z = zeros(length(clabels))
    theta = zeros(length(clabels))
    radius = zeros(length(clabels))
    radius_sph = zeros(length(clabels))
    theta_sph = zeros(length(clabels))
    phi_sph = zeros(length(clabels))

    idx2 = 1
    for idx1 in locs_l:(locs_l + locs_n - 1)
        l = elc_file[idx1]
        l[1] == ' ' && (l = l[2:end])
        x[idx2], y[idx2], z[idx2] = parse.(Float64, split(l, ' '))
        idx2 += 1
    end
    idx2 = 1
    for idx1 in (locs_l + 1 + locs_n):(locs_l + (2 * locs_n))
        clabels[idx2] = elc_file[idx1]
        idx2 += 1
    end
    x = normalize_minmax(x)
    y = normalize_minmax(y)
    z = normalize_minmax(z)

    locs = DataFrame(:channel=>1:length(clabels), :labels=>clabels, :loc_theta=>theta, :loc_radius=>radius, :loc_x=>x, :loc_y=>y, :loc_z=>z, :loc_radius_sph=>radius_sph, :loc_theta_sph=>theta_sph, :loc_phi_sph=>phi_sph)

    locs = _round_locs(locs)

    locs_cart2sph!(locs)
    locs_cart2pol!(locs)

    return locs

end

"""
    import_locs_tsv(file_name)

Load channel locations from TSV file.

# Arguments

- `file_name::String`

# Returns

- `locs::DataFrame`
"""
function import_locs_tsv(file_name::String)

    isfile(file_name) || throw(ArgumentError("$file_name not found."))
    splitext(file_name)[2] == ".tsv" || throw(ArgumentError("Not a TSV file."))

    locs = CSV.read(file_name, header=true, delim="\t", ignorerepeated=true, DataFrame)

    colnames = lowercase.(names(locs))
    DataFrames.rename!(locs, Symbol.(colnames))

    "labels" in colnames && (clabels = lstrip.(locs[!, "labels"]))
    "label" in colnames && (clabels = lstrip.(locs[!, "label"]))
    "site" in colnames && (clabels = lstrip.(locs[!, "site"]))

    x = zeros(length(clabels))
    y = zeros(length(clabels))
    z = zeros(length(clabels))
    radius = zeros(length(clabels))
    theta = zeros(length(clabels))
    radius_sph = zeros(length(clabels))
    theta_sph = zeros(length(clabels))
    phi_sph = zeros(length(clabels))
    
    "x" in colnames && (x = Float64.(locs[!, "x"]))
    "y" in colnames && (y = Float64.(locs[!, "y"]))
    "z" in colnames && (z = Float64.(locs[!, "z"]))
    "theta" in colnames && (theta = Float64.(locs[!, "theta"]))
    "radius" in colnames && (radius = Float64.(locs[!, "radius"]))
    "radius" in colnames && (radius_sph = Float64.(locs[!, "radius"]))
    "radius_sph" in colnames && (radius_sph = locs[!, "radius_sph"])
    "theta" in colnames && (theta_sph = Float64.(locs[!, "theta"]))
    "theta_sph" in colnames && (theta_sph = locs[!, "theta_sph"])
    "phi" in colnames && (phi_sph = locs[!, "phi"])
    "phi_sph" in colnames && (phi_sph = locs[!, "phi_sph"])

    locs = DataFrame(:channel=>1:length(clabels), :labels=>clabels, :loc_theta=>theta, :loc_radius=>radius, :loc_x=>x, :loc_y=>y, :loc_z=>z, :loc_radius_sph=>radius_sph, :loc_theta_sph=>theta_sph, :loc_phi_sph=>phi_sph)

    locs = _round_locs(locs)

    locs_cart2sph!(locs)
    locs_cart2pol!(locs)

    return locs

end

"""
    import_locs_sfp(file_name)

Load channel locations from SFP file.

# Arguments

- `file_name::String`

# Returns

- `locs::DataFrame`
"""
function import_locs_sfp(file_name::String)

    isfile(file_name) || throw(ArgumentError("$file_name not found."))
    splitext(file_name)[2] == ".sfp" || throw(ArgumentError("Not a SFP file."))
    
    locs = CSV.read(file_name, header=false, DataFrame)
    if size(locs, 2) != 4
        locs = CSV.read(file_name, header=false, delim="/t", ignorerepeated=true, DataFrame)
    end
    if size(locs, 2) != 4
        locs = CSV.read(file_name, header=false, delim=" ", ignorerepeated=true, DataFrame)
    end
    if size(locs, 2) != 4
        throw(ArgumentError("File $file_name cannot be opened, check delimeters."))
    end

    DataFrames.rename!(locs, [:label, :x, :y, :z])

    clabels = lstrip.(locs[!, "label"])

    x = Float64.(locs[!, :x])
    y = Float64.(locs[!, :y])
    z = Float64.(locs[!, :z])

    # x, y, z positions must be within -1..+1
    t = x[1]
    x, y, z = _locnorm(x, y, z)
    t -= x[1]
    # sometimes positions are shifted along x-axis, remove the shift
    x .+= abs(t)

    radius = zeros(length(clabels))
    theta = zeros(length(clabels))
    radius_sph = zeros(length(clabels))
    theta_sph = zeros(length(clabels))
    phi_sph = zeros(length(clabels))

    locs = DataFrame(:channel=>1:length(clabels), :labels=>clabels, :loc_theta=>theta, :loc_radius=>radius, :loc_x=>x, :loc_y=>y, :loc_z=>z, :loc_radius_sph=>radius_sph, :loc_theta_sph=>theta_sph, :loc_phi_sph=>phi_sph)

    locs = _round_locs(locs)

    locs_cart2sph!(locs)
    locs_cart2pol!(locs)

    return locs

end

"""
    import_locs_csd(file_name)

Load channel locations from CSD file.

# Arguments

- `file_name::String`

# Returns

- `locs::DataFrame`
"""
function import_locs_csd(file_name::String)

    isfile(file_name) || throw(ArgumentError("$file_name not found."))
    splitext(file_name)[2] == ".csd" || throw(ArgumentError("Not a CSD file."))

    locs = CSV.read(file_name, skipto=3, delim=' ', header=false, ignorerepeated=true, DataFrame)

    DataFrames.rename!(locs, [:labels, :theta_sph, :phi_sph, :radius_sph, :x, :y, :z, :surface])
    clabels = lstrip.(locs[!, "labels"])

    x = Float64.(locs[!, "x"])
    y = Float64.(locs[!, "y"])
    z = Float64.(locs[!, "z"])
    radius_sph = Float64.(locs[!, "radius_sph"])
    theta_sph = Float64.(locs[!, "theta_sph"])
    phi_sph = Float64.(locs[!, "phi_sph"])

    radius = zeros(length(x))
    theta = zeros(length(y))
    for idx in 1:length(x)
        radius[idx], theta[idx] = sph2pol(radius_sph[idx], theta_sph[idx], phi_sph[idx])
    end

    locs = DataFrame(:channel=>1:length(clabels), :labels=>clabels, :loc_theta=>theta, :loc_radius=>radius, :loc_x=>x, :loc_y=>y, :loc_z=>z, :loc_radius_sph=>radius_sph, :loc_theta_sph=>theta_sph, :loc_phi_sph=>phi_sph)

    locs = _round_locs(locs)

    return locs

end

"""
    import_locs_geo(file_name)

Load channel locations from GEO file.

# Arguments

- `file_name::String`

# Returns

- `locs::DataFrame`
"""
function import_locs_geo(file_name::String)

    isfile(file_name) || throw(ArgumentError("$file_name not found."))
    splitext(file_name)[2] == ".geo" || throw(ArgumentError("Not a GEO file."))

    f = open(file_name, "r")
    locs = readlines(f)
    close(f)

    l1 = 0
    l2 = 0
    for idx in 1:length(locs)
        locs[idx] == "View\"\"{" && (l1 = idx + 1)
        locs[idx] == "};" && (l2 = idx - 1)
    end
    locs = locs[l1+1:2:l2]

    clabels = repeat([""], length(locs))
    x = zeros(length(locs))
    y = zeros(length(locs))
    z = zeros(length(locs))

    p = r"(.+)(\(.+\)){(.+)}"
    for idx in 1:length(clabels)
        m = match(p, locs[idx])
        clabels[idx] = replace(m[3], "\"" => "")
        tmp = replace(m[2], "(" => "")
        tmp = replace(tmp, ")" => "")
        x[idx], y[idx], z[idx], = parse.(Float64, split(tmp, ", "))
    end

    x, y, z = _locnorm(x, y, z)

    # center x at 0
    x_adj = x[findfirst(isequal("Cz"), clabels)]
    x .-= x_adj
    x, y, z = _locnorm(x, y, z)

    radius = zeros(length(clabels))
    theta = zeros(length(clabels))
    radius_sph = zeros(length(clabels))
    theta_sph = zeros(length(clabels))
    phi_sph = zeros(length(clabels))

    locs = DataFrame(:channel=>1:length(clabels), :labels=>clabels, :loc_theta=>theta, :loc_radius=>radius, :loc_x=>x, :loc_y=>y, :loc_z=>z, :loc_radius_sph=>radius_sph, :loc_theta_sph=>theta_sph, :loc_phi_sph=>phi_sph)

    locs = _round_locs(locs)

    locs = locs_cart2sph(locs)
    locs = locs_cart2pol(locs)

    return locs

end

"""
    import_locs_mat(file_name)

Load channel locations from MAT file.

# Arguments

- `file_name::String`

# Returns

- `locs::DataFrame`
"""
function import_locs_mat(file_name::String)

    isfile(file_name) || throw(ArgumentError("$file_name not found."))
    splitext(file_name)[2] == ".mat" || throw(ArgumentError("Not a MAT file."))

    dataset = matread(file_name)
    x = dataset["Cpos"][1, :]
    y = dataset["Cpos"][2, :]
    r = dataset["Rxy"]
    ch_n = length(x)
    clabels = string.(vec(dataset["Cnames"]))

    # x, y, z positions must be within -1..+1
    x, y = _locnorm(x, y)

    z = zeros(ch_n)
    radius = zeros(length(clabels))
    theta = zeros(length(clabels))
    radius_sph = zeros(length(clabels))
    theta_sph = zeros(length(clabels))
    phi_sph = zeros(length(clabels))

    locs = DataFrame(:channel=>1:length(clabels), :labels=>clabels, :loc_theta=>theta, :loc_radius=>radius, :loc_x=>x, :loc_y=>y, :loc_z=>z, :loc_radius_sph=>radius_sph, :loc_theta_sph=>theta_sph, :loc_phi_sph=>phi_sph)

    locs = _round_locs(locs)
    locs_cart2sph!(locs)
    locs_cart2pol!(locs)

    return locs
    
end
