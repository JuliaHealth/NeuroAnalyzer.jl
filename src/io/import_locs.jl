export import_locs
export import_locs_ced
export import_locs_locs
export import_locs_csd
export import_locs_elc
export import_locs_tsv
export import_locs_sfp
export import_locs_geo
export import_locs_mat
export import_locs_txt
export import_locs_dat
export import_locs_asc

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
- TXT
- DAT
- ASC

This is a meta-function that triggers appropriate `import_locs_*()` function. File format is detected based on file extension.

# Arguments

- `file_name::String`: name of the file to load

# Returns

- `locs::DataFrame`
"""
function import_locs(file_name::String)

    @assert isfile(file_name) "File $file_name cannot be loaded."

    _info("Send standard locations for your channels to adam.wysokinski@neuroanalyzer.org")
    _info("Nose direction is set at '+Y'")

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
    elseif splitext(file_name)[2] == ".txt"
        locs = import_locs_txt(file_name)
    elseif splitext(file_name)[2] == ".dat"
        locs = import_locs_dat(file_name)
    elseif splitext(file_name)[2] == ".asc"
        locs = import_locs_asc(file_name)
    else
        @error "Unknown file format."
    end

    _locs_round!(locs)

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

    @assert isfile(file_name) "$file_name not found."
    @assert lowercase(splitext(file_name)[2]) == ".ced" "Not CED file."

    locs = CSV.read(file_name, delim="\t", stringtype=String, DataFrame)

    colnames = lowercase.(names(locs))
    DataFrames.rename!(locs, Symbol.(colnames))

    "labels" in colnames && (clabels = string.(lstrip.(locs[!, "labels"])))
    "label" in colnames && (clabels = string.(lstrip.(locs[!, "label"])))
    "name" in colnames && (clabels = string.(lstrip.(locs[!, "name"])))

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

    locs = DataFrame(:labels=>clabels, :loc_radius=>radius, :loc_theta=>theta, :loc_x=>x, :loc_y=>y, :loc_z=>z, :loc_radius_sph=>radius_sph, :loc_theta_sph=>theta_sph, :loc_phi_sph=>phi_sph)

    locs_swapxy!(locs, polar=true, cart=true, spherical=true)
    locs_flipx!(locs, polar=true, cart=false, spherical=false)

    locs_normalize!(locs)
    _locs_round!(locs)

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

    @assert isfile(file_name) "$file_name not found."
    @assert lowercase(splitext(file_name)[2]) == ".locs" "This is not LOCS file."

    locs = CSV.read(file_name, header=false, delim="\t", stringtype=String, DataFrame)

    DataFrames.rename!(locs, [:number, :theta, :radius, :labels])

    clabels = string.(lstrip.(locs[!, "labels"]))

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

    theta_sph = theta
    radius_sph = radius

    locs = DataFrame(:labels=>clabels, :loc_radius=>radius, :loc_theta=>theta, :loc_x=>x, :loc_y=>y, :loc_z=>z, :loc_radius_sph=>radius_sph, :loc_theta_sph=>theta_sph, :loc_phi_sph=>phi_sph)

    locs_swapxy!(locs, polar=true, cart=false, spherical=false)
    locs_flipx!(locs, polar=true, cart=false, spherical=false)

    locs[!, :loc_phi_sph] = zeros(nrow(locs))

    locs_normalize!(locs)
    _locs_round!(locs)

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

    @assert isfile(file_name) "$file_name not found."
    @assert lowercase(splitext(file_name)[2]) == ".elc" "This is not ELC file."

    f = open(file_name, "r")
    elc_file = readlines(f)
    close(f)

    locs_n = 0
    locs_l = 0
    for idx in eachindex(elc_file)
        if occursin("NumberPositions", elc_file[idx])
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

    locs = DataFrame(:labels=>clabels, :loc_radius=>radius, :loc_theta=>theta, :loc_x=>x, :loc_y=>y, :loc_z=>z, :loc_radius_sph=>radius_sph, :loc_theta_sph=>theta_sph, :loc_phi_sph=>phi_sph)

    locs_cart2sph!(locs)
    locs_cart2pol!(locs)

    locs_normalize!(locs)
    _locs_round!(locs)

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

    @assert isfile(file_name) "$file_name not found."
    @assert lowercase(splitext(file_name)[2]) == ".tsv" "This is not TSV file."

    locs = CSV.read(file_name, header=true, delim="\t", ignorerepeated=true, stringtype=String, DataFrame)

    colnames = lowercase.(names(locs))
    DataFrames.rename!(locs, Symbol.(colnames))

    "labels" in colnames && (clabels = string.(lstrip.(locs[!, "labels"])))
    "label" in colnames && (clabels = string.(lstrip.(locs[!, "label"])))
    "name" in colnames && (clabels = string.(lstrip.(locs[!, "name"])))
    "site" in colnames && (clabels = string.(lstrip.(locs[!, "site"])))

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

    locs = DataFrame(:labels=>clabels, :loc_radius=>radius, :loc_theta=>theta, :loc_x=>x, :loc_y=>y, :loc_z=>z, :loc_radius_sph=>radius_sph, :loc_theta_sph=>theta_sph, :loc_phi_sph=>phi_sph)

    locs_cart2sph!(locs)
    locs_cart2pol!(locs)

    locs_normalize!(locs)
    _locs_round!(locs)

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

    @assert isfile(file_name) "$file_name not found."
    @assert lowercase(splitext(file_name)[2]) == ".sfp" "This is not SFP file."

    locs = CSV.read(file_name, header=false, stringtype=String, DataFrame)
    _info("Checking TAB as delimeter")
    size(locs, 2) != 4 && (locs = CSV.read(file_name, header=false, delim="/t", ignorerepeated=true, stringtype=String, DataFrame))
    _info("Checking SPACE as delimeter")
    size(locs, 2) != 4 && (locs = CSV.read(file_name, header=false, delim=" ", ignorerepeated=true, stringtype=String, DataFrame))
    @assert size(locs, 2) == 4 "File $file_name cannot be opened, check delimeters."

    DataFrames.rename!(locs, [:label, :x, :y, :z])

    clabels = string.(lstrip.(locs[!, "label"]))

    x = Float64.(locs[!, :x])
    y = Float64.(locs[!, :y])
    z = Float64.(locs[!, :z])

    # x, y, z positions must be within -1..+1
    t = x[1]
    x, y, z = _locs_norm(x, y, z)
    t -= x[1]
    # sometimes positions are shifted along x-axis, remove the shift
    x .+= abs(t)

    radius = zeros(length(clabels))
    theta = zeros(length(clabels))
    radius_sph = zeros(length(clabels))
    theta_sph = zeros(length(clabels))
    phi_sph = zeros(length(clabels))

    locs = DataFrame(:labels=>clabels, :loc_radius=>radius, :loc_theta=>theta, :loc_x=>x, :loc_y=>y, :loc_z=>z, :loc_radius_sph=>radius_sph, :loc_theta_sph=>theta_sph, :loc_phi_sph=>phi_sph)

    # center X-axis at 0, 0, 0
    x_range = (maximum(x) + abs(minimum(x))) / 2
    locs[:, :loc_x] .-= x_range

    locs_cart2sph!(locs)
    locs_cart2pol!(locs)

    locs_normalize!(locs)
    _locs_round!(locs)

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

    @assert isfile(file_name) "$file_name not found."
    @assert lowercase(splitext(file_name)[2]) == ".csd" "This is not CSD file."

    locs = CSV.read(file_name, skipto=3, delim=' ', header=false, ignorerepeated=true, stringtype=String, DataFrame)

    DataFrames.rename!(locs, [:labels, :theta_sph, :phi_sph, :radius_sph, :x, :y, :z, :surface])
    clabels = string.(lstrip.(locs[!, "labels"]))

    x = Float64.(locs[!, "x"])
    y = Float64.(locs[!, "y"])
    z = Float64.(locs[!, "z"])
    radius_sph = Float64.(locs[!, "radius_sph"])
    theta_sph = Float64.(locs[!, "theta_sph"])
    phi_sph = Float64.(locs[!, "phi_sph"])

    radius = zeros(length(x))
    theta = zeros(length(y))
    for idx in eachindex(x)
        radius[idx], theta[idx] = sph2pol(radius_sph[idx], theta_sph[idx], phi_sph[idx])
    end

    locs = DataFrame(:labels=>clabels, :loc_radius=>radius, :loc_theta=>theta, :loc_x=>x, :loc_y=>y, :loc_z=>z, :loc_radius_sph=>radius_sph, :loc_theta_sph=>theta_sph, :loc_phi_sph=>phi_sph)

    locs_normalize!(locs)
    _locs_round!(locs)

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

    @assert isfile(file_name) "$file_name not found."
    @assert lowercase(splitext(file_name)[2]) == ".geo" "This is not GEO file."

    f = open(file_name, "r")
    locs = readlines(f)
    close(f)

    l1 = 0
    l2 = 0
    for idx in eachindex(locs)
        locs[idx] == "View\"\"{" && (l1 = idx + 1)
        locs[idx] == "};" && (l2 = idx - 1)
    end
    locs = locs[l1+1:2:l2]

    clabels = repeat([""], length(locs))
    x = zeros(length(locs))
    y = zeros(length(locs))
    z = zeros(length(locs))

    p = r"(.+)(\(.+\)){(.+)}"
    for idx in eachindex(clabels)
        m = match(p, locs[idx])
        clabels[idx] = replace(m[3], "\"" => "")
        tmp = replace(m[2], "(" => "")
        tmp = replace(tmp, ")" => "")
        x[idx], y[idx], z[idx], = parse.(Float64, split(tmp, ", "))
    end

    x, y, z = _locs_norm(x, y, z)

    # center x at 0
    x_adj = x[findfirst(isequal("Cz"), clabels)]
    x .-= x_adj
    x, y, z = _locs_norm(x, y, z)

    radius = zeros(length(clabels))
    theta = zeros(length(clabels))
    radius_sph = zeros(length(clabels))
    theta_sph = zeros(length(clabels))
    phi_sph = zeros(length(clabels))

    locs = DataFrame(:labels=>clabels, :loc_radius=>radius, :loc_theta=>theta, :loc_x=>x, :loc_y=>y, :loc_z=>z, :loc_radius_sph=>radius_sph, :loc_theta_sph=>theta_sph, :loc_phi_sph=>phi_sph)

    locs = locs_cart2sph(locs)
    locs = locs_cart2pol(locs)

    locs_normalize!(locs)
    _locs_round!(locs)

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

    @assert isfile(file_name) "$file_name not found."
    @assert lowercase(splitext(file_name)[2]) == ".mat" "This is not MAT file."

    dataset = matread(file_name)
    x = dataset["Cpos"][1, :]
    y = dataset["Cpos"][2, :]
    r = dataset["Rxy"]
    ch_n = length(x)
    clabels = string.(vec(dataset["Cnames"]))

    # x, y, z positions must be within -1..+1
    x, y = _locs_norm(x, y)

    z = zeros(ch_n)
    radius = zeros(length(clabels))
    theta = zeros(length(clabels))
    radius_sph = zeros(length(clabels))
    theta_sph = zeros(length(clabels))
    phi_sph = zeros(length(clabels))

    locs = DataFrame(:labels=>clabels, :loc_radius=>radius, :loc_theta=>theta, :loc_x=>x, :loc_y=>y, :loc_z=>z, :loc_radius_sph=>radius_sph, :loc_theta_sph=>theta_sph, :loc_phi_sph=>phi_sph)

    locs_cart2sph!(locs)
    locs_cart2pol!(locs)

    locs_normalize!(locs)
    _locs_round!(locs)

    return locs

end

"""
    import_locs_txt(file_name)

Load channel locations from TXT file.

# Arguments

- `file_name::String`

# Returns

- `locs::DataFrame`
"""
function import_locs_txt(file_name::String)

    @assert isfile(file_name) "$file_name not found."
    @assert lowercase(splitext(file_name)[2]) == ".txt" "This is not TXT file."

    locs = CSV.read(file_name, header=true, delim="\t", stringtype=String, DataFrame)

    DataFrames.rename!(locs, [:labels, :theta, :phi])
    clabels = lstrip.(locs[!, "labels"])

    x = zeros(length(clabels))
    y = zeros(length(clabels))
    z = zeros(length(clabels))

    radius_sph = ones(length(clabels))
    theta_sph = Float64.(locs[!, "theta"])
    phi_sph = Float64.(locs[!, "phi"])

    radius = radius_sph
    theta = theta_sph

    locs = DataFrame(:labels=>clabels, :loc_radius=>radius, :loc_theta=>theta, :loc_x=>x, :loc_y=>y, :loc_z=>z, :loc_radius_sph=>radius_sph, :loc_theta_sph=>theta_sph, :loc_phi_sph=>phi_sph)

    locs_sph2cart!(locs)
    locs_swapxy!(locs, polar=false, cart=true, spherical=false)
    locs_rotx!(locs, a=90, polar=false, cart=true, spherical=false)

    q1_idx = locs[!, :loc_x] .>= 0 .&& locs[!, :loc_y] .>= 0
    q2_idx = locs[!, :loc_x] .< 0 .&& locs[!, :loc_y] .>= 0
    q3_idx = locs[!, :loc_x] .< 0 .&& locs[!, :loc_y] .< 0
    q4_idx = locs[!, :loc_x] .>= 0 .&& locs[!, :loc_y] .< 0

    x = locs[!, :loc_x]
    y = locs[!, :loc_y]

    locs[q1_idx, :loc_x] .= -locs[q1_idx, :loc_x]
    locs[q2_idx, :loc_x] .= -locs[q2_idx, :loc_x]
    locs[q2_idx, :loc_y] .= -locs[q2_idx, :loc_y]
    locs[q3_idx, :loc_x] .= -locs[q3_idx, :loc_x]
    locs[q3_idx, :loc_y] .= -locs[q3_idx, :loc_y]
    locs[q4_idx, :loc_x] .= -locs[q4_idx, :loc_x]

    locs_cart2sph!(locs)
    locs_sph2pol!(locs)

    locs_normalize!(locs)
    _locs_round!(locs)

    return locs

end

"""
    import_locs_dat(file_name)

Load channel locations from DAT file.

# Arguments

- `file_name::String`

# Returns

- `locs::DataFrame`
"""
function import_locs_dat(file_name::String)

    @assert isfile(file_name) "$file_name not found."
    @assert lowercase(splitext(file_name)[2]) == ".dat" "Not DAT file."

    locs = CSV.read(file_name, ignorerepeated=true, delim=' ', stringtype=String, header=0, DataFrame)
    if ncol(locs) == 4
        if typeof(locs[!, 2]) == Vector{String}
            colnames = ["channel", "labels", "x", "y"]
        else
            colnames = ["channel", "x", "y", "z"]
        end
    elseif ncol(locs) == 5
        colnames = ["channel", "labels", "x", "y", "z"]
    end

    DataFrames.rename!(locs, colnames)
    if "labels" in colnames
        clabels = string.(lstrip.(locs[!, "labels"]))
    else
        clabels = string.(locs[!, "channel"])
    end

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

    locs = DataFrame(:labels=>clabels, :loc_radius=>radius, :loc_theta=>theta, :loc_x=>x, :loc_y=>y, :loc_z=>z, :loc_radius_sph=>radius_sph, :loc_theta_sph=>theta_sph, :loc_phi_sph=>phi_sph)

    locs_center!(locs, polar=false, spherical=false)
    locs_cart2pol!(locs)
    locs_cart2sph!(locs)
    locs_normalize!(locs)
    _locs_round!(locs)

    return locs

end

"""
    import_locs_asc(file_name)

Load channel locations from ASC file.

# Arguments

- `file_name::String`

# Returns

- `locs::DataFrame`
"""
function import_locs_asc(file_name::String)

    @assert isfile(file_name) "$file_name not found."
    @assert lowercase(splitext(file_name)[2]) == ".asc" "Not ASC file."

    buffer = readlines(file_name)
    # remove comments
    for idx in length(buffer):-1:1
        buffer[idx][1] == ';' && deleteat!(buffer, idx)
    end
    labels = String[]
    for idx in 1:length(buffer)
        buffer[idx][1] == '#' && push!(labels, buffer[idx])
    end
    labels_regexp = match.(r"\#.+ (.+)", labels)
    clabels = String[]
    for idx in 1:length(labels_regexp)
        push!(clabels, labels_regexp[idx][1])
    end
    buffer = buffer[length(labels) + 1:2 * length(labels)]
    locs = zeros(length(labels), 4)
    locs_regexp = match.(r"([0-9]+ +)([0-9]+ +)([0-9]+\.[0-9]+ +)([0-9]+\.[0-9]+ +)([0-9]+\.[0-9]+ +)([0-9]+\.[0-9]+)", buffer)
    for idx in 1:length(labels_regexp)
        locs[idx, 1] = parse(Float64, strip(locs_regexp[idx][3]))
        locs[idx, 2] = parse(Float64, strip(locs_regexp[idx][4]))
        locs[idx, 3] = parse(Float64, strip(locs_regexp[idx][5]))
        locs[idx, 4] = parse(Float64, strip(locs_regexp[idx][6]))
    end

    x = locs[:, 1]
    y = locs[:, 2]
    z = zeros(length(clabels))
    radius = zeros(length(clabels))
    theta = zeros(length(clabels))
    radius_sph = zeros(length(clabels))
    theta_sph = zeros(length(clabels))
    phi_sph = zeros(length(clabels))

    locs = DataFrame(:labels=>clabels, :loc_radius=>radius, :loc_theta=>theta, :loc_x=>x, :loc_y=>y, :loc_z=>z, :loc_radius_sph=>radius_sph, :loc_theta_sph=>theta_sph, :loc_phi_sph=>phi_sph)

    locs_center!(locs, polar=false, spherical=false)
    locs_flipy!(locs)
    locs_cart2pol!(locs)
    locs_cart2sph!(locs)
    locs_normalize!(locs)
    _locs_round!(locs)

    return locs

end