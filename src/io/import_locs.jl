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
export import_locs_csv

"""
    import_locs(file_name)

Load channel locations from a file. The format is detected from the extension.

Supported formats: CED, ELC, LOCS, TSV, SFP, CSD, GEO, MAT, TXT, DAT, ASC, CSV.

# Arguments

- `file_name::String`: path to the locations file

# Returns

- `DataFrame`

# Throws

- `ArgumentError` if the file does not exist or the extension is not recognized
"""
function import_locs(file_name::String)::DataFrame

    isfile(file_name) ||
        throw(ArgumentError("File $file_name cannot be loaded."))

    _info("Send standard locations for your channels to adam.wysokinski@neuroanalyzer.org")
    _info("Nose direction is set at '+Y'")

    # extract extension once; used for dispatch and error message
    ext = lowercase(splitext(file_name)[2])

    # dispatch table - each sub-function handles its own validation.
    locs = if ext == ".ced"
        import_locs_ced(file_name)
    elseif ext == ".elc"
        import_locs_elc(file_name)
    elseif ext == ".locs"
        import_locs_locs(file_name)
    elseif ext == ".tsv"
        import_locs_tsv(file_name)
    elseif ext == ".sfp"
        import_locs_sfp(file_name)
    elseif ext == ".csd"
        import_locs_csd(file_name)
    elseif ext == ".geo"
        import_locs_geo(file_name)
    elseif ext == ".mat"
        import_locs_mat(file_name)
    elseif ext == ".txt"
        import_locs_txt(file_name)
    elseif ext == ".dat"
        import_locs_dat(file_name)
    elseif ext == ".asc"
        import_locs_asc(file_name)
    elseif ext == ".csv"
        import_locs_csv(file_name)
    else
        throw(ArgumentError(
            "Unknown locations file format \"$ext\". " *
            "Supported: .ced, .elc, .locs, .tsv, .sfp, .csd, .geo, " *
            ".mat, .txt, .dat, .asc, .csv"))
    end

    _locs_round!(locs)

    return locs

end

"""
    import_locs_ced(file_name)

Load channel locations from a CED (EEGLAB) file.

# Arguments

- `file_name::String`

# Returns

- `DataFrame`
"""
function import_locs_ced(file_name::String)::DataFrame

    isfile(file_name) ||  throw(ArgumentError("$file_name not found."))
    lowercase(splitext(file_name)[2]) == ".ced" ||
        throw(ArgumentError("$file_name is not a CED file."))

    locs_raw = CSV.read(file_name; delim = "\t", stringtype = String, DataFrame)
    colnames  = lowercase.(names(locs_raw))
    DataFrames.rename!(locs_raw, Symbol.(colnames))

    clabels = if "labels" in colnames
        string.(lstrip.(locs_raw[!, "labels"]))
    elseif "label" in colnames
        string.(lstrip.(locs_raw[!, "label"]))
    elseif "name" in colnames
        string.(lstrip.(locs_raw[!, "name"]))
    else
        throw(ArgumentError("$file_name contains no recognized label column (labels/label/name)."))
    end

    n = length(clabels)
    x = "x" in colnames ? Float64.(locs_raw[!, "x"]) : zeros(n)
    y = "y" in colnames ? Float64.(locs_raw[!, "y"]) : zeros(n)
    z = "z" in colnames ? Float64.(locs_raw[!, "z"]) : zeros(n)
    theta = "theta" in colnames ? Float64.(locs_raw[!, "theta"]) : zeros(n)
    radius = "radius" in colnames ? Float64.(locs_raw[!, "radius"]) : zeros(n)
    radius_sph = "sph_radius" in colnames ? Float64.(locs_raw[!, "sph_radius"]) : zeros(n)
    theta_sph = "sph_theta" in colnames ? Float64.(locs_raw[!, "sph_theta"]) : zeros(n)
    phi_sph = "sph_phi" in colnames ? Float64.(locs_raw[!, "sph_phi"]) : zeros(n)

    locs = DataFrame(
        :label => clabels,
        :loc_radius => radius,
        :loc_theta => theta,
        :loc_x => x,
        :loc_y => y,
        :loc_z => z,
        :loc_radius_sph => radius_sph,
        :loc_theta_sph => theta_sph,
        :loc_phi_sph => phi_sph
    )

    locs_swapxy!(locs, polar = true, cart = true, spherical = true)
    locs_flipx!(locs, polar = true, cart = false, spherical = false)
    locs_normalize!(locs)
    _locs_round!(locs)

    return locs

end

"""
    import_locs_locs(file_name)

Load channel locations from an EEGLAB LOCS file.

# Arguments

- `file_name::String`

# Returns

- `DataFrame`
"""
function import_locs_locs(file_name::String)::DataFrame

    isfile(file_name) || throw(ArgumentError("$file_name not found."))
    lowercase(splitext(file_name)[2]) == ".locs" ||
        throw(ArgumentError("$file_name is not a LOCS file."))

    locs_raw = CSV.read(file_name, header = false, delim = "\t",
                        stringtype = String, DataFrame)
    DataFrames.rename!(locs_raw, [:number, :theta, :radius, :label])

    clabels = string.(lstrip.(locs_raw[!, :label]))
    n = length(clabels)
    theta = Float64.(locs_raw[!, "theta"])
    radius = Float64.(locs_raw[!, "radius"])

    locs = DataFrame(
        :label => clabels,
        :loc_radius => radius,
        :loc_theta => theta,
        :loc_x => zeros(n),
        :loc_y => zeros(n),
        :loc_z => zeros(n),
        :loc_radius_sph => copy(radius),
        :loc_theta_sph => copy(theta),
        :loc_phi_sph => zeros(n)
    )

    locs_swapxy!(locs, polar = true, cart = false, spherical = false)
    locs_flipx!(locs, polar = true, cart = false, spherical = false)
    locs[!, :loc_phi_sph] .= 0.0
    locs_normalize!(locs)
    _locs_round!(locs)

    return locs

end

"""
    import_locs_elc(file_name)

Load channel locations from an ELC file.

# Arguments

- `file_name::String`

# Returns

- `DataFrame`
"""
function import_locs_elc(file_name::String)::DataFrame

    isfile(file_name) || throw(ArgumentError("$file_name not found."))
    lowercase(splitext(file_name)[2]) == ".elc" ||
        throw(ArgumentError("$file_name is not an ELC file."))

    elc_file = readlines(file_name)

    locs_n = 0; locs_l = 0
    for idx in eachindex(elc_file)
        if occursin("NumberPositions", elc_file[idx])
            locs_n = parse(Int64, replace(elc_file[idx], "NumberPositions=" => ""))
            locs_l = idx + 2
        end
    end

    n = locs_n
    clabels = repeat([""], n)
    x = zeros(n); y = zeros(n); z = zeros(n)

    for (i, line_idx) in enumerate(locs_l:(locs_l + n - 1))
        l = lstrip(elc_file[line_idx])
        x[i], y[i], z[i] = parse.(Float64, split(l, ' '))
    end
    for (i, line_idx) in enumerate((locs_l + 1 + n):(locs_l + 2 * n))
        clabels[i] = elc_file[line_idx]
    end

    x = normalize_minmax(x)
    y = normalize_minmax(y)
    z = normalize_minmax(z)

    locs = DataFrame(
        :label => clabels,
        :loc_radius => zeros(n),
        :loc_theta => zeros(n),
        :loc_x => x,
        :loc_y => y,
        :loc_z => z,
        :loc_radius_sph => zeros(n),
        :loc_theta_sph => zeros(n),
        :loc_phi_sph => zeros(n)
    )

    locs_cart2sph!(locs); locs_cart2pol!(locs)
    locs_normalize!(locs); _locs_round!(locs)

    return locs

end

"""
    import_locs_tsv(file_name)

Load channel locations from a TSV (BIDS-style) file.

# Arguments

- `file_name::String`

# Returns

- `DataFrame`
"""
function import_locs_tsv(file_name::String)::DataFrame

    isfile(file_name) || throw(ArgumentError("$file_name not found."))
    lowercase(splitext(file_name)[2]) == ".tsv" ||
        throw(ArgumentError("$file_name is not a TSV file."))

    locs_raw = CSV.read(file_name, header = true, delim = "\t",
                        ignorerepeated = true, stringtype = String, DataFrame)
    colnames = lowercase.(names(locs_raw))
    DataFrames.rename!(locs_raw, Symbol.(colnames))

    clabels = if "labels" in colnames
        string.(lstrip.(locs_raw[!, "labels"]))
    elseif "label" in colnames
        string.(lstrip.(locs_raw[!, "label"]))
    elseif "name" in colnames
        string.(lstrip.(locs_raw[!, "name"]))
    elseif "site" in colnames
        string.(lstrip.(locs_raw[!, "site"]))
    else
        throw(ArgumentError("$file_name contains no recognized label column (labels/label/name/site)."))
    end

    n = length(clabels)
    x = "x" in colnames ? Float64.(locs_raw[!, "x"]) : zeros(n)
    y = "y" in colnames ? Float64.(locs_raw[!, "y"]) : zeros(n)
    z = "z" in colnames ? Float64.(locs_raw[!, "z"]) : zeros(n)
    theta = "theta" in colnames ? Float64.(locs_raw[!, "theta"]) : zeros(n)
    radius = "radius" in colnames ? Float64.(locs_raw[!, "radius"]) : zeros(n)
    radius_sph = "radius_sph" in colnames ? Float64.(locs_raw[!, "radius_sph"]) :
                 "radius" in colnames ? Float64.(locs_raw[!, "radius"]) : zeros(n)
    theta_sph  = "theta_sph" in colnames ? Float64.(locs_raw[!, "theta_sph"]) :
                 "theta" in colnames ? Float64.(locs_raw[!, "theta"]) : zeros(n)
    phi_sph    = "phi_sph" in colnames ? Float64.(locs_raw[!, "phi_sph"]) :
                 "phi" in colnames ? Float64.(locs_raw[!, "phi"]) : zeros(n)

    locs = DataFrame(
        :label => clabels,
        :loc_radius => radius,
        :loc_theta => theta,
        :loc_x => x,
        :loc_y => y,
        :loc_z => z,
        :loc_radius_sph => radius_sph,
        :loc_theta_sph => theta_sph,
        :loc_phi_sph => phi_sph
    )

    locs_cart2sph!(locs); locs_cart2pol!(locs)
    locs_normalize!(locs); _locs_round!(locs)

    return locs

end

"""
    import_locs_sfp(file_name)

Load channel locations from an SFP file.

# Arguments

- `file_name::String`

# Returns

- `DataFrame`
"""
function import_locs_sfp(file_name::String)::DataFrame

    isfile(file_name) || throw(ArgumentError("$file_name not found."))
    lowercase(splitext(file_name)[2]) == ".sfp" ||
        throw(ArgumentError("$file_name is not an SFP file."))

    # try common delimiters in order; SFP files are inconsistently delimited.
    locs_raw = CSV.read(file_name, header = false, stringtype = String, DataFrame)
    if size(locs_raw, 2) != 4
        _info("Checking TAB as delimiter")
        locs_raw = CSV.read(file_name, header = false, delim = "\t",
                            ignorerepeated = true, stringtype = String, DataFrame)
    end
    if size(locs_raw, 2) != 4
        _info("Checking SPACE as delimiter")
        locs_raw = CSV.read(file_name, header = false, delim = " ",
                            ignorerepeated = true, stringtype = String, DataFrame)
    end
    size(locs_raw, 2) == 4 ||
        throw(ArgumentError("$file_name could not be parsed — check delimiters."))

    DataFrames.rename!(locs_raw, [:label, :x, :y, :z])
    clabels = string.(lstrip.(locs_raw[!, "label"]))
    x = Float64.(locs_raw[!, :x])
    y = Float64.(locs_raw[!, :y])
    z = Float64.(locs_raw[!, :z])

    t = x[1]
    x, y, z = _locs_norm(x, y, z)
    t -= x[1]
    x .+= abs(t)

    n = length(clabels)
    locs = DataFrame(
        :label => clabels,
        :loc_radius => zeros(n),
        :loc_theta => zeros(n),
        :loc_x => x,
        :loc_y => y,
        :loc_z => z,
        :loc_radius_sph => zeros(n),
        :loc_theta_sph => zeros(n),
        :loc_phi_sph => zeros(n)
    )

    x_range = (maximum(x) + abs(minimum(x))) / 2
    locs[:, :loc_x] .-= x_range
    locs_cart2sph!(locs); locs_cart2pol!(locs)
    locs_normalize!(locs); _locs_round!(locs)

    return locs

end

"""
    import_locs_csd(file_name)

Load channel locations from a CSD file.

# Arguments

- `file_name::String`

# Returns

- `DataFrame`
"""
function import_locs_csd(file_name::String)::DataFrame

    isfile(file_name) || throw(ArgumentError("$file_name not found."))
    lowercase(splitext(file_name)[2]) == ".csd" ||
        throw(ArgumentError("$file_name is not a CSD file."))

    locs_raw = CSV.read(file_name, skipto = 3, delim = ' ', header = false,
                        ignorerepeated = true, stringtype = String, DataFrame)
    DataFrames.rename!(locs_raw, [:label, :theta_sph, :phi_sph, :radius_sph, :x, :y, :z, :surface])

    clabels = string.(lstrip.(locs_raw[!, :label]))
    x = Float64.(locs_raw[!, :x])
    y = Float64.(locs_raw[!, :y])
    z = Float64.(locs_raw[!, :z])
    radius_sph = Float64.(locs_raw[!, :radius_sph])
    theta_sph = Float64.(locs_raw[!, :theta_sph])
    phi_sph = Float64.(locs_raw[!, :phi_sph])
    n = length(x)

    radius = zeros(n); theta = zeros(n)
    for idx in 1:n
        radius[idx], theta[idx] = sph2pol(radius_sph[idx], theta_sph[idx], phi_sph[idx])
    end

    locs = DataFrame(
        :label => clabels,
        :loc_radius => radius,
        :loc_theta => theta,
        :loc_x => x,
        :loc_y => y,
        :loc_z => z,
        :loc_radius_sph => radius_sph,
        :loc_theta_sph => theta_sph,
        :loc_phi_sph => phi_sph
    )

    locs_normalize!(locs); _locs_round!(locs)

    return locs

end

"""
    import_locs_geo(file_name)

Load channel locations from a GEO file.

# Arguments

- `file_name::String`

# Returns

- `DataFrame`
"""
function import_locs_geo(file_name::String)::DataFrame

    isfile(file_name) || throw(ArgumentError("$file_name not found."))
    lowercase(splitext(file_name)[2]) == ".geo" ||
        throw(ArgumentError("$file_name is not a GEO file."))

    raw = readlines(file_name)
    l1 = 0; l2 = 0
    for idx in eachindex(raw)
        raw[idx] == "View\"\"{" && (l1 = idx + 1)
        raw[idx] == "};" && (l2 = idx - 1)
    end
    entries = raw[(l1 + 1):2:l2]
    n = length(entries)

    clabels = repeat([""], n)
    x = zeros(n); y = zeros(n); z = zeros(n)
    p = r"(.+)(\(.+\)){(.+)}"
    for idx in 1:n
        m = match(p, entries[idx])
        clabels[idx] = replace(m[3], "\"" => "")
        tmp = replace(replace(m[2], "(" => ""), ")" => "")
        x[idx], y[idx], z[idx] = parse.(Float64, split(tmp, ", "))
    end

    x, y, z = _locs_norm(x, y, z)
    cz_idx = findfirst(isequal("Cz"), clabels)
    if !isnothing(cz_idx)
        x .-= x[cz_idx]
        x, y, z = _locs_norm(x, y, z)
    end

    locs = DataFrame(
        :label => clabels,
        :loc_radius => zeros(n),
        :loc_theta => zeros(n),
        :loc_x => x,
        :loc_y => y,
        :loc_z => z,
        :loc_radius_sph => zeros(n),
        :loc_theta_sph => zeros(n),
        :loc_phi_sph => zeros(n)
    )

    locs_cart2sph!(locs); locs_cart2pol!(locs)
    locs_normalize!(locs); _locs_round!(locs)

    return locs

end

"""
    import_locs_mat(file_name)

Load channel locations from a MATLAB MAT file.

# Arguments

- `file_name::String`

# Returns

- `DataFrame`
"""
function import_locs_mat(file_name::String)::DataFrame

    isfile(file_name) || throw(ArgumentError("$file_name not found."))
    lowercase(splitext(file_name)[2]) == ".mat" ||
        throw(ArgumentError("$file_name is not a MAT file."))

    dataset = matread(file_name)
    x = dataset["Cpos"][1, :]
    y = dataset["Cpos"][2, :]
    r = dataset["Rxy"]
    ch_n = length(x)
    clabels = string.(vec(dataset["Cnames"]))

    x_norm, y_norm = _locs_norm(x, y)
    x = x_norm .* r
    y = y_norm .* r

    clabels = replace.(clabels, r"MEG\s*0+" => "MEG ")
    clabels = replace.(clabels, "  " => " ")

    n = length(clabels)
    locs = DataFrame(
        :label => clabels,
        :loc_radius => zeros(n),
        :loc_theta => zeros(n),
        :loc_x => x,
        :loc_y => y,
        :loc_z => zeros(n),
        :loc_radius_sph => zeros(n),
        :loc_theta_sph => zeros(n),
        :loc_phi_sph => zeros(n))

    locs_cart2sph!(locs)
    locs_cart2pol!(locs)
    _locs_round!(locs)

    return locs

end

"""
    import_locs_txt(file_name)

Load channel locations from a TXT file (spherical theta/phi format).

# Arguments

- `file_name::String`

# Returns

- `DataFrame`
"""
function import_locs_txt(file_name::String)::DataFrame

    isfile(file_name) || throw(ArgumentError("$file_name not found."))
    lowercase(splitext(file_name)[2]) == ".txt" ||
        throw(ArgumentError("$file_name is not a TXT file."))

    locs_raw = CSV.read(file_name, header = true, delim = "\t",
                        stringtype = String, DataFrame)
    DataFrames.rename!(locs_raw, [:label, :theta, :phi])

    clabels = string.(lstrip.(locs_raw[!, :label]))
    n = length(clabels)
    theta_sph = Float64.(locs_raw[!, "theta"])
    phi_sph = Float64.(locs_raw[!, "phi"])

    locs = DataFrame(
        :label => clabels,
        :loc_radius => ones(n),
        :loc_theta => copy(theta_sph),
        :loc_x => zeros(n),
        :loc_y => zeros(n),
        :loc_z => zeros(n),
        :loc_radius_sph => ones(n),
        :loc_theta_sph => theta_sph,
        :loc_phi_sph => phi_sph
    )

    locs_sph2cart!(locs)
    locs_swapxy!(locs; polar = false, cart = true, spherical = false)
    locs_rotx!(locs;   a = 90, polar = false, cart = true, spherical = false)

    q1 = locs[!, :loc_x] .>= 0 .&& locs[!, :loc_y] .>= 0
    q2 = locs[!, :loc_x] .< 0  .&& locs[!, :loc_y] .>= 0
    q3 = locs[!, :loc_x] .< 0  .&& locs[!, :loc_y] .< 0
    q4 = locs[!, :loc_x] .>= 0 .&& locs[!, :loc_y] .< 0

    locs[q1, :loc_x] .= -locs[q1, :loc_x]
    locs[q2, :loc_x] .= -locs[q2, :loc_x];  locs[q2, :loc_y] .= -locs[q2, :loc_y]
    locs[q3, :loc_x] .= -locs[q3, :loc_x];  locs[q3, :loc_y] .= -locs[q3, :loc_y]
    locs[q4, :loc_x] .= -locs[q4, :loc_x]

    locs_cart2sph!(locs)
    locs_sph2pol!(locs)
    locs_normalize!(locs)
    _locs_round!(locs)

    return locs

end

"""
    import_locs_dat(file_name)

Load channel locations from a DAT file.

# Arguments

- `file_name::String`

# Returns

- `DataFrame`
"""
function import_locs_dat(file_name::String)::DataFrame

    isfile(file_name) || throw(ArgumentError("$file_name not found."))
    lowercase(splitext(file_name)[2]) == ".dat" ||
        throw(ArgumentError("$file_name is not a DAT file."))

    locs_raw = CSV.read(file_name, ignorerepeated = true, delim = ' ',
                        stringtype = String, header = 0, DataFrame)

    # detect column layout from number of columns and type of column 2
    colnames = if ncol(locs_raw) == 4
        typeof(locs_raw[!, 2]) == Vector{String} ?
            ["channel", "labels", "x", "y"] : ["channel", "x", "y", "z"]
    elseif ncol(locs_raw) == 5
        ["channel", "labels", "x", "y", "z"]
    else
        throw(ArgumentError(
            "$file_name has $(ncol(locs_raw)) columns; expected 4 or 5."))
    end

    DataFrames.rename!(locs_raw, colnames)
    clabels = "labels" in colnames ?
        string.(lstrip.(locs_raw[!, "labels"])) : string.(locs_raw[!, "channel"])

    n = length(clabels)
    x = "x" in colnames ? Float64.(locs_raw[!, "x"]) : zeros(n)
    y = "y" in colnames ? Float64.(locs_raw[!, "y"]) : zeros(n)
    z = "z" in colnames ? Float64.(locs_raw[!, "z"]) : zeros(n)
    theta = "theta" in colnames ? Float64.(locs_raw[!, "theta"]) : zeros(n)
    radius = "radius" in colnames ? Float64.(locs_raw[!, "radius"]) : zeros(n)
    radius_sph = "sph_radius" in colnames ? Float64.(locs_raw[!, "sph_radius"]) : zeros(n)
    theta_sph = "sph_theta" in colnames ? Float64.(locs_raw[!, "sph_theta"]) : zeros(n)
    phi_sph = "sph_phi" in colnames ? Float64.(locs_raw[!, "sph_phi"]) : zeros(n)

    locs = DataFrame(
        :label => clabels,
        :loc_radius => radius,
        :loc_theta => theta,
        :loc_x => x,
        :loc_y => y,
        :loc_z => z,
        :loc_radius_sph => radius_sph,
        :loc_theta_sph => theta_sph,
        :loc_phi_sph => phi_sph
    )

    locs_center!(locs; polar = false, spherical = false)
    locs_cart2pol!(locs); locs_cart2sph!(locs)
    locs_normalize!(locs); _locs_round!(locs)

    return locs

end

"""
    import_locs_asc(file_name)

Load channel locations from an ASC file.

# Arguments

- `file_name::String`

# Returns

- `DataFrame`
"""
function import_locs_asc(file_name::String)::DataFrame

    isfile(file_name) || throw(ArgumentError("$file_name not found."))
    lowercase(splitext(file_name)[2]) == ".asc" ||
        throw(ArgumentError("$file_name is not an ASC file."))

    buffer = readlines(file_name)
    filter!(l -> !startswith(l, ';'), buffer) # remove comments

    label_lines = filter(l -> startswith(l, '#'), buffer)
    clabels = [m[1] for m in match.(r"\#.+ (.+)", label_lines)]
    n = length(clabels)

    data_lines = buffer[(n + 1):(2 * n)]
    locs_mat = zeros(n, 4)
    re = r"([0-9]+ +)([0-9]+ +)([0-9]+\.[0-9]+ +)([0-9]+\.[0-9]+ +)([0-9]+\.[0-9]+ +)([0-9]+\.[0-9]+)"
    for (idx, m) in enumerate(match.(re, data_lines))
        locs_mat[idx, 1] = parse(Float64, strip(m[3]))
        locs_mat[idx, 2] = parse(Float64, strip(m[4]))
        locs_mat[idx, 3] = parse(Float64, strip(m[5]))
        locs_mat[idx, 4] = parse(Float64, strip(m[6]))
    end

    locs = DataFrame(
        :label => clabels,
        :loc_radius => zeros(n),
        :loc_theta => zeros(n),
        :loc_x => locs_mat[:, 1],
        :loc_y => locs_mat[:, 2],
        :loc_z => zeros(n),
        :loc_radius_sph => zeros(n),
        :loc_theta_sph => zeros(n),
        :loc_phi_sph => zeros(n)
    )

    locs_center!(locs; polar = false, spherical = false)
    locs_flipy!(locs)
    locs_cart2pol!(locs)
    locs_cart2sph!(locs)
    locs_normalize!(locs)
    _locs_round!(locs)

    return locs

end

"""
    import_locs_csv(file_name)

Load channel locations from a NeuroAnalyzer standard CSV file.

# Arguments

- `file_name::String`

# Returns

- `DataFrame`
"""
function import_locs_csv(file_name::String)::DataFrame

    isfile(file_name) || throw(ArgumentError("$file_name not found."))
    lowercase(splitext(file_name)[2]) == ".csv" ||
        throw(ArgumentError("$file_name is not a CSV file."))

    locs = CSV.read(file_name, header = true, delim = ",",
                    stringtype = String, DataFrame)

    expected = ["label", "loc_radius", "loc_theta", "loc_x", "loc_y", "loc_z",
                "loc_radius_sph", "loc_theta_sph", "loc_phi_sph"]
    names(locs) == expected ||
        throw(ArgumentError(
            "$file_name is not a NeuroAnalyzer locs CSV file. " *
            "Expected columns: $(join(expected, ", "))."))

    return locs

end
