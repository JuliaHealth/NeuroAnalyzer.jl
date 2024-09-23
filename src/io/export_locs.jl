export export_locs

"""
    export_locs(obj; <keyword arguments>)

Export channel locations data, format is based on `file_name` extension (.csv, .ced, .locs or .tsv)

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `file_name::String`
- `overwrite::Bool=false`

# Returns

Nothing
"""
function export_locs(obj::NeuroAnalyzer.NEURO; file_name::String, overwrite::Bool=false)::Nothing

    @assert !(isfile(file_name) && !overwrite) "File $file_name cannot be saved, to overwrite use overwrite=true."

    channels = _find_bylabel(obj.locs, obj.locs[!, :label])
    clabels = obj.locs[!, :label]
    theta = obj.locs[!, :loc_theta]
    radius = obj.locs[!, :loc_radius]
    x = obj.locs[!, :loc_x]
    y = obj.locs[!, :loc_y]
    z = obj.locs[!, :loc_z]
    radius_sph = obj.locs[!, :loc_radius_sph]
    theta_sph = obj.locs[!, :loc_theta_sph]
    phi_sph = obj.locs[!, :loc_phi_sph]

    if splitext(file_name)[2] == ".ced"
        df = DataFrame(Number=channels, labels=clabels, theta=theta, radius=radius, X=x, Y=y, Z=z, sph_theta=theta_sph, sph_phi=phi_sph, sph_radius=radius_sph)
        CSV.write(file_name, df, delim="\t", header=true)
    elseif splitext(file_name)[2] == ".locs"
        df = DataFrame(Number=channels, theta=theta, radius=radius, labels=clabels)
        CSV.write(file_name, df, delim="\t", header=false)
    elseif splitext(file_name)[2] == ".tsv"
        df = DataFrame(labels=clabels, x=x, y=y, z=z, theta=theta, radius=radius, radius_sph=radius_sph, theta_sph=theta_sph, phi_sph=phi_sph)
        CSV.write(file_name, df, delim="\t", header=true)
    elseif splitext(file_name)[2] == ".csv"
        CSV.write(file_name, obj.locs, header=true)
    else
        @error "$file_name format must be .csv, .ced, .locs or .tsv."
    end

    return nothing

end

"""
    export_locs(locs; <keyword arguments>)

Export channel locations, format is based on `file_name` extension (.ced, .locs, .tsv)

# Arguments

- `locs::DataFrame`
- `file_name::String`
- `overwrite::Bool=false`

# Returns

Nothing
"""
function export_locs(locs::DataFrame; file_name::String, overwrite::Bool=false)::Nothing

    @assert !(isfile(file_name) && !overwrite) "File $file_name cannot be saved, to overwrite use overwrite=true."

    channels = _find_bylabel(obj.locs, obj.locs[!, :label])
    clabels = locs[!, :label]
    theta = locs[!, :loc_theta]
    radius = locs[!, :loc_radius]
    x = locs[!, :loc_x]
    y = locs[!, :loc_y]
    z = locs[!, :loc_z]
    radius_sph = locs[!, :loc_radius_sph]
    theta_sph = locs[!, :loc_theta_sph]
    phi_sph = locs[!, :loc_phi_sph]

    if splitext(file_name)[2] == ".ced"
        df = DataFrame(Number=channels, labels=clabels, theta=theta, radius=radius, X=x, Y=y, Z=z, sph_theta=theta_sph, sph_phi=phi_sph, sph_radius=radius_sph, head=true)
        CSV.write(file_name, df, delim="\t")
    elseif splitext(file_name)[2] == ".locs"
        df = DataFrame(Number=channels, theta=theta, radius=radius, labels=clabels)
        CSV.write(file_name, df, delim="\t", header=false)
    elseif splitext(file_name)[2] == ".tsv"
        df = DataFrame(labels=clabels, x=x, y=y, z=z, theta=theta, radius=radius, radius_sph=radius_sph, theta_sph=theta_sph, phi_sph=phi_sph)
        CSV.write(file_name, df, delim="\t", header=true)
    else
        @error "file_name format must be .ced, .locs or .tsv."
    end

    return nothing

end
