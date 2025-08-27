export zipratio

"""
    zipratio(obj)

Calculate zip ratio for the object data (ratio of zip-compressed to uncompressed object data).

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `zip_ratio::Array{Float64, 3}`

# Notes

Zip ratio is an intermediate marker of signal complexity (with lower values indicating lower complexity).

`zipratio()` requires `zip` (Linux/Mac) or `zip.exe` (Windows) command in the `PATH`.
"""
function zipratio(obj::NeuroAnalyzer.NEURO)

    tmp_name, _ = mktemp()
    zip_name = tmp_name * ".zip"
    tmp_name *= ".csv"
    export_csv(obj, file_name=tmp_name, names=false, header=false, epoch_time=false, components=false, markers=false, locs=false, history=false, overwrite=true)

    zip_cmd = ""
    if Sys.iswindows()
        zip_cmd = "zip.exe"
    elseif Sys.islinux()
        zip_cmd = "zip"
    elseif Sys.isapple()
        zip_cmd = "zip"
    end
    Sys.which(zip_cmd) === nothing && (@error "Unknown command: $zip_cmd")

    _info("Attempting to compress the object")
    run(`$zip_cmd -0 -q $zip_name $tmp_name`)
    zip_size_0 = filesize(zip_name)
    rm(zip_name)
    run(`$zip_cmd -9 -q $zip_name $tmp_name`)
    zip_size_9 = filesize(zip_name)
    rm(zip_name)

    zip_ratio = zip_size_9 / zip_size_0

    return zip_ratio

end
