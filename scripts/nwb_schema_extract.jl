# Extract the NWB schema cached in an NWB file written by pynwb (/specifications) into
# nwb_schema/<namespace>/<version>/<name>.json, used by export_nwb() to cache the schema in
# exported files. Run from the repository root:
#
#   julia --project=. scripts/nwb_schema_extract.jl <file.nwb>
#
# Update the licence files of nwb_schema/ from pynwb (nwb-schema/license.txt) and hdmf
# (common/hdmf-common-schema/license.txt) if they changed.
using NeuroAnalyzer
H = NeuroAnalyzer.HDF5
out = joinpath(pwd(), "nwb_schema")
H.h5open(ARGS[1], "r") do f
    sp = f["specifications"]
    for ns in keys(sp), v in keys(sp[ns])
        d = joinpath(out, ns, v)
        mkpath(d)
        for k in keys(sp[ns][v])
            write(joinpath(d, k * ".json"), read(sp[ns][v][k]))
        end
        println(ns, " ", v, ": ", length(keys(sp[ns][v])), " files")
    end
end
