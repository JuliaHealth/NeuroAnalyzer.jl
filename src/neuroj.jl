"""
    neuroj_version()

Shows NeuroJ and imported packages versions.
"""
function neuroj_version()
    m = Pkg.Operations.Context().env.manifest
    println("NeuroJ version: $(m[findfirst(v->v.name=="NeuroJ", m)].version)")
    println("Imported packages:")
    required_packages = ["CSV",
                         "DataFrames",
                         "Distances",
                         "DSP",
                         "FFTW",
                         "HypothesisTests",
                         "InformationMeasures",
                         "Interpolations",
                         "JLD2",
                         "LinearAlgebra",
                         "MultivariateStats",
                         "Plots",
                         "Polynomials",
                         "ScatteredInterpolation",
                         "Simpson",
                         "StatsKit",
                         "StatsPlots"]
    for idx in 1:length(required_packages)
        pkg = lpad(required_packages[idx], 25 - length(idx), " ")
        pkg_ver = m[findfirst(v->v.name==required_packages[idx], m)].version
        println("$pkg $pkg_ver")
    end
end

"""
    neuroj_reload_plugins()

Reload NeuroJ plugins. Plugins path is: ~/Documents/NeuroJ/plugins/
"""
function neuroj_reload_plugins()
    isdir(expanduser("~/Documents/NeuroJ/plugins/")) || mkpath(expanduser("~/Documents/NeuroJ/plugins/"))
    plugins_path=expanduser("~/Documents/NeuroJ/plugins/")
    cd(plugins_path)
    plugins = readdir(plugins_path)
    for f in plugins
        if splitext(f)[2] == ".jl"
            include(plugins_path * f)
        end
    end
end