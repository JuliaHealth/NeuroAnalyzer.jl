"""
    neuroj_version()

Show NeuroJ and imported packages versions.
"""
function neuroj_version()
    m = Pkg.Operations.Context().env.manifest
    println("   NeuroJ version: $(m[findfirst(v->v.name=="NeuroJ", m)].version)")
    println("          Threads: $(Threads.nthreads()) [set using using the `JULIA_NUM_THREADS` environment variable]")
    println("Imported packages:")
    required_packages = ["ColorSchemes",
                         "CSV",
                         "CubicSplines",
                         "DataFrames",
                         "Distances",
                         "DSP",
                         "FFTW",
                         "FindPeaks1D",
                         "HypothesisTests",
                         "InformationMeasures",
                         "Interpolations",
                         "JLD2",
                         "Loess",
                         "MultivariateStats",
                         "Plots",
                         "Polynomials",
                         "ScatteredInterpolation",
                         "Simpson",
                         "StatsFuns",
                         "StatsKit",
                         "StatsModels",
                         "StatsPlots",
                         "Wavelets"]
    for idx in 1:length(required_packages)
        pkg = lpad(required_packages[idx], 25 - length(idx), " ")
        pkg_ver = m[findfirst(v->v.name==required_packages[idx], m)].version
        println("$pkg $pkg_ver ")
    end
end

"""
    neuroj_plugins_reload()

Reload NeuroJ plugins. Plugins path is: `~/Documents/NeuroJ/plugins/`.
"""
function neuroj_plugins_reload()
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

"""
    neuroj_plugins_list()

List NeuroJ plugins. Plugins path is: `~/Documents/NeuroJ/plugins/`.
"""
function neuroj_plugins_list()
    isdir(expanduser("~/Documents/NeuroJ/plugins/")) || mkpath(expanduser("~/Documents/NeuroJ/plugins/"))
    plugins_path=expanduser("~/Documents/NeuroJ/plugins/")
    cd(plugins_path)
    plugins = readdir(plugins_path)
    for idx in 1:length(plugins)
        println("$idx. $(replace(plugins[idx], ".jl" => "()"))")
    end
end