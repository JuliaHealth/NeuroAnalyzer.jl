"""
    neuroj_version()

Show NeuroJ and imported packages versions.
"""
function neuroj_version()
    m = Pkg.Operations.Context().env.manifest
    println("    Julia version: $VERSION")
    println("   NeuroJ version: $(m[findfirst(v->v.name=="NeuroJ", m)].version)")
    if CUDA.functional()
        println("     CUDA version: $(CUDA.version())")
    else
        println("     CUDA version: CUDA not available")
    end
    println("     Plugins path: $(expanduser("~/Documents/NeuroJ/plugins/"))")
    println("          Threads: $(Threads.nthreads()) [set using using the `JULIA_NUM_THREADS` environment variable]")
    println("Imported packages:")
    required_packages = ["ColorSchemes",
                         "CSV",
                         "CubicSplines",
                         "CUDA",
                         "DataFrames",
                         "Distances",
                         "DSP",
                         "FFTW",
                         "FindPeaks1D",
                         "Git",
                         "GLM",
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
    for idx1 in 1:length(plugins)
        plugin = readdir(plugins[idx1] * "/src/")
        for idx2 in 1:length(plugin)
            if splitext(plugin[idx2])[2] == ".jl"
                include(plugins_path * plugins[idx1] * "/src/" * plugin[idx2])
            end
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
        println("$idx. $(replace(plugins[idx]))")
    end
end

"""
    neuroj_plugins_remove(plugin)

Remove NeuroJ plugin.

# Attributes

- `plugin::String`: plugin name
"""
function neuroj_plugins_remove(plugin::String)
    @warn "This will remove the whole $plugin directory, along with its file contents."
    isdir(expanduser("~/Documents/NeuroJ/plugins/")) || mkpath(expanduser("~/Documents/NeuroJ/plugins/"))
    plugins_path=expanduser("~/Documents/NeuroJ/plugins/")
    cd(plugins_path)
    plugins = readdir(plugins_path)
    plugin in plugins || throw(ArgumentError("Plugin $plugin does not exist."))
    try
        rm(plugin, recursive=true)
    catch err
        @warn "Cannot remove $plugin directory."
    end
    neuroj_plugins_reload()
end

"""
    neuroj_plugins_install(plugin)

Install NeuroJ plugin.

# Attributes

- `plugin::String`: plugin Git repository URL
"""
function neuroj_plugins_install(plugin::String)
    isdir(expanduser("~/Documents/NeuroJ/plugins/")) || mkpath(expanduser("~/Documents/NeuroJ/plugins/"))
    plugins_path=expanduser("~/Documents/NeuroJ/plugins/")
    cd(plugins_path)
    try
        run(`$(git()) clone $plugin`)
    catch err
        @warn "Cannot install $plugin."
    end
    neuroj_plugins_reload()
end

"""
    neuroj_plugins_update(plugin)

Install NeuroJ plugin.

# Attributes

- `plugin::String`: plugin to update; if empty, update all
"""
function neuroj_plugins_update(plugin::Union{String, Nothing}=nothing)
    isdir(expanduser("~/Documents/NeuroJ/plugins/")) || mkpath(expanduser("~/Documents/NeuroJ/plugins/"))
    plugins_path=expanduser("~/Documents/NeuroJ/plugins/")
    cd(plugins_path)
    plugins = readdir(plugins_path)
    if plugin === nothing
        for idx in 1:length(plugins)
            cd(plugins[idx])
            println(plugins[idx])
            try
                run(`$(git()) pull`)
            catch err
                @warn "Cannot update $(plugins[idx])."
            end
            cd(plugins_path)
        end
    else
        plugin in plugins || throw(ArgumentError("Plugin $plugin does not exist."))
        cd(plugin)
        try
            run(`$(git()) pull`)
        catch err
            @warn "Cannot update $plugin."
        end
        cd(plugins_path)
    end
    neuroj_plugins_reload()
end