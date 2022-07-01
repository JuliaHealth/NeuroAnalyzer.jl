"""
    neuroj_version()

Show NeuroJ and imported packages versions.
"""
function neuroj_version()
    m = Pkg.Operations.Context().env.manifest
    println("   NeuroJ version: $(m[findfirst(v->v.name=="NeuroJ", m)].version)")
    println("     Plugins path: $(expanduser("~/Documents/NeuroJ/plugins/"))")
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
    plugins_folders = readdir(plugins_path)
    for idx in 1:length(plugins_folders)
        for f in readdir(plugins_folders[idx])
            if splitext(f)[2] == ".jl"
                include(plugins_folders[idx] * "/" * f)
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
        println("$idx. $(replace(plugins[idx], ".jl" => "()"))")
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
end

"""
    neuroj_plugins_install(plugin)

Install NeuroJ plugin.

# Attributes

- `plugin::String`: plugin URL
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
end

"""
    neuroj_plugins_update(plugin)

Install NeuroJ plugin.

# Attributes

- `plugin::String`: plugin to update; if empty, update all
"""
function neuroj_plugins_update(plugin::Union{String, Nothing})
    isdir(expanduser("~/Documents/NeuroJ/plugins/")) || mkpath(expanduser("~/Documents/NeuroJ/plugins/"))
    plugins_path=expanduser("~/Documents/NeuroJ/plugins/")
    cd(plugins_path)
    plugins = readdir(plugins_path)
    plugin in plugins || throw(ArgumentError("Plugin $plugin does not exist."))
    if plugin === nothing
        for idx in 1:length(plugins_folders)
            cd(plugins[idx])
            try
                run(`$(git()) pull`)
            catch err
                @warn "Cannot update $(plugins[idx])."
            end
        end
    else
        cd(plugin)
        try
            run(`$(git()) pull`)
        catch err
            @warn "Cannot update $plugin."
        end
    end
end