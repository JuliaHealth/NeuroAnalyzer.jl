export na_info
export na_plugins_reload
export na_plugins_list
export na_plugins_remove
export na_plugins_install
export na_plugins_update
export na_set_use_cuda
export na_set_progress_bar
export na_set_prefs
export na_set_verbose
export na_version

"""
    na_info()

Show NeuroAnalyzer and imported packages versions.
"""
function na_info()

    println("     NeuroAnalyzer: $(NeuroAnalyzer.VER)")
    println("NeuroAnalyzer path: $(NeuroAnalyzer.PATH)")
    println("     Julia version: $VERSION")
    if CUDA.functional()
        println("              CUDA: $(CUDA.runtime_version()) (use_cuda = $use_cuda)")
    else
        println("              CUDA: not available (use_cuda = $use_cuda)")
    end
    println("      Plugins path: $plugins_path")
    println("    Resources path: $res_path")
    println(" Show progress bar: $progress_bar")
    println("           Verbose: $verbose")
    println("           Threads: $(Threads.nthreads()) [set using `JULIA_NUM_THREADS` environment variable or Julia --threads command-line option]")
    Threads.nthreads() < length(Sys.cpu_info()) || @info "For best performance, `JULIA_NUM_THREADS` ($(Threads.nthreads())) should be less than number of CPU threads ($(length(Sys.cpu_info())))."
    if "JULIA_COPY_STACKS" in keys(ENV) && ENV["JULIA_COPY_STACKS"] == "1"
        @info "Environment variable `JULIA_COPY_STACKS` is set to 1, multi-threading may not work correctly"
    end
    println()

    if isfile(joinpath(NeuroAnalyzer.PATH, "Manifest.toml"))
        println("Imported packages:")
        required_packages = [
            "Cairo",
            "ColorSchemes",
            "ContinuousWavelets",
            "CSV",
            "CubicSplines",
            "CUDA",
            "DataFrames",
            "Deconvolution",
            "DICOM",
            "Distances",
            "DSP",
            "FFTW",
            "FileIO",
            "FindPeaks1D",
            "FourierTools",
            "GeometryBasics",
            "Git",
            "GLM",
            "Gtk",
            "HypothesisTests",
            "Images",
            "ImageBinarization",
            "ImageMorphology",
            "InformationMeasures",
            "Interpolations",
            "Jacobi",
            "JLD2",
            "JSON",
            "Loess",
            "MAT",
            "MultivariateStats",
            "NPZ",
            "Plots",
            "Polynomials",
            "Preferences",
            "PrettyTables",
            "ProgressMeter",
            "SavitzkyGolay",
            "ScatteredInterpolation",
            "Simpson",
            "StatsFuns",
            "StatsKit",
            "StatsModels",
            "StatsPlots",
            "TimeZones",
            "TOML",
            "Wavelets",
            "WaveletsExt" ]
        versions = TOML.parsefile(joinpath(NeuroAnalyzer.PATH, "Manifest.toml"))["deps"]
        for idx in 1:length(required_packages)
            pkg = lpad(required_packages[idx], 25 - length(idx), " ")
            pkg_ver = versions[required_packages[idx]][1]["version"]
            println("$pkg $pkg_ver")
        end
    else
        @warn "Manifest.toml file could not be found in $(NeuroAnalyzer.PATH), cannot report versions of imported packages."
    end
end

"""
    na_plugins_reload()

Reload NeuroAnalyzer plugins.
"""
function na_plugins_reload()

    @assert isdir(plugins_path) "Folder $plugins_path cannot be opened."

    path_tmp = pwd()
    cd(plugins_path)
    plugins = readdir(plugins_path)
    for idx1 in 1:length(plugins)
        plugin = readdir(joinpath(plugins[idx1], "src"))
        for idx2 in 1:length(plugin)
            if splitext(plugin[idx2])[2] == ".jl"
                include(joinpath(plugins_path, plugins[idx1], "src", plugin[idx2]))
                _info(" Loaded: $(plugin[idx2])")
            end
        end
    end
    cd(path_tmp)

end

"""
    na_plugins_list()

List NeuroAnalyzer plugins.
"""
function na_plugins_list()

    @assert isdir(plugins_path) "Folder $plugins_path cannot be opened."

    path_tmp = pwd()
    cd(plugins_path)
    plugins = readdir(plugins_path)
    println("Available plugins:")
    for idx in 1:length(plugins)
        println("$idx. $(plugins[idx])")
    end
    cd(path_tmp)

end

"""
    na_plugins_remove(plugin)

Remove NeuroAnalyzer `plugin`.

# Arguments

- `plugin::String`: plugin name
"""
function na_plugins_remove(plugin::String)

    _warn("This will remove the whole $plugin directory, along with its file contents.")
    @assert isdir(plugins_path) "Folder $plugins_path cannot be opened."

    path_tmp = pwd()
    cd(plugins_path)
    plugins = readdir(plugins_path)
    @assert plugin in plugins "Plugin $plugin cannot be loaded."
    try
        rm(plugin, recursive=true)
    catch
        @error "Cannot remove $plugin directory."
    end
    na_plugins_reload()
    cd(path_tmp)

end

"""
    na_plugins_install(plugin)

Install NeuroAnalyzer `plugin`.

# Arguments

- `plugin::String`: plugin Git repository URL
"""
function na_plugins_install(plugin::String)

    @assert isdir(plugins_path) "Folder $plugins_path cannot be opened."

    path_tmp = pwd()
    cd(plugins_path)
    try
        run(`$(git()) clone $plugin`)
    catch
        @error "Cannot install $plugin."
    end
    na_plugins_reload()
    cd(path_tmp)

end

"""
    na_plugins_update(plugin)

Install NeuroAnalyzer `plugin`.

# Arguments

- `plugin::String`: plugin to update; if empty, update all
"""
function na_plugins_update(plugin::Union{String, Nothing}=nothing)

    @assert isdir(plugins_path) "Folder $plugins_path cannot be opened."

    path_tmp = pwd()
    cd(plugins_path)
    plugins = readdir(plugins_path)
    if plugin === nothing
        for idx in 1:length(plugins)
            cd(plugins[idx])
            @info "Updating: $(plugins[idx])"
            try
                run(`$(git()) pull`)
            catch
                @error "Cannot update $(plugins[idx])."
            end
            cd(plugins_path)
        end
    else
        @assert plugin in plugins "Plugin $plugin cannot be loaded."
        cd(plugin)
        try
            run(`$(git()) pull`)
        catch
            @error "Cannot update $plugin."
        end
        cd(plugins_path)
    end
    na_plugins_reload()
    cd(path_tmp)

end

"""
    na_set_use_cuda(use_cuda)

Change `use_cuda` preference.

# Arguments

- `use_cuda::Bool`: value
"""
function na_set_use_cuda(use_cuda::Bool)

    @set_preferences!("use_cuda" => use_cuda)
    _info("New option value set, restart your Julia session for this change to take effect")

end

"""
    na_set_progress_bar(progress_bar)

Change `progress_bar` preference.

# Arguments

- `progress_bar::Bool`: value
"""
function na_set_progress_bar(progress_bar::Bool)

    @set_preferences!("progress_bar" => progress_bar)
    _info("New option value set, restart your Julia session for this change to take effect")

end

"""
    na_set_prefs(use_cuda, plugins_path, progress_bar, verbose)

Save NeuroAnalyzer preferences.

# Arguments

- `use_cuda::Bool`
- `progress_bar::Bool`
- `verbose::Bool`
"""
function na_set_prefs(; use_cuda::Bool, progress_bar::Bool, verbose::Bool)

    @set_preferences!("use_cuda" => use_cuda)
    @set_preferences!("progress_bar" => progress_bar)
    @set_preferences!("verbose" => verbose)

end

"""
    na_set_verbose(verbose)

Change `verbose` preference.

# Arguments

- `verbose::Bool`: value
"""
function na_set_verbose(verbose::Bool)

    @set_preferences!("verbose" => verbose)
    _info("New option value set, restart your Julia session for this change to take effect")

end

"""
    na_version()

Convert NeuroAnalyzer version to string.

# Returns

- `VER::String`
"""
function na_version()

    return string(Int(NeuroAnalyzer.VER.major)) * "." * string(Int(NeuroAnalyzer.VER.minor)) * "." * string(Int(NeuroAnalyzer.VER.patch))

end