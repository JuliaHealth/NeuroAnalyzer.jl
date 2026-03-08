export na_info
export na_set_progress_bar
export na_set_prefs
export na_set_verbose
export na_set_colors
export na_set_exclude_bads
export na_version

"""
    na_info()

Show NeuroAnalyzer and imported packages versions.

# Arguments

Nothing

# Returns

  - `Nothing`
"""
function na_info()::Nothing

    println("     NeuroAnalyzer: $(NeuroAnalyzer.VER)")
    println("NeuroAnalyzer path: $(NeuroAnalyzer.PATH)")
    println("     Julia version: $VERSION")
    println("      Plugins path: $plugins_path")
    println("    Resources path: $res_path")
    println(" Show progress bar: $progress_bar")
    println("           Verbose: $(NeuroAnalyzer.verbose)")
    println("      Exclude bads: $(NeuroAnalyzer.exclude_bads)")
    println("            Colors: $(NeuroAnalyzer.colors)")
    println(
        "           Threads: $(Threads.nthreads()) [set using `JULIA_NUM_THREADS` environment variable or Julia --threads command-line option]",
    )
    println()
    Threads.nthreads() < length(Sys.cpu_info()) || println(
        "For best performance, environment variable `JULIA_NUM_THREADS` ($(Threads.nthreads())) should be less than number of CPU threads ($(length(Sys.cpu_info())))",
    )
    if "JULIA_COPY_STACKS" in keys(ENV) && ENV["JULIA_COPY_STACKS"] == "1"
        println("Environment variable `JULIA_COPY_STACKS` is set to 1, multi-threading may not work correctly")
    end
    println()

    na_pkg = dirname(dirname(Base.find_package("NeuroAnalyzer")))

    if isfile(joinpath(na_pkg, "Manifest.toml"))
        println("Imported packages:")
        required_packages = [
            "Cairo",
            "ColorSchemes",
            "ComplexityMeasures",
            "ContinuousWavelets",
            "Crayons",
            "CSV",
            "DataFrames",
            "Deconvolution",
            "DICOM",
            "Dierckx",
            "Distances",
            "DSP",
            "Einsum",
            "FFTW",
            "FileIO",
            "FindPeaks1D",
            "FourierTools",
            "FractalDimensions",
            "GeometryBasics",
            "Git",
            "GLM",
            "GLMakie",
            "GR",
            "Gtk4",
            "Hurst",
            "HypothesisTests",
            "Images",
            "ImageBinarization",
            "ImageMorphology",
            "InformationMeasures",
            "Interpolations",
            "Jacobi",
            "JLD2",
            "JSON",
            "LibSerialPort",
            "LinRegOutliers",
            "Loess",
            "MAT",
            "MLJ",
            "MultivariateStats",
            "NPZ",
            "PiGPIO",
            "Plots",
            "Polynomials",
            "Preferences",
            "PrettyTables",
            "ProgressMeter",
            "SavitzkyGolay",
            "ScatteredInterpolation",
            "StatsFuns",
            "StatsKit",
            "StatsModels",
            "StatsPlots",
            "Suppressor",
            "TimeZones",
            "TOML",
            "WAV",
            "Wavelets",
            "WaveletsExt",
            "XDF",
        ]
        versions = TOML.parsefile(joinpath(na_pkg, "Manifest.toml"))["deps"]
        for idx in eachindex(required_packages)
            pkg = lpad(required_packages[idx], 25 - length(idx), " ")
            pkg_ver = versions[required_packages[idx]][1]["version"]
            println("$pkg $pkg_ver")
        end
    else
        _warn("Manifest.toml file could not be found in $(na_pkg), cannot report versions of imported packages.")
    end

    return nothing

end

"""
    na_set_progress_bar(value)

Change `progress_bar` preference.

# Arguments

  - `value::Bool`: value

# Returns

  - `Nothing`
"""
function na_set_progress_bar(value::Bool)::Nothing

    NeuroAnalyzer.progress_bar = value
    @set_preferences!("progress_bar" => value)

    return nothing

end

"""
    na_set_verbose(value)

Change `verbose` preference.

# Arguments

  - `value::Bool`: value

# Returns

  - `Nothing`
"""
function na_set_verbose(value::Bool)::Nothing

    NeuroAnalyzer.verbose = value
    @set_preferences!("verbose" => value)

    return nothing

end

"""
    na_set_colors(value)

Change `colors` preference.

# Arguments

  - `value::Bool`: value

# Returns

  - `Nothing`
"""
function na_set_colors(value::Bool)::Nothing

    NeuroAnalyzer.colors = value
    @set_preferences!("colors" => value)

    return nothing

end

"""
    na_set_prefs(; <keyword arguments>)

Set and save NeuroAnalyzer preferences.

# Arguments

  - `progress_bar::Bool`
  - `verbose::Bool`
  - `exclude_bads::Bool`
  - `colors::Bool`

# Returns

  - `Nothing`
"""
function na_set_prefs(; progress_bar::Bool, verbose::Bool, exclude_bads::Bool, colors::Bool)::Nothing

    @set_preferences!("progress_bar" => progress_bar)
    @set_preferences!("verbose" => verbose)
    @set_preferences!("exclude_bads" => exclude_bads)
    @set_preferences!("colors" => colors)

    return nothing

end

"""
    na_version()

Convert NeuroAnalyzer version to string.

# Arguments

Nothing

# Returns

  - `VER::String`
"""
function na_version()::String

    VER =
        string(Int(NeuroAnalyzer.VER.major)) *
        "." *
        string(Int(NeuroAnalyzer.VER.minor)) *
        "." *
        string(Int(NeuroAnalyzer.VER.patch))

    return VER

end

"""
    na_set_exclude_bads(value)

Change `exclude_bads` preference.

# Arguments

  - `value::Bool`: value

# Returns

  - `Nothing`
"""
function na_set_exclude_bads(value::Bool)::Nothing

    NeuroAnalyzer.exclude_bads = value
    @set_preferences!("exclude_bads" => value)

    return nothing

end
