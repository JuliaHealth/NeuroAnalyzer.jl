export na_info
export na_set_use_cuda
export na_set_progress_bar
export na_set_prefs
export na_set_verbose
export na_set_exclude_bads
export na_version

"""
    na_info()

Show NeuroAnalyzer and imported packages versions.

# Arguments

Nothing

# Returns

Nothing
"""
function na_info()::Nothing

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
    println()
    Threads.nthreads() < length(Sys.cpu_info()) || println("For best performance, environment variable `JULIA_NUM_THREADS` ($(Threads.nthreads())) should be less than number of CPU threads ($(length(Sys.cpu_info())))")
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
            "ContinuousWavelets",
            "Crayons",
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
            "FIRLSFilterDesign",
            "FourierTools",
            "FractalDimensions",
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
            "LibSerialPort",
            "LinRegOutliers",
            "Loess",
            "MAT",
            "MLJ",
            "MultivariateStats",
            "NeuroStats",
            "NPZ",
            "PiGPIO",
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
            "WAV",
            "Wavelets",
            "WaveletsExt",
            "XDF" ]
        versions = TOML.parsefile(joinpath(na_pkg, "Manifest.toml"))["deps"]
        for idx in eachindex(required_packages)
            pkg = lpad(required_packages[idx], 25 - length(idx), " ")
            pkg_ver = versions[required_packages[idx]][1]["version"]
            println("$pkg $pkg_ver")
        end
    else
        @warn "Manifest.toml file could not be found in $(na_pkg), cannot report versions of imported packages."
    end

    return nothing

end

"""
    na_set_use_cuda(value)

Change `use_cuda` preference.

# Arguments

- `value::Bool`: value

# Returns

Nothing
"""
function na_set_use_cuda(value::Bool)::Nothing

    @set_preferences!("use_cuda" => value)
    _info("New option value set, restart your Julia session for this change to take effect")

    return nothing

end

"""
    na_set_progress_bar(value)

Change `progress_bar` preference.

# Arguments

- `value::Bool`: value

# Returns

Nothing
"""
function na_set_progress_bar(value::Bool)::Nothing

    @set_preferences!("progress_bar" => value)
    _info("New option value set, restart your Julia session for this change to take effect")

    return nothing

end

"""
    na_set_verbose(value)

Change `verbose` preference.

# Arguments

- `value::Bool`: value

# Returns

Nothing
"""
function na_set_verbose(value::Bool)::Nothing

    @set_preferences!("verbose" => value)
    _info("New option value set, restart your Julia session for this change to take effect")

    return nothing

end

"""
    na_set_prefs(; use_cuda, progress_bar, verbose, scheduler)

Save NeuroAnalyzer preferences.

# Arguments

- `use_cuda::Bool`
- `progress_bar::Bool`
- `verbose::Bool`
- `exclude_bads::Bool`

# Returns

Nothing
"""
function na_set_prefs(; use_cuda::Bool, progress_bar::Bool, verbose::Bool, exclude_bads::Bool)::Nothing

    @set_preferences!("use_cuda" => use_cuda)
    @set_preferences!("progress_bar" => progress_bar)
    @set_preferences!("verbose" => verbose)
    @set_preferences!("exclude_bads" => exclude_bads)

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

    VER = string(Int(NeuroAnalyzer.VER.major)) * "." * string(Int(NeuroAnalyzer.VER.minor)) * "." * string(Int(NeuroAnalyzer.VER.patch))

    return VER

end


"""
    na_set_exclude_bads(value)

Change `exclude_bads` preference.

# Arguments

- `value::Bool`: value

# Returns

Nothing
"""
function na_set_exclude_bads(value::Bool)::Nothing

    @set_preferences!("exclude_bads" => value)
    _info("New option value set, restart your Julia session for this change to take effect")

    return nothing

end