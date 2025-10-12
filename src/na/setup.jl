export na_info
export na_set_use_gpu
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

Nothing
"""
function na_info()::Nothing

    println("     NeuroAnalyzer: $(NeuroAnalyzer.VER)")
    println("NeuroAnalyzer path: $(NeuroAnalyzer.PATH)")
    println("     Julia version: $VERSION")
    println("      Plugins path: $plugins_path")
    println("    Resources path: $res_path")
    println(" Show progress bar: $progress_bar")
    println("           Use GPU: $use_gpu")
    if use_gpu
        if na_gpu === :cuda
#            println("          CUDA GPU: $(CUDA.runtime_version())")
            println("              CUDA: $(CUDA.runtime_version())")
        elseif na_gpu === :amdgpu
            println("           AMD GPU: $(AMDGPU.HIP.name(AMDGPU.device()))")
            println("        AMDGPU HIP: $(string(AMDGPU.HIP.runtime_version()))")
            println("     AMDGPU rocFFT: $(string(AMDGPU.rocFFT.version()))")
        end
    end
    println("           Verbose: $verbose")
    println("      Exclude bads: $exclude_bads")
    println("            Colors: $colors")
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
            "AMDGPU",
            "Cairo",
            "ColorSchemes",
            "ComplexityMeasures",
            "ContinuousWavelets",
            "Crayons",
            "CSV",
            "CUDA",
            "DataFrames",
            "Deconvolution",
            "DICOM",
            "Dierckx",
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
            "GLMakie",
            "GR",
            "Gtk",
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
        _warn("Manifest.toml file could not be found in $(na_pkg), cannot report versions of imported packages.")
    end

    return nothing

end

"""
    na_set_use_gpu(value)

Change `use_gpu` preference.

# Arguments

- `value::Bool`: value

# Returns

Nothing
"""
function na_set_use_gpu(value::Bool)::Nothing

    @set_preferences!("use_gpu" => value)
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
    na_set_colors(value)

Change `colors` preference.

# Arguments

- `value::Bool`: value

# Returns

Nothing
"""
function na_set_colors(value::Bool)::Nothing

    @set_preferences!("colors" => value)
    _info("New option value set, restart your Julia session for this change to take effect")

    return nothing

end

"""
    na_set_prefs(; use_gpu, progress_bar, verbose, exclude_bads, colors)

Set and save NeuroAnalyzer preferences.

# Arguments

- `use_gpu::Bool`
- `progress_bar::Bool`
- `verbose::Bool`
- `exclude_bads::Bool`
- `colors::Bool`

# Returns

Nothing
"""
function na_set_prefs(; use_gpu::Bool, progress_bar::Bool, verbose::Bool, exclude_bads::Bool, colors::Bool)::Nothing

    @set_preferences!("use_gpu" => use_gpu)
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