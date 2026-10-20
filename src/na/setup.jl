export info
export set_progress_bar
export set_prefs
export set_verbose
export set_colors
export set_exclude_bads
export version

"""Helper: print NA settings"""
function _print_settings()::Nothing
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
        "           Threads: $(Threads.nthreads()) [set using `JULIA_NUM_THREADS` " *
        "environment variable or Julia --threads command-line option]",
    )
    return nothing
end

"""Helper: print threads warning"""
function _print_thread_warnings()::Nothing
    cpu_threads = length(Sys.cpu_info())
    n_threads   = Threads.nthreads()
    if n_threads >= cpu_threads
        println(
            "For best performance, `JULIA_NUM_THREADS` ($n_threads) should be " *
            "less than the number of CPU threads ($cpu_threads).",
        )
    end
    if get(ENV, "JULIA_COPY_STACKS", "") == "1"
        println(
            "Environment variable `JULIA_COPY_STACKS` is set to 1; " *
            "multi-threading may not work correctly.",
        )
    end
    return nothing
end

"""Helper: print package versions"""
function _print_package_versions()::Nothing
    # versions of the packages actually loaded (read from their own Project.toml),
    # independent of the environment NeuroAnalyzer was installed in
    println("Imported packages:")
    for pkg in _NA_REQUIRED_PACKAGES
        m = getfield(NeuroAnalyzer, Symbol(pkg))
        ver = something(pkgversion(m), "stdlib")
        println("$(lpad(pkg, 25)) $ver")
    end
    return nothing
end

"""
    info()

Show NeuroAnalyzer and imported packages versions.

# Arguments

Nothing

# Returns

- `Nothing`
"""
function info()::Nothing
    _print_settings()
    _print_thread_warnings()
    println()
    _print_package_versions()
    return nothing
end

"""
    set_progress_bar(value)

Change `progress_bar` preference.

# Arguments

- `value::Bool`: value

# Returns

- `Nothing`
"""
function set_progress_bar(value::Bool)::Nothing
    NeuroAnalyzer.progress_bar = value
    @set_preferences!("progress_bar" => value)

    return nothing
end

"""
    set_verbose(value)

Change `verbose` preference.

# Arguments

- `value::Bool`: value

# Returns

- `Nothing`
"""
function set_verbose(value::Bool)::Nothing
    NeuroAnalyzer.verbose = value
    @set_preferences!("verbose" => value)

    return nothing
end

"""
    set_colors(value)

Change `colors` preference.

# Arguments

- `value::Bool`: value

# Returns

- `Nothing`
"""
function set_colors(value::Bool)::Nothing
    NeuroAnalyzer.colors = value
    @set_preferences!("colors" => value)

    return nothing
end

"""
    set_prefs(; <keyword arguments>)

Set and save NeuroAnalyzer preferences.

# Arguments

- `progress_bar::Bool`
- `verbose::Bool`
- `exclude_bads::Bool`
- `colors::Bool`

# Returns

- `Nothing`
"""
function set_prefs(;
    progress_bar::Bool,
    verbose::Bool,
    exclude_bads::Bool,
    colors::Bool,
)::Nothing
    set_progress_bar(progress_bar)
    set_verbose(verbose)
    set_exclude_bads(exclude_bads)
    set_colors(colors)

    return nothing
end

"""
    version()

Convert NeuroAnalyzer version to string.

# Arguments

Nothing

# Returns

- `VER::String`
"""
function version()::String
    VER =
        string(Int(NeuroAnalyzer.VER.major)) *
        "." *
        string(Int(NeuroAnalyzer.VER.minor)) *
        "." *
        string(Int(NeuroAnalyzer.VER.patch))

    return VER
end

"""
    set_exclude_bads(value)

Change `exclude_bads` preference.

# Arguments

- `value::Bool`: value

# Returns

- `Nothing`
"""
function set_exclude_bads(value::Bool)::Nothing
    NeuroAnalyzer.exclude_bads = value
    @set_preferences!("exclude_bads" => value)

    return nothing
end
