export plugins_reload
export plugins_list
export plugins_remove
export plugins_install
export plugins_update

"""
    plugins_reload()

Reload NeuroAnalyzer plugins.

# Arguments

Nothing

# Returns

- `Nothing`
"""
function plugins_reload()::Nothing
    # validate
    isdir(plugins_path) ||
        throw(ArgumentError("Folder $plugins_path cannot be opened."))

    for plugin_name in readdir(plugins_path)
        src_path = joinpath(plugins_path, plugin_name, "src")

        isdir(src_path) || continue  # skip if no src/ folder

        for filename in readdir(src_path)
            if splitext(filename)[2] == ".jl"
                include(joinpath(src_path, filename))
                _info("Loaded: $filename")
            end
        end
    end

    return nothing
end

"""
    plugins_list()

List NeuroAnalyzer plugins.

# Arguments

Nothing

# Returns

- `Nothing`
"""
function plugins_list()::Nothing
    # validate
    isdir(plugins_path) ||
        throw(ArgumentError("Folder $plugins_path cannot be opened."))

    plugins = Base.filter(isdir, readdir(plugins_path; join = true))
    println("Available plugins:")
    for (idx, plugin) in enumerate(basename.(plugins))
        println("$idx. $plugin")
    end

    return nothing
end

"""
    plugins_remove(plugin)

Remove NeuroAnalyzer plugin.

# Arguments

- `plugin::String`: plugin name

# Returns

- `Nothing`
"""
function plugins_remove(plugin::String)::Nothing
    # validate
    isdir(plugins_path) ||
        throw(ArgumentError("Folder $plugins_path cannot be opened."))

    plugins = Base.filter(isdir, readdir(plugins_path; join = true))
    plugin_path = joinpath(plugins_path, plugin)

    plugin_path in plugins ||
        throw(ArgumentError("Plugin $plugin does not exist."))

    _warn("This will remove the whole $plugin directory and all its contents.")

    try
        rm(plugin_path; recursive = true)
        _info("Removed plugin: $plugin")
    catch e
        @error "Cannot remove $plugin directory." exception=e
        return nothing  # abort - don't reload if removal failed
    end

    plugins_reload()
    return nothing
end

"""Helper: install plugin from URL"""
function _install_from_remote(plugin::String)::Nothing
    try
        run(Cmd(`$(git()) clone $plugin`; dir = plugins_path))
    catch e
        throw(ErrorException("Cannot clone $plugin: $e"))
    end
    return nothing
end

"""Helper: install plugin from local file"""
function _install_from_archive(plugin::String)::Nothing
    plugin = abspath(plugin)
    isfile(plugin) ||
        throw(ArgumentError("File $plugin cannot be opened."))
    ext  = lowercase(splitext(plugin)[2])
    ext2 = lowercase(splitext(splitext(plugin)[1])[2])
    if ext == ".zip"
        Sys.which("unzip") === nothing &&
            throw(ErrorException("Required command not found: unzip"))
        _info("Installing from .ZIP archive")
        run(Cmd(`unzip -oq $plugin`; dir = plugins_path))
    elseif ext == ".gz" && ext2 == ".tar"
        Sys.which("tar") === nothing &&
            throw(ErrorException("Required command not found: tar"))
        _info("Installing from .TAR.GZ archive")
        run(Cmd(`tar --overwrite -xzf $plugin`; dir = plugins_path))
    else
        throw(ArgumentError("Plugin must be a .zip or .tar.gz file, got: $plugin"))
    end
    return nothing
end

"""
    plugins_install(plugin)

Install NeuroAnalyzer plugin from remote Git repository or from local .TAR.GZ/.ZIP archive (requires unzip or tar command to be available).

# Arguments

- `plugin::String`: plugin Git repository URL or file name (with full path)

# Returns

- `Nothing`
"""
function plugins_install(plugin::String)::Nothing
    # validate
    isdir(plugins_path) ||
        throw(ArgumentError("Folder $plugins_path cannot be opened."))

    # archives are extracted / repositories cloned into plugins_path
    # (without changing the working directory of the session)
    try
        if startswith(plugin, r"https?://|git@")
            _install_from_remote(plugin)
        else
            _install_from_archive(plugin)
        end
        plugins_reload()
    catch e
        @error "Installation of $plugin failed." exception=e
    end

    return nothing
end

"""Helper: update plugin from a git repository"""
function _update_plugin(plugin_path::String)::Nothing
    # validate
    name = basename(plugin_path)
    isdir(joinpath(plugin_path, ".git")) ||
        (_warn("Skipping $name: not a git repository."); return nothing)
    _info("Updating: $name")
    try
        run(`$(git()) -C $plugin_path pull`)
    catch e
        @error "Cannot update $name." exception=e
    end
    return nothing
end

"""
    plugins_update(plugin)

Update NeuroAnalyzer plugin(s).

# Arguments

- `plugin::String=""`: plugin to update; if empty, update all

# Returns

- `Nothing`
"""
function plugins_update(plugin::String = "")::Nothing
    # validate
    isdir(plugins_path) ||
        throw(ArgumentError("Folder $plugins_path cannot be opened."))

    plugins = Base.filter(isdir, readdir(plugins_path; join = true))

    if isempty(plugin)
        for plugin_path in plugins
            _update_plugin(plugin_path)
        end
    else
        plugin_path = joinpath(plugins_path, plugin)
        plugin_path in plugins ||
            throw(ArgumentError("Plugin $plugin does not exist."))
        _update_plugin(plugin_path)
    end

    plugins_reload()
    return nothing
end
