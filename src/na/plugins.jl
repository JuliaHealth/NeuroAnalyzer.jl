export na_plugins_reload
export na_plugins_list
export na_plugins_remove
export na_plugins_install
export na_plugins_update

"""
    na_plugins_reload()

Reload NeuroAnalyzer plugins.

# Arguments

Nothing

# Returns

- `Nothing`
"""
function na_plugins_reload()::Nothing
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
    na_plugins_list()

List NeuroAnalyzer plugins.

# Arguments

Nothing

# Returns

- `Nothing`
"""
function na_plugins_list()::Nothing
    isdir(plugins_path) ||
        throw(ArgumentError("Folder $plugins_path cannot be opened."))

    plugins = Base.filter(isdir, readdir(plugins_path; join=true))
    println("Available plugins:")
    for (idx, plugin) in enumerate(basename(p))
        println("$idx. $plugin")
    end

    return nothing
end

"""
    na_plugins_remove(plugin)

Remove NeuroAnalyzer plugin.

# Arguments

- `plugin::String`: plugin name

# Returns

- `Nothing`
"""
function na_plugins_remove(plugin::String)::Nothing
    isdir(plugins_path) ||
        throw(ArgumentError("Folder $plugins_path cannot be opened."))

    plugins = Base.filter(isdir, readdir(plugins_path; join=true))
    plugin_path = joinpath(plugins_path, plugin)

    plugin_path in plugins ||
        throw(ArgumentError("Plugin $plugin does not exist."))

    _warn("This will remove the whole $plugin directory and all its contents.")

    try
        rm(plugin_path; recursive=true)
        _info("Removed plugin: $plugin")
    catch e
        @error "Cannot remove $plugin directory." exception=e
        return nothing  # abort - don't reload if removal failed
    end

    na_plugins_reload()
    return nothing
end

"""Helper: install plugin from URL"""
function _install_from_remote(plugin::String)::Nothing
    try
        run(`$(git()) clone $plugin`)
    catch e
        throw(ErrorException("Cannot clone $plugin: $e"))
    end
    return nothing
end

"""Helper: install plugin from local file"""
function _install_from_archive(plugin::String)::Nothing
    plugin = abspath(plugin)  # resolve before cd changes context
    isfile(plugin) ||
        throw(ArgumentError("File $plugin cannot be opened."))
    ext  = lowercase(splitext(plugin)[2])
    ext2 = lowercase(splitext(splitext(plugin)[1])[2])
    if ext == ".zip"
        Sys.which("unzip") === nothing &&
            throw(ErrorException("Required command not found: unzip"))
        _info("Installing from .ZIP archive")
        run(`unzip -oq $plugin`)
    elseif ext == ".gz" && ext2 == ".tar"
        Sys.which("tar") === nothing &&
            throw(ErrorException("Required command not found: tar"))
        _info("Installing from .TAR.GZ archive")
        run(`tar --overwrite -xzf $plugin`)
    else
        throw(ArgumentError("Plugin must be a .zip or .tar.gz file, got: $plugin"))
    end
    return nothing
end

"""
    na_plugins_install(plugin)

Install NeuroAnalyzer plugin from remote Git repository or from local .TAR.GZ/.ZIP archive (requires unzip or tar command to be available).

# Arguments

- `plugin::String`: plugin Git repository URL or file name (with full path)

# Returns

- `Nothing`
"""
function na_plugins_install(plugin::String)::Nothing
    isdir(plugins_path) ||
        throw(ArgumentError("Folder $plugins_path cannot be opened."))

    path_tmp = pwd()
    cd(plugins_path)

    try
        if startswith(plugin, "http")
            _install_from_remote(plugin)
        else
            _install_from_archive(plugin)
        end
        na_plugins_reload()
    catch e
        @error "Installation of $plugin failed." exception=e
    finally
        cd(path_tmp)
    end

    return nothing
end

"""Helper: update plugin from a git repository"""
function _update_plugin(plugin_path::String)::Nothing
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
    na_plugins_update(plugin)

Update NeuroAnalyzer plugin(s).

# Arguments

- `plugin::String`: plugin to update; if empty, update all

# Returns

- `Nothing`
"""
function na_plugins_update(plugin::String = "")::Nothing
    isdir(plugins_path) ||
        throw(ArgumentError("Folder $plugins_path cannot be opened."))

    plugins = Base.filter(isdir, readdir(plugins_path; join=true))

    if isnothing(plugin)
        for plugin_path in plugins
            _update_plugin(plugin_path)
        end
    else
        plugin_path = joinpath(plugins_path, plugin)
        plugin_path in plugins ||
            throw(ArgumentError("Plugin $plugin does not exist."))
        _update_plugin(plugin_path)
    end

    na_plugins_reload()
    return nothing
end