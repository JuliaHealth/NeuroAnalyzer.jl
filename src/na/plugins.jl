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

Nothing
"""
function na_plugins_reload()::Nothing

    @assert isdir(plugins_path) "Folder $plugins_path cannot be opened."

    path_tmp = pwd()
    cd(plugins_path)
    plugins = readdir(plugins_path)
    for idx1 in eachindex(plugins)
        plugin = readdir(joinpath(plugins[idx1], "src"))
        for idx2 in eachindex(plugin)
            if splitext(plugin[idx2])[2] == ".jl"
                include(joinpath(plugins_path, plugins[idx1], "src", plugin[idx2]))
                _info(" Loaded: $(plugin[idx2])")
            end
        end
    end

    cd(path_tmp)

    return nothing

end

"""
    na_plugins_list()

List NeuroAnalyzer plugins.

# Arguments

Nothing

# Returns

Nothing
"""
function na_plugins_list()::Nothing

    @assert isdir(plugins_path) "Folder $plugins_path cannot be opened."

    path_tmp = pwd()
    cd(plugins_path)
    plugins = readdir(plugins_path)
    println("Available plugins:")
    for idx in eachindex(plugins)
        println("$idx. $(plugins[idx])")
    end

    cd(path_tmp)

    return nothing

end

"""
    na_plugins_remove(plugin)

Remove NeuroAnalyzer plugin.

# Arguments

- `plugin::String`: plugin name

# Returns

Nothing
"""
function na_plugins_remove(plugin::String)::Nothing

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

    return nothing

end

"""
    na_plugins_install(plugin)

Install NeuroAnalyzer plugin from remote Git repository or from local .TAR.GZ/.ZIP archive (requires unzip or tar command to be available).

# Arguments

- `plugin::String`: plugin Git repository URL or file name (with full path)

# Returns

Nothing
"""
function na_plugins_install(plugin::String)::Nothing

    @assert isdir(plugins_path) "Folder $plugins_path cannot be opened."

    path_tmp = pwd()
    cd(plugins_path)
    if occursin(r"http.*", plugin)
        # install from remote repository
        try
            run(`$(git()) clone $plugin`)
        catch
            @error "Cannot install $plugin."
        end
    else
        # install from local archive
        @assert isfile(plugin) "File $plugin cannot be opened."
        @assert lowercase(splitext(plugin)[2]) in [".zip", ".gz"] "PLUGIN must specify .ZIP/.TAR.GZ file."
        if lowercase(splitext(plugin)[2]) == ".zip"
            Sys.which("unzip") === nothing && (@error "Unknown command: unzip")
            _info("Installing from .ZIP archive")
            try
                run(`unzip -oq $plugin`);
            catch
                @error "Cannot install $plugin."
            end
        elseif lowercase(splitext(plugin)[2]) == ".gz" && lowercase(splitext(splitext(plugin)[1])[2]) == ".tar"
            Sys.which("tar") === nothing && (@error "Unknown command: tar")
            _info("Installing from .TAR.GZ archive")
            try
                run(`tar --overwrite  -xzf $plugin`);
            catch
                @error "Cannot install $plugin."
            end
        end
    end

    na_plugins_reload()
    cd(path_tmp)

    return nothing

end

"""
    na_plugins_update(plugin)

Update NeuroAnalyzer plugin(s).

# Arguments

- `plugin::String`: plugin to update; if empty, update all

# Returns

Nothing
"""
function na_plugins_update(plugin::String="")::Nothing

    @assert isdir(plugins_path) "Folder $plugins_path cannot be opened."

    path_tmp = pwd()
    cd(plugins_path)
    plugins = readdir(plugins_path)
    if plugin == ""
        for idx in eachindex(plugins)
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

    return nothing

end
