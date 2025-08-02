# suppress console output
_log_off() = Logging.disable_logging(Logging.Warn)
# restore console output
_log_on() = Logging.disable_logging(Logging.Debug)

function _info(s::String)::Nothing
    if verbose
        if colors
            println(LIGHT_CYAN_FG("[ Info: "), s)
        else
            println("[ Info: $s")
        end
    end
    return nothing
end
function _warn(s::String)::Nothing
    if verbose
        if colors
            println(YELLOW_FG("[ Warning: "), s)
        else
            println("[ Warning: $s")
        end
    end
    return nothing
end
function _deprecated(s::String)::Nothing
    if verbose
        if colors
            println(RED_FG("[ Error: "), "Function $s() is deprecated.")
        else
            println("[ Error: Function $s() is deprecated.")
        end
    end
    return nothing
end
function _deprecated(s1::String, s2::String)::Nothing
    if verbose
        if colors
            println(RED_FG("[ Error: "), "Function $s1() is deprecated, please use $s2() instead.")
        else
            println("[ Error: Function $s1() is deprecated, please use $s2() instead.")
        end
    end
    return nothing
end
function _wip()::Nothing
    if verbose
        if allow_wip
            if colors
                println(YELLOW_FG("[ Warning: "), "This function has the WIP (Work In Progress) status and is not ready for production use.")
            else
                println("[ Warning: This function has the WIP (Work In Progress) status and is not ready for production use.")
            end
        else
            if colors
                println(RED_FG("[ Error: "), "This function has the WIP (Work In Progress) status and is not ready for production use.")
            else
                println("[ Error: This function has the WIP (Work In Progress) status and is not ready for production use.")
            end
        end
    end
    return nothing
end

function _load_functions(f::String)
    @assert isdir("src/$f") "Directory src/$f does not exist."
    files = readdir("src/$f")
    if length(files) > 0
        _info("Loading sub-module: $f")
        for idx in files
            include("$f/$idx")
        end
    end
    return nothing
end
