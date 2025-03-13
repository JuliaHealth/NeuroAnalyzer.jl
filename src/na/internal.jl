# suppress console output
_log_off() = Logging.disable_logging(Logging.Warn)
# restore console output
_log_on() = Logging.disable_logging(Logging.Debug)

_info(s::String) = verbose && println(LIGHT_CYAN_FG("[ Info: "), s)
_warn(s::String) = verbose && println(YELLOW_FG("[ Warning: "), s)
_deprecated(s::String) = verbose && println(RED_FG("[ Error: "), "Function $s() is deprecated.")
_deprecated(s1::String, s2::String) = verbose && println(RED_FG("[ Error: "), "Function $s1() is deprecated, please use $s2() instead.")
_wip() = allow_wip ? println(YELLOW_FG("[ Warning: "), "This function has the WIP (Work In Progress) status and is not ready for production use.") : println(RED_FG("[ Error: "), "This function has the WIP (Work In Progress) status and is not ready for production use.")

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
