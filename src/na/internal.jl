# suppress console output (warnings, info and debug messages; errors are shown)
_log_off() = Logging.disable_logging(Logging.Warn)
# restore console output
_log_on() = Logging.disable_logging(Logging.BelowMinLevel)

# print a message with a (colored) prefix, unless verbose is off
function _msg(prefix::String, color::Crayons.Crayon, s::String)::Nothing
    verbose || return nothing
    if colors
        println(color("[ $prefix: "), s)
    else
        println("[ $prefix: $s")
    end
    return nothing
end

"""Display info message."""
_info(s::String)::Nothing = _msg("Info", Crayons.Box.LIGHT_CYAN_FG, s)

"""Display warning message."""
_warn(s::String)::Nothing = _msg("Warning", Crayons.Box.YELLOW_FG, s)

"""Display error message."""
_error(s::String)::Nothing = _msg("Error", Crayons.Box.RED_FG, s)

"""Display deprecated function message."""
_deprecated(s::String)::Nothing = _warn("Function $s() is deprecated.")

"""Display deprecated function message. Provide the name of an alternative function."""
_deprecated(s1::String, s2::String)::Nothing =
    _warn("Function $s1() is deprecated, please use $s2() instead.")

"""Display WIP (Work In Progress) message."""
function _wip()::Nothing
    m = "This function has the WIP (Work In Progress) status and is not ready for production use."
    # devel branch: warning; stable branch: error
    return allow_wip ? _warn(m) : _error(m)
end
