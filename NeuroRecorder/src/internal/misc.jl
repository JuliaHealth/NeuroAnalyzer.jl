_wip() = allow_wip == true ? (@warn "This function has the WIP (Work In Progress) status and is not ready for production use.") : (@error "This function has the WIP (Work In Progress) status and is not ready for production use.")

function _beep()
    beep, fs = wavread(joinpath(res_path, "beep.wav"))
    wavplay(beep, fs)
end

function _checkrpi()
    if Sys.which("pigpiod") === nothing
        return false
    else
        return true
    end
end

function _getch(c::Channel)
    t = REPL.TerminalMenus.terminal
    while true
        REPL.Terminals.raw!(t, true) || error("unable to switch to raw mode")
        keypress = Char(REPL.TerminalMenus.readkey(t.in_stream))
        REPL.Terminals.raw!(t, false) || error("unable to switch back from raw mode")
        put!(c, keypress)
    end
end
