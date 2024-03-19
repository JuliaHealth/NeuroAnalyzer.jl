_wip() = allow_wip == true ? (@warn "This function has the WIP (Work In Progress) status and is not ready for production use.") : (@error "This function has the WIP (Work In Progress) status and is not ready for production use.")

function _beep()
    beep, fs = wavread(joinpath(res_path, "beep.wav"))
    wavplay(beep, fs)
end

function _check_rpi()
    if Sys.which("pigpiod") === nothing
        return false
    else
        return true
    end
end

function kbd_listener(c::Channel)
    # based on https://discourse.julialang.org/t/how-to-detect-key-down-events/95011/2
    # run listener as separate task using channels, put keypresses in channel for main loop
    t = REPL.TerminalMenus.terminal
    while true
        REPL.Terminals.raw!(t, true) || error("Unable to switch to raw mode.")
        keypress = Char(REPL.TerminalMenus.readkey(t.in_stream))
        REPL.Terminals.raw!(t, false) || error("Unable to switch back from raw mode.")
        put!(c, keypress)
    end
end
