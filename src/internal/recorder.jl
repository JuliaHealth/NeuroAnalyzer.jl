function _beep()
    beep, fs = wavread(joinpath(res_path, "beep.wav"))
    wavplay(beep, fs)
end

function _check_rpi()
    if Sys.which("pigpiod") === nothing
        rpi = false
    else
        rpi = Pi()
    end
    return rpi
end

function _kbd_listener(c::Channel)
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

function _serial_listener(port_name::String="/dev/ttyUSB0"; baudrate::Int64=115200, mode=SP_MODE_READ, blocks::Int64=256, n::Int64=1, t::Real=0)

    # `blocks`: number of data blocks to record
    # `n`: number of records per block
    # `t`: recording time in seconds; if t > 0, blocks ignored and calculated based on recorded data

    @assert port_name in LibSerialPort.get_port_list() "$port_name does not exist."
    
    if Sys.isunix()
        @assert "dialout" in split(readchomp(`groups`), ' ') "User $(readchomp(`sh -c 'echo $USER'`)) does not belong to the dialout group."
    end

    sp = nothing
    try
        sp = LibSerialPort.open(port_name, baudrate, mode=mode)
        sleep(1)
    catch
        error("Serial port $port_name cannot be opened.")
    end

    @assert isopen(sp) "Serial port $port_name is not open."

    tp = Float64[]
    tmp_data = String[]

    if t == 0
        _info("Recording $blocks data-blocks from $port_name")
        _info("$n records per data-block")
        _beep()
        for idx in 1:(blocks * n)
            push!(tp, time())
            push!(tmp_data, String(readline(sp)))
            sleep(0.01)
        end
    else
        _info("Recording for $t seconds from $port_name")
        _info("$n records per data-block")
        _beep()
        t_start = time()
        while time() < t_start + t
            for idx in 1:n
                push!(tp, time())
                push!(tmp_data, String(readline(sp)))
                sleep(0.01)
            end
        end
        blocks = (length(tp) ÷ n)
    end
    _beep()
    _info("Recording finished")
    close(sp)

    # recalculate time points per block
    tp .-= tp[1]
    tp = tp[1:n:end]

    # calculate sampling rate
    sr = round(Int64, 1 / (tp[end] - tp[end - 1]))
    _info("Sampling rate: $sr Hz")
    tp = round.(tp, digits=3)

    # create data frame
    names = String[]
    push!(names, "time")
    for idx1 in 1:n
        push!(names, split.(tmp_data[1:n], ':')[idx1][1])
    end
    data = zeros(blocks, n + 1)
    idx = 1
    for idx1 in 1:n:length(tmp_data)
        block = split.(tmp_data[idx1:idx1 + (n - 1)], ':')
        for idx2 in 1:n
            data[idx, idx2 + 1] = parse(Float64, block[idx2][2])
        end
        idx += 1
    end
    data[:, 1] = tp

    df = DataFrame(data, names)

    return df

end
