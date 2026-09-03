"""
    _beep()

Play the NeuroAnalyzer notification sound.
"""
function _beep()::Nothing
    beep, fs = wavread(joinpath(res_path, "beep.wav"))
    wavplay(beep, fs)
    return nothing
end

"""
    _check_rpi()

Check whether a Raspberry Pi GPIO daemon (`pigpiod`) is available.

# Returns

- `Pi` object if `pigpiod` is found on PATH
- `false` otherwise
"""
function _check_rpi()::Union{Pi, Bool}
    Sys.which("pigpiod") === nothing && return false
    return Pi()
end

"""
    _kbd_listener(c)

Run a blocking keyboard listener that puts each keypress into channel `c`.

Must be launched as a separate `Task`. Based on:

https://discourse.julialang.org/t/how-to-detect-key-down-events/95011/2

# Arguments

- `c::Channel`: target channel for `Char` keypresses
"""
function _kbd_listener(c::Channel)::Nothing
    t = REPL.TerminalMenus.terminal
    while true
        REPL.Terminals.raw!(t, true) || error("Unable to switch to raw mode.")
        keypress = Char(REPL.TerminalMenus.readkey(t.in_stream))
        REPL.Terminals.raw!(t, false) || error("Unable to switch back from raw mode.")
        put!(c, keypress)
    end
    return nothing
end

"""
    _check_dialout()

On Unix systems, verify that the current user belongs to the `dialout` group.

Throws `ArgumentError` if not. No-op on non-Unix platforms.
"""
function _check_dialout()::Nothing
    Sys.isunix() || return nothing
    user   = readchomp(`sh -c 'echo $USER'`)
    groups = split(readchomp(`groups`), ' ')
    "dialout" in groups || throw(
        ArgumentError("User $user does not belong to the dialout group."),
    )
    return nothing
end

"""
    _serial_open(port_name; <keyword arguments>)

Open a serial port and return the `SerialPort` handle.

# Arguments

- `port_name::String="/dev/ttyACM0"`: serial device path
- `baudrate::Int64=115200`: baud rate
- `m`: access mode (default: `LibSerialPort.SP_MODE_READ`)

# Returns

- `SerialPort`
"""
function _serial_open(
    port_name::String = "/dev/ttyACM0";
    baudrate::Int64   = 115200,
    m                 = LibSerialPort.SP_MODE_READ,
)::SerialPort
    # validate
    port_name in LibSerialPort.get_port_list() ||
        throw(ArgumentError("Serial port $port_name does not exist."))
    _check_dialout()

    try
        sp = LibSerialPort.open(port_name, baudrate; mode = m)
        sleep(1)
        return sp
    catch err
        throw(ErrorException("Serial port $port_name cannot be opened: $err"))
    end
end

"""
    _serial_listener(sp)

Read one line from `sp` if bytes are available, otherwise return `nothing`.

# Arguments

- `sp::LibSerialPort.SerialPort`: open serial port

# Returns

- `String` if data is available, `nothing` otherwise
"""
function _serial_listener(sp::LibSerialPort.SerialPort)::Union{String, Nothing}
    !isopen(sp) && return nothing
    return bytesavailable(sp) > 0 ? String(readline(sp)) : nothing
end

"""
    _serial_close(sp)

Close a serial port. Accepts an open `SerialPort` handle.

No-op if the port is already closed.
"""
function _serial_close(sp::LibSerialPort.SerialPort)::Nothing
    isopen(sp) && close(sp)
    return nothing
end

"""
    _serial_close(port_name)

Close a serial port. Accepts a port name string.

No-op if the port is already closed.
"""
function _serial_close(port_name::String)::Nothing
    sp = LibSerialPort.open(port_name)
    isopen(sp) && close(sp)
    return nothing
end

"""
    _serial_recorder(port_name; <keyword arguments>)

Record data from a serial port and return it as a `DataFrame`.

Each "block" consists of `n` colon-separated `"key:value"` records. Recording stops after `blocks` blocks have been collected, or after `t` seconds (if `t > 0`).

# Arguments

- `port_name::String="/dev/ttyUSB0"`: serial device path
- `baudrate::Int64=115200`: baud rate
- `m`: access mode (default: `LibSerialPort.SP_MODE_READ`)
- `blocks::Int64=256`: number of data blocks to record (ignored when `t > 0`)
- `n::Int64=1`: number of records per data block
- `t::Real=0`: recording duration in seconds; if `> 0`, overrides `blocks`

# Returns

- `DataFrame` with a `"time"` column and one column per signal channel
"""
function _serial_recorder(
    port_name::String = "/dev/ttyUSB0";
    baudrate::Int64   = 115200,
    m                 = LibSerialPort.SP_MODE_READ,
    blocks::Int64     = 256,
    n::Int64          = 1,
    t::Real           = 0,
)::DataFrame
    port_name in LibSerialPort.get_port_list() ||
        throw(ArgumentError("Serial port $port_name does not exist."))
    _check_dialout()
    n >= 1 || throw(ArgumentError("n must be ≥ 1."))
    t >= 0 || throw(ArgumentError("t must be ≥ 0."))
    blocks >= 1 || throw(ArgumentError("blocks must be ≥ 1."))

    local sp
    try
        sp = LibSerialPort.open(port_name, baudrate; mode = m)
        sleep(1)
    catch err
        throw(ErrorException("Serial port $port_name cannot be opened: $err"))
    end
    isopen(sp) || throw(ArgumentError("Serial port $port_name is not open."))

    tp = Float64[]
    tmp_data = String[]

    _beep()
    if t == 0
        _info("Recording $blocks data-blocks from $port_name ($n record(s) per block)")
        for _ = 1:(blocks * n)
            push!(tp, time())
            push!(tmp_data, String(readline(sp)))
            sleep(0.01)
        end
    else
        _info("Recording for $t seconds from $port_name ($n record(s) per block)")
        t_start = time()
        while time() < t_start + t
            for _ = 1:n
                push!(tp, time())
                push!(tmp_data, String(readline(sp)))
                sleep(0.01)
            end
        end
        blocks = length(tp) ÷ n
    end
    _beep()
    _info("Recording finished")
    close(sp)

    # recalculate timestamps relative to recording start, one per block
    tp .-= tp[1]
    tp = tp[1:n:end]

    # estimate sampling rate from the last inter-block interval
    fs = round(Int64, (blocks - 1) / tp[end])
    _info("Sampling rate: $fs Hz")
    tp = round.(tp; digits = 4)

    # parse channel names from the first block
    col_names = ["time"]
    for idx = 1:n
        push!(col_names, split(tmp_data[idx], ':')[1])
    end

    # parse values into a (blocks × n+1) matrix
    data = zeros(Float64, blocks, n + 1)
    for (blk_idx, start) in enumerate(1:n:length(tmp_data))
        stop  = min(start + n - 1, length(tmp_data))
        block = split.(tmp_data[start:stop], ':')
        for rec_idx in eachindex(block)
            data[blk_idx, rec_idx + 1] = parse(Float64, last(block[rec_idx]))
        end
    end
    data[:, 1] = tp

    return DataFrame(data, col_names)
end
