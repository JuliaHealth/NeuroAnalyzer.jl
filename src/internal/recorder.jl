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

Must be launched as a separate `Task`. Stops when `c` is closed (after the next keypress). Based on:

https://discourse.julialang.org/t/how-to-detect-key-down-events/95011/2

# Arguments

- `c::Channel`: target channel for `Char` keypresses
- `timestamp::Bool=false`: if `true`, put `(keypress, time())` tuples, timed when the key is read
"""
function _kbd_listener(c::Channel; timestamp::Bool = false)::Nothing
    t = REPL.TerminalMenus.terminal
    while isopen(c)
        REPL.Terminals.raw!(t, true) || error("Unable to switch to raw mode.")
        keypress = try
            Char(REPL.TerminalMenus.readkey(t.in_stream))
        finally
            REPL.Terminals.raw!(t, false)
        end
        tk = time()
        isopen(c) || break
        put!(c, timestamp ? (keypress, tk) : keypress)
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

# serial ports opened by `_serial_open`, by port name, so that `_serial_close(port_name)`
# can close the handle that is actually in use (LibSerialPort cannot look it up by name)
const _SERIAL_PORTS = Dict{String, LibSerialPort.SerialPort}()
const _SERIAL_PORTS_LOCK = ReentrantLock()

"""
    _serial_open(port_name; <keyword arguments>)

Open a serial port and return the `SerialPort` handle.

The handle is registered under `port_name`, so it can also be closed with `_serial_close(port_name)`. Throws if the port is already open through `_serial_open`.

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

    lock(_SERIAL_PORTS_LOCK) do
        haskey(_SERIAL_PORTS, port_name) && isopen(_SERIAL_PORTS[port_name]) &&
            throw(ArgumentError("Serial port $port_name is already open."))
    end

    sp = try
        LibSerialPort.open(port_name, baudrate; mode = m)
    catch err
        throw(ErrorException("Serial port $port_name cannot be opened: $err"))
    end
    lock(() -> (_SERIAL_PORTS[port_name] = sp), _SERIAL_PORTS_LOCK)
    sleep(1)

    return sp
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

No-op if the port is already closed. The handle is removed from the registry of ports opened by `_serial_open`.
"""
function _serial_close(sp::LibSerialPort.SerialPort)::Nothing
    isopen(sp) && close(sp)
    lock(_SERIAL_PORTS_LOCK) do
        Base.filter!(p -> p.second !== sp, _SERIAL_PORTS)
    end
    return nothing
end

"""
    _serial_close(port_name)

Close a serial port opened by `_serial_open`. Accepts a port name string.

No-op if the port is not open (or was not opened by `_serial_open`).
"""
function _serial_close(port_name::String)::Nothing
    sp = lock(() -> pop!(_SERIAL_PORTS, port_name, nothing), _SERIAL_PORTS_LOCK)
    !isnothing(sp) && isopen(sp) && close(sp)
    return nothing
end

"""
    _serial_try_open(port_name; <keyword arguments>)

Open a serial port like `_serial_open`, but return `nothing` (with an info message) instead of throwing when the port does not exist or cannot be opened. Keyword arguments are passed to `_serial_open`.

# Returns

- `Union{SerialPort, Nothing}`
"""
function _serial_try_open(port_name::String; kwargs...)::Union{SerialPort, Nothing}
    return try
        _serial_open(port_name; kwargs...)
    catch err
        _info("Serial port $port_name is not available: $(err isa Exception ? sprint(showerror, err) : err)")
        nothing
    end
end

"""
    _serial_button!(sp, gpio, t0, t_press, t_release, pressed)

Read all pending button events (`"<gpio>:1"` = pressed, `"<gpio>:0"` = released) from `sp`. The time of each state change (seconds since `t0`) is appended to `t_press` or `t_release`.

# Returns

- `Bool`: button state after the events
"""
function _serial_button!(
    sp::LibSerialPort.SerialPort,
    gpio::Int64,
    t0::Float64,
    t_press::Vector{Float64},
    t_release::Vector{Float64},
    pressed::Bool,
)::Bool
    while true
        line = _serial_listener(sp)
        isnothing(line) && break
        # lines may end with "\r" (CRLF)
        ev = strip(line)
        if ev == "$gpio:1" && !pressed
            push!(t_press, time() - t0)
            pressed = true
        elseif ev == "$gpio:0" && pressed
            push!(t_release, time() - t0)
            pressed = false
        end
    end
    return pressed
end

"""
    _serial_flush(sp, seconds)

Drain the serial port `sp` for `seconds` seconds, discarding all incoming data. Used to discard stale buffered samples before recording begins.
"""
function _serial_flush(sp, seconds::Real)::Nothing
    t0 = time()
    while time() - t0 <= seconds
        _serial_listener(sp)
    end
    return nothing
end

"""
    _wait_spacebar()

Block until the SPACEBAR is pressed in the terminal (raw mode; the terminal mode is always restored).
"""
function _wait_spacebar()::Nothing
    while true
        ret = ccall(:jl_tty_set_mode, Int32, (Ptr{Cvoid}, Int32), stdin.handle, true)
        ret == 0 || error("Unable to switch the terminal to raw mode.")
        kbd_key = try
            read(stdin, Char)
        finally
            ccall(:jl_tty_set_mode, Int32, (Ptr{Cvoid}, Int32), stdin.handle, false)
        end
        kbd_key == ' ' && return nothing
    end
end

"""
    _serial_recorder(port_name; <keyword arguments>)

Record data from a serial port and return it as a `DataFrame`.

Each "block" consists of `n` colon-separated `"key:value"` records. Recording stops after `blocks` blocks have been collected, or after `t` seconds (if `t > 0`).

# Arguments

- `port_name::String="/dev/ttyUSB0"`: serial device path
- `baudrate::Int64=115200`: baud rate
- `m`: access mode (default: `LibSerialPort.SP_MODE_READ`)
- `blocks::Int64=256`: number of data blocks to record (ignored when `t > 0`); must be ≥ 2
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
    n >= 1 || throw(ArgumentError("n must be ≥ 1."))
    t >= 0 || throw(ArgumentError("t must be ≥ 0."))
    # ≥ 2 blocks are needed to estimate the sampling rate
    blocks >= 2 || throw(ArgumentError("blocks must be ≥ 2."))

    # validates the port, opens it and registers the handle
    sp = _serial_open(port_name; baudrate = baudrate, m = m)
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
    _serial_close(sp)

    # recalculate timestamps relative to recording start, one per block
    tp .-= tp[1]
    tp = tp[1:n:end]

    # estimate sampling rate as the mean block rate over the whole recording
    blocks >= 2 && tp[end] > 0 ||
        throw(ArgumentError("At least 2 data blocks are needed to estimate the sampling rate."))
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
