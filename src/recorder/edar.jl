export iedar
export edar

# Grove GSR Arduino sketch: 50 Hz, lines "gsr:<value>"
const _EDA_FS = 50
const _EDA_BAUDRATE = 19200
const _EDA_REGEX = r"gsr:([0-9]+\.[0-9]+)"

# value of one EDA serial line, or nothing if the line is not an EDA sample
function _parse_eda_sample(line::String)::Union{Float64, Nothing}
    m = match(_EDA_REGEX, line)
    return isnothing(m) ? nothing : parse(Float64, m.captures[1])
end

# NEURO object with one EDA channel from the recorded samples
function _build_eda_object(eda_signal::Vector{Float64})::NeuroAnalyzer.NEURO
    isempty(eda_signal) && throw(ArgumentError("No EDA samples were recorded."))
    obj = create_object(; data_type = "eda")
    # append_channel!() needs the sampling rate of the (empty) object
    obj.header.recording[:sampling_rate] = _EDA_FS
    append_channel!(
        obj;
        data = reshape(eda_signal, 1, :, 1),
        label = ["eda1"],
        type = ["eda"],
        unit = ["μS"],
    )
    create_time!(obj; fs = _EDA_FS)
    return obj
end

# window of iedar(): live plot of the recorded EDA, RECORD button and status
function _iedar_build(duration::Int64)
    n = duration * _EDA_FS
    t = collect(0:(n - 1)) ./ _EDA_FS

    fig = GLMakie.Figure(; size = (900, 500))
    ax = GLMakie.Axis(
        fig[1, 1:3];
        title = "Electrodermal activity",
        xlabel = "Time [s]",
        ylabel = "Skin conductance [μS]",
        _AXIS_LOCK_KWARGS...,
    )
    GLMakie.xlims!(ax, 0, duration)
    GLMakie.ylims!(ax, 0, 10)
    pts = Observable(GLMakie.Point2f[])
    GLMakie.lines!(ax, pts; color = :black, linewidth = 1)
    bt_record = GLMakie.Button(fig[2, 1]; label = "RECORD", width = 120)
    GLMakie.Label(fig[2, 2], "Status:"; halign = :right)
    lb_status = GLMakie.Label(fig[2, 3], _REC_STATUS[1]; halign = :left, tellwidth = false)

    # show the first n_rec samples; the y-axis follows the signal
    function show!(eda::Vector{Float64}, n_rec::Int64)
        pts[] = GLMakie.Point2f.(t[1:n_rec], eda[1:n_rec])
        if n_rec > 0
            lo, hi = extrema(eda[1:n_rec])
            m = max(0.1 * (hi - lo), 0.5)
            GLMakie.ylims!(ax, lo - m, hi + m)
        end
        return nothing
    end

    return (; fig, ax, bt_record, lb_status, pts, show!)
end

"""
    iedar(; <keyword arguments>)

Record electrodermal activity (EDA) in GUI mode.

EDA (also known as Galvanic Skin Response / skin conductance) is recorded using a Grove GSR sensor via an Arduino connected over a virtual serial port. The recorded signal is plotted live.

Nominal sampling rate is 50 Hz. If the window is closed before the end of the recording, the samples recorded so far are returned.

# Arguments

- `duration::Int64=20`: recording duration in seconds; must be > 0
- `port_name::String="/dev/ttyUSB0"`: serial port to which the Arduino is connected

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object (EDA) with one channel at 50 Hz
"""
function iedar(;
    duration::Int64 = 20,
    port_name::String = "/dev/ttyUSB0",
)::NeuroAnalyzer.NEURO
    duration > 0 || throw(ArgumentError("duration must be > 0."))
    sp = _serial_open(port_name; baudrate = _EDA_BAUDRATE)

    n = duration * _EDA_FS
    eda_signal = fill(NaN, n)
    # shared with the recording task
    n_rec = Threads.Atomic{Int64}(0)
    phase = Threads.Atomic{Int64}(_REC_READY)
    aborted = Threads.Atomic{Bool}(false)
    rec_task = Ref{Union{Nothing, Task}}(nothing)

    w = _iedar_build(duration)

    on(w.bt_record.clicks) do _
        phase[] == _REC_READY || return nothing
        phase[] = _REC_PREPARING
        rec_task[] = Threads.@spawn begin
            try
                # discard data buffered while waiting for the user
                _serial_flush(sp, 2)
                _beep()
                phase[] = _REC_RECORDING
                while n_rec[] < n && !aborted[]
                    line = _serial_listener(sp)
                    isnothing(line) && (sleep(0.001); continue)
                    v = _parse_eda_sample(line)
                    if !isnothing(v)
                        eda_signal[n_rec[] + 1] = v
                        n_rec[] += 1
                    end
                end
                aborted[] && return nothing
                phase[] = _REC_FINISHED
                _beep()
            catch err
                phase[] = _REC_ERROR
                @error "EDA recording failed" exception = (err, catch_backtrace())
            end
        end
        return nothing
    end

    shown = Ref(-1)
    function refresh!()
        st = _REC_STATUS[phase[] + 1]
        if w.lb_status.text[] != st
            w.lb_status.text[] = st
            # RECORD works once: grayed out after it has been pressed
            phase[] == _REC_READY || _button_disable!(w.bt_record)
        end
        k = n_rec[]
        k != shown[] && (w.show!(eda_signal, k); shown[] = k)
        return nothing
    end
    _recorder_window_run(w.fig, "NeuroRecorder: iedar()", phase, refresh!)

    # the window is closed: stop the recording task (if any) and release the port
    aborted[] = true
    isnothing(rec_task[]) || wait(rec_task[])
    _serial_close(sp)

    return _build_eda_object(eda_signal[1:n_rec[]])
end

"""
    edar(; <keyword arguments>)

Record electrodermal activity (EDA) in CLI mode.

EDA (also known as Galvanic Skin Response / skin conductance) is recorded using a Grove GSR sensor via an Arduino connected over a virtual serial port.

Nominal sampling rate is 50 Hz.

# Arguments

- `duration::Int64=20`: recording duration in seconds; must be > 0
- `port_name::String="/dev/ttyUSB0"`: serial port to which the Arduino is connected

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object (EDA) with one channel at 50 Hz
"""
function edar(;
    duration::Int64 = 20,
    port_name::String = "/dev/ttyUSB0",
)::NeuroAnalyzer.NEURO
    duration > 0 || throw(ArgumentError("duration must be > 0."))
    sp = _serial_open(port_name; baudrate = _EDA_BAUDRATE)

    n = duration * _EDA_FS
    eda_signal = zeros(n)

    try
        println("NeuroRecorder: EDA")
        println("==================")
        println("   Duration: $duration [seconds]")
        println("Serial port: $port_name")
        println()
        println("Ready to start, press SPACEBAR to begin recording")
        println()

        _wait_spacebar()

        println("The recording will start after a beep")
        # discard data buffered while waiting for the user
        _serial_flush(sp, 2)
        println()
        _beep()
        print("Recording .")

        idx = 1
        ts = time()
        while idx <= n
            line = _serial_listener(sp)
            if isnothing(line)
                sleep(0.001)
            else
                v = _parse_eda_sample(line)
                if !isnothing(v)
                    eda_signal[idx] = v
                    idx += 1
                end
            end
            if time() - ts >= 1.0
                print(".")
                ts = time()
            end
        end
    finally
        _serial_close(sp)
    end
    _beep()
    println()
    println()
    println("Recording finished.")

    return _build_eda_object(eda_signal)
end
