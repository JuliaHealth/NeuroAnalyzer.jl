export itpt
export tpt

# =============================================================================
# shared helpers used by both itpt() and tpt()
# =============================================================================

# Serial communication constants
const _TPT_BAUDRATE = 19200 # MMA7660 Arduino firmware baud rate
const _TPT_FS = 50 # accelerometer sampling rate [Hz]
const _TPT_WARMUP_SECS = 2 # seconds to flush stale serial data before recording

# Expected serial line format from the Arduino firmware:
#   "tpt: X Y Z AX AY AZ"
# where X/Y/Z are raw orientation ints and AX/AY/AZ are computed acceleration floats.
const _TPT_SERIAL_REGEX =
    r"(tpt\: )(\-*[0-9]+) (\-*[0-9]+) (\-*[0-9]+) (\-*[0-9]+\.[0-9]+) (\-*[0-9]+\.[0-9]+) (\-*[0-9]+\.[0-9]+)"

"""
    _parse_tpt_sample(line) -> NamedTuple | nothing

Attempt to parse one TPT serial data line into its six sensor values.

Returns a named tuple `(x, y, z, accx, accy, accz)` on success, or `nothing` if the line does not match the expected format.

Expected format: `"tpt: <x> <y> <z> <accx> <accy> <accz>"`
- `x`, `y`, `z` — raw MMA7660 orientation register values (integers)
- `accx`, `accy`, `accz` — calibrated accelerations in m/s² (floats)
"""
function _parse_tpt_sample(line::String)
    m = match(_TPT_SERIAL_REGEX, line)
    isnothing(m) && return nothing
    return (
        x = parse(Float64, m.captures[2]),
        y = parse(Float64, m.captures[3]),
        z = parse(Float64, m.captures[4]),
        accx = parse(Float64, m.captures[5]),
        accy = parse(Float64, m.captures[6]),
        accz = parse(Float64, m.captures[7]),
    )
end

"""
    _build_tpt_object(ch_x, ch_y, ch_z, ch_accx, ch_accy, ch_accz) -> NeuroAnalyzer.NEURO

Assemble the six TPT channel arrays into a `NeuroAnalyzer.NEURO` object.

Channels:
- `pos_x / pos_y / pos_z` — raw MMA7660 orientation (type `"orient"`, unitless)
- `acc_x / acc_y / acc_z` — calibrated acceleration (type `"accel"`, unit `"m/s²"`)
"""
function _build_tpt_object(
    ch_x::Vector{Float64}, ch_y::Vector{Float64}, ch_z::Vector{Float64},
    ch_accx::Vector{Float64}, ch_accy::Vector{Float64}, ch_accz::Vector{Float64},
)::NeuroAnalyzer.NEURO
    # stack channels as rows, add a singleton epoch dimension: (6, n_samples, 1)
    signal = reshape(Matrix([ch_x ch_y ch_z ch_accx ch_accy ch_accz]'), 6, :, 1)

    obj = create_object(; data_type = "tpt")
    # add_channel!() needs the sampling rate of the (empty) object
    obj.header.recording[:sampling_rate] = _TPT_FS
    add_channel!(
        obj;
        data  = signal,
        label = ["pos_x", "pos_y", "pos_z", "acc_x", "acc_y", "acc_z"],
        type  = ["orient", "orient", "orient", "accel", "accel", "accel"],
        unit  = ["", "", "", "m/s²", "m/s²", "m/s²"],
    )
    create_time!(obj; fs = _TPT_FS)
    return obj
end

# =============================================================================

# window of itpt(): finger image (idle / pinch), live plot of the acceleration, RECORD
# button and status
function _itpt_build(duration::Int64)
    n = duration * _TPT_FS
    t = collect(0:(n - 1)) ./ _TPT_FS
    img_idle = rotr90(FileIO.load(joinpath(res_path, "finger_nopinch.png")))
    img_pinch = rotr90(FileIO.load(joinpath(res_path, "finger_pinch.png")))

    fig = GLMakie.Figure(; size = (1300, 650))
    ax_img = GLMakie.Axis(
        fig[1, 1];
        aspect = GLMakie.DataAspect(),
        width = 450,
        _AXIS_LOCK_KWARGS...,
    )
    GLMakie.hidedecorations!(ax_img)
    GLMakie.hidespines!(ax_img)
    img = Observable(img_idle)
    GLMakie.image!(ax_img, img)

    ax = GLMakie.Axis(
        fig[1, 2:4];
        title = "Acceleration",
        xlabel = "Time [s]",
        ylabel = "Acceleration [m/s²]",
        _AXIS_LOCK_KWARGS...,
    )
    GLMakie.xlims!(ax, 0, duration)
    GLMakie.ylims!(ax, -1, 1)
    pts = [Observable(GLMakie.Point2f[]) for _ = 1:3]
    for (p, c, l) in zip(pts, (:dodgerblue3, :darkorange2, :forestgreen), ("X", "Y", "Z"))
        GLMakie.lines!(ax, p; color = c, linewidth = 1, label = l)
    end
    GLMakie.axislegend(ax; position = :rt, orientation = :horizontal, framevisible = false)

    bt_record = GLMakie.Button(fig[2, 2]; label = "RECORD", width = 120)
    GLMakie.Label(fig[2, 3], "Status:"; halign = :right)
    lb_status = GLMakie.Label(fig[2, 4], _REC_STATUS[1]; halign = :left, tellwidth = false)

    # show the first n_rec samples of the acceleration (x, y, z); the y-axis follows the signal
    function show!(acc::NTuple{3, Vector{Float64}}, n_rec::Int64)
        for idx = 1:3
            pts[idx][] = GLMakie.Point2f.(t[1:n_rec], acc[idx][1:n_rec])
        end
        if n_rec > 0
            lo = minimum(minimum(a[1:n_rec]) for a in acc)
            hi = maximum(maximum(a[1:n_rec]) for a in acc)
            m = max(0.1 * (hi - lo), 0.5)
            GLMakie.ylims!(ax, lo - m, hi + m)
        end
        return nothing
    end
    # finger image: pinching while recording
    pinch!(p::Bool) = (img[] = p ? img_pinch : img_idle; nothing)

    return (; fig, ax, bt_record, lb_status, pts, img, show!, pinch!)
end

"""
    itpt(; <keyword arguments>)

Perform the **Two-Point Pinch Test (TPT)** in GUI mode.

Orientation and acceleration are sampled at $_TPT_FS Hz from an **MMA7660** accelerometer connected via an Arduino on a USB virtual serial port.

The test window shows a finger graphic, a live plot of the acceleration and a status label. Pressing **RECORD** flushes stale serial data for $_TPT_WARMUP_SECS seconds, then captures exactly `duration × $_TPT_FS` samples before closing automatically. If the window is closed before the end of the test, the samples recorded so far are returned.

# Arguments

- `duration::Int64=20`: recording duration in seconds
- `port_name::String="/dev/ttyUSB0"`: serial port path for the Arduino

# Returns

- `NeuroAnalyzer.NEURO` with six channels sampled at $_TPT_FS Hz:

| Label   | Type     | Unit  | Description                  |
|---------|----------|-------|------------------------------|
| `pos_x` | `orient` |       | raw MMA7660 X orientation    |
| `pos_y` | `orient` |       | raw MMA7660 Y orientation    |
| `pos_z` | `orient` |       | raw MMA7660 Z orientation    |
| `acc_x` | `accel`  | m/s²  | calibrated X acceleration    |
| `acc_y` | `accel`  | m/s²  | calibrated Y acceleration    |
| `acc_z` | `accel`  | m/s²  | calibrated Z acceleration    |
"""
function itpt(;
    duration::Int64 = 20,
    port_name::String = "/dev/ttyUSB0",
)::NeuroAnalyzer.NEURO
    duration > 0 || throw(ArgumentError("duration must be > 0."))
    sp = _serial_open(port_name; baudrate = _TPT_BAUDRATE)

    n_samples = duration * _TPT_FS
    tpt_ch_x = zeros(n_samples)
    tpt_ch_y = zeros(n_samples)
    tpt_ch_z = zeros(n_samples)
    tpt_ch_accx = zeros(n_samples)
    tpt_ch_accy = zeros(n_samples)
    tpt_ch_accz = zeros(n_samples)
    # shared with the recording task
    n_rec = Threads.Atomic{Int64}(0)
    phase = Threads.Atomic{Int64}(_REC_READY)
    aborted = Threads.Atomic{Bool}(false)
    rec_task = Ref{Union{Nothing, Task}}(nothing)

    w = _itpt_build(duration)

    on(w.bt_record.clicks) do _
        phase[] == _REC_READY || return nothing
        phase[] = _REC_PREPARING
        rec_task[] = Threads.@spawn begin
            try
                # flush stale serial data before recording begins
                _serial_flush(sp, _TPT_WARMUP_SECS)
                _beep()
                phase[] = _REC_RECORDING
                # collect n_samples valid samples (or until the window is closed)
                while n_rec[] < n_samples && !aborted[]
                    line = _serial_listener(sp)
                    isnothing(line) && (sleep(0.001); continue)
                    sample = _parse_tpt_sample(line)
                    if !isnothing(sample)
                        k = n_rec[] + 1
                        tpt_ch_x[k] = sample.x
                        tpt_ch_y[k] = sample.y
                        tpt_ch_z[k] = sample.z
                        tpt_ch_accx[k] = sample.accx
                        tpt_ch_accy[k] = sample.accy
                        tpt_ch_accz[k] = sample.accz
                        n_rec[] = k
                    end
                end
                aborted[] && return nothing
                phase[] = _REC_FINISHED
                _beep()
            catch err
                phase[] = _REC_ERROR
                @error "TPT recording failed" exception = (err, catch_backtrace())
            end
        end
        return nothing
    end

    shown = Ref(-1)
    pinching = Ref(false)
    function refresh!()
        st = _REC_STATUS[phase[] + 1]
        if w.lb_status.text[] != st
            w.lb_status.text[] = st
            # RECORD works once: grayed out after it has been pressed
            phase[] == _REC_READY || _button_disable!(w.bt_record)
        end
        recording = phase[] == _REC_RECORDING
        recording == pinching[] || (w.pinch!(recording); pinching[] = recording)
        k = n_rec[]
        if k != shown[]
            w.show!((tpt_ch_accx, tpt_ch_accy, tpt_ch_accz), k)
            shown[] = k
        end
        return nothing
    end
    _recorder_window_run(w.fig, "NeuroTester: itpt()", phase, refresh!)

    # the window is closed: stop the recording task (if any) and release the port
    aborted[] = true
    isnothing(rec_task[]) || wait(rec_task[])
    _serial_close(sp)

    r = 1:n_rec[]
    return _build_tpt_object(
        tpt_ch_x[r],
        tpt_ch_y[r],
        tpt_ch_z[r],
        tpt_ch_accx[r],
        tpt_ch_accy[r],
        tpt_ch_accz[r],
    )
end

"""
    tpt(; <keyword arguments>)

Perform the **Two-Point Pinch Test (TPT)** in CLI (terminal) mode.

Functionally identical to [`itpt`](@ref) but runs entirely in the terminal without a graphical window. Press `SPACEBAR` to begin; the test starts after a beep and a $_TPT_WARMUP_SECS-second serial flush.

# Arguments

- `duration::Int64=20`: recording duration in seconds
- `port_name::String="/dev/ttyUSB0"`: serial port path for the Arduino

# Returns

- `NeuroAnalyzer.NEURO` with six channels sampled at $_TPT_FS Hz:

| Label   | Type     | Unit  | Description                  |
|---------|----------|-------|------------------------------|
| `pos_x` | `orient` |       | raw MMA7660 X orientation    |
| `pos_y` | `orient` |       | raw MMA7660 Y orientation    |
| `pos_z` | `orient` |       | raw MMA7660 Z orientation    |
| `acc_x` | `accel`  | m/s²  | calibrated X acceleration    |
| `acc_y` | `accel`  | m/s²  | calibrated Y acceleration    |
| `acc_z` | `accel`  | m/s²  | calibrated Z acceleration    |
"""
function tpt(;
    duration::Int64   = 20,
    port_name::String = "/dev/ttyUSB0",
)::NeuroAnalyzer.NEURO
    duration > 0 || throw(ArgumentError("duration must be > 0."))
    sp = _serial_open(port_name; baudrate = _TPT_BAUDRATE)

    n_samples = duration * _TPT_FS
    tpt_ch_x = zeros(n_samples)
    tpt_ch_y = zeros(n_samples)
    tpt_ch_z = zeros(n_samples)
    tpt_ch_accx = zeros(n_samples)
    tpt_ch_accy = zeros(n_samples)
    tpt_ch_accz = zeros(n_samples)

    try
        println("NeuroTester: TPT")
        println("================")
        println("   Duration: $duration [seconds]")
        println("Serial port: $port_name")
        println()
        println("Ready to start — press SPACEBAR to begin the test")
        println()

        _wait_spacebar()

        println("The test will start after a beep")

        # flush stale serial data accumulated while waiting for the user
        _serial_flush(sp, _TPT_WARMUP_SECS)

        println()
        _beep()
        print("   Pinch the thumb and the index finger as quickly as possible")

        # collect exactly n_samples valid samples
        idx = 1
        while idx <= n_samples
            line = _serial_listener(sp)
            if isnothing(line)
                sleep(0.001)
            else
                sample = _parse_tpt_sample(line)
                if !isnothing(sample)
                    tpt_ch_x[idx] = sample.x
                    tpt_ch_y[idx] = sample.y
                    tpt_ch_z[idx] = sample.z
                    tpt_ch_accx[idx] = sample.accx
                    tpt_ch_accy[idx] = sample.accy
                    tpt_ch_accz[idx] = sample.accz
                    idx += 1
                end
            end
        end
    finally
        _serial_close(sp)
    end
    _beep()
    println()
    println()
    println("Testing completed")

    return _build_tpt_object(
        tpt_ch_x,
        tpt_ch_y,
        tpt_ch_z,
        tpt_ch_accx,
        tpt_ch_accy,
        tpt_ch_accz,
    )
end
