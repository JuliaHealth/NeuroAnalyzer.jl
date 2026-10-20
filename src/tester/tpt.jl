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

"""
    itpt(; duration, port_name)

Perform the **Two-Point Pinch Test (TPT)** in **GUI mode**.

Orientation and acceleration are sampled at $_TPT_FS Hz from an **MMA7660** accelerometer connected via an Arduino on a USB virtual serial port.

The test window shows a finger graphic and a status label. Pressing **RECORD** flushes stale serial data for $_TPT_WARMUP_SECS seconds, then captures exactly `duration × $_TPT_FS` samples before closing automatically.

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

    img_idle = read_from_png(joinpath(res_path, "finger_nopinch.png"))
    img_pinch = read_from_png(joinpath(res_path, "finger_pinch.png"))

    n_samples = duration * _TPT_FS
    tpt_ch_x = zeros(n_samples)
    tpt_ch_y = zeros(n_samples)
    tpt_ch_z = zeros(n_samples)
    tpt_ch_accx = zeros(n_samples)
    tpt_ch_accy = zeros(n_samples)
    tpt_ch_accz = zeros(n_samples)
    n_rec = 0 # number of recorded samples
    aborted = false # set when the window is closed
    rec_task = nothing # recording task (started by RECORD)

    # =========================================================================
    function _activate(app)
        win = GtkApplicationWindow(app, "NeuroAnalyzer: itpt()")
        Gtk4.default_size(win, Int64(img_idle.width), Int64(img_idle.height) + 100)

        can = GtkCanvas()
        can.content_width = Int64(img_idle.width)
        can.content_height = Int64(img_idle.height)

        g = GtkGrid()
        g.column_homogeneous = false
        g.column_spacing = 5
        g.row_spacing = 5

        bt_record = GtkButton("RECORD")
        bt_record.tooltip_text = "Start recording"

        lb_status1 = GtkLabel("Status:")
        lb_status1.halign = 2
        lb_status2 = GtkLabel("READY TO START")
        lb_status2.halign = 1

        g[1:2, 1] = can
        g[1:2, 2] = bt_record
        g[1, 3]   = lb_status1
        g[2, 3]   = lb_status2

        vbox = GtkBox(:v)
        push!(vbox, g)
        push!(win, vbox)
        Gtk4.show(win)

        # closing the window stops the recording
        signal_connect(win, "close-request") do widget
            aborted = true
            return false
        end

        # initial canvas draw — idle (unpinched) finger image
        @guarded draw(can) do widget
            ctx = getgc(can)
            Cairo.set_source_surface(ctx, img_idle, 0, 0)
            return Cairo.paint(ctx)
        end

        # --- RECORD button ---
        @guarded signal_connect(bt_record, "clicked") do widget
            bt_record.sensitive = false
            rec_task = Threads.@spawn begin
                try
                    # flush stale serial data before recording begins
                    @idle_add lb_status2.label = "PREPARING"
                    _serial_flush(sp, _TPT_WARMUP_SECS)

                    # signal start and switch to pinch image
                    _beep()
                    @idle_add @guarded draw(can) do widget
                        ctx = getgc(can)
                        Cairo.set_source_surface(ctx, img_pinch, 0, 0)
                        return Cairo.paint(ctx)
                    end
                    @idle_add lb_status2.label = "RECORDING"

                    # collect n_samples valid samples (or until the window is closed)
                    while n_rec < n_samples && !aborted
                        line = _serial_listener(sp)
                        isnothing(line) && (sleep(0.001); continue)
                        sample = _parse_tpt_sample(line)
                        if !isnothing(sample)
                            n_rec += 1
                            tpt_ch_x[n_rec] = sample.x
                            tpt_ch_y[n_rec] = sample.y
                            tpt_ch_z[n_rec] = sample.z
                            tpt_ch_accx[n_rec] = sample.accx
                            tpt_ch_accy[n_rec] = sample.accy
                            tpt_ch_accz[n_rec] = sample.accz
                        end
                    end
                    aborted && return nothing

                    # finalize
                    _beep()
                    @idle_add lb_status2.label = "FINISHED"
                    @idle_add @guarded draw(can) do widget
                        ctx = getgc(can)
                        Cairo.set_source_surface(ctx, img_idle, 0, 0)
                        return Cairo.paint(ctx)
                    end
                    sleep(2)
                    @idle_add close(win)
                catch err
                    @error "TPT recording failed" exception = (err, catch_backtrace())
                finally
                    # ensure the port is always closed, even if the task errors
                    _serial_close(sp)
                end
            end
        end
    end
    # =========================================================================

    app = GtkApplication("org.neuroanalyzer.itpt")
    Gtk4.signal_connect(_activate, app, :activate)
    Gtk4.GLib.stop_main_loop()
    Gtk4.run(app)

    # the window is closed: stop the recording task (if any) and release the port
    aborted = true
    isnothing(rec_task) || wait(rec_task)
    _serial_close(sp)

    r = 1:n_rec
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
    tpt(; duration, port_name)

Perform the **Two-Point Pinch Test (TPT)** in **CLI (terminal) mode**.

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
