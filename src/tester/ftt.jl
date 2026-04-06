export iftt
export ftt
# =============================================================================
# shared post-processing helpers used by both iftt() and ftt()
# =============================================================================

"""
    _pack_taps(flat_t, flat_d, counts) -> (t_per_trial, d_per_trial)

Restructure flat tap-time and tap-duration vectors into per-trial nested vectors, using `counts` to know how many taps belong to each trial.

Taps are popped from the END of the flat vectors (which accumulate in chronological order), so the inner loop runs in reverse-trial order and the result is reversed before returning.

Returns two `Vector{Vector{Float64}}` of the same length as `counts`.
"""
function _pack_taps(
    flat_t::Vector{Float64},
    flat_d::Vector{Float64},
    counts::Vector{Int64},
)::Tuple{Vector{Vector{Float64}}, Vector{Vector{Float64}}}
    n = length(counts)
    t_packed = Vector{Vector{Float64}}(undef, n)
    d_packed = Vector{Vector{Float64}}(undef, n)
    for idx in n:-1:1
        # Pop `counts[idx]` entries; popping from the back gives reverse order,
        # so reverse before storing.
        tk = [pop!(flat_t) for _ in 1:counts[idx]]
        td = [pop!(flat_d) for _ in 1:counts[idx]]
        t_packed[idx] = round.(reverse!(tk); digits = 1)
        d_packed[idx] = round.(reverse!(td); digits = 1)
    end
    return t_packed, d_packed
end

"""
    _trim_taps!(t_vec, d_vec, counts, max_ms)

Remove taps that fall outside the time window `[0, max_ms]` from each per-trial vector in-place (iterating backwards to allow safe deletion). `counts[i]` is decremented for every removed tap.
"""
function _trim_taps!(
    t_vec::Vector{Vector{Float64}},
    d_vec::Vector{Vector{Float64}},
    counts::Vector{Int64},
    max_ms::Float64,
)::Nothing
    for idx1 in eachindex(t_vec)
        for idx2 in length(t_vec[idx1]):-1:1
            if t_vec[idx1][idx2] > max_ms
                deleteat!(t_vec[idx1], idx2)
                deleteat!(d_vec[idx1], idx2)
                counts[idx1] -= 1
            end
        end
    end
    return nothing
end

"""
    _dedup_taps!(t_vec, d_vec, counts)

Remove duplicate time-point entries from each per-trial vector in-place. `counts[i]` is updated to reflect the number of unique taps remaining. Duplicates can arise from serial-port bounce or rapid repeated events.
"""
function _dedup_taps!(
    t_vec::Vector{Vector{Float64}},
    d_vec::Vector{Vector{Float64}},
    counts::Vector{Int64},
)::Nothing
    for idx in eachindex(t_vec)
        keep = if length(unique(t_vec[idx])) != length(t_vec[idx])
            # Keep only the first occurrence of each time point
            unique(i -> t_vec[idx][i], eachindex(t_vec[idx]))
        else
            collect(eachindex(t_vec[idx]))  # no duplicates - keep all
        end
        counts[idx]  = length(keep)
        t_vec[idx]   = t_vec[idx][keep]
        d_vec[idx]   = d_vec[idx][keep]
    end
    return nothing
end

# =============================================================================

"""
    iftt(; <keyword arguments>)

Perform the **Finger Tapping Test (FTT)** in **GUI mode**.

The subject taps a key as fast as possible during each *trial* window, and intentionally suppresses tapping during each *interval* window. Both active and suppressed taps are recorded.

**Input devices (in priority order):**

1. Serial port - a push-button wired to a Raspberry Pi GPIO and connected to the PC via USB (virtual serial port).
2. Keyboard - the `SPACEBAR` key (fallback when no serial port is configured or the port cannot be opened).

# Arguments

- `duration::Int64=20`: single trial duration in seconds
- `trials::Int64=2`: number of trials
- `interval::Int64=2`: rest interval between trials in seconds (set to `0` to omit rest intervals)
- `gpio::Int64=-1`: Raspberry Pi GPIO pin the switch is wired to (e.g. `gpio=23` for board pin 16); **required** when `port_name` is set
- `port_name::String=""`: serial port path (e.g. `"/dev/ttyACM0"`); leave empty (`""`) to use the keyboard

# Returns

A `NamedTuple` with fields:

| Field        | Type                      | Description                            |
|--------------|---------------------------|----------------------------------------|
| `taps`       | `Vector{Int64}`           | tap count per trial                    |
| `tap_t`      | `Vector{Vector{Float64}}` | tap onset times [ms] per trial         |
| `tap_d`      | `Vector{Vector{Float64}}` | tap durations [ms] per trial           |
| `taps_int`   | `Vector{Int64}`           | tap count per interval                 |
| `tap_t_int`  | `Vector{Vector{Float64}}` | tap onset times [ms] per interval      |
| `tap_d_int`  | `Vector{Vector{Float64}}` | tap durations [ms] per interval        |

All times are relative to the **start of their respective trial or interval**.
"""
function iftt(;
    duration::Int64   = 20,
    trials::Int64     = 2,
    interval::Int64   = 2,
    gpio::Int64       = -1,
    port_name::String = "",
)::@NamedTuple{
    taps::Vector{Int64},
    tap_t::Vector{Vector{Float64}},
    tap_d::Vector{Vector{Float64}},
    taps_int::Vector{Int64},
    tap_t_int::Vector{Vector{Float64}},
    tap_d_int::Vector{Vector{Float64}},
}
    # only throw when a port is given but no GPIO is specified
    port_name != "" && gpio == -1 &&
        throw(ArgumentError("gpio must be specified when port_name is set."))

    # probe serial port availability; fall back to keyboard if unavailable
    sp = nothing
    if port_name != ""
        sp = _serial_open(port_name)
        if sp === nothing
            _info("Serial port $port_name is not available - keyboard SPACEBAR will be used")
            port_name = ""
        else
            _serial_close(sp)
        end
    end

    img_idle  = read_from_png(joinpath(res_path, "finger_noclick.png"))
    img_press = read_from_png(joinpath(res_path, "finger_click.png"))

    # per-trial accumulators
    result       = zeros(Int64, trials)   # tap count per trial
    int_result   = zeros(Int64, trials)   # tap count per interval

    # keyboard-mode: nested per-trial vectors (pushed as each trial completes)
    t_kp         = Vector{Vector{Float64}}()
    d_kp         = Vector{Vector{Float64}}()
    int_t_kp     = Vector{Vector{Float64}}()
    int_d_kp     = Vector{Vector{Float64}}()

    # serial-mode: flat chronological vectors (restructured in post-processing)
    t_kp_flat    = Vector{Float64}()
    d_kp_flat    = Vector{Float64}()
    int_t_kp_flat = Vector{Float64}()
    int_d_kp_flat = Vector{Float64}()

    # trial/interval start times (keyboard mode - used for relative time offsets)
    t_trial_start = Vector{Float64}()
    t_int_start   = Vector{Float64}()

    # intra-trial buffers (keyboard mode - swapped out at the start of each trial)
    t_kp_tmp     = Vector{Float64}()
    d_kp_tmp     = Vector{Float64}()
    int_t_kp_tmp = Vector{Float64}()
    int_d_kp_tmp = Vector{Float64}()

    key_pressed  = false  # debounce flag shared across key handlers

    # =========================================================================
    function _activate(app)
        win = GtkApplicationWindow(app, "NeuroAnalyzer: iftt()")
        Gtk4.default_size(win, Int64(img_idle.width), Int64(img_idle.height) + 100)

        # canvas that shows the finger graphic (idle or pressed)
        can = GtkCanvas()
        can.content_width  = Int64(img_idle.width)
        can.content_height = Int64(img_idle.height)

        g1 = GtkGrid()
        g1.column_homogeneous = true
        g1.column_spacing = 20
        g1.row_spacing    = 20
        g1.margin_start   = 5
        g1.margin_end     = 5
        g1.margin_top     = 5
        g1.margin_bottom  = 5

        bt_start      = GtkButton("START")
        bt_start.tooltip_text = "Start the test"

        lb_status1    = GtkLabel("Status:");        lb_status1.halign = 2
        lb_status2    = GtkLabel("READY TO START"); lb_status2.halign = 1
        lb_trial1     = GtkLabel("Trial #:");       lb_trial1.halign  = 2
        lb_trial2     = GtkLabel("-");              lb_trial2.halign  = 1
        lb_interval1  = GtkLabel("Interval #:");    lb_interval1.halign = 2
        lb_interval2  = GtkLabel("-");              lb_interval2.halign = 1

        g1[1:3, 1] = can
        g1[1, 2]   = lb_status1;  g1[3, 2] = lb_status2
        g1[1, 3]   = lb_trial1;   g1[3, 3] = lb_trial2
        g1[1, 4]   = lb_interval1; g1[3, 4] = lb_interval2
        g1[1:3, 5] = GtkLabel("")
        g1[1:3, 6] = bt_start

        vbox = GtkBox(:v)
        push!(vbox, g1)
        push!(win, vbox)
        Gtk4.show(win)

        # initial canvas draw - idle finger image
        @guarded draw(can) do widget
            ctx = getgc(can)
            Cairo.set_source_surface(ctx, img_idle, 0, 0)
            Cairo.paint(ctx)
        end

        # --- keyboard event handlers (active only during TEST / INTERVAL phases) ---
        win_key = Gtk4.GtkEventControllerKey(win)

        signal_connect(win_key, "key-pressed") do widget, keyval, keycode, state
            if keyval == 32  # SPACEBAR
                if lb_status2.label == "TEST" && !key_pressed
                    push!(t_kp_tmp, time())
                    key_pressed = true
                elseif lb_status2.label == "INTERVAL" && !key_pressed
                    push!(int_t_kp_tmp, time())
                    key_pressed = true
                end
            end
            return sleep(0.1)  # simple debounce
        end

        signal_connect(win_key, "key-released") do widget, keyval, keycode, state
            if keyval == 32  # SPACEBAR
                if lb_status2.label == "TEST" && key_pressed
                    push!(d_kp_tmp, time())
                    key_pressed = false
                elseif lb_status2.label == "INTERVAL" && key_pressed
                    push!(int_d_kp_tmp, time())
                    key_pressed = false
                end
            end
            return sleep(0.1)
        end

        # --- START button ---
        return signal_connect(bt_start, "clicked") do widget
            bt_start.sensitive = false

            if port_name == ""
                # ---- keyboard input path ----
                Threads.@spawn begin
                    for idx in 1:trials
                        # reset intra-trial buffers for this trial
                        t_kp_tmp     = Vector{Float64}()
                        d_kp_tmp     = Vector{Float64}()
                        int_t_kp_tmp = Vector{Float64}()
                        int_d_kp_tmp = Vector{Float64}()

                        # begin trial
                        _beep()
                        @idle_add @guarded draw(can) do widget
                            ctx = getgc(can)
                            Cairo.set_source_surface(ctx, img_press, 0, 0)
                            Cairo.paint(ctx)
                        end
                        @idle_add lb_status2.label  = "TEST"
                        @idle_add lb_trial2.label   = strip("$idx of $trials")
                        @idle_add lb_interval2.label = "-"

                        push!(t_trial_start, time())
                        sleep(duration)
                        _beep()

                        # commit this trial's taps
                        push!(t_kp, t_kp_tmp)
                        push!(d_kp, d_kp_tmp)
                        result[idx] = length(t_kp_tmp)

                        # rest interval (skip if interval == 0)
                        if interval > 0
                            @idle_add @guarded draw(can) do widget
                                ctx = getgc(can)
                                Cairo.set_source_surface(ctx, img_idle, 0, 0)
                                Cairo.paint(ctx)
                            end
                            @idle_add lb_status2.label   = "INTERVAL"
                            @idle_add lb_trial2.label    = "-"
                            @idle_add lb_interval2.label = strip("$idx of $trials")

                            push!(t_int_start, time())
                            sleep(interval)

                            push!(int_t_kp, int_t_kp_tmp)
                            push!(int_d_kp, int_d_kp_tmp)
                            int_result[idx] = length(int_t_kp_tmp)
                        end
                    end
                    @idle_add close(win)
                end

            else
                # ---- serial port input path ----
                Threads.@spawn begin
                    for idx in 1:trials
                        _beep()
                        @idle_add @guarded draw(can) do widget
                            ctx = getgc(can)
                            Cairo.set_source_surface(ctx, img_press, 0, 0)
                            Cairo.paint(ctx)
                        end
                        @idle_add lb_status2.label   = "TEST"
                        @idle_add lb_trial2.label    = strip("$idx of $trials")
                        @idle_add lb_interval2.label = "-"

                        key_pressed = false
                        sp = _serial_open(port_name)
                        t_trial = time()
                        while time() - t_trial <= duration
                            t = time() - t_trial
                            serial_key = _serial_listener(sp)
                            if serial_key == "$gpio:1" && !key_pressed
                                push!(t_kp_flat, t)
                                result[idx] += 1
                                key_pressed = true
                            elseif serial_key == "$gpio:0" && key_pressed
                                push!(d_kp_flat, t)
                                key_pressed = false
                            end
                            sleep(0.1)
                        end
                        _serial_close(sp)
                        _beep()

                        # drop an unmatched press at the very end of the window
                        if length(d_kp_flat) < sum(result)
                            pop!(t_kp_flat)
                            result[idx] -= 1
                        end

                        if interval > 0
                            @idle_add @guarded draw(can) do widget
                                ctx = getgc(can)
                                Cairo.set_source_surface(ctx, img_idle, 0, 0)
                                Cairo.paint(ctx)
                            end
                            @idle_add lb_status2.label   = "INTERVAL"
                            @idle_add lb_interval2.label = strip("$idx of $trials")
                            @idle_add lb_trial2.label    = "-"

                            key_pressed = false
                            sp = _serial_open(port_name)
                            t_int = time()
                            while time() - t_int <= interval
                                t = time() - t_int
                                serial_key = _serial_listener(sp)
                                if serial_key == "$gpio:1" && !key_pressed
                                    push!(int_t_kp_flat, t)
                                    int_result[idx] += 1
                                    key_pressed = true
                                elseif serial_key == "$gpio:0" && key_pressed
                                    push!(int_d_kp_flat, t)
                                    key_pressed = false
                                end
                                sleep(0.1)
                            end
                            _serial_close(sp)

                            if length(int_d_kp_flat) < sum(int_result)
                                pop!(int_t_kp_flat)
                                int_result[idx] -= 1
                            end
                        end
                    end
                    @idle_add close(win)
                end
            end
        end
    end # _activate
    # =========================================================================

    app = GtkApplication("org.neuroanalyzer.iftt")
    Gtk4.signal_connect(_activate, app, :activate)
    Gtk4.GLib.stop_main_loop()
    Gtk4.run(app)

    # =========================================================================
    # post-processing (runs after the GTK window closes)
    # =========================================================================

    if port_name == ""
        # ---- keyboard path: t_kp/d_kp are already nested per-trial ----

        # align press/release counts (a tap released after window end has no release)
        for idx in 1:trials
            if length(t_kp[idx]) != length(d_kp[idx])
                l = min(length(t_kp[idx]), length(d_kp[idx]))
                t_kp[idx]  = t_kp[idx][1:l]
                d_kp[idx]  = d_kp[idx][1:l]
                result[idx] = l
            end
            if length(int_t_kp[idx]) != length(int_d_kp[idx])
                l = min(length(int_t_kp[idx]), length(int_d_kp[idx]))
                int_t_kp[idx]  = int_t_kp[idx][1:l]
                int_d_kp[idx]  = int_d_kp[idx][1:l]
                int_result[idx] = l
            end

            # convert absolute epoch times → duration since trial/interval start [ms]
            d_kp[idx]     = round.((d_kp[idx]     .- t_kp[idx])          .* 1000; digits = 1)
            t_kp[idx]     = round.((t_kp[idx]     .- t_trial_start[idx]) .* 1000; digits = 1)
            int_d_kp[idx] = round.((int_d_kp[idx] .- int_t_kp[idx])      .* 1000; digits = 1)
            int_t_kp[idx] = round.((int_t_kp[idx] .- t_int_start[idx])   .* 1000; digits = 1)
        end

        # remove taps that slipped past the trial/interval end boundary
        _trim_taps!(t_kp, d_kp, result, Float64(duration * 1000))
        _trim_taps!(int_t_kp, int_d_kp, int_result, Float64(interval * 1000))

        return (
            taps = result, tap_t = t_kp, tap_d = d_kp,
            taps_int = int_result, tap_t_int = int_t_kp, tap_d_int = int_d_kp,
        )

    else
        # ---- serial path: flat vectors need restructuring into per-trial ----

        # convert raw seconds → milliseconds and compute durations from press times
        d_kp_flat     = round.((d_kp_flat     .- t_kp_flat)     .* 1000; digits = 1)
        int_d_kp_flat = round.((int_d_kp_flat .- int_t_kp_flat) .* 1000; digits = 1)
        t_kp_flat     = round.(t_kp_flat     .* 1000; digits = 1)
        int_t_kp_flat = round.(int_t_kp_flat .* 1000; digits = 1)

        # pack flat vectors into per-trial nested vectors
        t_keypressed, d_keypressed = _pack_taps(t_kp_flat,     d_kp_flat,     result)
        int_t_keypressed, int_d_keypressed = _pack_taps(int_t_kp_flat, int_d_kp_flat, int_result)

        _trim_taps!(t_keypressed, d_keypressed, result, Float64(duration * 1000))
        _trim_taps!(int_t_keypressed, int_d_keypressed, int_result, Float64(interval * 1000))

        # remove any duplicate time points introduced by serial-port bounce
        _dedup_taps!(t_keypressed, d_keypressed, result)
        _dedup_taps!(int_t_keypressed, int_d_keypressed, int_result)

        return (
            taps     = result,     tap_t     = t_keypressed,     tap_d     = d_keypressed,
            taps_int = int_result, tap_t_int = int_t_keypressed, tap_d_int = int_d_keypressed,
        )
    end
end

"""
    ftt(; duration, trials, interval, gpio, port_name)

Perform the **Finger Tapping Test (FTT)** in **CLI (terminal) mode**.

Functionally identical to [`iftt`](@ref) but without a graphical window - the test runs entirely in the terminal. Suitable for headless or scripted use.

**Input devices (in priority order):**

1. Raspberry Pi GPIO - direct hardware button via `pigpiod` daemon.
2. Serial port - button wired to RPi GPIO, connected via USB.
3. Keyboard - `SPACEBAR` fallback.

!!! note "Keyboard tap durations"

When using the keyboard, only tap *count* and *onset time* are reliably captured. Durations are fixed at **100 ms** as a placeholder.

# Arguments

- `duration::Int64=20`: single trial duration in seconds
- `trials::Int64=2`: number of trials
- `interval::Int64=2`: rest interval between trials in seconds (set to `0` to omit rest intervals)
- `gpio::Int64=-1`: Raspberry Pi GPIO pin the switch is wired to (e.g. `gpio=23` for board pin 16); **required** when `port_name` is set
- `port_name::String=""`: serial port path (e.g. `"/dev/ttyACM0"`); leave empty (`""`) to use the keyboard

# Returns

A `NamedTuple` with fields:

| Field        | Type                      | Description                            |
|--------------|---------------------------|----------------------------------------|
| `taps`       | `Vector{Int64}`           | tap count per trial                    |
| `tap_t`      | `Vector{Vector{Float64}}` | tap onset times [ms] per trial         |
| `tap_d`      | `Vector{Vector{Float64}}` | tap durations [ms] per trial           |
| `taps_int`   | `Vector{Int64}`           | tap count per interval                 |
| `tap_t_int`  | `Vector{Vector{Float64}}` | tap onset times [ms] per interval      |
| `tap_d_int`  | `Vector{Vector{Float64}}` | tap durations [ms] per interval        |

All times are relative to the **start of their respective trial or interval**.
"""
function ftt(;
    duration::Int64   = 20,
    trials::Int64     = 2,
    interval::Int64   = 2,
    gpio::Int64       = -1,
    port_name::String = "",
)::@NamedTuple{
    taps::Vector{Int64},
    tap_t::Vector{Vector{Float64}},
    tap_d::Vector{Vector{Float64}},
    taps_int::Vector{Int64},
    tap_t_int::Vector{Vector{Float64}},
    tap_d_int::Vector{Vector{Float64}},
}
    # only throw when a port is given but no GPIO is specified
    port_name != "" && gpio == -1 &&
        throw(ArgumentError("gpio must be specified when port_name is set."))

    sp  = nothing
    rpi = false

    if gpio != -1 && port_name == ""
        # direct Raspberry Pi GPIO mode via pigpiod
        rpi = _check_rpi()
        if rpi != false
            set_mode(rpi, gpio, PiGPIO.INPUT)
        else
            _info("Could not detect pigpiod daemon - keyboard SPACEBAR will be used")
        end
    elseif port_name != ""
        # serial port mode
        sp = _serial_open(port_name)
        if sp === nothing
            _info("Serial port $port_name is not available - keyboard SPACEBAR will be used")
            port_name = ""
        else
            _serial_close(sp)
        end
    end

    # print session header
    println("NeuroTester: FTT")
    println("================")
    println("  Trials: $trials")
    println("Duration: $duration [seconds]")
    println("Interval: $interval [seconds]")
    if rpi isa PiGPIO.Pi
        println("  Button: RPi GPIO $gpio")
    elseif port_name != ""
        println("  Button: serial port $port_name GPIO $gpio")
    else
        println("  Button: SPACEBAR")
    end
    println()

    # wait for user to start
    if rpi isa PiGPIO.Pi || !isnothing(sp)
        println("Ready to start - press the BUTTON to begin")
    else
        println("Ready to start - press SPACEBAR to begin")
    end
    println()

    if !(rpi isa PiGPIO.Pi) && isnothing(sp)
        # keyboard: block until SPACEBAR
        while true
            ret = ccall(:jl_tty_set_mode, Int32, (Ptr{Cvoid}, Int32), stdin.handle, true)
            ret == 0 || error("Unable to switch terminal to raw mode.")
            kbd_key = read(stdin, Char)
            ccall(:jl_tty_set_mode, Int32, (Ptr{Cvoid}, Int32), stdin.handle, false)
            kbd_key == ' ' && break
        end
    elseif rpi isa PiGPIO.Pi
        # RPi GPIO: block until button pressed
        while true
            PiGPIO.read(rpi, gpio) != false && break
            sleep(0.001)  # avoid pegging the CPU
        end
    elseif !isnothing(sp)
        # serial: block until button event arrives
        sp = _serial_open(port_name)
        while true
            a = _serial_listener(sp)
            (a !== nothing && a == "$gpio:1") && break
        end
        _serial_close(sp)
    end

    print("The test will start after a beep")
    sleep(1)

    # =========================================================================
    # per-run accumulators
    # =========================================================================
    result     = zeros(Int64, trials)
    int_result = zeros(Int64, trials)

    t_kp      = Vector{Float64}()
    d_kp      = Vector{Float64}()
    int_t_kp  = Vector{Float64}()
    int_d_kp  = Vector{Float64}()

    key_pressed = false

    # =========================================================================
    # main test loop - branched by input device
    # =========================================================================

    if !(rpi isa PiGPIO.Pi) && isnothing(sp)
        # ---- keyboard input path ----

        # build a timeline of segment boundaries [ms] for all trials + intervals
        # odd segments (1, 3, 5, …) are trials; even segments (2, 4, 6, …) are intervals
        l_seg      = duration + interval                            # length of one trial+interval block
        n_segs     = 2 * trials + 1
        t_segments = zeros(n_segs)
        for idx in 1:trials
            t_segments[(idx * 2) - 1] = l_seg * (idx - 1)           # trial start
            t_segments[idx * 2]       = (l_seg * idx) - interval    # trial end / interval start
        end
        t_segments .*= 1000
        t_segments[end] = ((trials * duration) + (trials * interval)) * 1000  # total end

        channel = Channel(_kbd_listener, 1024)   # async keyboard event producer
        stop    = false
        r       = 0
        t_raw   = Float64[]    # raw key-press times relative to t_s [ms]
        seg_idx = 1            # pointer into t_segments for the next segment boundary
        trial_n = 1            # trial counter for status printing
        int_n   = 1            # interval counter for status printing

        t_s = time()
        while !stop
            sleep(0.1)

            # drain all key presses that arrived since last iteration
            while !isempty(channel)
                c = take!(channel)
                if c == ' '
                    r += 1
                    push!(t_raw, time() - t_s)
                end
            end

            # check whether a new segment boundary has been crossed
            if (time() - t_s) * 1000 >= t_segments[seg_idx]
                _beep()
                println()
                if iseven(seg_idx)
                    println()
                    print("Interval $int_n: DO NOT press the SPACEBAR button")
                    int_n += 1
                else
                    println()
                    print("   Trial $trial_n: press the SPACEBAR as quickly as possible")
                    trial_n += 1
                end
                seg_idx += 1
            end

            # stop when we reach the total end time
            if (time() - t_s) * 1000 >= t_segments[end]
                close(channel)
                stop = true
            end
        end

        # assign each key press to a trial or interval based on its timestamp
        t_raw = round.(t_raw .* 1000; digits = 3)
        for press_ms in t_raw
            for seg in 1:(2 * trials)
                if t_segments[seg] <= press_ms <= t_segments[seg + 1]
                    if iseven(seg + 1)       # odd segment index → trial
                        idx3 = (seg + 1) ÷ 2
                        result[idx3] += 1
                        push!(t_kp, press_ms)
                        push!(d_kp, 100.0)   # keyboard: duration fixed at 100 ms
                    else                     # even segment index → interval
                        idx3 = seg ÷ 2
                        int_result[idx3] += 1
                        push!(int_t_kp, press_ms)
                        push!(int_d_kp, 100.0)
                    end
                end
            end
        end
        println()

    elseif !isnothing(sp)
        # ---- serial input path ----
        println()
        for idx in 1:trials
            _beep()
            println()
            print("   Trial $idx: press the BUTTON as quickly as possible")

            key_pressed = false
            sp = _serial_open(port_name)
            t_trial = time()
            while time() - t_trial <= duration
                t = time() - t_trial
                serial_key = _serial_listener(sp)
                if serial_key == "$gpio:1" && !key_pressed
                    push!(t_kp, t);  result[idx] += 1;  key_pressed = true
                elseif serial_key == "$gpio:0" && key_pressed
                    push!(d_kp, t);  key_pressed = false
                end
                sleep(0.1)
            end
            _serial_close(sp)
            _beep()
            if length(d_kp) < sum(result)
                pop!(t_kp);  result[idx] -= 1
            end

            println(); println()
            print("Interval $idx: DO NOT press the BUTTON")

            key_pressed = false
            sp = _serial_open(port_name)
            t_int = time()
            while time() - t_int <= interval
                t = time() - t_int
                serial_key = _serial_listener(sp)
                if serial_key == "$gpio:1" && !key_pressed
                    push!(int_t_kp, t);  int_result[idx] += 1;  key_pressed = true
                elseif serial_key == "$gpio:0" && key_pressed
                    push!(int_d_kp, t);  key_pressed = false
                end
                sleep(0.1)
            end
            _serial_close(sp)
            if length(int_d_kp) < sum(int_result)
                pop!(int_t_kp);  int_result[idx] -= 1
            end
            println()
        end

    elseif rpi isa PiGPIO.Pi
        # ---- Raspberry Pi direct GPIO path ----
        debounce_ms = 50   # minimum time between state changes to count as new event [ms]
        println()
        for idx in 1:trials
            _beep()
            println()
            print("   Trial $idx: press the BUTTON as quickly as possible")

            key_state        = 0
            key_last_state   = 0
            last_debounce_ms = 0.0

            t_trial = time()
            while time() <= t_trial + duration
                t       = time() - t_trial
                rpi_key = PiGPIO.read(rpi, gpio)

                rpi_key != key_last_state && (last_debounce_ms = time() * 1000)

                if (time() * 1000 - last_debounce_ms) > debounce_ms && rpi_key != key_state
                    key_state = rpi_key
                    if key_state == 1
                        push!(t_kp, t);  result[idx] += 1
                    else
                        push!(d_kp, t)
                    end
                end
                key_last_state = rpi_key
                sleep(0.001)
            end

            _beep()
            if length(d_kp) < sum(result)
                pop!(t_kp);  result[idx] -= 1
            end

            println(); println()
            print("Interval $idx: DO NOT press the BUTTON")

            key_state        = 0
            key_last_state   = 0
            last_debounce_ms = 0.0

            t_int = time()
            while time() <= t_int + interval
                rpi_key = PiGPIO.read(rpi, gpio)
                rpi_key != key_last_state && (last_debounce_ms = time() * 1000)
                if (time() * 1000 - last_debounce_ms) > debounce_ms && rpi_key != key_state
                    t = time() - t_int
                    key_state = rpi_key
                    if key_state == 1
                        push!(int_t_kp, t);  int_result[idx] += 1
                    else
                        push!(int_d_kp, t)
                    end
                end
                key_last_state = rpi_key
                sleep(0.001)
            end

            if length(int_d_kp) < sum(int_result)
                pop!(int_t_kp);  int_result[idx] -= 1
            end
            println()
        end
    end

    println()
    println("Testing completed")

    # =========================================================================
    # post-processing
    # =========================================================================

    # convert raw seconds → ms and compute durations (RPi / serial paths only;
    # keyboard path already stores ms and fixed 100 ms durations)
    if rpi isa PiGPIO.Pi || !isnothing(sp)
        d_kp     = round.((d_kp     .- t_kp)     .* 1000; digits = 1)
        int_d_kp = round.((int_d_kp .- int_t_kp) .* 1000; digits = 1)
        t_kp     = round.(t_kp     .* 1000; digits = 1)
        int_t_kp = round.(int_t_kp .* 1000; digits = 1)
    end

    # pack flat vectors into per-trial nested vectors
    t_keypressed, d_keypressed = _pack_taps(t_kp, d_kp, result)
    int_t_keypressed, int_d_keypressed = _pack_taps(int_t_kp, int_d_kp, int_result)

    # for the keyboard path, times are global (relative to t_s); subtract per-trial
    # offsets to make them relative to each trial / interval start
    if isnothing(sp) && !(rpi isa PiGPIO.Pi)
        for idx in eachindex(t_keypressed)
            offset = (idx - 1) * (duration + interval) * 1000
            t_keypressed[idx] = round.(t_keypressed[idx] .- offset; digits = 3)
        end
        for idx in eachindex(int_t_keypressed)
            offset = (idx * duration + (idx - 1) * interval) * 1000
            int_t_keypressed[idx] = round.(int_t_keypressed[idx] .- offset; digits = 1)
        end
    end

    _trim_taps!(t_keypressed,     d_keypressed,     result,     Float64(duration * 1000))
    _trim_taps!(int_t_keypressed, int_d_keypressed, int_result, Float64(interval * 1000))

    # Remove duplicate time points (can arise from serial bounce)
    _dedup_taps!(t_keypressed, d_keypressed, result)
    _dedup_taps!(int_t_keypressed, int_d_keypressed, int_result)

    return (
        taps = result, tap_t = t_keypressed, tap_d = d_keypressed,
        taps_int = int_result, tap_t_int = int_t_keypressed, tap_d_int = int_d_keypressed,
    )
end
