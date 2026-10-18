export iftt
export ftt

# =============================================================================
# shared helpers used by both iftt() and ftt()
# =============================================================================

# fixed tap duration [s] used for the terminal keyboard, which reports key presses only
const _FTT_KBD_TAP_DURATION = 0.1

"""
    _ftt_taps(press, release)

Tap onsets and durations from press / release times (seconds since the start of the trial or interval). A final press without a release (still pressed at the end of the window) is dropped.

# Returns

- `Int64`: number of taps
- `Vector{Float64}`: tap onsets [ms]
- `Vector{Float64}`: tap durations [ms]
"""
function _ftt_taps(
    press::Vector{Float64},
    release::Vector{Float64},
)::Tuple{Int64, Vector{Float64}, Vector{Float64}}
    n = min(length(press), length(release))
    return n, press[1:n] .* 1000, (release[1:n] .- press[1:n]) .* 1000
end

# results of the completed trials / intervals as returned by iftt() and ftt()
function _ftt_result(
    tr_press::Vector{Vector{Float64}},
    tr_rel::Vector{Vector{Float64}},
    int_press::Vector{Vector{Float64}},
    int_rel::Vector{Vector{Float64}},
    n::Int64,
)::@NamedTuple{
    taps::Vector{Int64},
    tap_t::Vector{Vector{Float64}},
    tap_d::Vector{Vector{Float64}},
    taps_int::Vector{Int64},
    tap_t_int::Vector{Vector{Float64}},
    tap_d_int::Vector{Vector{Float64}},
}
    tr = [_ftt_taps(tr_press[idx], tr_rel[idx]) for idx = 1:n]
    it = [_ftt_taps(int_press[idx], int_rel[idx]) for idx = 1:n]
    return (
        taps = Int64[r[1] for r in tr],
        tap_t = Vector{Float64}[r[2] for r in tr],
        tap_d = Vector{Float64}[r[3] for r in tr],
        taps_int = Int64[r[1] for r in it],
        tap_t_int = Vector{Float64}[r[2] for r in it],
        tap_d_int = Vector{Float64}[r[3] for r in it],
    )
end

function _ftt_check(duration::Int64, trials::Int64, interval::Int64, gpio::Int64, port_name::String)::Nothing
    duration > 0 || throw(ArgumentError("duration must be > 0."))
    trials > 0 || throw(ArgumentError("trials must be > 0."))
    interval >= 0 || throw(ArgumentError("interval must be ≥ 0."))
    port_name != "" && gpio == -1 &&
        throw(ArgumentError("gpio must be specified when port_name is set."))
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

- `duration::Int64=20`: single trial duration in seconds; must be > 0
- `trials::Int64=2`: number of trials; must be > 0
- `interval::Int64=2`: rest interval after each trial in seconds (set to `0` to omit rest intervals)
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

All times are relative to the **start of their respective trial or interval**. If the window is closed before the end of the test, only the completed trials are returned.
"""
function iftt(;
    duration::Int64 = 20,
    trials::Int64 = 2,
    interval::Int64 = 2,
    gpio::Int64 = -1,
    port_name::String = "",
)::@NamedTuple{
    taps::Vector{Int64},
    tap_t::Vector{Vector{Float64}},
    tap_d::Vector{Vector{Float64}},
    taps_int::Vector{Int64},
    tap_t_int::Vector{Vector{Float64}},
    tap_d_int::Vector{Vector{Float64}},
}
    _ftt_check(duration, trials, interval, gpio, port_name)

    # the serial port stays open for the whole test; keyboard if it is not available
    sp = port_name == "" ? nothing : _serial_try_open(port_name)
    port_name != "" && isnothing(sp) && _info("Keyboard SPACEBAR will be used")

    img_idle = read_from_png(joinpath(res_path, "finger_noclick.png"))
    img_press = read_from_png(joinpath(res_path, "finger_click.png"))

    # press / release times [s] since the start of each trial and interval
    tr_press = [Float64[] for _ = 1:trials]
    tr_rel = [Float64[] for _ = 1:trials]
    int_press = [Float64[] for _ = 1:trials]
    int_rel = [Float64[] for _ = 1:trials]

    # state shared by the test task and the keyboard handlers
    phase = :idle # :trial or :interval while taps are recorded
    cur = 0 # current trial
    t0 = 0.0 # start of the current trial / interval
    key_pressed = false
    n_done = 0 # completed trials (trial + interval)
    aborted = false # set when the window is closed
    test_task = nothing

    # record one trial or interval of `len` seconds
    function _ftt_phase(ph::Symbol, idx::Int64, len::Int64)
        press, rel = ph === :trial ? (tr_press[idx], tr_rel[idx]) : (int_press[idx], int_rel[idx])
        key_pressed = false
        cur = idx
        t0 = time()
        phase = ph
        while time() - t0 < len && !aborted
            if isnothing(sp)
                # keyboard events are recorded by the key handlers
                sleep(0.01)
            else
                key_pressed = _serial_button!(sp, gpio, t0, press, rel, key_pressed)
                sleep(0.001)
            end
        end
        phase = :idle
        return nothing
    end

    # =========================================================================
    function _activate(app)
        win = GtkApplicationWindow(app, "NeuroAnalyzer: iftt()")
        Gtk4.default_size(win, Int64(img_idle.width), Int64(img_idle.height) + 100)

        # canvas that shows the finger graphic (idle or pressed)
        can = GtkCanvas()
        can.content_width = Int64(img_idle.width)
        can.content_height = Int64(img_idle.height)

        g1 = GtkGrid()
        g1.column_homogeneous = true
        g1.column_spacing = 20
        g1.row_spacing = 20
        g1.margin_start = 5
        g1.margin_end = 5
        g1.margin_top = 5
        g1.margin_bottom = 5

        bt_start = GtkButton("START")
        bt_start.tooltip_text = "Start the test"

        lb_status1 = GtkLabel("Status:")
        lb_status1.halign = 2
        lb_status2 = GtkLabel("READY TO START")
        lb_status2.halign = 1
        lb_trial1 = GtkLabel("Trial #:")
        lb_trial1.halign = 2
        lb_trial2 = GtkLabel("-")
        lb_trial2.halign = 1
        lb_interval1 = GtkLabel("Interval #:")
        lb_interval1.halign = 2
        lb_interval2 = GtkLabel("-")
        lb_interval2.halign = 1

        g1[1:3, 1] = can
        g1[1, 2] = lb_status1
        g1[3, 2] = lb_status2
        g1[1, 3] = lb_trial1
        g1[3, 3] = lb_trial2
        g1[1, 4] = lb_interval1
        g1[3, 4] = lb_interval2
        g1[1:3, 5] = GtkLabel("")
        g1[1:3, 6] = bt_start

        vbox = GtkBox(:v)
        push!(vbox, g1)
        push!(win, vbox)
        Gtk4.show(win)

        # closing the window stops the test
        signal_connect(win, "close-request") do widget
            aborted = true
            return false
        end

        function _show(img)
            @guarded draw(can) do widget
                ctx = getgc(can)
                Cairo.set_source_surface(ctx, img, 0, 0)
                return Cairo.paint(ctx)
            end
            return nothing
        end
        _show(img_idle)

        # --- keyboard event handlers (record only during a trial / interval) ---
        win_key = Gtk4.GtkEventControllerKey(win)

        signal_connect(win_key, "key-pressed") do widget, keyval, keycode, state
            # SPACEBAR; key auto-repeat is ignored while the key is held
            if keyval == 32 && isnothing(sp) && phase !== :idle && !key_pressed
                push!(phase === :trial ? tr_press[cur] : int_press[cur], time() - t0)
                key_pressed = true
            end
            return nothing
        end

        signal_connect(win_key, "key-released") do widget, keyval, keycode, state
            if keyval == 32 && isnothing(sp) && phase !== :idle && key_pressed
                push!(phase === :trial ? tr_rel[cur] : int_rel[cur], time() - t0)
                key_pressed = false
            end
            return nothing
        end

        # --- START button ---
        return signal_connect(bt_start, "clicked") do widget
            bt_start.sensitive = false
            test_task = Threads.@spawn begin
                try
                    for idx = 1:trials
                        _beep()
                        @idle_add _show(img_press)
                        @idle_add lb_status2.label = "TEST"
                        @idle_add lb_trial2.label = "$idx of $trials"
                        @idle_add lb_interval2.label = "-"
                        _ftt_phase(:trial, idx, duration)
                        aborted && break
                        _beep()

                        # rest interval (skip if interval == 0)
                        if interval > 0
                            @idle_add _show(img_idle)
                            @idle_add lb_status2.label = "INTERVAL"
                            @idle_add lb_trial2.label = "-"
                            @idle_add lb_interval2.label = "$idx of $trials"
                            _ftt_phase(:interval, idx, interval)
                            aborted && break
                        end
                        n_done = idx
                    end
                    aborted || @idle_add close(win)
                catch err
                    @error "FTT failed" exception = (err, catch_backtrace())
                end
            end
        end
    end # _activate
    # =========================================================================

    app = GtkApplication("org.neuroanalyzer.iftt")
    Gtk4.signal_connect(_activate, app, :activate)
    Gtk4.GLib.stop_main_loop()
    Gtk4.run(app)

    # the window is closed: stop the test task (if any) and release the port
    aborted = true
    isnothing(test_task) || wait(test_task)
    isnothing(sp) || _serial_close(sp)

    n_done < trials && _info("Test interrupted: $n_done of $trials trials completed")

    return _ftt_result(tr_press, tr_rel, int_press, int_rel, n_done)
end

"""
    ftt(; <keyword arguments>)

Perform the **Finger Tapping Test (FTT)** in **CLI (terminal) mode**.

Functionally identical to [`iftt`](@ref) but without a graphical window - the test runs entirely in the terminal. Suitable for headless or scripted use.

**Input devices (in priority order):**

1. Raspberry Pi GPIO - direct hardware button via `pigpiod` daemon (`gpio` set, `port_name` empty).
2. Serial port - button wired to RPi GPIO, connected via USB.
3. Keyboard - `SPACEBAR` fallback.

!!! note "Keyboard tap durations"

When using the keyboard, only tap *count* and *onset time* are captured (the terminal reports key presses, not releases). Durations are fixed at **100 ms** as a placeholder.

# Arguments

- `duration::Int64=20`: single trial duration in seconds; must be > 0
- `trials::Int64=2`: number of trials; must be > 0
- `interval::Int64=2`: rest interval after each trial in seconds (set to `0` to omit rest intervals)
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
    duration::Int64 = 20,
    trials::Int64 = 2,
    interval::Int64 = 2,
    gpio::Int64 = -1,
    port_name::String = "",
)::@NamedTuple{
    taps::Vector{Int64},
    tap_t::Vector{Vector{Float64}},
    tap_d::Vector{Vector{Float64}},
    taps_int::Vector{Int64},
    tap_t_int::Vector{Vector{Float64}},
    tap_d_int::Vector{Vector{Float64}},
}
    _ftt_check(duration, trials, interval, gpio, port_name)

    # input device: :rpi (pigpiod), :serial (kept open for the whole test) or :keyboard
    sp = nothing
    rpi = false
    if gpio != -1 && port_name == ""
        rpi = _check_rpi()
        if rpi isa PiGPIO.Pi
            set_mode(rpi, gpio, PiGPIO.INPUT)
        else
            _info("Could not detect pigpiod daemon - keyboard SPACEBAR will be used")
        end
    elseif port_name != ""
        sp = _serial_try_open(port_name)
        isnothing(sp) && _info("Keyboard SPACEBAR will be used")
    end
    device = rpi isa PiGPIO.Pi ? :rpi : (isnothing(sp) ? :keyboard : :serial)

    # press / release times [s] since the start of each trial and interval
    tr_press = [Float64[] for _ = 1:trials]
    tr_rel = [Float64[] for _ = 1:trials]
    int_press = [Float64[] for _ = 1:trials]
    int_rel = [Float64[] for _ = 1:trials]

    # keyboard: key presses with the time they were read
    channel = device === :keyboard ? Channel{Tuple{Char, Float64}}(c -> _kbd_listener(c; timestamp = true), 1024) : nothing

    # record one trial or interval of `len` seconds
    function _ftt_phase(press::Vector{Float64}, rel::Vector{Float64}, len::Int64)
        t0 = time()
        if device === :serial
            pressed = false
            while time() - t0 < len
                pressed = _serial_button!(sp, gpio, t0, press, rel, pressed)
                sleep(0.001)
            end
        elseif device === :rpi
            debounce_ms = 50 # minimum time between state changes to count as a new event [ms]
            key_state = 0
            key_last_state = 0
            last_debounce_ms = 0.0
            while time() - t0 < len
                rpi_key = PiGPIO.read(rpi, gpio)
                rpi_key != key_last_state && (last_debounce_ms = time() * 1000)
                if (time() * 1000 - last_debounce_ms) > debounce_ms && rpi_key != key_state
                    key_state = rpi_key
                    push!(key_state == 1 ? press : rel, time() - t0)
                end
                key_last_state = rpi_key
                sleep(0.001)
            end
        else
            # key presses made before this phase (e.g. during the beep) are discarded
            while time() - t0 < len
                while isready(channel)
                    c, tk = take!(channel)
                    (c == ' ' && t0 <= tk < t0 + len) || continue
                    push!(press, tk - t0)
                    push!(rel, tk - t0 + _FTT_KBD_TAP_DURATION)
                end
                sleep(0.005)
            end
        end
        return nothing
    end

    btn = device === :keyboard ? "SPACEBAR" : "BUTTON"
    try
        # print session header
        println("NeuroTester: FTT")
        println("================")
        println("  Trials: $trials")
        println("Duration: $duration [seconds]")
        println("Interval: $interval [seconds]")
        if device === :rpi
            println("  Button: RPi GPIO $gpio")
        elseif device === :serial
            println("  Button: serial port $port_name GPIO $gpio")
        else
            println("  Button: SPACEBAR")
        end
        println()
        println("Ready to start - press the $btn to begin")
        println()

        # wait for the user to start
        if device === :keyboard
            while true
                c, _ = take!(channel)
                c == ' ' && break
            end
        elseif device === :rpi
            while PiGPIO.read(rpi, gpio) == false
                sleep(0.001)
            end
        else
            pressed = false
            t_start = Float64[]
            while isempty(t_start)
                pressed = _serial_button!(sp, gpio, time(), t_start, Float64[], pressed)
                sleep(0.001)
            end
        end

        print("The test will start after a beep")
        sleep(1)
        println()

        for idx = 1:trials
            _beep()
            println()
            print("   Trial $idx: press the $btn as quickly as possible")
            _ftt_phase(tr_press[idx], tr_rel[idx], duration)
            _beep()
            println()
            if interval > 0
                println()
                print("Interval $idx: DO NOT press the $btn")
                _ftt_phase(int_press[idx], int_rel[idx], interval)
                println()
            end
        end
    finally
        isnothing(sp) || _serial_close(sp)
        isnothing(channel) || close(channel)
    end

    println()
    println("Testing completed")

    return _ftt_result(tr_press, tr_rel, int_press, int_rel, trials)
end
