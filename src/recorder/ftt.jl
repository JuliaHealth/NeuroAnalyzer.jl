export iftt
export ftt

"""
    iftt(; duration, trials, interval)

Perform Finger Tapping Test (FTT) in GUI mode. Use computer keyboard (SPACEBAR key) or switch panel attached to Raspberry Pi via a GPIO pin. Number of taps, time points and durations of taps are recorded. Also, taps during intervals (when the study subject should suppress tapping) are recorded.

# Arguments

- `duration::Int64=5`: single trial duration in seconds
- `trials::Int64=2`: number of trials
- `interval::Int64=2`: interval between trials in seconds

# Returns

Named tuple containing:
- `taps::Vector{Int64}`: number of taps per trial
- `tap_t::Vector{Vector{Float64}}`: taps time point [s]
- `tap_d::Vector{Vector{Float64}}`: taps duration [s]
- `taps_int::Vector{Int64}`: number of taps per trial during intervals
- `tap_t_int::Vector{Vector{Float64}}`: taps time point [s] during intervals
- `tap_d_int::Vector{Vector{Float64}}`: taps duration [s] during intervals
"""
function iftt(; duration::Int64=5, trials::Int64=2, interval::Int64=2)

    img1 = read_from_png(joinpath(res_path, "finger_noclick.png"))
    img2 = read_from_png(joinpath(res_path, "finger_click.png"))

    win = GtkWindow("NeuroRecorder: iftt()", img1.width, img1.height + 100)
    set_gtk_property!(win, :border_width, 20)
    set_gtk_property!(win, :resizable, false)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    set_gtk_property!(win, :startup_id, "org.neuroanalyzer")

    can = GtkCanvas(img1.width, img1.height)

    g1 = GtkGrid()
    set_gtk_property!(g1, :column_homogeneous, true)
    set_gtk_property!(g1, :column_spacing, 20)
    set_gtk_property!(g1, :row_spacing, 20)
    g2 = GtkGrid()
    set_gtk_property!(g2, :column_homogeneous, true)
    set_gtk_property!(g2, :column_spacing, 20)
    set_gtk_property!(g2, :row_spacing, 20)

    bt_start = GtkButton("START")
    set_gtk_property!(bt_start, :tooltip_text, "Start the test")
    
    lb_status1 = GtkLabel("Status:")
    lb_status2 = GtkLabel("READY TO START")
    set_gtk_property!(lb_status1, :halign, 2)
    set_gtk_property!(lb_status2, :halign, 1)
    lb_trial1 = GtkLabel("Trial #:")
    lb_trial2 = GtkLabel("-")
    set_gtk_property!(lb_trial1, :halign, 2)
    set_gtk_property!(lb_trial2, :halign, 1)
    lb_interval1 = GtkLabel("Interval #:")
    lb_interval2 = GtkLabel("-")
    set_gtk_property!(lb_interval1, :halign, 2)
    set_gtk_property!(lb_interval2, :halign, 1)

    g1[1:3, 1] = can
    g1[1, 2] = lb_status1
    g1[3, 2] = lb_status2
    g1[1, 3] = lb_trial1
    g1[3, 3] = lb_trial2
    g1[1, 4] = lb_interval1
    g1[3, 4] = lb_interval2
    g1[1:3, 5] = GtkLabel("")
    g1[1:3, 6] = bt_start
    vbox1 = GtkBox(:v)
    push!(vbox1, g1)

    push!(win, vbox1)
    showall(win)

    @guarded draw(can) do widget
        ctx = getgc(can)
        Cairo.set_source_surface(ctx, img1, 0, 0)
        Cairo.paint(ctx)
    end

    result = zeros(Int64, trials)
    int_result = zeros(Int64, trials)
    t_s = nothing
    t_kp = Vector{Float64}()
    d_kp = Vector{Float64}()
    int_t_kp = Vector{Float64}()
    int_d_kp = Vector{Float64}()
    key_pressed = false

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        if k == 32 && get_gtk_property(lb_status2, :label, String) == "TEST"
            t = parse(Int64, split(get_gtk_property(lb_trial2, :label, String), ' ')[1])
            kp = key_pressed
            if kp == false
                push!(t_kp, time())
                result[t] += 1
                key_pressed = true
            end
        elseif k == 32 && get_gtk_property(lb_status2, :label, String) == "INTERVAL"
            t = parse(Int64, split(get_gtk_property(lb_interval2, :label, String), ' ')[1])
            kp = key_pressed
            if kp == false
                push!(int_t_kp, time())
                int_result[t] += 1
                key_pressed = true
            end
        end
    end

    signal_connect(win, "key-release-event") do widget, event
        k = event.keyval
        if k == 32 && get_gtk_property(lb_status2, :label, String) == "TEST"
            t = parse(Int64, split(get_gtk_property(lb_trial2, :label, String), ' ')[1])
            kp = key_pressed
            if kp == true
                push!(d_kp, time() - t_kp[end])
                key_pressed = false
            end
        elseif k == 32 && get_gtk_property(lb_status2, :label, String) == "INTERVAL"
            t = parse(Int64, split(get_gtk_property(lb_interval2, :label, String), ' ')[1])
            kp = key_pressed
            if kp == true
                push!(int_d_kp, time() - int_t_kp[end])
                key_pressed = false
            end
        end
    end

    signal_connect(bt_start, "clicked") do widget
        set_gtk_property!(bt_start, :sensitive, false)
        t_s = time()
        Threads.@spawn begin
            for idx in 1:trials
                _beep()
                @guarded draw(can) do widget
                    ctx = getgc(can)
                    Cairo.set_source_surface(ctx, img2, 0, 0)
                    Cairo.paint(ctx)
                end
                set_gtk_property!(lb_status2, :label, "TEST")
                l = strip(string(idx) * " of " * string(trials))
                set_gtk_property!(lb_trial2, :label, l)
                set_gtk_property!(lb_interval2, :label, "-")
                t1 = time()
                while time() <= t1 + duration
                end
                _beep()
                if length(d_kp) < sum(result)
                    pop!(t_kp)
                    result[idx] -= 1
                end
                @guarded draw(can) do widget
                    ctx = getgc(can)
                    Cairo.set_source_surface(ctx, img1, 0, 0)
                    Cairo.paint(ctx)
                end
                set_gtk_property!(lb_status2, :label, "INTERVAL")
                set_gtk_property!(lb_trial2, :label, "-")
                l = strip(string(idx) * " of " * string(trials))
                set_gtk_property!(lb_interval2, :label, l)
                t1 = time()
                while time() <= t1 + interval
                end
                if length(int_d_kp) < sum(int_result)
                    pop!(int_t_kp)
                    int_result[idx] -= 1
                end
            end

            # Interacting with GTK from a thread other than the main thread is
            # generally not allowed, so we register an idle callback instead.
            Gtk.GLib.g_idle_add(nothing) do user_data
                Gtk.destroy(win)
            end
        end
    end

    cnd = Condition()
    signal_connect(win, :destroy) do widget
        notify(cnd)
    end
    @async Gtk.gtk_main()
    wait(cnd)

    # format time points
    t_kp = t_kp .- t_s
    t_keypressed = Vector{Vector{Float64}}()
    d_keypressed = Vector{Vector{Float64}}()
    for idx1 in trials:-1:1
        tk = zeros(result[idx1])
        td = zeros(result[idx1])
        for idx2 in 1:result[idx1]
            tk[idx2] = pop!(t_kp)
            td[idx2] = pop!(d_kp)
        end
        reverse!(tk)
        reverse!(td)
        tk = tk .- ((idx1 - 1) * (duration + interval))
        push!(t_keypressed, tk)
        push!(d_keypressed, td)
    end
    reverse!(t_keypressed)
    reverse!(d_keypressed)

    # format time points
    int_t_kp = int_t_kp .- t_s
    int_t_keypressed = Vector{Vector{Float64}}()
    int_d_keypressed = Vector{Vector{Float64}}()
    for idx1 in trials:-1:1
        tk = zeros(int_result[idx1])
        td = zeros(int_result[idx1])
        for idx2 in 1:int_result[idx1]
            tk[idx2] = pop!(int_t_kp)
            td[idx2] = pop!(int_d_kp)
        end
        reverse!(tk)
        reverse!(td)
        tk = tk .- ((idx1 * duration) + ((idx1 - 1) * interval))
        push!(int_t_keypressed, tk)
        push!(int_d_keypressed, td)
    end
    reverse!(int_t_keypressed)
    reverse!(int_d_keypressed)

    # remove out of time boundary taps
    for idx1 in 1:trials
        for idx2 in length(t_keypressed[idx1]):-1:1
            if t_keypressed[idx1][idx2] > duration
                deleteat!(t_keypressed[idx1], idx2)
                deleteat!(d_keypressed[idx1], idx2)
                result[idx1] -= 1
            end
        end
    end
    for idx1 in 1:trials
        for idx2 in length(int_t_keypressed[idx1]):-1:1
            if int_t_keypressed[idx1][idx2] > interval
                deleteat!(int_t_keypressed[idx1], idx2)
                deleteat!(int_d_keypressed[idx1], idx2)
                int_result[idx1] -= 1
            end
        end
    end

    return (taps=result, tap_t=t_keypressed, tap_d=d_keypressed, taps_int=int_result, tap_t_int=int_t_keypressed, tap_d_int=int_d_keypressed)

end

"""
    ftt(; duration, trials, interval, gpio)

Perform Finger Tapping Test (FTT) in CLI mode. Use computer keyboard (SPACEBAR key) or switch panel attached to Raspberry Pi via a GPIO pin. Number of taps, time points and durations of taps are recorded. Also, taps during intervals (when the study subject should suppress tapping) are recorded. When using computer keyboard, only the number of taps and their time points are recorded; tap durations are set to -1.

# Arguments

- `duration::Int64=5`: single trial duration in seconds
- `trials::Int64=2`: number of trials
- `interval::Int64=2`: interval between trials in seconds
- `gpio::Int64=-1`: Raspberry Pi GPIO to which the switch is connected (e.g. `gpio=23` is Pi board pin 16); set to -1 to disable using GPIO
- `port_name::String=""`: serial port to which the switch is connected (e.g. `/dev/ttyACM0`); set to "" to disable using serial port

# Returns

Named tuple containing:
- `taps::Vector{Int64}`: number of taps per trial
- `tap_t::Vector{Vector{Float64}}`: taps time point [s]
- `tap_d::Vector{Vector{Float64}}`: taps duration [s]
- `taps_int::Vector{Int64}`: number of taps per trial during intervals
- `tap_t_int::Vector{Vector{Float64}}`: taps time point [s] during intervals
- `tap_d_int::Vector{Vector{Float64}}`: taps duration [s] during intervals
"""
function ftt(; duration::Int64=5, trials::Int64=2, interval::Int64=2, gpio::Int64=-1, port_name::String="")

    @assert !(port_name != "" && gpio == -1) "If serial port is used, GPIO must be specified."

    sp = nothing

    if gpio != -1 && port_name == ""
        # check if running on RPi
        rpi = _check_rpi()
        if rpi != false
            set_mode(rpi, gpio, PiGPIO.INPUT)
        else
            @info "Could not detect pigpiod daemon, keyboard SPACEBAR key will be used"
        end
    elseif port_name != ""
        sp = _serial_open(port_name)
        if sp === nothing
            @info "Serial port $port_name is not available, keyboard SPACEBAR key will be used"
            port_name = ""
        else
            _serial_close(sp)
        end
        rpi = false
    else
        @info "GPIO disabled, keyboard SPACEBAR key will be used"
        rpi = false
    end
    
    println("NeuroRecorder: FTT")
    println("==================")
    println("  Trials: $trials")
    println("Duration: $duration [seconds]")
    println("Interval: $interval [seconds]")
    if rpi != false
        println("  Button: RPi GPIO $gpio")
    elseif port_name != ""
        println("  Button: serial port $port_name GPIO $gpio")
    else
        println("  Button: SPACEBAR")
    end

    println()
    if rpi != false || sp !== nothing
        println("Ready to start, press the BUTTON to begin the test")
    else
        println("Ready to start, press SPACEBAR to begin the test")
    end
    println()

    if rpi == false && sp === nothing
        while true
            ret = ccall(:jl_tty_set_mode, Int32, (Ptr{Cvoid},Int32), stdin.handle, true)
            ret == 0 || error("Unable to switch to raw mode.")
            kbd_key = read(stdin, Char)
            ccall(:jl_tty_set_mode, Int32, (Ptr{Cvoid},Int32), stdin.handle, false)
            if kbd_key == ' '
                break
            else
                continue
            end
        end
    elseif rpi != false
        while true
            PiGPIO.read(rpi, gpio) != false && break
        end
    elseif sp !== nothing
        sp = _serial_open(port_name)
        while true
            a = nothing
            a = _serial_listener(sp)
            (a !== nothing && a == "$gpio:1") && break
        end
        _serial_close(sp)
    end
    print("The test will start after a beep")
    sleep(1)

    result = zeros(Int64, trials)
    int_result = zeros(Int64, trials)
    t_s = nothing
    t_kp = Vector{Float64}()
    d_kp = Vector{Float64}()
    int_t_kp = Vector{Float64}()
    int_d_kp = Vector{Float64}()

    key_pressed = false
    rpi_key = false
    kbd_key = nothing

    if rpi == false && sp === nothing
        # use computer keyboard
        t_s = time()
        # calculate segments time points
        t_segments = zeros(1 + 2 * trials)
        l_seg = duration + interval
        for idx in 1:trials
            t_segments[(idx * 2) - 1] = l_seg * (idx - 1)
        end
        for idx in 1:trials
            t_segments[idx * 2] = (l_seg * idx) - interval
        end
        t_segments .+= t_s
        t_e = (trials * duration) + (trials * interval)
        t_segments[end] = t_s + t_e
        channel = Channel(_kbd_listener, 1024) # Start task, 10 is buffer size for channel
        stop = false
        r = 0
        t = Float64[]
        idx = 1
        idx1 = 1
        idx2 = 1
        while !stop
            sleep(0.001)
            if time() >= t_segments[end]
                stop = true
                close(channel)
                break
            end
            while !isempty(channel) # process all keypresses
                c = take!(channel)
                if c == ' '
                    r += 1
                    push!(t, time())
                end
            end
            if time() >= t_segments[idx]
                _beep()
                println()
                if iseven(idx)
                    println()
                    print("Interval $idx2: DO NOT press the SPACEBAR button")
                    idx2 += 1
                else
                    println()
                    print("   Trial $idx1: press the SPACEBAR button as quickly as possible")
                    idx1 += 1
                end
                idx += 1
            end
        end
        # format time points
        t .-= t_segments[1]
        t = round.(t, digits=4)
        t_segments .-= t_segments[1]
        if r > 0
            for idx1 in 1:r
                for idx2 in 1:(2 * trials)
                    if t[idx1] in t_segments[idx2]:0.001:t_segments[idx2 + 1]
                        if iseven(idx2 + 1)
                            idx3 = (idx2 + 1) ÷ 2
                            result[idx3] += 1
                            push!(t_kp, t[idx1])
                            push!(d_kp, -1)
                        end
                        if isodd(idx2 + 1)
                            idx3 = idx2 ÷ 2
                            int_result[idx3] += 1
                            push!(int_t_kp, t[idx1])
                            push!(int_d_kp, -1)
                        end
                    end
                end
            end
        end
        println()
    elseif sp !== nothing
        # use serial port
        println()
        t_s = time()
        for idx in 1:trials
            _beep()
            println()
            print("   Trial $idx: press the BUTTON button as quickly as possible")
            t1 = time()
            key_pressed = false
            sp = _serial_open(port_name)
            while time() <= t1 + duration
                serial_key = _serial_listener(sp)
                t2 = time()
                if serial_key == "$gpio:1"
                    if key_pressed == false
                        # key is pressed
                        push!(t_kp, t2)
                        result[idx] += 1
                        key_pressed = true
                        sleep(0.05)
                        continue
                    end
                elseif serial_key == "$gpio:0"
                    if key_pressed == true
                        # key is released
                        push!(d_kp, t2 - t_kp[end])
                        key_pressed = false
                        sleep(0.05)
                        continue
                    end
                end
            end
            _serial_close(sp)
            _beep()
            if length(d_kp) < sum(result)
                pop!(t_kp)
                result[idx] -= 1
            end
            println()
            println()
            print("Interval $idx: DO NOT press the BUTTON button")
            t1 = time()
            key_pressed = false
            sp = _serial_open(port_name)
            while time() <= t1 + interval
                serial_key = _serial_listener(sp)
                t2 = time()
                if serial_key == "$gpio:1"
                    if key_pressed == false
                        # key is pressed
                        push!(int_t_kp, t2)
                        int_result[idx] += 1
                        key_pressed = true
                        sleep(0.01)
                        continue
                    end
                elseif serial_key == "$gpio:0"
                    if key_pressed == true
                        # key is released
                        push!(int_d_kp, t2 - int_t_kp[end])
                        key_pressed = false
                        sleep(0.01)
                        continue
                    end
                end
            end
            _serial_close(sp)
            if length(int_d_kp) < sum(int_result)
                pop!(int_t_kp)
                int_result[idx] -= 1
            end
            println()
        end
        # format time points
        t_kp = t_kp .- t_s
        int_t_kp = int_t_kp .- t_s
    elseif rpi !== false
        # use RPi
        debounce_delay = 50 # ms
        println()
        t_s = time()
        for idx in 1:trials
            _beep()
            println()
            print("   Trial $idx: press the BUTTON button as quickly as possible")
            t1 = time()
            key_state = false
            key_last_state = false
            key_pressed = false
            last_debounce_time = 0
            while time() <= t1 + duration
                rpi_key = PiGPIO.read(rpi, gpio)
                rpi_key != key_last_state && (last_debounce_time = time())
                if ((time() - last_debounce_time) > debounce_delay)
                    if rpi_key != key_state
                        key_state = rpi_key
                        if key_state == true
                            # key is pressed
                            print("_")
                            push!(t_kp, time())
                            result[idx] += 1
                        else
                            # key is released
                            print("-")
                            push!(d_kp, time() - t_kp[end])
                        end
                    end
                end
                key_last_state = rpi_key;
    #=
                t2 = time()
                if rpi_key == true
                    if key_pressed == false
                        # key is pressed
                        push!(t_kp, t2)
                        result[idx] += 1
                        key_pressed = true
                        sleep(0.15)
                        continue
                    end
                else
                    if key_pressed == true
                        # key is released
                        push!(d_kp, t2 - t_kp[end])
                        key_pressed = false
                        sleep(0.15)
                        continue
                    end
                end
=#
            end
            @show result
            @show t_kp
            @show d_kp
            return nothing
            _beep()
            if length(d_kp) < sum(result)
                pop!(t_kp)
                result[idx] -= 1
            end
            println()
            println()
            print("Interval $idx: DO NOT press the BUTTON button")
            t1 = time()
            key_pressed = false
            while time() <= t1 + interval
                rpi_key = PiGPIO.read(rpi, gpio)
                t2 = time()
                if rpi_key == true
                    if key_pressed == false
                        # key is pressed
                        push!(int_t_kp, t2)
                        int_result[idx] += 1
                        key_pressed = true
                        sleep(0.15)
                        continue
                    end
                else
                    if key_pressed == true
                        # key is released
                        push!(int_d_kp, t2 - int_t_kp[end])
                        key_pressed = false
                        sleep(0.15)
                        continue
                    end
                end
            end
            if length(int_d_kp) < sum(int_result)
                pop!(int_t_kp)
                int_result[idx] -= 1
            end
            println()
        end
        # format time points
        t_kp = t_kp .- t_s
        int_t_kp = int_t_kp .- t_s
    end

    println()
    println("Testing completed.")

    t_keypressed = Vector{Vector{Float64}}()
    d_keypressed = Vector{Vector{Float64}}()
    for idx1 in trials:-1:1
        tk = zeros(result[idx1])
        td = zeros(result[idx1])
        for idx2 in 1:result[idx1]
            tk[idx2] = pop!(t_kp)
            td[idx2] = pop!(d_kp)
        end
        reverse!(tk)
        reverse!(td)
        tk = tk .- ((idx1 - 1) * (duration + interval))
        push!(t_keypressed, round.(tk, digits=4))
        push!(d_keypressed, round.(td, digits=4))
    end
    reverse!(t_keypressed)
    reverse!(d_keypressed)

    # format time points
    int_t_keypressed = Vector{Vector{Float64}}()
    int_d_keypressed = Vector{Vector{Float64}}()
    for idx1 in trials:-1:1
        tk = zeros(int_result[idx1])
        td = zeros(int_result[idx1])
        for idx2 in 1:int_result[idx1]
            tk[idx2] = pop!(int_t_kp)
            td[idx2] = pop!(int_d_kp)
        end
        reverse!(tk)
        reverse!(td)
        tk = tk .- ((idx1 * duration) + ((idx1 - 1) * interval))
        push!(int_t_keypressed, round.(tk, digits=4))
        push!(int_d_keypressed, round.(td, digits=4))
    end
    reverse!(int_t_keypressed)
    reverse!(int_d_keypressed)

    # remove out of time boundary taps
    for idx1 in 1:trials
        for idx2 in length(t_keypressed[idx1]):-1:1
            if t_keypressed[idx1][idx2] > duration
                deleteat!(t_keypressed[idx1], idx2)
                deleteat!(d_keypressed[idx1], idx2)
                result[idx1] -= 1
            end
        end
    end
    for idx1 in 1:trials
        for idx2 in length(int_t_keypressed[idx1]):-1:1
            if int_t_keypressed[idx1][idx2] > interval
                deleteat!(int_t_keypressed[idx1], idx2)
                deleteat!(int_d_keypressed[idx1], idx2)
                int_result[idx1] -= 1
            end
        end
    end

    # remove duplicates
    d_idx = Vector{Vector{Int64}}()
    for idx in 1:length(t_keypressed)
        if length(unique(t_keypressed[idx])) != length(t_keypressed[idx])
            push!(d_idx, unique(i -> t_keypressed[idx][i], eachindex(t_keypressed[idx])))
        else
            push!(d_idx, eachindex(t_keypressed[idx]))
        end
        result[idx] = length(d_idx[idx])
        t_keypressed[idx] = t_keypressed[idx][d_idx[idx]]
        d_keypressed[idx] = d_keypressed[idx][d_idx[idx]]
    end
    d_idx = Vector{Vector{Int64}}()
    for idx in 1:length(int_t_keypressed)
        if length(unique(int_t_keypressed[idx])) != length(int_t_keypressed[idx])
            push!(d_idx, unique(i -> int_t_keypressed[idx][i], eachindex(int_t_keypressed[idx])))
        else
            push!(d_idx, eachindex(int_t_keypressed[idx]))
        end
        int_result[idx] = length(d_idx[idx])
        int_t_keypressed[idx] = int_t_keypressed[idx][d_idx[idx]]
        int_d_keypressed[idx] = int_d_keypressed[idx][d_idx[idx]]
    end

    return (taps=result, tap_t=t_keypressed, tap_d=d_keypressed, taps_int=int_result, tap_t_int=int_t_keypressed, tap_d_int=int_d_keypressed)

end
