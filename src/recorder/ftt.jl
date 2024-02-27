export iftt
export ftt

"""
    iftt(; duration, trials, interval)

Perform Finger Tapping Test (FTT) in GUI mode. Use computer keyboard (SPACE key) or switch panel attached to Raspberry Pi via a GPIO pin. Number of taps, time points and durations of taps are recorded. Also, taps during intervals (when the study subject should suppress tapping) are recorded.

# Arguments

- `duration::Int64=10`: single trial duration in seconds
- `trials::Int64=6`: number of trials
- `interval::Int64=10`: interval between trials in seconds

# Returns

Named tuple containing:
- `taps::Vector{Int64}`: number of taps per trial
- `tap_t::Vector{Vector{Float64}}`: taps time point [s]
- `tap_d::Vector{Vector{Float64}}`: taps duration [s]
- `taps_int::Vector{Int64}`: number of taps per trial during intervals
- `tap_t_int::Vector{Vector{Float64}}`: taps time point [s] during intervals
- `tap_d_int::Vector{Vector{Float64}}`: taps duration [s] during intervals
"""
function iftt(; duration::Int64=10, trials::Int64=6, interval::Int64=10)

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
                # sleep(duration)
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

Perform Finger Tapping Test (FTT) in CLI mode. Switch panel attached to Raspberry Pi via a GPIO pin. Number of taps, time points and durations of taps are recorded. Also, taps during intervals (when the study subject should suppress tapping) are recorded.

# Arguments

- `duration::Int64=10`: single trial duration in seconds
- `trials::Int64=6`: number of trials
- `interval::Int64=10`: interval between trials in seconds
- `gpio::Int64=23`: Raspberry Pi GPIO to which the switch is connected (default is GPIO 23 = BOARD 16 pin)

# Returns

Named tuple containing:
- `taps::Vector{Int64}`: number of taps per trial
- `tap_t::Vector{Vector{Float64}}`: taps time point [s]
- `tap_d::Vector{Vector{Float64}}`: taps duration [s]
- `taps_int::Vector{Int64}`: number of taps per trial during intervals
- `tap_t_int::Vector{Vector{Float64}}`: taps time point [s] during intervals
- `tap_d_int::Vector{Vector{Float64}}`: taps duration [s] during intervals
"""
function ftt(; duration::Int64=10, trials::Int64=6, interval::Int64=10, gpio::Int64=23)

    _wip()

    # check if running on RPi
    rpi = _checkrpi()
    if rpi == true
        rpi = Pi()
        set_mode(rpi, gpio, PiGPIO.INPUT)
    else
        @info "Could not detect pigpiod daemon, will use keyboard spacebar key"
    end

    println("NeuroRecorder: FTT")
    println("==================")
    println("  Trials: $trials")
    println("Duration: $duration [seconds]")
    println("Interval: $interval [seconds]")
    if rpi != false
        println("  Button: RPi GPIO $gpio")
    else
        println("  Button: SPACEBAR")
    end

    println("")
    println("Ready to start, press the BUTTON to begin the test")
    
    # while (rpi != false && PiGPIO.read(rpi, gpio) == false) && (_getch() != 32)
    rpi_key = false
    kbd_key = Channel(_getch, 10) # Start task, 10 is buffer size for channel
    c = Char
    while c != ' '
        if rpi != false
            rpi_key = PiGPIO.read(rpi, gpio)
        end
        rpi_key != false && break
        while !isempty(kbd_key) # Process all keypresses
            c = take!(kbd_key)
            println(c)
            if c == ' '
                println("quitting")
                close(channel)
                break
            end
        end
    end
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

    t_s = time()
    for idx in 1:trials
        _beep()
        @info "   Trial: $idx - press the button as quickly as possible"
        t1 = time()
        # sleep(duration)
        while time() <= t1 + duration
            sleep(0.1)
            rpi != false && (rpi_key = PiGPIO.read(rpi, gpio))
            kbd_key = _getch()
            t2 = time()
            if (rpi_key == true) || (kbd_key == 32)
                if key_pressed == false
                    print("_")
                    # key is pressed                        
                    push!(t_kp, t2)
                    result[idx] += 1
                    key_pressed = true
                    continue
                end
            else
                if key_pressed == true
                    print("-")
                    # key is released
                    push!(d_kp, t2 - t_kp[end])
                    key_pressed = false
                    continue
                end
            end
        end
        _beep()
        if length(d_kp) < sum(result)
            pop!(t_kp)
            result[idx] -= 1
        end
        @info "Interval: $idx - DO NOT press the button"
        t1 = time()
        while time() <= t1 + interval
            rpi != false && (rpi_key = PiGPIO.read(rpi, gpio))
            kbd_key = _getch()
            t2 = time()
            if (rpi_key == true) || (kbd_key == 32)
                if key_pressed == false
                    # key is pressed
                    print("_")
                    push!(int_t_kp, t2)
                    int_result[idx] += 1
                    key_pressed = true
                    sleep(0.05)
                    continue
                end
            else
                if key_pressed == true
                    # key is released
                    print("-")
                    push!(int_d_kp, t2 - int_t_kp[end])
                    key_pressed = false
                    sleep(0.05)
                    continue
                end
            end
        end
        if length(int_d_kp) < sum(int_result)
            pop!(int_t_kp)
            int_result[idx] -= 1
        end
    end

    println("Testing completed.")

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
