export itpt
export tpt

"""
    itpt(; <keyword arguments>)

Perform Two-point Pinch Test (TPT) in GUI mode. TPT is recorded using MMA7660 accelerometer via Arduino attached to the PC via USB cable (virtual serial port). Sampling rate is 50 Hz.

# Arguments

- `duration::Int64=20`: recording duration in seconds
- `trials::Int64=2`: number of trials
- `interval::Int64=2`: interval between trials in seconds
- `port_name::String="/dev/ttyUSB0"`: serial port to which the Arduino is connected

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function itpt(; duration::Int64=20, trials::Int64=2, interval::Int64=2, port_name::String="/dev/ttyUSB0")::NeuroAnalyzer.NEURO

    sp = _serial_open(port_name, baudrate=19200)
    @assert !isnothing(sp) _info("Serial port $port_name is not available")

    img1 = read_from_png(joinpath(res_path, "finger_nopinch.png"))
    img2 = read_from_png(joinpath(res_path, "finger_pinch.png"))

    win = GtkWindow("NeuroRecorder: itpt()", img1.width, img1.height + 100)
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

    # sampling rate is 50 Hz = 20 ms per loop
    fs = 50
    t = round.(collect(0:1/fs:((duration * trials) + (interval * (trials - 1)))), digits=3)
    tpt_ch_x = zeros(length(t))
    tpt_ch_y = zeros(length(t))
    tpt_ch_z = zeros(length(t))
    tpt_ch_accx = zeros(length(t))
    tpt_ch_accy = zeros(length(t))
    tpt_ch_accz = zeros(length(t))

    signal_connect(bt_start, "clicked") do widget
        set_gtk_property!(bt_start, :sensitive, false)
        Threads.@spawn begin
            for trial_idx in 1:trials
                _beep()
                @guarded draw(can) do widget
                    ctx = getgc(can)
                    Cairo.set_source_surface(ctx, img2, 0, 0)
                    Cairo.paint(ctx)
                end
                set_gtk_property!(lb_status2, :label, "TEST")
                l = strip(string(trial_idx) * " of " * string(trials))
                set_gtk_property!(lb_trial2, :label, l)
                set_gtk_property!(lb_interval2, :label, "-")
                idx = 1
                while idx <= length(duration * fs)
                    t1 = time() * 1000
                    sp_signal = _serial_listener(sp)
                    if !isnothing(sp_signal)
                        m = match(r"(tpt\: )(\-*[0-9]+) (\-*[0-9]+) (\-*[0-9]+) (\-*[0-9]+\.[0-9]+) (\-*[0-9]+\.[0-9]+) (\-*[0-9]+\.[0-9]+)", sp_signal)
                        if !isnothing(m)
                            if length(m.captures) == 7
                                tpt_ch_x[idx] = parse(Float64, m.captures[2])
                                tpt_ch_y[idx] = parse(Float64, m.captures[3])
                                tpt_ch_z[idx] = parse(Float64, m.captures[4])
                                tpt_ch_accx[idx] = parse(Float64, m.captures[5])
                                tpt_ch_accy[idx] = parse(Float64, m.captures[6])
                                tpt_ch_accz[idx] = parse(Float64, m.captures[7])
                                idx += 1
                            end
                        end
                    end
                    t2 = time() * 1000
                    sleep((20 - (t2 - t1)) / 1000)
                end
                _beep()
                @guarded draw(can) do widget
                    ctx = getgc(can)
                    Cairo.set_source_surface(ctx, img1, 0, 0)
                    Cairo.paint(ctx)
                end
                set_gtk_property!(lb_status2, :label, "INTERVAL")
                l = strip(string(trial_idx) * " of " * string(trials))
                set_gtk_property!(lb_interval2, :label, l)
                set_gtk_property!(lb_trial2, :label, "-")
                idx = 1
                while idx <= length(interval * fs)
                    sp_signal = _serial_listener(sp)
                    if !isnothing(sp_signal)
                        m = match(r"(tpt\: )(\-*[0-9]+) (\-*[0-9]+) (\-*[0-9]+) (\-*[0-9]+\.[0-9]+) (\-*[0-9]+\.[0-9]+) (\-*[0-9]+\.[0-9]+)", sp_signal)
                        if !isnothing(m)
                            if length(m.captures) == 7
                                tpt_ch_x[idx] = parse(Float64, m.captures[2])
                                tpt_ch_y[idx] = parse(Float64, m.captures[3])
                                tpt_ch_z[idx] = parse(Float64, m.captures[4])
                                tpt_ch_accx[idx] = parse(Float64, m.captures[5])
                                tpt_ch_accy[idx] = parse(Float64, m.captures[6])
                                tpt_ch_accz[idx] = parse(Float64, m.captures[7])
                                idx += 1
                            end
                        end
                    end
                end
            end
            _serial_close(sp)
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

    # tpt_ch_x = tpt_ch_x[1:(end - 1)]
    # tpt_ch_y = tpt_ch_y[1:(end - 1)]
    # tpt_ch_z = tpt_ch_z[1:(end - 1)]
    # tpt_ch_accx = tpt_ch_accx[1:(end - 1)]
    # tpt_ch_accy = tpt_ch_accy[1:(end - 1)]
    # tpt_ch_accz = tpt_ch_accz[1:(end - 1)]
    # t = round.(t[1:(end - 1)], digits=3)
    tpt_signal = Matrix([tpt_ch_x tpt_ch_y tpt_ch_z tpt_ch_accx tpt_ch_accy tpt_ch_accz]')
    tpt_signal = reshape(tpt_signal, 6, :, 1)

    obj = create_object(data_type="tpt")
    add_channel!(obj, data=tpt_signal, label=["pos_x", "pos_y", "pos_z", "acc_x", "acc_y", "acc_z"], type=["orient", "orient", "orient", "accel", "accel", "accel"], unit=["", "", "", "m/s²", "m/s²", "m/s²"])
    create_time!(obj, fs=fs)

    return obj

end

"""
    tpt(; <keyword arguments>)

Perform Two-point Pinch Test (TPT) in CLI mode. TPT is recorded using MMA7660 accelerometer via Arduino attached to the PC via USB cable (virtual serial port). Sampling rate is 50 Hz.

# Arguments

- `duration::Int64=20`: recording duration in seconds
- `port_name::String="/dev/ttyUSB0"`: serial port to which the Arduino is connected

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function tpt(; duration::Int64=20, port_name::String="/dev/ttyUSB0")::NeuroAnalyzer.NEURO

    sp = _serial_open(port_name, baudrate=19200)
    @assert !isnothing(sp) "Serial port $port_name is not available"

    println("NeuroRecorder: TPT")
    println("==================")
    println("   Duration: $duration [seconds]")
    println("Serial port: $port_name")
    println()
    println("Ready to start, press SPACEBAR to begin recording")
    println()

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

    println("The recording will start after a beep")
    ts = time()
    while time() - ts <= 2
        _serial_listener(sp)
    end
    println()
    _beep()
    print("Recording .")

    # sampling rate is 50 Hz = 20 ms per loop
    fs = 50
    t = collect(0:1/fs:duration)
    tpt_ch_x = zeros(length(t))
    tpt_ch_y = zeros(length(t))
    tpt_ch_z = zeros(length(t))
    tpt_ch_accx = zeros(length(t))
    tpt_ch_accy = zeros(length(t))
    tpt_ch_accz = zeros(length(t))

    idx = 1
    ts = time()
    while idx <= length(tpt_ch_x)
        sp_signal = _serial_listener(sp)
        if !isnothing(sp_signal)
            m = match(r"(tpt\: )(\-*[0-9]+) (\-*[0-9]+) (\-*[0-9]+) (\-*[0-9]+\.[0-9]+) (\-*[0-9]+\.[0-9]+) (\-*[0-9]+\.[0-9]+)", sp_signal)
            if !isnothing(m)
                if length(m.captures) == 7
                    tpt_ch_x[idx] = parse(Float64, m.captures[2])
                    tpt_ch_y[idx] = parse(Float64, m.captures[3])
                    tpt_ch_z[idx] = parse(Float64, m.captures[4])
                    tpt_ch_accx[idx] = parse(Float64, m.captures[5])
                    tpt_ch_accy[idx] = parse(Float64, m.captures[6])
                    tpt_ch_accz[idx] = parse(Float64, m.captures[7])
                    idx += 1
                end
            end
        end
        if time() - ts >= 1.0
            print(".")
            ts = time()
        end
    end
    _serial_close(sp)
    _beep()
    println()
    println()
    println("Recording finished.")

    tpt_ch_x = tpt_ch_x[1:(end - 1)]
    tpt_ch_y = tpt_ch_y[1:(end - 1)]
    tpt_ch_z = tpt_ch_z[1:(end - 1)]
    tpt_ch_accx = tpt_ch_accx[1:(end - 1)]
    tpt_ch_accy = tpt_ch_accy[1:(end - 1)]
    tpt_ch_accz = tpt_ch_accz[1:(end - 1)]
    t = round.(t[1:(end - 1)], digits=3)
    tpt_signal = Matrix([tpt_ch_x tpt_ch_y tpt_ch_z tpt_ch_accx tpt_ch_accy tpt_ch_accz]')
    tpt_signal = reshape(tpt_signal, 6, :, 1)

    obj = create_object(data_type="tpt")
    add_channel!(obj, data=tpt_signal, label=["pos_x", "pos_y", "pos_z", "acc_x", "acc_y", "acc_z"], type=["orient", "orient", "orient", "accel", "accel", "accel"], unit=["", "", "", "m/s²", "m/s²", "m/s²"])
    create_time!(obj, fs=fs)

    return obj

end
