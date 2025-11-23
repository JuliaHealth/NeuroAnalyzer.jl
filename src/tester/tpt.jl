export itpt
export tpt

"""
    itpt(; <keyword arguments>)

Perform Two-point Pinch Test (TPT) in GUI mode. TPT is recorded using MMA7660 accelerometer via Arduino attached to the PC via USB cable (virtual serial port). Sampling rate is 50 Hz.

# Arguments

- `duration::Int64=20`: recording duration in seconds
- `port_name::String="/dev/ttyUSB0"`: serial port to which the Arduino is connected

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function itpt(; duration::Int64=20, port_name::String="/dev/ttyUSB0")::NeuroAnalyzer.NEURO

    sp = _serial_open(port_name, baudrate=19200)
    @assert !isnothing(sp) _info("Serial port $port_name is not available")

    img1 = read_from_png(joinpath(res_path, "finger_nopinch.png"))
    img2 = read_from_png(joinpath(res_path, "finger_pinch.png"))

    # sampling rate is 50 Hz = 20 ms per loop
    fs = 50
    t = collect(0:1/fs:duration)
    tpt_ch_x = zeros(length(t))
    tpt_ch_y = zeros(length(t))
    tpt_ch_z = zeros(length(t))
    tpt_ch_accx = zeros(length(t))
    tpt_ch_accy = zeros(length(t))
    tpt_ch_accz = zeros(length(t))

    function _activate(app)

        win = GtkApplicationWindow(app, "NeuroAnalyzer: itpt()")
        win.width_request = Int64(img1.width)
        win.height_request = Int64(img1.height) + 100

        can = GtkCanvas(Int64(img1.width), Int64(img1.height))

        g = GtkGrid()
        g.column_homogeneous = false
        g.column_spacing = 5
        g.row_spacing = 5

        bt_record = GtkButton("RECORD")
        bt_record.tooltip_text = "Start recording"

        lb_status1 = GtkLabel("Status:")
        lb_status2 = GtkLabel("READY TO START")
        lb_status1.halign = 2
        lb_status2.halign = 1

        g[1:2, 1] = can
        g[1:2, 2] = bt_record
        g[1, 3] = lb_status1
        g[2, 3] = lb_status2
        vbox = GtkBox(:v)
        push!(vbox, g)

        push!(win, vbox)

        Gtk4.show(win)

        @guarded draw(can) do widget
            ctx = getgc(can)
            Cairo.set_source_surface(ctx, img1, 0, 0)
            Cairo.paint(ctx)
        end

        @guarded signal_connect(bt_record, "clicked") do widget
            bt_record.sensitive = false
            Threads.@spawn begin
                @idle_add lb_status2.label = "PREPARING"
                ts = time()
                while time() - ts <= 2
                    _serial_listener(sp)
                end
                _beep()
                @idle_add @guarded draw(can) do widget
                    ctx = getgc(can)
                    Cairo.set_source_surface(ctx, img2, 0, 0)
                    Cairo.paint(ctx)
                end
                @idle_add lb_status2.label = "RECORDING"
                idx = 1
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
                end
                _serial_close(sp)
                _beep()
                @idle_add lb_status2.label = "FINISHED"
                @idle_add @guarded draw(can) do widget
                    ctx = getgc(can)
                    Cairo.set_source_surface(ctx, img1, 0, 0)
                    Cairo.paint(ctx)
                end
                sleep(2)
                @idle_add close(win)
            end
        end
    end

    app = GtkApplication("org.neuroanalyzer.itpt")
    Gtk4.signal_connect(_activate, app, :activate)
    Gtk4.GLib.stop_main_loop()
    Gtk4.run(app)

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

    println("NeuroTester: TPT")
    println("================")
    println("   Duration: $duration [seconds]")
    println("Serial port: $port_name")
    println()
    println("Ready to start, press SPACEBAR to begin the test")
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

    println("The test will start after a beep")
    ts = time()
    while time() - ts <= 2
        _serial_listener(sp)
    end
    println()
    _beep()
    print("   Pinch the thumb and the index finger as quickly as possible")

    # sampling rate is 50 Hz = 20 ms per loop
    fs = 50
    t = collect(0:1/fs:duration)
    t = t[1:(end - 1)]
    tpt_ch_x = zeros(length(t))
    tpt_ch_y = zeros(length(t))
    tpt_ch_z = zeros(length(t))
    tpt_ch_accx = zeros(length(t))
    tpt_ch_accy = zeros(length(t))
    tpt_ch_accz = zeros(length(t))

    idx = 1
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
    end
    _serial_close(sp)
    _beep()
    println()
    println()
    println("Testing completed")
    tpt_signal = Matrix([tpt_ch_x tpt_ch_y tpt_ch_z tpt_ch_accx tpt_ch_accy tpt_ch_accz]')
    tpt_signal = reshape(tpt_signal, 6, :, 1)

    obj = create_object(data_type="tpt")
    add_channel!(obj, data=tpt_signal, label=["pos_x", "pos_y", "pos_z", "acc_x", "acc_y", "acc_z"], type=["orient", "orient", "orient", "accel", "accel", "accel"], unit=["", "", "", "m/s²", "m/s²", "m/s²"])
    create_time!(obj, fs=fs)

    return obj

end
