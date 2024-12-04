export itpt
export tpt

"""
    itpt(; <keyword arguments>)

Record Two-point Pinch Test (TPT) in GUI mode. TPT is recorded using MMA7660 accelerometer via Arduino attached to the PC via USB cable (virtual serial port). Sampling rate is 50 Hz.

# Arguments

- `duration::Int64=20`: single trial duration in seconds
- `trials::Int64=2`: number of trials
- `interval::Int64=2`: interval between trials in seconds
- `port_name::String="/dev/ttyUSB0"`: serial port to which the Arduino is connected

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
- `Nothing` if recording fails
"""
function itpt(; duration::Int64=20, trials::Int64=2, interval::Int64=2, port_name::String="/dev/ttyUSB0")::Union{NeuroAnalyzer.NEURO, Nothing}

    sp = _serial_open(port_name, baudrate=19200)
    @assert !isnothing(sp) @info "Serial port $port_name is not available"

    # sampling rate is 50 Hz = 20 ms per loop
    fs = 50
    t = collect(0:1/fs:duration)
    tpt_ch_x = repeat([NaN], length(t))
    tpt_ch_y = repeat([NaN], length(t))
    tpt_ch_z = repeat([NaN], length(t))
    tpt_ch_accx = repeat([NaN], length(t))
    tpt_ch_accy = repeat([NaN], length(t))
    tpt_ch_accz = repeat([NaN], length(t))

    p = Plots.plot(ylims=(0, 10),
                   xlims=(t[1], t[end]),
                   legend=false,
                   palette=:darktest,
                   size=(800, 400),
                   margins=20Plots.px,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=8,
                   ytickfontsize=8)

    win = GtkWindow("NeuroRecorder: itpt()", p.attr[:size][1], p.attr[:size][2] + 40)
    set_gtk_property!(win, :border_width, 20)
    set_gtk_property!(win, :resizable, false)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    set_gtk_property!(win, :startup_id, "org.neuroanalyzer")
    can = GtkCanvas(Int32(p.attr[:size][1]), Int32(p.attr[:size][2]))
    g = GtkGrid()
    set_gtk_property!(g, :column_homogeneous, false)
    set_gtk_property!(g, :column_spacing, 5)
    set_gtk_property!(g, :row_spacing, 5)

    bt_record = GtkButton("RECORD")
    set_gtk_property!(bt_record, :tooltip_text, "Start recording")

    lb_status1 = GtkLabel("Status:")
    lb_status2 = GtkLabel("READY TO START")
    set_gtk_property!(lb_status1, :halign, 2)
    set_gtk_property!(lb_status2, :halign, 1)

    g[1:2, 1] = can
    g[1:2, 2] = bt_record
    g[1, 3] = lb_status1
    g[2, 3] = lb_status2
    vbox = GtkBox(:v)
    push!(vbox, g)

    push!(win, vbox)
    showall(win)

    @guarded draw(can) do widget
        p = Plots.plot(t,
                       tpt_ch_x,
                       ms=0.5,
                       lw=0.5,
                       ylims=(0, 20),
                       xlims=(t[1], t[end]),
                       legend=false,
                       palette=:darktest,
                       size=(800, 400),
                       margins=20Plots.px,
                       xlabelfontsize=8,
                       ylabelfontsize=8,
                       xtickfontsize=8,
                       ytickfontsize=8)
        p = Plots.plot!(t,
                        tpt_ch_y)
        p = Plots.plot!(t,
                        tpt_ch_z)
        p = Plots.plot!(t,
                        tpt_ch_accx)
        p = Plots.plot!(t,
                        tpt_ch_accy)
        p = Plots.plot!(t,
                        tpt_ch_accz)
        ctx = getgc(can)
        show(io, MIME("image/png"), p)
        img = read_from_png(io)
        set_source_surface(ctx, img, 0, 0)
        Cairo.paint(ctx)
    end

    @guarded signal_connect(bt_record, "clicked") do widget
        set_gtk_property!(bt_record, :sensitive, false)
        Threads.@spawn begin
            set_gtk_property!(lb_status2, :label, "PREPARING")
            ts = time()
            while time() - ts <= 2
                _serial_listener(sp)
            end
            _beep()
            set_gtk_property!(lb_status2, :label, "RECORDING")
            t_refresh = time()
            idx = 1
            while idx <= length(tpt_ch_x)
                if time() - t_refresh >= 0.1
                    draw(can)
                    t_refresh = time()
                end
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
            draw(can)
            _serial_close(sp)
            set_gtk_property!(lb_status2, :label, "FINISHED")
            _beep()
            sleep(2)

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

    if length(tpt_ch_x) > 0
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
    else
        return nothing
    end

end

"""
    tpt(; <keyword arguments>)

Record electrodermal activity (EDA), also called Galvanic Skin Response (GSR) or skin conductance, in CLI mode. EDA is recorded using Groove GSR sensor via Arduino attached to the PC via USB cable (virtual serial port).

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
    eda_signal = zeros(length(t))

    idx = 1
    ts = time()
    while idx <= length(eda_signal)
        sp_signal = _serial_listener(sp)
        if !isnothing(sp_signal)
            m = match(r"(gsr\:)([0-9]+\.[0-9]+)", sp_signal)
            if !isnothing(m)
                if length(m.captures) == 2
                    eda_signal[idx] = parse(Float64, m.captures[2])
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

    if length(eda_signal) > 0
        eda_signal = eda_signal[1:(end - 1)]
        t = round.(t[1:(end - 1)], digits=3)
        eda_signal = reshape(eda_signal, 1, :, 1)
        obj = create_object(data_type="eda")
        add_channel!(obj, data=eda_signal, label=["eda1"], type=["eda"], unit=["µS"])
        create_time!(obj, fs=fs)
        return obj
    else
        return nothing
    end

end
