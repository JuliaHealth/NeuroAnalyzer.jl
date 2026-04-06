export iedar
export edar

"""
    iedar(; <keyword arguments>)

Record electrodermal activity (EDA) in GUI mode.

EDA (also known as Galvanic Skin Response / skin conductance) is recorded using a Grove GSR sensor via an Arduino connected over a virtual serial port.

Nominal sampling rate is 50 Hz.

# Arguments

- `duration::Int64=20`: recording duration in seconds
- `port_name::String="/dev/ttyUSB0"`: serial port to which the Arduino is connected

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object (EDA) with one channel at 50 Hz
"""
function iedar(;
    duration::Int64 = 20,
    port_name::String = "/dev/ttyUSB0",
)::NeuroAnalyzer.NEURO
    sp = _serial_open(port_name; baudrate = 19200)
    isnothing(sp) && throw(ArgumentError("Serial port $port_name is not available."))

    fs = 50   # Grove GSR Arduino sketch runs at 50 Hz
    t = collect(0:(1 / fs):duration)
    eda_signal = fill(NaN, length(t))

    p = Plots.plot(;
        ylims          = (0, 10),
        xlims          = (t[1], t[end]),
        legend         = false,
        palette        = :darktest,
        size           = (800, 400),
        margins        = 20Plots.px,
        xlabelfontsize = 8,
        ylabelfontsize = 8,
        xtickfontsize  = 8,
        ytickfontsize  = 8,
    )

    function _activate(app)
        win = GtkApplicationWindow(app, "NeuroRecorder: iedar()")
        Gtk4.default_size(win, p.attr[:size][1], p.attr[:size][2] + 40)

        can              = GtkCanvas()
        can.content_width  = p.attr[:size][1]
        can.content_height = p.attr[:size][2]

        g                    = GtkGrid()
        g.column_homogeneous = false
        g.column_spacing     = 5
        g.row_spacing        = 5
        g.margin_start       = 5
        g.margin_end         = 5
        g.margin_top         = 5
        g.margin_bottom      = 5

        bt_record              = GtkButton("RECORD")
        bt_record.tooltip_text = "Start recording"

        lb_status1        = GtkLabel("Status:")
        lb_status2        = GtkLabel("READY TO START")
        lb_status1.halign = 2
        lb_status2.halign = 1

        g[1:2, 1] = can
        g[1:2, 2] = bt_record
        g[1, 3]   = lb_status1
        g[2, 3]   = lb_status2
        vbox = GtkBox(:v)
        push!(vbox, g)
        push!(win, vbox)
        Gtk4.show(win)

        @guarded draw(can) do widget
            p_draw = Plots.plot(
                t, eda_signal;
                mc             = :black,
                ms             = 0.5,
                lw             = 0.5,
                lc             = :black,
                ylims          = (0, 10),
                xlims          = (t[1], t[end]),
                legend         = false,
                palette        = :darktest,
                size           = (800, 400),
                margins        = 20Plots.px,
                xlabelfontsize = 8,
                ylabelfontsize = 8,
                xtickfontsize  = 8,
                ytickfontsize  = 8,
            )
            io  = IOBuffer()
            ctx = getgc(can)
            withenv("GKSwstype" => "100") do
                png(p_draw, io)
            end
            seek(io, 0)
            img = read_from_png(io)
            set_source_surface(ctx, img, 0, 0)
            Cairo.paint(ctx)
        end

        return @guarded signal_connect(bt_record, "clicked") do widget
            bt_record.sensitive = false
            Threads.@spawn begin
                @idle_add lb_status2.label = "PREPARING"
                ts = time()
                while time() - ts <= 2
                    _serial_listener(sp)
                end
                _beep()
                @idle_add lb_status2.label = "RECORDING"
                t_refresh = time()
                idx = 1
                while idx <= length(eda_signal)
                    if time() - t_refresh >= 0.1
                        @idle_add draw(can)
                        t_refresh = time()
                    end
                    sp_signal = _serial_listener(sp)
                    if !isnothing(sp_signal)
                        m = match(r"gsr:([0-9]+\.[0-9]+)", sp_signal)
                        if !isnothing(m)
                            eda_signal[idx] = parse(Float64, m.captures[1])
                            idx += 1
                        end
                    end
                end
                @idle_add draw(can)
                _serial_close(sp)
                @idle_add lb_status2.label = "FINISHED"
                _beep()
                sleep(2)
                @idle_add close(win)
            end
        end
    end

    app = GtkApplication("org.neuroanalyzer.iedar")
    Gtk4.signal_connect(_activate, app, :activate)
    Gtk4.GLib.stop_main_loop()
    Gtk4.run(app)

    # trim the pre-allocated trailing sample
    eda_signal = reshape(eda_signal[1:(end - 1)], 1, :, 1)

    obj = create_object(; data_type = "eda")
    add_channel!(obj; data = eda_signal, label = ["eda1"], type = ["eda"], unit = ["µS"])
    create_time!(obj; fs = fs)

    return obj
end

"""
    edar(; <keyword arguments>)

Record electrodermal activity (EDA) in CLI mode.

EDA (also known as Galvanic Skin Response / skin conductance) is recorded using a Grove GSR sensor via an Arduino connected over a virtual serial port.

Nominal sampling rate is 50 Hz.

# Arguments

- `duration::Int64=20`: recording duration in seconds
- `port_name::String="/dev/ttyUSB0"`: serial port to which the Arduino is connected

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object (EDA) with one channel at 50 Hz
"""
function edar(;
    duration::Int64   = 20,
    port_name::String = "/dev/ttyUSB0",
)::NeuroAnalyzer.NEURO
    sp = _serial_open(port_name; baudrate = 19200)
    isnothing(sp) && throw(ArgumentError("Serial port $port_name is not available."))

    println("NeuroRecorder: EDA")
    println("==================")
    println("   Duration: $duration [seconds]")
    println("Serial port: $port_name")
    println()
    println("Ready to start, press SPACEBAR to begin recording")
    println()

    # wait for spacebar
    while true
        ret = ccall(:jl_tty_set_mode, Int32, (Ptr{Cvoid}, Int32), stdin.handle, true)
        ret == 0 || error("Unable to switch to raw mode.")
        kbd_key = read(stdin, Char)
        ccall(:jl_tty_set_mode, Int32, (Ptr{Cvoid}, Int32), stdin.handle, false)
        kbd_key == ' ' && break
    end

    println("The recording will start after a beep")
    ts = time()
    while time() - ts <= 2
        _serial_listener(sp)
    end
    println()
    _beep()
    print("Recording .")

    fs = 50
    t = collect(0:(1 / fs):duration)
    eda_signal = zeros(length(t))

    idx = 1
    ts  = time()
    while idx <= length(eda_signal)
        sp_signal = _serial_listener(sp)
        if !isnothing(sp_signal)
            m = match(r"gsr:([0-9]+\.[0-9]+)", sp_signal)
            if !isnothing(m)
                eda_signal[idx] = parse(Float64, m.captures[1])
                idx += 1
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

    # trim the pre-allocated trailing sample and build the NEURO object
    eda_signal = reshape(eda_signal[1:(end - 1)], 1, :, 1)

    obj = create_object(; data_type = "eda")
    add_channel!(obj; data = eda_signal, label = ["eda1"], type = ["eda"], unit = ["µS"])
    create_time!(obj; fs = fs)

    return obj
end