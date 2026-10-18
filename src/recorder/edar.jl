export iedar
export edar

# Grove GSR Arduino sketch: 50 Hz, lines "gsr:<value>"
const _EDA_FS = 50
const _EDA_BAUDRATE = 19200
const _EDA_REGEX = r"gsr:([0-9]+\.[0-9]+)"

# value of one EDA serial line, or nothing if the line is not an EDA sample
function _parse_eda_sample(line::String)::Union{Float64, Nothing}
    m = match(_EDA_REGEX, line)
    return isnothing(m) ? nothing : parse(Float64, m.captures[1])
end

# NEURO object with one EDA channel from the recorded samples
function _build_eda_object(eda_signal::Vector{Float64})::NeuroAnalyzer.NEURO
    isempty(eda_signal) && throw(ArgumentError("No EDA samples were recorded."))
    obj = create_object(; data_type = "eda")
    add_channel!(
        obj;
        data = reshape(eda_signal, 1, :, 1),
        label = ["eda1"],
        type = ["eda"],
        unit = ["μS"],
    )
    create_time!(obj; fs = _EDA_FS)
    return obj
end

"""
    iedar(; <keyword arguments>)

Record electrodermal activity (EDA) in GUI mode.

EDA (also known as Galvanic Skin Response / skin conductance) is recorded using a Grove GSR sensor via an Arduino connected over a virtual serial port.

Nominal sampling rate is 50 Hz. If the window is closed before the end of the recording, the samples recorded so far are returned.

# Arguments

- `duration::Int64=20`: recording duration in seconds; must be > 0
- `port_name::String="/dev/ttyUSB0"`: serial port to which the Arduino is connected

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object (EDA) with one channel at 50 Hz
"""
function iedar(;
    duration::Int64 = 20,
    port_name::String = "/dev/ttyUSB0",
)::NeuroAnalyzer.NEURO
    duration > 0 || throw(ArgumentError("duration must be > 0."))
    sp = _serial_open(port_name; baudrate = _EDA_BAUDRATE)

    n = duration * _EDA_FS
    t = collect(0:(n - 1)) ./ _EDA_FS
    eda_signal = fill(NaN, n)
    n_rec = 0             # number of recorded samples
    aborted = false       # set when the window is closed
    rec_task = nothing    # recording task (started by RECORD)

    plot_kw = (
        ylims = (0, 10),
        xlims = (t[1], t[end]),
        legend = false,
        palette = :darktest,
        size = (800, 400),
        margins = 20Plots.PlotMeasures.px,
        xlabelfontsize = 8,
        ylabelfontsize = 8,
        xtickfontsize = 8,
        ytickfontsize = 8,
    )

    function _activate(app)
        win = GtkApplicationWindow(app, "NeuroRecorder: iedar()")
        Gtk4.default_size(win, plot_kw.size[1], plot_kw.size[2] + 40)

        can = GtkCanvas()
        can.content_width = plot_kw.size[1]
        can.content_height = plot_kw.size[2]

        g = GtkGrid()
        g.column_homogeneous = false
        g.column_spacing = 5
        g.row_spacing = 5
        g.margin_start = 5
        g.margin_end = 5
        g.margin_top = 5
        g.margin_bottom = 5

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

        # closing the window stops the recording
        signal_connect(win, "close-request") do widget
            aborted = true
            return false
        end

        @guarded draw(can) do widget
            p_draw = Plots.plot(
                t, eda_signal;
                mc = :black,
                ms = 0.5,
                lw = 0.5,
                lc = :black,
                plot_kw...,
            )
            io = IOBuffer()
            ctx = getgc(can)
            withenv("GKSwstype" => "100") do
                return png(p_draw, io)
            end
            seek(io, 0)
            img = read_from_png(io)
            set_source_surface(ctx, img, 0, 0)
            return Cairo.paint(ctx)
        end

        return @guarded signal_connect(bt_record, "clicked") do widget
            bt_record.sensitive = false
            rec_task = Threads.@spawn begin
                try
                    @idle_add lb_status2.label = "PREPARING"
                    _serial_flush(sp, 2)
                    _beep()
                    @idle_add lb_status2.label = "RECORDING"
                    t_refresh = time()
                    while n_rec < n && !aborted
                        if time() - t_refresh >= 0.1
                            @idle_add draw(can)
                            t_refresh = time()
                        end
                        line = _serial_listener(sp)
                        isnothing(line) && (sleep(0.001); continue)
                        v = _parse_eda_sample(line)
                        if !isnothing(v)
                            n_rec += 1
                            eda_signal[n_rec] = v
                        end
                    end
                    aborted && return nothing
                    @idle_add draw(can)
                    @idle_add lb_status2.label = "FINISHED"
                    _beep()
                    sleep(2)
                    @idle_add close(win)
                catch err
                    @error "EDA recording failed" exception = (err, catch_backtrace())
                finally
                    _serial_close(sp)
                end
            end
        end
    end

    app = GtkApplication("org.neuroanalyzer.iedar")
    Gtk4.signal_connect(_activate, app, :activate)
    Gtk4.GLib.stop_main_loop()
    Gtk4.run(app)

    # the window is closed: stop the recording task (if any) and release the port
    aborted = true
    isnothing(rec_task) || wait(rec_task)
    _serial_close(sp)

    return _build_eda_object(eda_signal[1:n_rec])
end

"""
    edar(; <keyword arguments>)

Record electrodermal activity (EDA) in CLI mode.

EDA (also known as Galvanic Skin Response / skin conductance) is recorded using a Grove GSR sensor via an Arduino connected over a virtual serial port.

Nominal sampling rate is 50 Hz.

# Arguments

- `duration::Int64=20`: recording duration in seconds; must be > 0
- `port_name::String="/dev/ttyUSB0"`: serial port to which the Arduino is connected

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object (EDA) with one channel at 50 Hz
"""
function edar(;
    duration::Int64 = 20,
    port_name::String = "/dev/ttyUSB0",
)::NeuroAnalyzer.NEURO
    duration > 0 || throw(ArgumentError("duration must be > 0."))
    sp = _serial_open(port_name; baudrate = _EDA_BAUDRATE)

    n = duration * _EDA_FS
    eda_signal = zeros(n)

    try
        println("NeuroRecorder: EDA")
        println("==================")
        println("   Duration: $duration [seconds]")
        println("Serial port: $port_name")
        println()
        println("Ready to start, press SPACEBAR to begin recording")
        println()

        _wait_spacebar()

        println("The recording will start after a beep")
        # discard data buffered while waiting for the user
        _serial_flush(sp, 2)
        println()
        _beep()
        print("Recording .")

        idx = 1
        ts = time()
        while idx <= n
            line = _serial_listener(sp)
            if isnothing(line)
                sleep(0.001)
            else
                v = _parse_eda_sample(line)
                if !isnothing(v)
                    eda_signal[idx] = v
                    idx += 1
                end
            end
            if time() - ts >= 1.0
                print(".")
                ts = time()
            end
        end
    finally
        _serial_close(sp)
    end
    _beep()
    println()
    println()
    println("Recording finished.")

    return _build_eda_object(eda_signal)
end
