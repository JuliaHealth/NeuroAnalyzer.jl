export ieda
export eda

"""
    ieda(; duration, port_name)

Record electrodermal activity (EDA), also called Galvanic Skin Response (GSR) or skin conductance, in GUI mode. EDA is recorded using Groove GSR sensor via Arduino attached to the PC via USB cable (virtual serial port).

# Arguments

- `duration::Int64=20`: recording duration in seconds
- `port_name::String="/dev/ttyUSB0"`: serial port to which the Arduino is connected

# Returns

Named tuple containing:
- `eda::Vector{Int64}`: EDA signal, conductance in μS
- `time::Vector{Int64}`: time points
- `f::Int64`: sampling frequency
"""
function ieda(; duration::Int64=20, port_name::String="/dev/ttyUSB0")

    sp = _serial_open(port_name, baudrate=9600)
    @assert !isnothing(sp) @info "Serial port $port_name is not available"

    p = Plots.plot(yticks=(0, 5),
                   xticks=nothing,
                   legend=false,
                   palette=:darktest,
                   size=(800, 400),
                   margins=20Plots.px,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=8,
                   ytickfontsize=8)

    win = GtkWindow("NeuroRecorder: ieda()", p.attr[:size][1], p.attr[:size][2] + 40)
    set_gtk_property!(win, :border_width, 20)
    set_gtk_property!(win, :resizable, false)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    set_gtk_property!(win, :startup_id, "org.neuroanalyzer")
    can = GtkCanvas(Int32(p.attr[:size][1]), Int32(p.attr[:size][2]))
    g = GtkGrid()
    set_gtk_property!(g, :column_homogeneous, false)
    set_gtk_property!(g, :column_spacing, 10)
    set_gtk_property!(g, :row_spacing, 10)

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

    eda_signal = Vector{Float64}()

    @guarded draw(can) do widget
        p = Plots.plot(eda_signal,
                       ylims=(0, 5),
                       xticks=nothing,
                       legend=false,
                       palette=:darktest,
                       lc=:blue,
                       lw=0.5,
                       size=(800, 400),
                       margins=20Plots.px,
                       xlabelfontsize=8,
                       ylabelfontsize=8,
                       xtickfontsize=8,
                       ytickfontsize=8)
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
            t_1 = time()
            while time() - t_1 <= 2
                _serial_listener(sp)
            end
            _beep()
            set_gtk_property!(lb_status2, :label, "RECORDING")
            t_1 = time()
            t_refresh = time()
            while time() - t_1 <= duration
                sp_signal = _serial_listener(sp)
                if !isnothing(sp_signal)
                    m = match(r"(gsr\:)([0-9]+\.[0-9]+)", sp_signal)
                    if !isnothing(m)
                        if length(m.captures) == 2
                            push!(eda_signal, parse(Float64, m.captures[2]))
                        end
                    end
                end
                sleep(0.01)
                if time() - t_refresh >= 0.1
                    draw(can)
                    t_refresh = time()
                end
            end
            _beep()
            set_gtk_property!(lb_status2, :label, "FINISHED")
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

    _serial_close(sp)

    if length(eda_signal) > 0
        t = round.(linspace(0, duration, length(eda_signal)), digits=3)
        f = round(Int64, 1 / (t[2] - t[1]))
        eda_signal = reshape(eda_signal, 1, :, 1)
        obj = create_object(data_type="eda")
        add_channel!(obj, data=eda_signal, label=["eda1"], type=["eda"], unit=["µS"])
        create_ti1
        return obj
    else
        return nothing
    end

end

"""
    eda(; duration, port_name)

Record electrodermal activity (EDA), also called Galvanic Skin Response (GSR) or skin conductance, in CLI mode. EDA is recorded using Groove GSR sensor via Arduino attached to the PC via USB cable (virtual serial port).

# Arguments

- `duration::Int64=20`: recording duration in seconds
- `port_name::String="/dev/ttyUSB0"`: serial port to which the switch is connected

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function eda(; duration::Int64=20, port_name::String="/dev/ttyUSB0")

    sp = _serial_open(port_name, baudrate=9600)
    @assert !isnothing(sp) "Serial port $port_name is not available"

    println("NeuroRecorder: EDA")
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
    sleep(1)
    println()
    _beep()
    print("Recording ")

    eda_signal = Vector{Float64}()

    t_1 = time()
    counter = 0
    while time() - t_1 <= duration
        sp_signal = _serial_listener(sp)
        if !isnothing(sp_signal)
            m = match(r"(gsr\:)([0-9]+\.[0-9]+)", sp_signal)
            if !isnothing(m)
                if length(m.captures) == 2
                    push!(eda_signal, parse(Float64, m.captures[2]))
                end
            end
        end
        sleep(0.01)
        counter += 1
        if counter == 100
            print(".")
            counter = 0
        end
    end
    _serial_close(sp)
    _beep()
    println()
    println()
    println("Recording finished.")

    t = round.(linspace(0, duration, length(eda_signal)), digits=3)
    f = round(Int64, 1 / (t[2] - t[1]))

    eda_signal = reshape(eda_signal, 1, :, 1)
    obj = create_object(data_type="eda")
    add_channel!(obj, data=eda_signal, label=["eda1"], type=["eda"], unit=["µS"])
    create_time!(obj, fs=f)

    return obj

end
