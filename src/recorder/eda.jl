export ieda
export eda

"""
    ieda(; port_name)

Record electrodermal activity (EDA), also called Galvanic Skin Response (GSR) or skin conductance, in GUI mode. EDA is recorded using Groove GSR sensor via Arduino attached to the PC via USB cable (virtual serial port).

# Arguments

- `port_name::String="/dev/ttyUSB0"`: serial port to which the Arduino is connected

# Returns

Named tuple containing:
- `eda::Vector{Int64}`: EDA signal, conductance in μS
- `time::Vector{Int64}`: time points
- `f::Int64`: sampling frequency
"""
function ieda(; port_name::String="/dev/ttyUSB0")

    sp = nothing
    sp = _serial_open(port_name)
    @assert !isnothing(sp) @info "Serial port $port_name is not available"

    t = collect(0:1:10)
    y = nothing

    p = Plots.plot(ylims=(-10, 10),
                   xticks=(t[1]:t[end]),
                   legend=false,
                   palette=:darktest,
                   size=(800, 400),
                   margins=20Plots.px,
                   titlefontsize=10,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=8,
                   ytickfontsize=8)

    win = GtkWindow("NeuroRecorder: ieda()", 800, 400)
    set_gtk_property!(win, :border_width, 20)
    set_gtk_property!(win, :resizable, false)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    set_gtk_property!(win, :startup_id, "org.neuroanalyzer")

    can = GtkCanvas(800, 600)

    g1 = GtkGrid()
    set_gtk_property!(g1, :column_homogeneous, true)
    set_gtk_property!(g1, :column_spacing, 20)
    set_gtk_property!(g1, :row_spacing, 20)
    g2 = GtkGrid()
    set_gtk_property!(g2, :column_homogeneous, true)
    set_gtk_property!(g2, :column_spacing, 20)
    set_gtk_property!(g2, :row_spacing, 20)

    bt_start = GtkButton("RECORD")
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
        show(io, MIME("image/png"), Plots.plot(obj,
                                                  ch=ch[ch_first]:ch[ch_last],
                                                  seg=(time1, time2),
                                                  mono=true,
                                                  title=""))
        img = read_from_png(io)
        set_source_surface(ctx, img, 0, 0)
        Cairo.paint(ctx)
    end

    eda_signal = Vector{Float64}()
    time_pts = Vector{Float64}()
    f = nothing

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        if k == 32 && get_gtk_property(lb_status2, :label, String) == "TEST"
            #t = parse(Int64, split(get_gtk_property(lb_trial2, :label, String), ' ')[1])
            if !key_pressed
                push!(t_kp_tmp, time())
                key_pressed = true
            end
        elseif k == 32 && get_gtk_property(lb_status2, :label, String) == "INTERVAL"
            #t = parse(Int64, split(get_gtk_property(lb_interval2, :label, String), ' ')[1])
            if !key_pressed
                push!(int_t_kp_tmp, time())
                key_pressed = true
            end
        end
        sleep(0.1)
    end

    signal_connect(win, "key-release-event") do widget, event
        k = event.keyval
        if k == 32 && get_gtk_property(lb_status2, :label, String) == "TEST"
            #t = parse(Int64, split(get_gtk_property(lb_trial2, :label, String), ' ')[1])
            if key_pressed
                push!(d_kp_tmp, time())
                key_pressed = false
            end
        elseif k == 32 && get_gtk_property(lb_status2, :label, String) == "INTERVAL"
            #t = parse(Int64, split(get_gtk_property(lb_interval2, :label, String), ' ')[1])
            if key_pressed
                push!(int_d_kp_tmp, time())
                key_pressed = false
            end
        end
        sleep(0.1)
    end

    signal_connect(bt_start, "clicked") do widget
        if port_name == ""
            set_gtk_property!(bt_start, :sensitive, false)
            Threads.@spawn begin
                for idx in 1:trials
                    t_kp_tmp = Vector{Float64}()
                    d_kp_tmp = Vector{Float64}()
                    int_t_kp_tmp = Vector{Float64}()
                    int_d_kp_tmp = Vector{Float64}()
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
                    push!(t_idx, t1)
                    while time() <= t1 + duration
                    end
                    _beep()
                    push!(t_kp, t_kp_tmp)
                    push!(d_kp, d_kp_tmp)
                    result[idx] = length(t_kp_tmp)
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
                    push!(int_idx, t1)
                    while time() <= t1 + interval
                    end
                    push!(int_t_kp, int_t_kp_tmp)
                    push!(int_d_kp, int_d_kp_tmp)
                    int_result[idx] = length(int_t_kp_tmp)
                end

                # Interacting with GTK from a thread other than the main thread is
                # generally not allowed, so we register an idle callback instead.
                Gtk.GLib.g_idle_add(nothing) do user_data
                    Gtk.destroy(win)
                end
            end
        elseif !isnothing(sp)
            # use serial port
            set_gtk_property!(bt_start, :sensitive, false)
            t_s = nothing
            t_kp = Vector{Float64}()
            d_kp = Vector{Float64}()
            int_t_kp = Vector{Float64}()
            int_d_kp = Vector{Float64}()
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
                    key_pressed = false
                    sp = _serial_open(port_name)
                    t_1 = time()
                    while time() - t_1 <= duration
                        t = time() - t_1
                        serial_key = _serial_listener(sp)
                        if serial_key == "$gpio:1"
                            if !key_pressed
                                # key is pressed
                                push!(t_kp, t)
                                result[idx] += 1
                                key_pressed = true
                            end
                        elseif serial_key == "$gpio:0"
                            if key_pressed
                                # key is released
                                push!(d_kp, t)
                                key_pressed = false
                            end
                        end
                        sleep(0.1)
                    end
                    _serial_close(sp)
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
                    key_pressed = false
                    sp = _serial_open(port_name)
                    t_2 = time()
                    while time() - t_2 <= interval
                        t = time() - t_2
                        serial_key = _serial_listener(sp)
                        if serial_key == "$gpio:1"
                            if !key_pressed
                                # key is pressed
                                push!(int_t_kp, t)
                                int_result[idx] += 1
                                key_pressed = true
                            end
                        elseif serial_key == "$gpio:0"
                            if key_pressed
                                # key is released
                                push!(int_d_kp, t)
                                key_pressed = false
                            end
                        end
                        sleep(0.1)
                    end
                    _serial_close(sp)
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
    end

    cnd = Condition()
    signal_connect(win, :destroy) do widget
        notify(cnd)
    end
    @async Gtk.gtk_main()
    wait(cnd)

    if port_name == ""
        for idx in 1:trials
            if length(t_kp[idx]) != length(d_kp[idx])
                l = minimum([length(t_kp[idx]), length(d_kp[idx])])
                t_kp[idx] = t_kp[idx][1:l]
                d_kp[idx] = d_kp[idx][1:l]
                result[idx] = length(t_kp[idx])
            end
            if length(int_t_kp[idx]) != length(int_d_kp[idx])
                l = minimum([length(int_t_kp[idx]), length(int_d_kp[idx])])
                int_t_kp[idx] = int_t_kp[idx][1:l]
                int_d_kp[idx] = int_d_kp[idx][1:l]
                int_result[idx] = length(int_t_kp[idx])
            end
            d_kp[idx] = round.((d_kp[idx] .- t_kp[idx]) .* 1000, digits=1)
            t_kp[idx] = round.((t_kp[idx] .- t_idx[idx]) .* 1000, digits=1)
            int_d_kp[idx] = round.((int_d_kp[idx] .- int_t_kp[idx]) .* 1000, digits=1)
            int_t_kp[idx] = round.((int_t_kp[idx] .- int_idx[idx]) .* 1000, digits=1)
        end

        # remove out of time boundary taps
        for idx1 in 1:trials
            for idx2 in length(t_kp[idx1]):-1:1
                if t_kp[idx1][idx2] > duration * 1000
                    deleteat!(t_kp[idx1], idx2)
                    deleteat!(d_kp[idx1], idx2)
                    result[idx1] -= 1
                end
            end
        end
        for idx1 in 1:trials
            for idx2 in length(int_t_kp[idx1]):-1:1
                if int_t_kp[idx1][idx2] > interval * 1000
                    deleteat!(int_t_kp[idx1], idx2)
                    deleteat!(int_d_kp[idx1], idx2)
                    int_result[idx1] -= 1
                end
            end
        end

        return (taps=result, tap_t=t_kp, tap_d=d_kp, taps_int=int_result, tap_t_int=int_t_kp, tap_d_int=int_d_kp)
    elseif !isnothing(sp)
        # format time points
        d_kp = round.((d_kp .- t_kp) .* 1000, digits=1)
        int_d_kp = round.((int_d_kp .- int_t_kp) .* 1000, digits=1)
        t_kp = round.(t_kp .* 1000, digits=1)
        int_t_kp = round.(int_t_kp .* 1000, digits=1)

        # format time points
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
            push!(t_keypressed, round.(tk, digits=1))
            push!(d_keypressed, round.(td, digits=1))
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
            push!(int_t_keypressed, round.(tk, digits=1))
            push!(int_d_keypressed, round.(td, digits=1))
        end
        reverse!(int_t_keypressed)
        reverse!(int_d_keypressed)

        # remove out of time boundary taps
        for idx1 in 1:trials
            for idx2 in length(t_keypressed[idx1]):-1:1
                if t_keypressed[idx1][idx2] > (duration * idx1 * 1000)
                    deleteat!(t_keypressed[idx1], idx2)
                    deleteat!(d_keypressed[idx1], idx2)
                    result[idx1] -= 1
                end
            end
        end
        for idx1 in 1:trials
            for idx2 in length(int_t_keypressed[idx1]):-1:1
                if int_t_keypressed[idx1][idx2] > (interval * idx1 * 1000)
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
end

"""
    eda(; duration, port_name)

Record electrodermal activity (EDA), also called Galvanic Skin Response (GSR) or skin conductance, in CLI mode. EDA is recorded using Groove GSR sensor via Arduino attached to the PC via USB cable (virtual serial port).

# Arguments

- `duration::Int64=5`: recording duration in seconds
- `port_name::String="/dev/ttyUSB0"`: serial port to which the switch is connected

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function eda(; duration::Int64=5, port_name::String="/dev/ttyUSB0")

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
