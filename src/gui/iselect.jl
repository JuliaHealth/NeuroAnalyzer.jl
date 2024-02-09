export iselect

"""
    iselect(m)

Interactive selection of matrix area.

# Arguments

- `m::AbstractMatrix`

# Returns

- `r::Vector{Int64}`: list of row indices
- `c::Vector{Int64}`: list of column indices
"""
function iselect(m::AbstractMatrix)

    p = heatmap(m,
                framestyle=:none,
                cb=false;
                fill=:darktest,
                margins=-100Plots.px,
                legend=false)
    
    dim_x = size(m, 2)
    dim_y = size(m, 1)
    size_x = p.attr[:size][1] ÷ dim_x
    size_y = p.attr[:size][2] ÷ dim_y

    win = GtkWindow("NeuroAnalyzer: iselect()", p.attr[:size][1] + 2, p.attr[:size][2] + 2)
    set_gtk_property!(win, :border_width, 0)
    set_gtk_property!(win, :resizable, false)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    set_gtk_property!(win, :startup_id, "org.neuroanalyzer")
    can = GtkCanvas(p.attr[:size][1] + 2, p.attr[:size][2] + 2)
    push!(win, can)
    showall(win)

    x = Int64[]
    y = Int64[]

    @guarded draw(can) do widget
        show(io, MIME("image/png"), p)
        img = read_from_png(io)
        ctx = getgc(can)
        Cairo.set_source_surface(ctx, img, 1, 1)
        Cairo.paint(ctx)
    end

    can.mouse.button1press = @guarded (widget, event) -> begin
        x_pos = round(event.x)
        y_pos = round(event.y)
        push!(x, x_pos)
        push!(y, y_pos)
        ctx = getgc(widget)
        if length(x) == 1
            Gtk.arc(ctx, x[end], y[end], 2.5, 0, 2*pi)
            Gtk.set_source_rgb(ctx, 1, 0, 0)
            Gtk.fill(ctx)
        else
            Gtk.set_line_cap(ctx, Cairo.CAIRO_LINE_CAP_ROUND)
            Gtk.move_to(ctx, x[end - 1], y[end - 1])
            Gtk.line_to(ctx, x[end], y[end])
            Gtk.set_source_rgb(ctx, 1, 0, 0)
            Gtk.set_line_width(ctx, 5.0);
        end
        Gtk.stroke(ctx)
        Gtk.reveal(widget)
    end

    can.mouse.button3press = @guarded (widget, event) -> begin
        if length(x) > 2
            push!(x, x[1])
            push!(y, y[1])
            ctx = getgc(widget)
            Gtk.set_line_cap(ctx, Cairo.CAIRO_LINE_CAP_ROUND)
            Gtk.move_to(ctx, x[end - 1], y[end - 1])
            Gtk.line_to(ctx, x[end], y[end])
            Gtk.set_source_rgb(ctx, 1, 0, 0)
            Gtk.set_line_width(ctx, 5.0);
            Gtk.stroke(ctx)
            Gtk.reveal(widget)
        end
    end

    can.mouse.button2press = @guarded (widget, event) -> begin
        if length(x) > 0
            show(io, MIME("image/png"), p)
            img = read_from_png(io)
            ctx = getgc(can)
            Cairo.set_source_surface(ctx, img, 1, 1)
            Cairo.paint(ctx)
            pop!(x)
            pop!(y)
            if length(x) > 0
                Gtk.arc(ctx, x[1], y[1], 2.5, 0, 2*pi)
                Gtk.set_source_rgb(ctx, 1, 0, 0)
                Gtk.fill(ctx)
            end
            if length(x) > 1
                Gtk.set_line_cap(ctx, Cairo.CAIRO_LINE_CAP_ROUND)
                Gtk.set_source_rgb(ctx, 1, 0, 0)
                Gtk.set_line_width(ctx, 5.0);
                for idx in 2:length(x)
                    Gtk.move_to(ctx, x[idx - 1], y[idx - 1])
                    Gtk.line_to(ctx, x[idx], y[idx])
                end
            end
            Gtk.stroke(ctx)
            Gtk.reveal(widget)
        end
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        s = event.state
        if s == 20
            if k == 115 # s
                file_name = save_dialog("Pick image file", GtkNullContainer(), (GtkFileFilter("*.png", name="All supported formats"), "*.png"))
                    if file_name != ""
                        if splitext(file_name)[2] in [".png"]
                            surface_buf = Gtk.cairo_surface(can)
                            Cairo.write_to_png(surface_buf, file_name)
                        else
                            warn_dialog("Incorrect file name!")
                        end
                    end
            elseif k == 113 # q
                Gtk.destroy(win)
                x = nothing
                y = nothing
            end
        end
        if k == 65293
            Gtk.destroy(win)
        end
    end

    cnd = Condition()
    signal_connect(win, :destroy) do widget
        notify(cnd)
    end
    @async Gtk.gtk_main()
    wait(cnd)

    if x !== nothing && y !== nothing && length(x) > 0 && length(y) > 0 
        if length(x) > 1 && length(y) > 1 && x[end] == x[1] && y[end] == y[1]
            pop!(x)
            pop!(y)
        end
        r = div.(x, size_x) .+ 1
        c = div.(y, size_y) .+ 1
        return r, c
    else
        return Int64[], Int64[]
    end

end
