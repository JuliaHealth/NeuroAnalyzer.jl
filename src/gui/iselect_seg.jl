export iselect_seg

"""
    iselect_seg(m; c, extract, v)

Interactive selection of matrix area.

# Arguments

- `m::AbstractMatrix`
- `c::Bool=false`: if true, select circular segment
- `extract::Bool=false`
- `v::Bool=false`

# Returns

- `r1::Int64`: upper-left corner
- `r2::Int64`: bottom-right corner
- `c1::Int64`: upper-left corner
- `c2::Int64`: bottom-right corner

or

- `seg::Union{AbstractMatrix, AbstractVector}`: extracted segment
"""
function iselect_seg(m::AbstractMatrix; c::Bool=false, extract::Bool=false, v::Bool=false)

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

    win = GtkWindow("NeuroAnalyzer: iselect_seg()", p.attr[:size][1] + 2, p.attr[:size][2] + 2)
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
        ctx = getgc(widget)
        if length(x) >= 0 && length(x) < 2
            push!(x, x_pos)
            push!(y, y_pos)
        end
        if length(x) == 1
            Gtk.arc(ctx, x[1], y[1], 2, 0, 2*pi)
            Gtk.set_source_rgb(ctx, 1, 0, 0)
            Gtk.fill(ctx)
            Gtk.stroke(ctx)
            Gtk.reveal(widget)
        else length(x) == 2
            if c == true
                Gtk.arc(ctx, x[1], y[1], distance((x[1], y[1]), (x[2], y[2])), 0, 2*pi)
            else
                Gtk.rectangle(ctx, x[1], y[1], x[2] - x[1], y[2] - y[1])
            end
            Gtk.set_source_rgb(ctx, 1, 0, 0)
            Gtk.set_line_width(ctx, 4.0);
            Gtk.stroke(ctx)
            Gtk.reveal(widget)
        end
    end

    can.mouse.button2press = @guarded (widget, event) -> begin
        show(io, MIME("image/png"), p)
        img = read_from_png(io)
        ctx = getgc(can)
        Cairo.set_source_surface(ctx, img, 1, 1)
        Cairo.paint(ctx)
        if length(x) > 0
            pop!(x)
            pop!(y)
        end
        if length(x) == 1
            Gtk.arc(ctx, x[1], y[1], 2, 0, 2*pi)
            Gtk.set_source_rgb(ctx, 1, 0, 0)
            Gtk.fill(ctx)
            Gtk.stroke(ctx)
        end
        Gtk.reveal(widget)
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
        c1 = div(x[1], size_x) .+ 1
        c2 = div(x[2], size_x) .+ 1
        r1 = div(y[1], size_y) .+ 1
        r2 = div(y[2], size_y) .+ 1
        if c == false
            r1 > r2 && ((r1, r2) = _swap(r1, r2))
            c1 > c2 && ((c1, c2) = _swap(c1, c2))
        end
        return extract == false ? (r1, c1, r2, c2) : seg_extract(m, (r1, c1, r2, c2), v=v, c=c)
    else
        return nothing
    end

end
