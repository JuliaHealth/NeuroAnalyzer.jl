export iselect_seg

"""
    iselect_seg(m; <keyword arguments>)

Interactive selection of matrix area.

# Arguments

- `m::AbstractMatrix`
- `shape::Symbol=:r`: selection shape:
    - `:r`: rectangular
    - `:c`: circular
    - `:p`: point
- `extract::Bool=false`: if true, return values of the matrix
- `v::Bool=false`: if true, return as vector (matrix m by rows over columns)

# Returns

- `r1::Int64`: upper-left corner
- `r2::Int64`: bottom-right corner
- `c1::Int64`: upper-left corner
- `c2::Int64`: bottom-right corner

or

- `seg::Union{Nothing, <:Real, Tuple{Int64, Int64}, Tuple{Int64, Int64, Int64, Int64}, Union{AbstractMatrix, AbstractVector, Tuple{AbstractVector, AbstractVector}}}`: extracted segment or its coordinates
"""
function iselect_seg(m::AbstractMatrix; shape::Symbol=:r, extract::Bool=false, v::Bool=false)::Union{Nothing, <:Real, Tuple{Int64, Int64}, Tuple{Int64, Int64, Int64, Int64}, Union{AbstractMatrix, AbstractVector, Tuple{AbstractVector, AbstractVector}}}

    _check_var(shape, [:r, :c, :p], "shape")

    size_x = size(m, 2)
    size_y = size(m, 1)

    p = Plots.heatmap(m,
                      framestyle=:none,
                      cb=false;
                      fill=:darktest,
                      margins=-100Plots.px,
                      legend=false,
                      size=(size_x, size_y))

    x_pos = Int64[]
    y_pos = Int64[]

    function _activate(app)

        win = GtkApplicationWindow(app, "NeuroAnalyzer: iview()")
        Gtk4.default_size(win, p.attr[:size][1] + 4, p.attr[:size][2] + 4)

        can = GtkCanvas()
        can.content_width = p.attr[:size][1]
        can.content_height = p.attr[:size][2]
        can.margin_start = 2
        can.margin_end = 2
        can.margin_top = 2
        can.margin_bottom = 2
        push!(win, can)

        Gtk4.show(win)

        @guarded draw(can) do widget
            show(io, MIME("image/png"), p)
            img = read_from_png(io)
            ctx = getgc(can)
            Cairo.set_source_surface(ctx, img, 1, 1)
            Cairo.paint(ctx)
        end

        function _lmb_click(_, _, x, y)
            x = round(x)
            y = round(y)
            ctx = getgc(can)
            if length(x) >= 0 && length(x_pos) < 2 && shape in [:r, :c]
                push!(x_pos, x)
                push!(y_pos, y)
            elseif shape === :p && length(x_pos) == 0
                push!(x_pos, x)
                push!(y_pos, y)
            end
            if length(x_pos) == 1
                Gtk4.arc(ctx, x_pos[1], y_pos[1], 2, 0, 2*pi)
                Gtk4.set_source_rgb(ctx, 0, 0, 0)
                Gtk4.fill(ctx)
                Gtk4.stroke(ctx)
                Gtk4.reveal(can)
            else length(x) == 2
                if shape === :c
                    Gtk4.arc(ctx, x_pos[1], y_pos[1], distance((x_pos[1], y_pos[1]), (x_pos[2], y_pos[2])), 0, 2*pi)
                elseif shape === :r
                    Gtk4.rectangle(ctx, x_pos[1], y_pos[1], x_pos[2] - x_pos[1], y_pos[2] - y_pos[1])
                end
                Gtk4.set_source_rgb(ctx, 0, 0, 0)
                Gtk4.set_line_width(ctx, 4.0);
                Gtk4.stroke(ctx)
                Gtk4.reveal(can)
            end
        end
        ggc_l = GtkGestureClick()
        ggc_l.button = 1
        push!(can, ggc_l)
        signal_connect(_lmb_click, ggc_l, "pressed")

        function _rmb_click(_, _, x, y)
            show(io, MIME("image/png"), p)
            img = read_from_png(io)
            ctx = getgc(can)
            Cairo.set_source_surface(ctx, img, 1, 1)
            Cairo.paint(ctx)
            if length(x_pos) > 0
                pop!(x_pos)
                pop!(y_pos)
            end
            if length(x_pos) == 1
                Gtk4.arc(ctx, x_pos[1], y_pos[1], 2, 0, 2*pi)
                Gtk4.set_source_rgb(ctx, 0, 0, 0)
                Gtk4.fill(ctx)
                Gtk4.stroke(ctx)
            end
            Gtk4.reveal(can)
        end
        ggc_r = GtkGestureClick()
        ggc_r.button = 3
        push!(can, ggc_r)
        signal_connect(_rmb_click, ggc_r, "pressed")

        win_key = Gtk4.GtkEventControllerKey(win)

        help = "Keyboard shortcuts:\n\nCtrl + Enter\t\tReturn selected segment\nCtrl + s\t\t\tSave selected segment as PNG\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tClose\n"

        signal_connect(win_key, "key-pressed") do widget, keyval, keycode, state
            if ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('s'))
                save_dialog("Pick an image file", win, ["*.png"]) do file_name
                    if file_name != ""
                        surface_buf = Gtk4.cairo_surface(can)
                        if Cairo.write_to_png(surface_buf, file_name) == Cairo.STATUS_SUCCESS
                            _info("Plot saved as: $file_name")
                        else
                            warn_dialog(_nill, "File $file_name cannot be written!", win)
                        end
                    end
                end
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('h'))
                info_dialog(_nill, help, win)
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == 0x0000ff0d) # Enter
                close(win)
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('q'))
                x_pos = nothing
                y_pos = nothing
                close(win)
            end
        end
    end

    app = GtkApplication("org.neuroanalyzer.iselect_seg")
    Gtk4.signal_connect(_activate, app, :activate)
    Gtk4.GLib.stop_main_loop()
    Gtk4.run(app)

    if x_pos === nothing && y_pos === nothing
        return nothing
    end

    if shape in [:r, :c]
        if x_pos !== nothing && y_pos !== nothing && length(x_pos) > 0 && length(y_pos) > 0
            if length(x_pos) > 1 && length(y_pos) > 1 && x_pos[end] == x_pos[1] && y_pos[end] == y_pos[1]
                pop!(x_pos)
                pop!(y_pos)
            end
            x_pos[1] > size_x && (x_pos[1] = size_x)
            x_pos[2] > size_x && (x_pos[2] = size_x)
            y_pos[1] > size_y && (y_pos[1] = size_y)
            y_pos[2] > size_y && (y_pos[2] = size_y)
            c1 = x_pos[1]
            c2 = x_pos[2]
            r1 = y_pos[1]
            r2 = y_pos[2]
            if shape === :r
                r1 > r2 && ((r1, r2) = _swap(r1, r2))
                c1 > c2 && ((c1, c2) = _swap(c1, c2))
            end
            c = shape == :c
            return !extract ? (r1, c1, r2, c2) : seg_extract(m, (r1, c1, r2, c2), v=v, c=c)
        end
    else
        x_pos[1] > size_x && (x_pos[1] = size_x)
        y_pos[1] > size_y && (y_pos[1] = size_y)
        c = x_pos[1]
        r = y_pos[1]
        return !extract ? (r, c) : m[r, c]
    end

end
