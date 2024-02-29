export iview_plot

"""
    iview_plot(p)

View plot object.

# Arguments

- `p::Plots.Plot{Plots.GRBackend}`
"""
function iview_plot(p::Plots.Plot{Plots.GRBackend})

    win = GtkWindow("NeuroAnalyzer: iview_plot()", p.attr[:size][1] + 2, p.attr[:size][2] + 2)
    set_gtk_property!(win, :border_width, 0)
    set_gtk_property!(win, :resizable, false)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    set_gtk_property!(win, :startup_id, "org.neuroanalyzer")
    can = GtkCanvas(p.attr[:size][1] + 2, p.attr[:size][2] + 2)
    push!(win, can)
    showall(win)

    @guarded draw(can) do widget
        show(io, MIME("image/png"), p)
        img = read_from_png(io)
        w = img.width
        h = img.height
        ctx = getgc(can)
        # Cairo.scale(ctx, 242/w, 241/h)
        Cairo.set_source_surface(ctx, img, 1, 2)
        Cairo.paint(ctx)
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
                            warn_dialog("Incorrect filename!")
                        end
                    end
            elseif k == 113 # q
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

    return nothing

end

"""
    iview_plot(file_name)

View PNG image.

# Arguments

- `file_name::String`
"""
function iview_plot(file_name::String)

    @assert isfile(file_name) "File $file_name cannot be opened."
    if splitext(file_name)[2] != ".png"
        @error "Incorrect filename!"
        return nothing
    end

    img = read_from_png(file_name)

    win = GtkWindow("NeuroAnalyzer: iview_plot()", img.width, img.height)
    set_gtk_property!(win, :border_width, 0)
    set_gtk_property!(win, :resizable, false)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    set_gtk_property!(win, :startup_id, "org.neuroanalyzer")
    can = GtkCanvas(img.width, img.height)
    push!(win, can)
    showall(win)

    @guarded draw(can) do widget
        ctx = getgc(can)
        Cairo.set_source_surface(ctx, img, 0, 0)
        Cairo.paint(ctx)
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
                            warn_dialog("Incorrect filename!")
                        end
                    end
            elseif k == 111 # o
                info_dialog("This feature has not been implemented yet.")
            elseif k == 113 # q
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

    return nothing

end

"""
    iview_plot(c)

View Cairo surface object.

# Arguments

- `c::Cairo.CairoSurfaceBase{UInt32}`
"""
function iview_plot(c::Cairo.CairoSurfaceBase{UInt32})

    win = GtkWindow("NeuroAnalyzer: iview_plot()", c.width, c.height)
    set_gtk_property!(win, :border_width, 0)
    set_gtk_property!(win, :resizable, false)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    set_gtk_property!(win, :startup_id, "org.neuroanalyzer")
    can = GtkCanvas(c.width, c.height)
    push!(win, can)
    showall(win)

    @guarded draw(can) do widget
        ctx = getgc(can)
        Cairo.set_source_rgb(ctx, 255, 255, 255)
        Cairo.set_source_surface(ctx, c, 1, 1)
        Cairo.paint(ctx)
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
                            warn_dialog("Incorrect filename!")
                        end
                    end
            elseif k == 113 # q
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

    return nothing

end