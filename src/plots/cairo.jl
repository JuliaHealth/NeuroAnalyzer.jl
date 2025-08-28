export plot2canvas
export resize_canvas
export add_topmargin_canvas
export add_to_canvas

"""
    plot2canvas(p)

Convert Plots.Plot to CairoSurfaceBase.

# Arguments

- `p::Plots.Plot{Plots.GRBackend}`

# Returns

- `c::Cairo.CairoSurfaceBase{UInt32}`
"""
function plot2canvas(p::Plots.Plot{Plots.GRBackend})::Cairo.CairoSurfaceBase{UInt32}

    p_size = p.attr[:size]
    c = CairoRGBSurface(p_size[1], p_size[2])
    cr = CairoContext(c)
    show(io, MIME("image/png"), p)
    img = read_from_png(io)
    Cairo.set_source_surface(cr, img, 0, 0)
    Cairo.paint(cr)

    return c

end

"""
    resize_canvas(c; <keyword arguments>)

Resize CairoSurfaceBase by a factor.

# Arguments

- `c::Cairo.CairoSurfaceBase{UInt32}`
- `r::Real`: resizing factor

# Returns

- `c_new::Cairo.CairoSurfaceBase{UInt32}`
"""
function resize_canvas(c::Cairo.CairoSurfaceBase{UInt32}; r::Real)::Cairo.CairoSurfaceBase{UInt32}

    c_new = CairoRGBSurface(ceil(Int64, c.width * r) - 1, round(Int64, c.height * r) - 1)
    cr = CairoContext(c_new)
    Cairo.scale(cr, r, r)
    # Cairo.set_source_rgb(cr, 256, 256, 256)
    # Cairo.rectangle(cr, 0, 0, ceil(Int64, c.width * r) + 10, ceil(Int64, c.height * r) + 10)
    # Cairo.fill(cr)
    Cairo.set_source_surface(cr, c, 0, 0)
    Cairo.paint(cr)

    return c_new

end

"""
    add_topmargin_canvas(c1, c2)

Resize CairoSurfaceBase to make space for another canvas

# Arguments

- `c1::Cairo.CairoSurfaceBase{UInt32}`
- `c2::Cairo.CairoSurfaceBase{UInt32}`

# Returns

- `c::Cairo.CairoSurfaceBase{UInt32}`
"""
function add_topmargin_canvas(c1::Cairo.CairoSurfaceBase{UInt32}, c2::Cairo.CairoSurfaceBase{UInt32})::Cairo.CairoSurfaceBase{UInt32}

    c = CairoRGBSurface(c1.width, c1.height + c2.height)
    cr = CairoContext(c)
    Cairo.set_source_rgb(cr, 256, 256, 256)
    Cairo.rectangle(cr, 0.0, 0.0, c1.width, c1.height + c2.height)
    Cairo.fill(cr)
    Cairo.set_source_surface(cr, c1, 0, c2.height)
    Cairo.paint(cr)

    return c

end

"""
    add_to_canvas(c1, c2; <keyword arguments>)

Place CairoSurfaceBase at another canvas at `x, y`. If `file_name` is provided, the plot is saved as PNG file.

# Arguments

- `c1::Cairo.CairoSurfaceBase{UInt32}`
- `c2::Cairo.CairoSurfaceBase{UInt32}`
- `x::Int64`
- `y::Int64`
- `title::String=""`: title of the subplot
- `view::Bool=true`: view the output image
- `file_name::String=""`: output image filename

# Returns

- `c::Cairo.CairoSurfaceBase{UInt32}`
"""
function add_to_canvas(c1::Cairo.CairoSurfaceBase{UInt32}, c2::Cairo.CairoSurfaceBase{UInt32}; x::Int64, y::Int64, title::String="", view::Bool=true, file_name::String="")::Cairo.CairoSurfaceBase{UInt32}

    if file_name != ""
        ext = lowercase(splitext(file_name)[2])
        @assert ext == ".png" "Filename extension must be .png"
        isfile(file_name) && _warn("File $file_name will be overwritten.")
    end

    c = CairoRGBSurface(c1.width, c1.height)
    cr = CairoContext(c)
    Cairo.set_source_surface(cr, c1, 0, 0)
    Cairo.paint(cr)
    Cairo.set_source_surface(cr, c2, x, y)
    Cairo.paint(cr)
    if title != ""
        Cairo.set_font_size(cr, 10.0)
        Cairo.set_source_rgb(cr, 0, 0, 0)
        extents = Cairo.text_extents(cr, title)
        x = x + div(c2.width, 2) + 6 - (extents[3] / 2 + extents[1])
        # y = y + div(c2.height, 2) + 6 - (extents[4] / 2 + extents[2])
        y = y + c2.height + 12
        # - (extents[4] / 2 + extents[2])
        Cairo.move_to(cr, x, y)
        Cairo.show_text(cr, title)
    end

    file_name != "" && Cairo.write_to_png(c, file_name)
    view && iview_plot(c)

    return c

end
