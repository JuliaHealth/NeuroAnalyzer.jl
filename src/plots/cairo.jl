export plot2canvas
export resize_canvas
export add_topmargin_canvas
export add_to_canvas

"""
    plot2canvas(p)

Render a `Plots.Plot` (GR backend) to an RGB Cairo surface.

# Arguments

- `p::Plots.Plot{Plots.GRBackend}`

# Returns

- `Cairo.CairoSurfaceBase{UInt32}`
"""
function plot2canvas(p::Plots.Plot{Plots.GRBackend})::Cairo.CairoSurfaceBase{UInt32}
    p_size = p.attr[:size]
    c = CairoRGBSurface(p_size[1], p_size[2])
    cr = CairoContext(c)

    # create an in-memory buffer, render the plot as PNG into it, then seek back to the start before reading
    io = IOBuffer()
    withenv("GKSwstype" => "100") do
        return png(p, io)
    end
    seekstart(io)
    img = read_from_png(io)

    Cairo.set_source_surface(cr, img, 0, 0)
    Cairo.paint(cr)

    return c
end

"""
    plot2canvas(fig)

Render a `GLMakie.Figure` to an RGB Cairo surface via a temporary PNG file.

# Arguments

- `fig::GLMakie.Figure`

# Returns

- `Cairo.CairoSurfaceBase{UInt32}`
"""
function plot2canvas(fig::GLMakie.Figure)::Cairo.CairoSurfaceBase{UInt32}
    p_size = size(fig.scene)
    c = CairoRGBSurface(p_size[1], p_size[2])
    cr = CairoContext(c)

    # GLMakie can only save to a file path, not a stream, so a temporary file is used
    # it is removed after reading to avoid leaking disk space
    fname = tempname() * ".png"
    try
        GLMakie.save(fname, fig)
        img = read_from_png(fname)
        Cairo.set_source_surface(cr, img, 0, 0)
        Cairo.paint(cr)
    finally
        isfile(fname) && rm(fname)
    end

    return c
end

"""
    resize_canvas(c; <keyword arguments>)

Resize a Cairo surface by a uniform scale factor.

# Arguments

- `c::Cairo.CairoSurfaceBase{UInt32}`
- `r::Real`: scale factor (e.g. `0.5` halves both dimensions)

# Returns

- `Cairo.CairoSurfaceBase{UInt32}`: new surface of size `(w*r, h*r)`
"""
function resize_canvas(
    c::Cairo.CairoSurfaceBase{UInt32};
    r::Real,
)::Cairo.CairoSurfaceBase{UInt32}
    # use round consistently for both axes to get the nearest integer size.
    new_w = round(Int64, c.width * r)
    new_h = round(Int64, c.height * r)

    c_new = CairoRGBSurface(new_w, new_h)
    cr = CairoContext(c_new)

    # fill background white before scaling so unpainted areas are not black
    Cairo.set_source_rgb(cr, 1.0, 1.0, 1.0)
    Cairo.rectangle(cr, 0.0, 0.0, new_w, new_h)
    Cairo.fill(cr)

    # apply the scale transform before painting the source surface so that
    # Cairo maps source pixels to the scaled coordinate space
    Cairo.scale(cr, r, r)
    Cairo.set_source_surface(cr, c, 0, 0)
    Cairo.paint(cr)

    return c_new
end

"""
    add_topmargin_canvas(c1, c2)

Create a new canvas with `c2` placed at the top and `c1` below it, effectively adding a top margin of `c2.height` pixels above `c1`.

# Arguments

- `c1::Cairo.CairoSurfaceBase{UInt32}`: main canvas
- `c2::Cairo.CairoSurfaceBase{UInt32}`: canvas to place above `c1`

# Returns

- `Cairo.CairoSurfaceBase{UInt32}`: combined canvas of height `c1.height + c2.height`
"""
function add_topmargin_canvas(
    c1::Cairo.CairoSurfaceBase{UInt32},
    c2::Cairo.CairoSurfaceBase{UInt32},
)::Cairo.CairoSurfaceBase{UInt32}
    total_h = c1.height + c2.height
    c = CairoRGBSurface(c1.width, total_h)
    cr = CairoContext(c)

    Cairo.set_source_rgb(cr, 1.0, 1.0, 1.0)
    Cairo.rectangle(cr, 0.0, 0.0, c1.width, total_h)
    Cairo.fill(cr)

    # paint c2 at the very top of the combined canvas
    Cairo.set_source_surface(cr, c2, 0, 0)
    Cairo.paint(cr)

    # paint c1 directly below c2 (offset by c2's height)
    Cairo.set_source_surface(cr, c1, 0, c2.height)
    Cairo.paint(cr)

    return c
end

"""
    add_to_canvas(c1, c2; <keyword arguments>)

Composite `c2` onto `c1` at position `(x, y)`, optionally adding a title label below `c2` and/or saving the result to a PNG file.

# Arguments

- `c1::Cairo.CairoSurfaceBase{UInt32}`: destination canvas
- `c2::Cairo.CairoSurfaceBase{UInt32}`: source canvas to place on `c1`
- `x::Int64`: left edge of `c2` in `c1` coordinates
- `y::Int64`: top edge of `c2` in `c1` coordinates
- `title::String=""`: label placed centred below `c2`
- `file_name::String=""`: if non-empty, save as PNG at this path

# Returns

- `c::Cairo.CairoSurfaceBase{UInt32}`
"""
function add_to_canvas(
    c1::Cairo.CairoSurfaceBase{UInt32},
    c2::Cairo.CairoSurfaceBase{UInt32};
    x::Int64,
    y::Int64,
    title::String = "",
    file_name::String = "",
)::Cairo.CairoSurfaceBase{UInt32}
    # create output canvas matching c1's size, paint c1 as background, then composite c2 at the requested position
    c = CairoRGBSurface(c1.width, c1.height)
    cr = CairoContext(c)
    Cairo.set_source_surface(cr, c1, 0, 0)
    Cairo.paint(cr)
    Cairo.set_source_surface(cr, c2, x, y)
    Cairo.paint(cr)

    # optionally render a centered title label below c2
    if title != ""
        Cairo.set_font_size(cr, 10.0)
        Cairo.set_source_rgb(cr, 0.0, 0.0, 0.0)
        extents = Cairo.text_extents(cr, title)

        title_x = x + div(c2.width, 2) + 6 - (extents[3] / 2 + extents[1])
        title_y = y + c2.height + 12

        Cairo.move_to(cr, title_x, title_y)
        Cairo.show_text(cr, title)
    end

    # optionally save the composed canvas to a PNG file
    if file_name != ""
        ext = lowercase(splitext(file_name)[2])
        ext == ".png" ||
            throw(ArgumentError("file_name extension must be .png, got \"$ext\"."))
        isfile(file_name) && _warn("File $file_name will be overwritten.")
        Cairo.write_to_png(c, file_name)
    end

    return c
end
