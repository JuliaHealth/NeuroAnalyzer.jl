# Interactive channel / location editor (GLMakie)
# opened by edit_locs(obj) and edit_locs!(obj)

# characters accepted by the numeric text boxes
_el_numchar(c::Char)::Bool = isdigit(c) || c in ('.', '-', '+', 'e', 'E')

# format a coordinate for a text box
_el_fmt(v::Real)::String = string(round(Float64(v); digits = 4))

# Cartesian positions of all location rows: from loc_x/y/z, or from spherical coordinates
function _el_xyz(locs::DataFrame, cart::Bool)
    cart && return Float64.(locs.loc_x), Float64.(locs.loc_y), Float64.(locs.loc_z)
    p = [
        sph2cart(locs.loc_radius_sph[i], locs.loc_theta_sph[i], locs.loc_phi_sph[i]) for
        i = 1:DataFrames.nrow(locs)
    ]
    return Float64.(first.(p)), Float64.(getindex.(p, 2)), Float64.(last.(p))
end

# 2-D (top view) positions: from loc_x/y, or from polar coordinates
function _el_top(locs::DataFrame, cart::Bool)
    cart && return Float64.(locs.loc_x), Float64.(locs.loc_y)
    p = [pol2cart(locs.loc_radius[i], locs.loc_theta[i]) for i = 1:DataFrames.nrow(locs)]
    return Float64.(first.(p)), Float64.(last.(p))
end

"""
    _edit_locs_build(obj)

Build the interactive location editor for `obj` without displaying it. All edits are made
on a copy of `obj`.

# Returns

Named tuple:

- `fig::GLMakie.Figure`: editor window
- `obj::NeuroAnalyzer.NEURO`: edited copy of `obj`
- `applied::Base.RefValue{Bool}`: `true` after the *Apply* button was clicked
- `ctl::Dict{Symbol, Any}`: widgets (used by the tests)
- `cur::Observable{Int64}`: current channel number
"""
function _edit_locs_build(obj::NeuroAnalyzer.NEURO)
    _check_chtype(obj, ["eeg", "meg"])
    nchannels(obj) > 0 || throw(ArgumentError("OBJ must contain ≥ 1 channel."))

    # working copy: everything is edited here, the input object is not touched
    w = deepcopy(obj)
    DataFrames.nrow(w.locs) == 0 && (w.locs = _initialize_locs(w))

    cur = Observable(1)
    applied = Ref(false)
    undo = DataFrame[]
    syncing = Ref(false)
    pending = Ref(:none)
    ctl = Dict{Symbol, Any}()

    # ---------------------------------------------------------------- #
    # layout                                                            #
    # ---------------------------------------------------------------- #
    fig = GLMakie.Figure(; size = (1500, 950))
    ctrl = fig[1, 1] = GLMakie.GridLayout(; valign = :top, tellheight = false)
    views = fig[1, 2] = GLMakie.GridLayout()
    GLMakie.colsize!(fig.layout, 1, GLMakie.Fixed(440))

    ax_xy = GLMakie.Axis(views[1, 1]; title = "Top (X-Y)", xlabel = "X", ylabel = "Y")
    ax_xz = GLMakie.Axis(views[1, 2]; title = "Front (X-Z)", xlabel = "X", ylabel = "Z")
    ax_yz = GLMakie.Axis(views[2, 1]; title = "Side (Y-Z)", xlabel = "Y", ylabel = "Z")
    ax_3d = GLMakie.Axis3(
        views[2, 2];
        title = "3-D",
        xlabel = "X",
        ylabel = "Y",
        zlabel = "Z",
    )
    for ax in (ax_xy, ax_xz, ax_yz)
        ax.aspect = DataAspect()
    end
    ax_3d.aspect = :data

    row = Ref(0)
    nextrow() = (row[] += 1)
    section(txt) = GLMakie.Label(
        ctrl[nextrow(), 1],
        txt;
        font = :bold,
        halign = :left,
        fontsize = 15,
    )
    # numeric text box: only number characters can be typed
    numbox(pos, width, txt = "") = GLMakie.Textbox(
        pos;
        width = width,
        validator = Float64,
        restriction = _el_numchar,
        stored_string = txt,
    )
    hgrid() = (ctrl[nextrow(), 1] = GLMakie.GridLayout(; halign = :left))

    status = Observable("Click an electrode in a 2-D view to select its channel.")

    # ---------------------------------------------------------------- #
    # helpers                                                           #
    # ---------------------------------------------------------------- #
    msg!(s) = (status[] = s; nothing)
    loc_row(ch::Int64) =
        (r = _find_bylabel(w.locs, labels(w)[ch]); r isa Int64 ? r : nothing)

    # rows modified by the geometric operations: reference and EOG electrodes are kept as
    # they are
    function op_rows()::Vector{Int64}
        skip = lowercase.(String.(get_channel(w; type = ["ref", "eog"])))
        return findall(l -> !(lowercase(l) in skip), w.locs[!, :label])
    end

    function snapshot!()
        push!(undo, copy(w.locs))
        length(undo) > 50 && popfirst!(undo)
        return nothing
    end

    # apply `f!` (modifying a DataFrame in-place) to the selected location rows
    function apply_op!(f!::Function, what::String; rows::Vector{Int64} = op_rows())
        isempty(rows) && return msg!("No locations to modify.")
        snapshot!()
        sub = w.locs[rows, :]
        try
            f!(sub)
        catch err
            pop!(undo)
            return msg!("$what failed: $(sprint(showerror, err))")
        end
        w.locs[rows, :] = sub
        redraw!()
        return msg!("$what applied.")
    end

    flags() = (
        polar = ctl[:cb_polar].checked[],
        cart = ctl[:cb_cart].checked[],
        spherical = ctl[:cb_sph].checked[],
    )

    # value shown in a numeric text box: typed text counts even if Enter was not pressed
    num(tb) = tryparse(Float64, something(tb.displayed_string[], ""))

    # ---------------------------------------------------------------- #
    # channel                                                           #
    # ---------------------------------------------------------------- #
    section("Channel")
    g = hgrid()
    bt_prev = GLMakie.Button(g[1, 1]; label = "◀", width = 36)
    menu_ch = GLMakie.Menu(g[1, 2]; options = labels(w), width = 220)
    bt_next = GLMakie.Button(g[1, 3]; label = "▶", width = 36)
    g = hgrid()
    GLMakie.Label(g[1, 1], "Label")
    tb_label = GLMakie.Textbox(g[1, 2]; width = 90, stored_string = labels(w)[1])
    GLMakie.Label(g[1, 3], "Type")
    menu_type = GLMakie.Menu(g[1, 4]; options = channel_types[2:end], width = 90)
    GLMakie.Label(g[1, 5], "Unit")
    unit_opts = [u == "" ? "none" : u for u in channel_units]
    menu_unit = GLMakie.Menu(g[1, 6]; options = unit_opts, width = 70)

    # ---------------------------------------------------------------- #
    # coordinates of the current channel                                #
    # ---------------------------------------------------------------- #
    section("Coordinates (stored values, original units)")
    coords = [
        ("Cartesian", [("X", :loc_x), ("Y", :loc_y), ("Z", :loc_z)]),
        ("Spherical", [("r", :loc_radius_sph), ("θ", :loc_theta_sph), ("φ", :loc_phi_sph)]),
        ("Polar", [("r", :loc_radius), ("θ", :loc_theta)]),
    ]
    tb_coord = Dict{Symbol, Any}()
    for (name, fields) in coords
        g = hgrid()
        GLMakie.Label(g[1, 1], name; width = 70, halign = :left)
        for (i, (lab, col)) in enumerate(fields)
            GLMakie.Label(g[1, 2i], lab)
            tb_coord[col] = numbox(g[1, 2i + 1], 80)
        end
    end

    # ---------------------------------------------------------------- #
    # operations on all locations                                       #
    # ---------------------------------------------------------------- #
    section("Transform all locations (REF and EOG excluded)")
    g = hgrid()
    GLMakie.Label(g[1, 1], "Apply to:")
    ctl[:cb_polar] = GLMakie.Checkbox(g[1, 2]; checked = true)
    GLMakie.Label(g[1, 3], "polar")
    ctl[:cb_cart] = GLMakie.Checkbox(g[1, 4]; checked = true)
    GLMakie.Label(g[1, 5], "Cartesian")
    ctl[:cb_sph] = GLMakie.Checkbox(g[1, 6]; checked = true)
    GLMakie.Label(g[1, 7], "spherical")

    axes3 = ["X", "Y", "Z"]
    g = hgrid()
    bt_flip = GLMakie.Button(g[1, 1]; label = "Flip", width = 90)
    menu_flip = GLMakie.Menu(g[1, 2]; options = axes3, width = 60)
    bt_swap = GLMakie.Button(g[1, 3]; label = "Swap X/Y", width = 90)
    bt_center = GLMakie.Button(g[1, 4]; label = "Center", width = 90)

    g = hgrid()
    bt_rot = GLMakie.Button(g[1, 1]; label = "Rotate", width = 90)
    menu_rot = GLMakie.Menu(g[1, 2]; options = axes3, width = 60)
    tb_rot = numbox(g[1, 3], 70, "0")

    g = hgrid()
    bt_scale = GLMakie.Button(g[1, 1]; label = "Scale", width = 90)
    tb_scale = numbox(g[1, 2], 70, "1.0")
    bt_norm = GLMakie.Button(g[1, 3]; label = "Normalize", width = 90)

    g = hgrid()
    bt_origin = GLMakie.Button(g[1, 1]; label = "Move origin", width = 90)
    menu_origin = GLMakie.Menu(g[1, 2]; options = axes3, width = 60)
    tb_origin = numbox(g[1, 3], 70, "0")

    g = hgrid()
    conversions = [
        "Cartesian → polar" => locs_cart2pol!,
        "Cartesian → spherical" => locs_cart2sph!,
        "spherical → Cartesian" => locs_sph2cart!,
        "spherical → polar" => locs_sph2pol!,
        "polar → Cartesian" => locs_pol2cart!,
        "polar → spherical" => locs_pol2sph!,
    ]
    bt_conv = GLMakie.Button(g[1, 1]; label = "Convert", width = 90)
    menu_conv = GLMakie.Menu(g[1, 2]; options = first.(conversions), width = 190)

    g = hgrid()
    bt_gen = GLMakie.Button(g[1, 1]; label = "Generate 10-10", width = 130)
    bt_undo = GLMakie.Button(g[1, 2]; label = "Undo", width = 90)

    # ---------------------------------------------------------------- #
    # file, view and window                                             #
    # ---------------------------------------------------------------- #
    section("File")
    g = hgrid()
    tb_file = GLMakie.Textbox(g[1, 1]; width = 220, placeholder = "locations file path")
    bt_load = GLMakie.Button(g[1, 2]; label = "Load", width = 70)
    bt_save = GLMakie.Button(g[1, 3]; label = "Save", width = 70)

    section("View")
    g = hgrid()
    cb_vcart = GLMakie.Checkbox(g[1, 1]; checked = true)
    GLMakie.Label(g[1, 2], "Cartesian coordinates")
    cb_labels = GLMakie.Checkbox(g[1, 3]; checked = false)
    GLMakie.Label(g[1, 4], "labels")

    g = hgrid()
    bt_delete = GLMakie.Button(g[1, 1]; label = "Delete channel", width = 120)
    bt_apply = GLMakie.Button(g[1, 2]; label = "Apply", width = 90)
    bt_cancel = GLMakie.Button(g[1, 3]; label = "Cancel", width = 90)

    GLMakie.Label(
        ctrl[nextrow(), 1],
        status;
        halign = :left,
        word_wrap = true,
        width = 420,
        color = :gray25,
    )

    # ---------------------------------------------------------------- #
    # drawing                                                           #
    # ---------------------------------------------------------------- #
    # show `txt` as the (stored) content of a text box; callbacks are muted by `syncing`
    function settext!(tb, txt::String)
        tb.stored_string[] = txt
        tb.displayed_string[] = txt
        return nothing
    end

    function sync_widgets!()
        syncing[] = true
        try
            ch = cur[]
            menu_ch.options[] = labels(w)
            menu_ch.i_selected[] = ch
            settext!(tb_label, labels(w)[ch])
            t = findfirst(==(w.header.recording[:channel_type][ch]), channel_types[2:end])
            !isnothing(t) && (menu_type.i_selected[] = t)
            u = findfirst(==(w.header.recording[:unit][ch]), channel_units)
            !isnothing(u) && (menu_unit.i_selected[] = u)
            r = loc_row(ch)
            for (col, tb) in tb_coord
                settext!(tb, isnothing(r) ? "0.0" : _el_fmt(w.locs[r, col]))
            end
        finally
            syncing[] = false
        end
        return nothing
    end

    function redraw!()
        cart = cb_vcart.checked[]
        lbl = String.(w.locs[!, :label])
        r = loc_row(cur[])
        tx, ty = _el_top(w.locs, cart)
        x, y, z = _el_xyz(w.locs, cart)
        for (ax, a, b) in ((ax_xy, tx, ty), (ax_xz, x, z), (ax_yz, y, z))
            empty!(ax)
            GLMakie.hlines!(ax, 0; color = (:gray, 0.4), linestyle = :dash)
            GLMakie.vlines!(ax, 0; color = (:gray, 0.4), linestyle = :dash)
            isempty(a) && continue
            GLMakie.scatter!(ax, a, b; color = :steelblue, markersize = 9)
            !isnothing(r) &&
                GLMakie.scatter!(ax, [a[r]], [b[r]]; color = :red, markersize = 16)
            cb_labels.checked[] &&
                GLMakie.text!(ax, a, b; text = lbl, fontsize = 10, offset = (4, 4))
            GLMakie.autolimits!(ax)
        end
        empty!(ax_3d)
        if !isempty(x)
            GLMakie.scatter!(ax_3d, x, y, z; color = :steelblue, markersize = 9)
            !isnothing(r) &&
                GLMakie.scatter!(
                    ax_3d,
                    [x[r]],
                    [y[r]],
                    [z[r]];
                    color = :red,
                    markersize = 16,
                )
            cb_labels.checked[] && GLMakie.text!(ax_3d, x, y, z; text = lbl, fontsize = 10)
        end
        sync_widgets!()
        return nothing
    end

    # ---------------------------------------------------------------- #
    # callbacks                                                         #
    # ---------------------------------------------------------------- #
    on(cur) do _
        pending[] = :none
        return redraw!()
    end
    on(_ -> (cur[] = cur[] == 1 ? nchannels(w) : cur[] - 1), bt_prev.clicks)
    on(_ -> (cur[] = cur[] == nchannels(w) ? 1 : cur[] + 1), bt_next.clicks)
    on(menu_ch.selection) do s
        syncing[] && return nothing
        ch = findfirst(==(s), labels(w))
        !isnothing(ch) && ch != cur[] && (cur[] = ch)
        return nothing
    end
    on(_ -> redraw!(), cb_vcart.checked)
    on(_ -> redraw!(), cb_labels.checked)

    # select a channel by clicking its electrode in a 2-D view
    for (ax, view) in ((ax_xy, :xy), (ax_xz, :xz), (ax_yz, :yz))
        on(GLMakie.events(ax.scene).mousebutton; priority = 2) do ev
            (ev.button == GLMakie.Mouse.left && ev.action == GLMakie.Mouse.press) ||
                return GLMakie.Consume(false)
            GLMakie.is_mouseinside(ax.scene) || return GLMakie.Consume(false)
            DataFrames.nrow(w.locs) == 0 && return GLMakie.Consume(false)
            cart = cb_vcart.checked[]
            x, y, z = _el_xyz(w.locs, cart)
            a, b =
                view === :xy ? _el_top(w.locs, cart) : view === :xz ? (x, z) : (y, z)
            p = GLMakie.mouseposition(ax)
            d = hypot.(a .- p[1], b .- p[2])
            i = argmin(d)
            span = max(maximum(a) - minimum(a), maximum(b) - minimum(b), eps())
            d[i] > 0.05 * span && return GLMakie.Consume(false)
            ch = findfirst(==(lowercase(w.locs[i, :label])), lowercase.(labels(w)))
            if isnothing(ch)
                msg!("$(w.locs[i, :label]) is not a channel of OBJ.")
                return GLMakie.Consume(false)
            end
            cur[] = ch
            return GLMakie.Consume(true)
        end
    end

    # text boxes without a button: commit typed text on Enter and when the box loses focus
    for tb in (tb_label, values(tb_coord)...)
        on(tb.focused) do f
            f && return nothing
            txt = tb.displayed_string[]
            (isnothing(txt) || txt == tb.stored_string[]) && return nothing
            tb === tb_label || !isnothing(tryparse(Float64, txt)) || return nothing
            tb.stored_string[] = txt
            return nothing
        end
    end

    on(tb_label.stored_string) do s
        syncing[] && return nothing
        old = labels(w)[cur[]]
        (isnothing(s) || s == old) && return nothing
        if isempty(strip(s)) || s in labels(w)
            msg!("Label \"$s\" is empty or already used.")
        else
            rename_channel!(w; ch = old, name = s)
            msg!("Channel $old renamed to $s.")
        end
        return redraw!()
    end

    on(menu_type.selection) do t
        (syncing[] || isnothing(t)) && return nothing
        t == w.header.recording[:channel_type][cur[]] && return nothing
        set_channel_type!(w; ch = labels(w)[cur[]], type = t)
        u = _ch_units(t)
        u in channel_units && (w.header.recording[:unit][cur[]] = u)
        msg!("Channel type set to $t.")
        return redraw!()
    end

    on(menu_unit.selection) do u
        (syncing[] || isnothing(u)) && return nothing
        w.header.recording[:unit][cur[]] = u == "none" ? "" : u
        return nothing
    end

    for (col, tb) in tb_coord
        on(tb.stored_string) do s
            syncing[] && return nothing
            v = tryparse(Float64, something(s, ""))
            isnothing(v) && return nothing
            snapshot!()
            r = loc_row(cur[])
            if isnothing(r)
                # channel without a location: add a row
                push!(
                    w.locs,
                    (labels(w)[cur[]], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
                    promote = true,
                )
                r = DataFrames.nrow(w.locs)
            end
            w.locs[r, col] = v
            return redraw!()
        end
    end

    on(bt_flip.clicks) do _
        f! = (locs_flipx!, locs_flipy!, locs_flipz!)[menu_flip.i_selected[]]
        return apply_op!(l -> f!(l; flags()...), "Flip $(menu_flip.selection[])")
    end
    on(_ -> apply_op!(l -> locs_swapxy!(l; flags()...), "Swap X/Y"), bt_swap.clicks)
    on(_ -> apply_op!(l -> locs_center!(l; flags()...), "Center"), bt_center.clicks)
    on(bt_rot.clicks) do _
        a = num(tb_rot)
        isnothing(a) && return msg!("Rotation angle must be a number.")
        f! = (locs_rotx!, locs_roty!, locs_rotz!)[menu_rot.i_selected[]]
        return apply_op!(l -> f!(l; a = a, flags()...), "Rotation by $a degrees")
    end
    on(bt_scale.clicks) do _
        s = num(tb_scale)
        (isnothing(s) || s <= 0) && return msg!("Scaling factor must be > 0.")
        return apply_op!(l -> locs_scale!(l; r = s, flags()...), "Scaling by $s")
    end
    on(_ -> apply_op!(l -> locs_normalize!(l; flags()...), "Normalization"), bt_norm.clicks)
    on(bt_origin.clicks) do _
        v = num(tb_origin)
        isnothing(v) && return msg!("Origin shift must be a number.")
        ax = (:x, :y, :z)[menu_origin.i_selected[]]
        return apply_op!(l -> locs_origin!(l; NamedTuple{(ax,)}((v,))...), "Origin shift")
    end
    on(bt_conv.clicks) do _
        name, f! = conversions[menu_conv.i_selected[]]
        return apply_op!(f!, name; rows = collect(1:DataFrames.nrow(w.locs)))
    end
    on(bt_gen.clicks) do _
        return apply_op!(
            locs_generate!,
            "10-10 positions";
            rows = collect(1:DataFrames.nrow(w.locs)),
        )
    end
    on(bt_undo.clicks) do _
        isempty(undo) && return msg!("Nothing to undo.")
        w.locs = pop!(undo)
        redraw!()
        return msg!("Last change undone.")
    end

    on(bt_load.clicks) do _
        f = String(strip(something(tb_file.displayed_string[], "")))
        isfile(f) || return msg!("File not found: $f")
        snapshot!()
        try
            load_locs!(w; file_name = f)
        catch err
            pop!(undo)
            return msg!("Loading failed: $(sprint(showerror, err))")
        end
        redraw!()
        return msg!("Locations loaded from $f.")
    end
    on(bt_save.clicks) do _
        f = String(strip(something(tb_file.displayed_string[], "")))
        isempty(f) && return msg!("Enter a file name.")
        if isfile(f) && pending[] !== :save
            pending[] = :save
            return msg!("$f exists: click Save again to overwrite it.")
        end
        pending[] = :none
        try
            export_locs(w; file_name = f, overwrite = true)
        catch err
            return msg!("Saving failed: $(sprint(showerror, err))")
        end
        return msg!("Locations saved to $f.")
    end

    on(bt_delete.clicks) do _
        ch = labels(w)[cur[]]
        nchannels(w) == 1 && return msg!("The last channel cannot be deleted.")
        if pending[] !== :delete
            pending[] = :delete
            return msg!("Click Delete channel again to delete $ch.")
        end
        pending[] = :none
        delete_channel!(w; ch = ch)
        # the undo history refers to the deleted channel's location rows
        empty!(undo)
        cur[] = min(cur[], nchannels(w))
        redraw!()
        return msg!("Channel $ch deleted.")
    end

    ctl[:apply] = bt_apply
    ctl[:cancel] = bt_cancel
    on(_ -> (applied[] = true), bt_apply.clicks)

    merge!(
        ctl,
        Dict(
            :prev => bt_prev,
            :next => bt_next,
            :menu_ch => menu_ch,
            :label => tb_label,
            :type => menu_type,
            :unit => menu_unit,
            :coord => tb_coord,
            :flip => bt_flip,
            :menu_flip => menu_flip,
            :swap => bt_swap,
            :center => bt_center,
            :rot => bt_rot,
            :tb_rot => tb_rot,
            :scale => bt_scale,
            :tb_scale => tb_scale,
            :norm => bt_norm,
            :origin => bt_origin,
            :tb_origin => tb_origin,
            :conv => bt_conv,
            :menu_conv => menu_conv,
            :gen => bt_gen,
            :undo => bt_undo,
            :file => tb_file,
            :load => bt_load,
            :save => bt_save,
            :delete => bt_delete,
            :status => status,
        ),
    )

    redraw!()

    return (; fig, obj = w, applied, ctl, cur)
end

"""
    _edit_locs_gui(obj)

Open the interactive location editor and wait until its window is closed.

# Returns

- `Union{Nothing, NeuroAnalyzer.NEURO}`: edited copy of `obj` if the changes were applied,
  `nothing` if cancelled
"""
function _edit_locs_gui(obj::NeuroAnalyzer.NEURO)::Union{Nothing, NeuroAnalyzer.NEURO}
    ed = _edit_locs_build(obj)
    screen = GLMakie.Screen(; title = "NeuroAnalyzer: edit_locs()")
    on(_ -> close(screen), ed.ctl[:apply].clicks)
    on(_ -> close(screen), ed.ctl[:cancel].clicks)
    display(screen, ed.fig)
    GLMakie.GLFW.MaximizeWindow(screen.glscreen)
    wait(screen)
    ed.applied[] || return nothing
    _history!(ed.obj, "edit_locs", NamedTuple(); note = "<not stored: interactive edits of the locations>")
    return ed.obj
end
