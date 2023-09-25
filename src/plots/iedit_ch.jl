export iedit_ch

"""
    iedit_ch(obj)

Interactive edit of signal channels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
"""
function iedit_ch(obj::NeuroAnalyzer.NEURO)

    # TO DO
    # - select channel by clicking its location
    # 

    _wip()

    @assert nchannels(obj) >= 1 "OBJ must contain ≥ 1 channel."

    if _has_locs(obj)
        locs = obj.locs
        locs_ch = locs[!, :channel]
    else
        _info("Creating empty channel locations data")
        locs = DataFrame(:channel=>signal_channels(obj), :labels=>labels(obj)[signal_channels(obj)], :loc_theta=>zeros(length(signal_channels(obj))), :loc_radius=>zeros(length(signal_channels(obj))), :loc_x=>zeros(length(signal_channels(obj))), :loc_y=>zeros(length(signal_channels(obj))), :loc_z=>zeros(length(signal_channels(obj))), :loc_radius_sph=>zeros(length(signal_channels(obj))), :loc_theta_sph=>zeros(length(signal_channels(obj))), :loc_phi_sph=>zeros(length(signal_channels(obj))))
        locs_ch = Int64[]
    end

    current_channel = 1
    channel_types = obj.header.recording[:channel_type]
    channel_labels = labels(obj)
    already_scaled1 = false
    already_scaled2 = false
    already_scaled3 = false
    already_scaled4 = false

    win = GtkWindow("NeuroAnalyzer: iedit_ch()", 1400, 900)
    set_gtk_property!(win, :border_width, 20)
    set_gtk_property!(win, :resizable, true)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)

    g = GtkGrid()
    set_gtk_property!(g, :column_homogeneous, false)
    set_gtk_property!(g, :column_spacing, 10)
    set_gtk_property!(g, :row_spacing, 10)

    g_opts = GtkGrid()
    set_gtk_property!(g_opts, :column_homogeneous, true)
    set_gtk_property!(g_opts, :column_spacing, 10)
    set_gtk_property!(g_opts, :row_spacing, 10)

    can1 = GtkCanvas(450, 450)
    can2 = GtkCanvas(450, 450)
    can3 = GtkCanvas(450, 450)
    can4 = GtkCanvas(450, 450)

    lab_chn = GtkLabel("Channel number:")
    set_gtk_property!(lab_chn, :halign, 2)
    entry_ch = GtkSpinButton(1, nchannels(obj), 1)
    set_gtk_property!(entry_ch, :value, current_channel)
    set_gtk_property!(entry_ch, :tooltip_text, "Channel number")
    bt_start = GtkButton("⇤")
    set_gtk_property!(bt_start, :tooltip_text, "Go to the first channel")
    bt_end = GtkButton("⇥")
    set_gtk_property!(bt_end, :tooltip_text, "Go to the last channel")

    lab_chtype = GtkLabel("Channel type:")
    set_gtk_property!(lab_chtype, :halign, 2)
    combo_chtype = GtkComboBoxText()
    for idx in string.(NeuroAnalyzer.channel_types[2:end])
        push!(combo_chtype, idx)
    end
    set_gtk_property!(combo_chtype, :active, findfirst(isequal(Symbol(channel_types[current_channel])), NeuroAnalyzer.channel_types) - 2)
    set_gtk_property!(combo_chtype, :tooltip_text, "Channel type")

    lab_chunits = GtkLabel("Channel units:")
    set_gtk_property!(lab_chunits, :halign, 2)
    combo_chunits = GtkComboBoxText()
    for idx in NeuroAnalyzer.channel_units
        push!(combo_chunits, idx)
    end
    set_gtk_property!(combo_chunits, :active, 0)
    set_gtk_property!(combo_chunits, :tooltip_text, "Channel units")

    lab_chlabel = GtkLabel("Channel label:")
    set_gtk_property!(lab_chlabel, :halign, 2)
    entry_label = GtkEntry()
    set_gtk_property!(entry_label, :text, channel_labels[current_channel])
    set_gtk_property!(entry_label, :tooltip_text, "Channel label")

    bt_delete = GtkButton("Delete")
    nchannels(obj) < 2 && set_gtk_property!(bt_delete, :sensitive, false)

    lab_loc_x = GtkLabel("Cartesian X:")
    set_gtk_property!(lab_loc_x, :halign, 2)
    set_gtk_property!(lab_loc_x, :tooltip_text, "Cartesian X (:loc_x)")
    lab_loc_y = GtkLabel("Cartesian Y:")
    set_gtk_property!(lab_loc_y, :halign, 2)
    set_gtk_property!(lab_loc_y, :tooltip_text, "Cartesian Y (:loc_y)")
    lab_loc_z = GtkLabel("Cartesian Z:")
    set_gtk_property!(lab_loc_z, :halign, 2)
    set_gtk_property!(lab_loc_z, :tooltip_text, "Cartesian Z (:loc_z)")
    lab_loc_z = GtkLabel("Cartesian Z:")
    set_gtk_property!(lab_loc_z, :halign, 2)
    set_gtk_property!(lab_loc_z, :tooltip_text, "Cartesian Z (:loc_z)")
    lab_loc_z = GtkLabel("Cartesian Z:")
    set_gtk_property!(lab_loc_z, :halign, 2)
    set_gtk_property!(lab_loc_z, :tooltip_text, "Cartesian Z (:loc_z)")
    lab_loc_theta_sph = GtkLabel("Spherical theta:")
    set_gtk_property!(lab_loc_theta_sph, :halign, 2)
    set_gtk_property!(lab_loc_theta_sph, :tooltip_text, "Spherical horiz. angle (:loc_theta_sph)")
    lab_loc_phi_sph = GtkLabel("Spherical phi:")
    set_gtk_property!(lab_loc_phi_sph, :halign, 2)
    set_gtk_property!(lab_loc_phi_sph, :tooltip_text, "Spherical azimuth angle (:loc_phi_sph)")
    lab_loc_radius_sph = GtkLabel("Spherical radius:")
    set_gtk_property!(lab_loc_radius_sph, :halign, 2)
    set_gtk_property!(lab_loc_radius_sph, :tooltip_text, "Spherical radius (:loc_radius_sph)")

    lab_loc_theta = GtkLabel("Polar theta:")
    set_gtk_property!(lab_loc_theta, :halign, 2)
    entry_loc_theta = GtkSpinButton(0, 360.0, 0.05)
    set_gtk_property!(entry_loc_theta, :tooltip_text, "Polar theta angle (:loc_theta)")
    set_gtk_property!(entry_loc_theta, :digits, 2)

    lab_loc_radius = GtkLabel("Polar radius:")
    set_gtk_property!(lab_loc_radius, :halign, 2)
    entry_loc_radius = GtkSpinButton(-1.5, 1.5, 0.01)
    set_gtk_property!(entry_loc_radius, :tooltip_text, "Polar radius (:loc_radius)")
    set_gtk_property!(entry_loc_radius, :digits, 2)
    if current_channel in signal_channels(obj)
        set_gtk_property!(entry_loc_radius, :sensitive, true)
        set_gtk_property!(entry_loc_radius, :value, locs[!, :loc_radius][current_channel])
    else
        set_gtk_property!(entry_loc_radius, :sensitive, false)
        set_gtk_property!(entry_loc_radius, :value, -1)
    end

    lab_ax_rot = GtkLabel("Rotate axis:")
    set_gtk_property!(lab_ax_rot, :halign, 2)
    lab_scale = GtkLabel("Scale:")
    set_gtk_property!(lab_scale, :halign, 2)
    lab_transform = GtkLabel("Transform:")
    set_gtk_property!(lab_transform, :halign, 2)

    bt_apply = GtkButton("Apply")
    bt_load = GtkButton("Load")
    bt_save = GtkButton("Save")
    bt_generate = GtkButton("Generate")

    if current_channel in locs_ch
        set_gtk_property!(entry_loc_theta, :value, locs[!, :loc_theta][current_channel])
        set_gtk_property!(entry_loc_radius, :value, locs[!, :loc_radius][current_channel])
    else
        set_gtk_property!(entry_loc_theta, :value, -1)
        set_gtk_property!(entry_loc_radius, :value, -1)
    end

    if current_channel in signal_channels(obj)
        set_gtk_property!(entry_loc_theta, :sensitive, true)
        set_gtk_property!(entry_loc_radius, :sensitive, true)
    else
        set_gtk_property!(entry_loc_theta, :sensitive, false)
        set_gtk_property!(entry_loc_radius, :sensitive, false)
    end

    lab_cart = GtkLabel("Plot using Cartesian coordinates:")
    set_gtk_property!(lab_cart, :halign, 2)
    cb_cart = GtkCheckButton("")
    set_gtk_property!(cb_cart, :tooltip_text, "Plot using Cartesian coordinates")
    set_gtk_property!(cb_cart, :active, false)

    g_opts[1:4, 1] = GtkLabel("Channel number")
    g_opts[1, 2] = bt_start
    g_opts[2, 2] = entry_ch
    g_opts[3, 2] = bt_end
    g_opts[4, 2] = bt_delete
    g_opts[1:4, 3] = GtkLabel("Channel properties")
    g_opts[1, 4] = lab_chlabel
    g_opts[2:4, 4] = entry_label
    g_opts[1, 5] = lab_chtype
    g_opts[2:4, 5] = combo_chtype
    g_opts[1, 6] = lab_chunits
    g_opts[2:4, 6] = combo_chunits
    g_opts[1:4, 7] = GtkLabel("Edit coordinates")
    g_opts[1, 8] = lab_loc_theta
    g_opts[2:4, 8] = entry_loc_theta
    g_opts[1, 9] = lab_loc_radius
    g_opts[2:4, 9] = entry_loc_radius
    g_opts[1, 10] = lab_loc_x
    g_opts[1, 11] = lab_loc_y
    g_opts[1, 12] = lab_loc_z
    g_opts[1, 13] = lab_loc_theta_sph
    g_opts[1, 14] = lab_loc_phi_sph
    g_opts[1, 15] = lab_loc_radius_sph
    g_opts[1:4, 16] = GtkLabel("Edit locs")
    g_opts[1, 17] = lab_ax_rot
    g_opts[1, 18] = lab_scale
    g_opts[1, 19] = lab_transform
    g_opts[1:4, 20] = GtkLabel("Locs operations")
    g_opts[1, 21] = bt_apply
    g_opts[2, 21] = bt_generate
    g_opts[3, 21] = bt_load
    g_opts[4, 21] = bt_save
    g_opts[1:3, 22] = lab_cart
    g_opts[4, 22] = cb_cart
    vbox = GtkBox(:v)
    push!(vbox, g_opts)

    g[1, 1:2] = vbox
    g[2, 1] = can1
    g[3, 1] = can2
    g[2, 2] = can3
    g[3, 2] = can4

    push!(win, g)

    showall(win)

    @guarded draw(can1) do widget
        cart = get_gtk_property(cb_cart, :active, Bool)
        ch = signal_channels(Symbol(obj.header.recording[:data_type]), channel_types)
        if current_channel in ch
            selected = current_channel
        else
            selected = 0
        end
            p = NeuroAnalyzer.plot_locs(locs, ch=ch, selected=selected, ch_labels=false, head_labels=false, cart=cart, plane=:xy, grid=true)
            img = read_from_png(io)
            ctx = getgc(can1)
            if already_scaled1 == false
                Cairo.scale(ctx, 0.70, 0.70)
                already_scaled1 = true
            end
            rectangle(ctx, 0, 0, 1200, 1200)
            set_source_rgb(ctx, 1, 1, 1)
            fill(ctx)
            show(io, MIME("image/png"), p)
            img = read_from_png(io)
            set_source_surface(ctx, img, 0, 0)
            paint(ctx)
    end

    @guarded draw(can2) do widget
        cart = get_gtk_property(cb_cart, :active, Bool)
        ch = signal_channels(Symbol(obj.header.recording[:data_type]), channel_types)
        if current_channel in ch
            selected = current_channel
        else
            selected = 0
        end
            p = NeuroAnalyzer.plot_locs(locs, ch=ch, selected=selected, ch_labels=false, head_labels=false, cart=cart, plane=:xz, grid=true)
            img = read_from_png(io)
            ctx = getgc(can2)
            if already_scaled2 == false
                Cairo.scale(ctx, 0.70, 0.70)
                already_scaled2 = true
            end
            rectangle(ctx, 0, 0, 1200, 1200)
            set_source_rgb(ctx, 1, 1, 1)
            fill(ctx)
            show(io, MIME("image/png"), p)
            img = read_from_png(io)
            set_source_surface(ctx, img, 0, 0)
            paint(ctx)
    end

    @guarded draw(can3) do widget
        cart = get_gtk_property(cb_cart, :active, Bool)
        ch = signal_channels(Symbol(obj.header.recording[:data_type]), channel_types)
        if current_channel in ch
            selected = current_channel
        else
            selected = 0
        end
            p = NeuroAnalyzer.plot_locs(locs, ch=ch, selected=selected, ch_labels=false, head_labels=false, cart=cart, plane=:yz, grid=true)
            img = read_from_png(io)
            ctx = getgc(can3)
            if already_scaled3 == false
                Cairo.scale(ctx, 0.70, 0.70)
                already_scaled3 = true
            end
            rectangle(ctx, 0, 0, 1200, 1200)
            set_source_rgb(ctx, 1, 1, 1)
            fill(ctx)
            show(io, MIME("image/png"), p)
            img = read_from_png(io)
            set_source_surface(ctx, img, 0, 0)
            paint(ctx)
    end

    @guarded draw(can4) do widget
        cart = get_gtk_property(cb_cart, :active, Bool)
        ch = signal_channels(Symbol(obj.header.recording[:data_type]), channel_types)
        if current_channel in ch
            selected = current_channel
        else
            selected = 0
        end
            p = NeuroAnalyzer.plot_locs3d(locs, ch=ch, selected=selected, ch_labels=false, head_labels=false, cart=cart);
            img = read_from_png(io)
            ctx = getgc(can4)
            if already_scaled4 == false
                Cairo.scale(ctx, 0.70, 0.70)
                already_scaled4 = true
            end
            rectangle(ctx, 0, 0, 1200, 1200)
            set_source_rgb(ctx, 1, 1, 1)
            fill(ctx)
            show(io, MIME("image/png"), p)
            img = read_from_png(io)
            set_source_surface(ctx, img, 0, 0)
            paint(ctx)
    end

    signal_connect(bt_start, "clicked") do widget
        current_channel = 1
        Gtk.@sigatom begin
            set_gtk_property!(entry_ch, :value, current_channel)
        end
    end

    signal_connect(bt_end, "clicked") do widget
        current_channel = nchannels(obj)
        Gtk.@sigatom begin
            set_gtk_property!(entry_ch, :value, current_channel)
        end
    end

    signal_connect(combo_chtype, "changed") do widget
        Gtk.@sigatom begin
            if current_channel in signal_channels(obj)
                set_gtk_property!(entry_loc_theta, :sensitive, true)
                set_gtk_property!(entry_loc_radius, :sensitive, true)
            else
                set_gtk_property!(entry_loc_theta, :sensitive, false)
                set_gtk_property!(entry_loc_radius, :sensitive, false)
            end
            if current_channel in locs_ch
                set_gtk_property!(entry_loc_theta, :value, locs[!, :loc_theta][current_channel])
                set_gtk_property!(entry_loc_radius, :value, locs[!, :loc_radius][current_channel])
            else
                set_gtk_property!(entry_loc_theta, :value, -1)
                set_gtk_property!(entry_loc_radius, :value, -1)
            end
        end
    end

    signal_connect(entry_ch, "value-changed") do widget
        current_channel = get_gtk_property(entry_ch, :value, Int64)
        Gtk.@sigatom begin
            set_gtk_property!(entry_label, :text, channel_labels[current_channel])
            set_gtk_property!(combo_chtype, :active, findfirst(isequal(Symbol(channel_types[current_channel])), NeuroAnalyzer.channel_types) - 2)
            if current_channel in signal_channels(obj)
                set_gtk_property!(entry_loc_theta, :sensitive, true)
                set_gtk_property!(entry_loc_radius, :sensitive, true)
            else
                set_gtk_property!(entry_loc_theta, :sensitive, false)
                set_gtk_property!(entry_loc_radius, :sensitive, false)
            end
            if current_channel in locs_ch
                set_gtk_property!(entry_loc_theta, :value, locs[!, :loc_theta][current_channel])
                set_gtk_property!(entry_loc_radius, :value, locs[!, :loc_radius][current_channel])
            else
                set_gtk_property!(entry_loc_theta, :value, -1)
                set_gtk_property!(entry_loc_radius, :value, -1)
            end
        end
        draw(can1)
        draw(can2)
        draw(can3)
        draw(can4)
    end

    return nothing

end