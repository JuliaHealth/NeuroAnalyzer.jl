export iedit_ch

"""
    iedit_ch(obj)

Interactive edit signal channels properties and locations.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
"""
function iedit_ch(obj::NeuroAnalyzer.NEURO)

    @assert obj.header.recording[:data_type] == "eeg" "Currently this function only works for EEG data."

    # TO DO: select channel by clicking its location
    # TO DO: generate locations
    # TO DO: other recording types

    obj_new = deepcopy(obj)

    function _refresh_plots()
        Gtk.@sigatom begin
            draw(can1)
            draw(can2)
            draw(can3)
            draw(can4)
        end
    end
    function _refresh_locs()
        Gtk.@sigatom begin
            if current_channel in ch_signal
                set_gtk_property!(entry_loc_theta, :sensitive, true)
                set_gtk_property!(entry_loc_radius, :sensitive, true)
                set_gtk_property!(entry_loc_x, :sensitive, true)
                set_gtk_property!(entry_loc_y, :sensitive, true)
                set_gtk_property!(entry_loc_z, :sensitive, true)
                set_gtk_property!(entry_loc_theta_sph, :sensitive, true)
                set_gtk_property!(entry_loc_radius_sph, :sensitive, true)
                set_gtk_property!(entry_loc_phi_sph, :sensitive, true)
            else
                set_gtk_property!(entry_loc_theta, :sensitive, false)
                set_gtk_property!(entry_loc_radius, :sensitive, false)
                set_gtk_property!(entry_loc_x, :sensitive, false)
                set_gtk_property!(entry_loc_y, :sensitive, false)
                set_gtk_property!(entry_loc_z, :sensitive, false)
                set_gtk_property!(entry_loc_theta_sph, :sensitive, false)
                set_gtk_property!(entry_loc_radius_sph, :sensitive, false)
                set_gtk_property!(entry_loc_phi_sph, :sensitive, false)
            end
            if current_channel in locs_ch
                set_gtk_property!(entry_loc_theta, :value, locs[!, :loc_theta][current_channel])
                set_gtk_property!(entry_loc_radius, :value, locs[!, :loc_radius][current_channel])
                set_gtk_property!(entry_loc_x, :value, locs[!, :loc_x][current_channel])
                set_gtk_property!(entry_loc_y, :value, locs[!, :loc_y][current_channel])
                set_gtk_property!(entry_loc_z, :value, locs[!, :loc_z][current_channel])
                set_gtk_property!(entry_loc_theta_sph, :value, locs[!, :loc_theta_sph][current_channel])
                set_gtk_property!(entry_loc_radius_sph, :value, locs[!, :loc_radius_sph][current_channel])
                set_gtk_property!(entry_loc_phi_sph, :value, locs[!, :loc_phi_sph][current_channel])
            else
                set_gtk_property!(entry_loc_theta, :value, 0)
                set_gtk_property!(entry_loc_radius, :value, 0)
                set_gtk_property!(entry_loc_x, :value, 0)
                set_gtk_property!(entry_loc_y, :value, 0)
                set_gtk_property!(entry_loc_z, :value, 0)
                set_gtk_property!(entry_loc_theta_sph, :value, 0)
                set_gtk_property!(entry_loc_radius_sph, :value, 0)
                set_gtk_property!(entry_loc_phi_sph, :value, 0)
            end
        end
    end

    if nchannels(obj) < 1
        _warn("OBJ must contain ≥ 1 channel.")
        warn_dialog("OBJ must contain ≥ 1 channel.")
        Gtk.destroy(win)
        return nothing
    end

    current_channel = 1
    ch_types = obj_new.header.recording[:channel_type]
    ch_units = obj_new.header.recording[:units]
    ch_labels = obj_new.header.recording[:labels]
    ch_signal = signal_channels(obj_new)

    if _has_locs(obj_new)
        locs = obj_new.locs
        locs_ch = _find_bylabel(obj_new.locs, ch_labels[ch_signal])
    else
        _initialize_locs!(obj_new)
        locs = obj_new.locs
        locs_ch = _find_bylabel(obj_new.locs, ch_labels[ch_signal])
    end

    # Gtk canvas / Cairo context should be scaled only once
    already_scaled1 = false
    already_scaled2 = false
    already_scaled3 = false
    already_scaled4 = false
    refresh = true

    win = GtkWindow("NeuroAnalyzer: iedit_ch()", 1000, 900)
    set_gtk_property!(win, :border_width, 20)
    set_gtk_property!(win, :resizable, false)
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
    entry_ch = GtkSpinButton(1, nchannels(obj_new), 1)
    set_gtk_property!(entry_ch, :value, current_channel)
    set_gtk_property!(entry_ch, :climb_rate, 0.1)
    set_gtk_property!(entry_ch, :tooltip_text, "Channel number")
    bt_start = GtkButton("⇤")
    set_gtk_property!(bt_start, :tooltip_text, "Go to the first channel")
    bt_end = GtkButton("⇥")
    set_gtk_property!(bt_end, :tooltip_text, "Go to the last channel")

    lab_chtype = GtkLabel("Type")
    #set_gtk_property!(lab_chtype, :halign, 2)
    combo_chtype = GtkComboBoxText()
    for idx in NeuroAnalyzer.channel_types[2:end]
        push!(combo_chtype, idx)
    end
    set_gtk_property!(combo_chtype, :active, findfirst(isequal(ch_types[current_channel]), NeuroAnalyzer.channel_types) - 2)

    lab_chunits = GtkLabel("Units")
    #set_gtk_property!(lab_chunits, :halign, 2)
    combo_chunits = GtkComboBoxText()
    for idx in NeuroAnalyzer.channel_units
        push!(combo_chunits, idx)
    end
    set_gtk_property!(combo_chunits, :active, findfirst(isequal(NeuroAnalyzer._ch_units(ch_types[current_channel])), NeuroAnalyzer.channel_units) - 1)

    lab_chlabel = GtkLabel("Label")
    #set_gtk_property!(lab_chlabel, :halign, 2)
    entry_label = GtkEntry()
    set_gtk_property!(entry_label, :text, ch_labels[current_channel])

    bt_delete = GtkButton("Delete channel")
    nchannels(obj_new) == 0 && set_gtk_property!(bt_delete, :sensitive, false)

    lab_loc_x = GtkLabel("Cartesian X")
    #set_gtk_property!(lab_loc_x, :halign, 2)
    entry_loc_x = GtkSpinButton(-1.5, 1.5, 0.01)
    set_gtk_property!(entry_loc_x, :tooltip_text, "Cartesian X coordinate (:loc_x)")
    set_gtk_property!(entry_loc_x, :digits, 2)
    lab_loc_y = GtkLabel("Cartesian Y")
    #set_gtk_property!(lab_loc_y, :halign, 2)
    entry_loc_y = GtkSpinButton(-1.5, 1.5, 0.01)
    set_gtk_property!(entry_loc_y, :tooltip_text, "Cartesian Y coordinate (:loc_y)")
    set_gtk_property!(entry_loc_y, :digits, 2)
    lab_loc_z = GtkLabel("Cartesian Z")
    #set_gtk_property!(lab_loc_z, :halign, 2)
    entry_loc_z = GtkSpinButton(-1.5, 1.5, 0.01)
    set_gtk_property!(entry_loc_z, :tooltip_text, "Cartesian Z coordinate (:loc_z)")
    set_gtk_property!(entry_loc_z, :digits, 2)
    lab_loc_theta_sph = GtkLabel("Spherical theta")
    #set_gtk_property!(lab_loc_theta_sph, :halign, 2)
    entry_loc_theta_sph = GtkSpinButton(-360.0, 360.0, 0.05)
    set_gtk_property!(entry_loc_theta_sph, :tooltip_text, "Spherical horizontal angle (:loc_theta_sph)")
    set_gtk_property!(entry_loc_theta_sph, :digits, 2)
    lab_loc_phi_sph = GtkLabel("Spherical phi")
    #set_gtk_property!(lab_loc_phi_sph, :halign, 2)
    entry_loc_phi_sph = GtkSpinButton(-360.0, 360.0, 0.05)
    set_gtk_property!(entry_loc_phi_sph, :tooltip_text, "Spherical azimuth angle (:loc_phi_sph)")
    set_gtk_property!(entry_loc_phi_sph, :digits, 2)
    lab_loc_radius_sph = GtkLabel("Spherical radius")
    #set_gtk_property!(lab_loc_radius_sph, :halign, 2)
    entry_loc_radius_sph = GtkSpinButton(-1.5, 1.5, 0.01)
    set_gtk_property!(entry_loc_radius_sph, :tooltip_text, "Spherical radius (:loc_radius_sph)")
    set_gtk_property!(entry_loc_radius_sph, :digits, 2)

    lab_loc_theta = GtkLabel("Polar theta")
    #set_gtk_property!(lab_loc_theta, :halign, 2)
    entry_loc_theta = GtkSpinButton(-360.0, 360.0, 0.05)
    set_gtk_property!(entry_loc_theta, :tooltip_text, "Polar theta angle (:loc_theta)")
    set_gtk_property!(entry_loc_theta, :digits, 2)
    lab_loc_radius = GtkLabel("Polar radius")
    #set_gtk_property!(lab_loc_radius, :halign, 2)
    entry_loc_radius = GtkSpinButton(-1.5, 1.5, 0.01)
    set_gtk_property!(entry_loc_radius, :tooltip_text, "Polar radius (:loc_radius)")
    set_gtk_property!(entry_loc_radius, :digits, 2)

    bt_swapxy = GtkButton("Swap XY")
    set_gtk_property!(bt_swapxy, :tooltip_text, "Swap X and Y axes")
    bt_flip = GtkButton("Flip")
    set_gtk_property!(bt_flip, :tooltip_text, "Flip over X/Y/Z axis")
    combo_flip = GtkComboBoxText()
    set_gtk_property!(combo_flip, :tooltip_text, "Axis to flip over")
    axes = ["X", "Y", "Z"]
    for idx in axes
        push!(combo_flip, idx)
    end
    set_gtk_property!(combo_flip, :active, 0)

    bt_ax_rot = GtkButton("Rotate around axis")
    combo_ax_rot = GtkComboBoxText()
    axes = ["X", "Y", "Z"]
    for idx in axes
        push!(combo_ax_rot, idx)
    end
    set_gtk_property!(combo_ax_rot, :active, 0)
    entry_ax_rot_degree = GtkSpinButton(-360, 360, 1.0)
    set_gtk_property!(entry_ax_rot_degree, :digits, 1)
    set_gtk_property!(entry_ax_rot_degree, :tooltip_text, "Rotation angle in degrees\nPositive angle rotates anti-clockwise for X and Z axes and clockwise for Y axis")
    set_gtk_property!(entry_ax_rot_degree, :value, 0)
    bt_scale = GtkButton("Scale")
    entry_scale = GtkSpinButton(0.1, 10.00, 0.1)
    set_gtk_property!(entry_scale, :tooltip_text, "Scaling factor")
    set_gtk_property!(entry_scale, :value, 1.0)
    bt_normalize = GtkButton("Normalize")
    set_gtk_property!(bt_normalize, :tooltip_text, "Normalize channel locations to fit the unit sphere")
    bt_transform = GtkButton("Transform")
    set_gtk_property!(bt_transform, :tooltip_text, "Transform coordinates from one set to another")
    combo_transform = GtkComboBoxText()
    transformations = ["Cartesian → polar", "Cartesian → spherical", "polar → Cartesian", "polar → spherical", "spherical → Cartesian", "spherical → polar"]
    for idx in transformations
        push!(combo_transform, idx)
    end
    set_gtk_property!(combo_transform, :active, 0)

    bt_load = GtkButton("Load")
    set_gtk_property!(bt_load, :tooltip_text, "Load location coordinates")
    bt_save = GtkButton("Save")
    set_gtk_property!(bt_save, :tooltip_text, "Save location coordinates")
    bt_generate = GtkButton("Generate")
    set_gtk_property!(bt_generate, :tooltip_text, "Generate location coordinates using 10-5 system")
    bt_apply = GtkButton("Apply")
    set_gtk_property!(bt_apply, :tooltip_text, "Apply changes and close this window")
    bt_cancel = GtkButton("Cancel")
    set_gtk_property!(bt_cancel, :tooltip_text, "Close this window and abandon changes")

    lab_cart = GtkLabel("Plot using Cartesian coordinates:")
    set_gtk_property!(lab_cart, :halign, 2)
    cb_plot_cart = GtkCheckButton("")
    set_gtk_property!(cb_plot_cart, :active, false)

    lab_hdlab = GtkLabel("Plot head labels:")
    set_gtk_property!(lab_hdlab, :halign, 2)
    cb_hdlab = GtkCheckButton("")
    set_gtk_property!(cb_hdlab, :active, false)

    cb_polar = GtkCheckButton("")
    set_gtk_property!(cb_polar, :halign, 3)
    set_gtk_property!(cb_polar, :active, true)
    set_gtk_property!(cb_polar, :tooltip_text, "Apply operations to polar coordinates")
    cb_cartesian = GtkCheckButton("")
    set_gtk_property!(cb_cartesian, :halign, 3)
    set_gtk_property!(cb_cartesian, :active, true)
    set_gtk_property!(cb_cartesian, :tooltip_text, "Apply operations to Cartesian coordinates")
    cb_spherical = GtkCheckButton("")
    set_gtk_property!(cb_spherical, :halign, 3)
    set_gtk_property!(cb_spherical, :active, true)
    set_gtk_property!(cb_spherical, :tooltip_text, "Apply operations to spherical coordinates")

    g_opts[1:3, 1] = GtkLabel("Channel number")
    g_opts[1, 2] = bt_start
    g_opts[2, 2] = entry_ch
    g_opts[3, 2] = bt_end
    g_opts[2, 3] = bt_delete
    g_opts[1:3, 4] = GtkLabel("Channel properties")
    g_opts[1, 5] = lab_chlabel
    g_opts[1, 6] = entry_label
    g_opts[2, 5] = lab_chtype
    g_opts[2, 6] = combo_chtype
    g_opts[3, 5] = lab_chunits
    g_opts[3, 6] = combo_chunits
    g_opts[1:3, 7] = GtkLabel("Channel coordinates")
    g_opts[1, 8] = lab_loc_radius
    g_opts[1, 9] = entry_loc_radius
    g_opts[2, 8] = lab_loc_theta
    g_opts[2, 9] = entry_loc_theta
    g_opts[1, 10] = lab_loc_x
    g_opts[1, 11] = entry_loc_x
    g_opts[2, 10] = lab_loc_y
    g_opts[2, 11] = entry_loc_y
    g_opts[3, 10] = lab_loc_z
    g_opts[3, 11] = entry_loc_z
    g_opts[1, 12] = lab_loc_radius_sph
    g_opts[1, 13] = entry_loc_radius_sph
    g_opts[2, 12] = lab_loc_theta_sph
    g_opts[2, 13] = entry_loc_theta_sph
    g_opts[3, 12] = lab_loc_phi_sph
    g_opts[3, 13] = entry_loc_phi_sph
    g_opts[1:3, 14] = GtkLabel("Edit locs")
    g_opts[1, 15] = GtkLabel("Apply to polar")
    g_opts[1, 16] = cb_polar
    g_opts[2, 15] = GtkLabel("Apply to Cartesian")
    g_opts[2, 16] = cb_cartesian
    g_opts[3, 15] = GtkLabel("Apply to spherical")
    g_opts[3, 16] = cb_spherical
    g_opts[1, 17] = bt_flip
    g_opts[2, 17] = combo_flip
    g_opts[3, 17] = bt_swapxy
    g_opts[1, 18] = bt_ax_rot
    g_opts[2, 18] = combo_ax_rot
    g_opts[3, 18] = entry_ax_rot_degree
    g_opts[1, 19] = bt_scale
    g_opts[2, 19] = entry_scale
    g_opts[3, 19] = bt_normalize
    g_opts[1, 20] = bt_transform
    g_opts[2, 20] = combo_transform
    g_opts[1:3, 21] = GtkLabel("Locs operations")
    g_opts[1, 22] = bt_generate
    g_opts[2, 22] = bt_load
    g_opts[3, 22] = bt_save
    g_opts[1, 23] = lab_hdlab
    g_opts[2, 23] = cb_hdlab
    g_opts[1, 24] = lab_cart
    g_opts[2, 24] = cb_plot_cart
    g_opts[1, 25] = bt_apply
    g_opts[2, 25] = bt_cancel
    vbox = GtkBox(:v)
    push!(vbox, g_opts)

    g[1, 1:2] = vbox
    g[2, 1] = can1
    g[3, 1] = can2
    g[2, 2] = can3
    g[3, 2] = can4

    push!(win, g)

    showall(win)

    _refresh_locs()

    @guarded draw(can1) do widget
        if refresh
            cart = get_gtk_property(cb_plot_cart, :active, Bool)
            hdlab = get_gtk_property(cb_hdlab, :active, Bool)
            if current_channel in ch_signal
                selected = current_channel
            else
                selected = 0
            end
            p = NeuroAnalyzer.plot_locs(locs, ch=ch_signal, selected=selected, ch_labels=false, head_labels=hdlab, cart=cart, plane=:xy, grid=true)
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
    end

    @guarded draw(can2) do widget
        if refresh
            cart = get_gtk_property(cb_plot_cart, :active, Bool)
            hdlab = get_gtk_property(cb_hdlab, :active, Bool)
            if current_channel in ch_signal
                selected = current_channel
            else
                selected = 0
            end
            p = NeuroAnalyzer.plot_locs(locs, ch=ch_signal, selected=selected, ch_labels=false, head_labels=hdlab, cart=cart, plane=:xz, grid=true)
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
    end

    @guarded draw(can3) do widget
        if refresh
            cart = get_gtk_property(cb_plot_cart, :active, Bool)
            hdlab = get_gtk_property(cb_hdlab, :active, Bool)
            if current_channel in ch_signal
                selected = current_channel
            else
                selected = 0
            end
            p = NeuroAnalyzer.plot_locs(locs, ch=ch_signal, selected=selected, ch_labels=false, head_labels=hdlab, cart=cart, plane=:yz, grid=true)
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
    end

    @guarded draw(can4) do widget
        if refresh
            cart = get_gtk_property(cb_plot_cart, :active, Bool)
            hdlab = get_gtk_property(cb_hdlab, :active, Bool)
            if current_channel in ch_signal
                selected = current_channel
            else
                selected = 0
            end
            p = NeuroAnalyzer.plot_locs3d(locs, ch=ch_signal, selected=selected, ch_labels=false, head_labels=hdlab, cart=cart);
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
    end

    signal_connect(bt_start, "clicked") do widget
        current_channel = 1
        Gtk.@sigatom begin
            set_gtk_property!(entry_ch, :value, current_channel)
        end
    end

    signal_connect(bt_end, "clicked") do widget
        current_channel = nchannels(obj_new)
        Gtk.@sigatom begin
            set_gtk_property!(entry_ch, :value, current_channel)
        end
    end

    signal_connect(entry_ch, "value-changed") do widget
        current_channel = get_gtk_property(entry_ch, :value, Int64)
        Gtk.@sigatom begin
            set_gtk_property!(entry_label, :text, ch_labels[current_channel])
            set_gtk_property!(combo_chtype, :active, findfirst(isequal(ch_types[current_channel]), NeuroAnalyzer.channel_types) - 2)
            set_gtk_property!(combo_chunits, :active, findfirst(isequal(ch_units[current_channel]), NeuroAnalyzer.channel_units) - 1)
        end
        refresh = false
        _refresh_locs()
        refresh = true
        _refresh_plots()
    end

    signal_connect(bt_delete, "clicked") do widget
        if ask_dialog("Delete channel $current_channel ?", "No", "Yes")
            delete_channel!(obj_new, ch=current_channel)
            current_channel > nchannels(obj_new) && (current_channel = nchannels(obj_new))
            ch_types = obj_new.header.recording[:channel_type]
            ch_units = obj_new.header.recording[:units]
            ch_labels = obj_new.header.recording[:labels]
            locs = obj_new.locs
            locs_ch = _find_bylabel(obj_new.locs, ch_labels[ch_signal])
            Gtk.@sigatom begin
                set_gtk_property!(entry_ch, :value, current_channel)
                set_gtk_property!(entry_label, :text, ch_labels[current_channel])
                set_gtk_property!(combo_chtype, :active, findfirst(isequal(ch_types[current_channel]), NeuroAnalyzer.channel_types) - 2)
                set_gtk_property!(combo_chunits, :active, findfirst(isequal(ch_units[current_channel]), NeuroAnalyzer.channel_units) - 1)
            end
            refresh = false
            _refresh_locs()
            refresh = true
            _refresh_plots()
        end
    end 

    signal_connect(entry_label, "changed") do widget
        ch_labels[current_channel] = get_gtk_property(entry_label, :text, String)
    end

    signal_connect(combo_chtype, "changed") do widget
        ch_types[current_channel] = string(NeuroAnalyzer.channel_types[get_gtk_property(combo_chtype, :active, Int64) + 2])
        ch_signal = signal_channels(obj_new)
        Gtk.@sigatom begin
            set_gtk_property!(combo_chunits, :active, findfirst(isequal(ch_units[current_channel]), NeuroAnalyzer.channel_units) - 1)
        end
        refresh = false
        _refresh_locs()
        refresh = true
        _refresh_plots()
    end

    signal_connect(combo_chunits, "changed") do widget
        ch_units[current_channel] = NeuroAnalyzer.channel_units[get_gtk_property(combo_chunits, :active, Int64) + 1]
    end

    signal_connect(cb_plot_cart, "clicked") do widget
        _refresh_plots()
    end

    signal_connect(cb_hdlab, "clicked") do widget
        _refresh_plots()
    end

    signal_connect(entry_loc_radius, "value-changed") do widget
        current_channel in locs_ch && (locs[current_channel, :loc_radius] = get_gtk_property(entry_loc_radius, :value, Float64))
        _refresh_plots()
    end

    signal_connect(entry_loc_theta, "value-changed") do widget
        current_channel in locs_ch && (locs[current_channel, :loc_theta] = get_gtk_property(entry_loc_theta, :value, Float64))
        _refresh_plots()
    end

    signal_connect(entry_loc_x, "value-changed") do widget
        current_channel in locs_ch && (locs[current_channel, :loc_x] = get_gtk_property(entry_loc_x, :value, Float64))
        _refresh_plots()
    end

    signal_connect(entry_loc_y, "value-changed") do widget
        current_channel in locs_ch && (locs[current_channel, :loc_y] = get_gtk_property(entry_loc_y, :value, Float64))
        _refresh_plots()
    end

    signal_connect(entry_loc_z, "value-changed") do widget
        current_channel in locs_ch && (locs[current_channel, :loc_z] = get_gtk_property(entry_loc_z, :value, Float64))
        _refresh_plots()
    end

    signal_connect(entry_loc_radius_sph, "value-changed") do widget
        current_channel in locs_ch && (locs[current_channel, :loc_radius_sph] = get_gtk_property(entry_loc_radius_sph, :value, Float64))
        _refresh_plots()
    end

    signal_connect(entry_loc_theta_sph, "value-changed") do widget
        current_channel in locs_ch && (locs[current_channel, :loc_theta_sph] = get_gtk_property(entry_loc_theta_sph, :value, Float64))
        _refresh_plots()
    end

    signal_connect(entry_loc_phi_sph, "value-changed") do widget
        current_channel in locs_ch && (locs[current_channel, :loc_phi_sph] = get_gtk_property(entry_loc_phi_sph, :value, Float64))
        _refresh_plots()
    end

    signal_connect(bt_flip, "clicked") do widget
        # do not modify "ref" and "eog" channels
        obj_tmp = deepcopy(obj_new)
        delete_channel!(obj_tmp, ch=get_channel_bytype(obj_tmp, type="eog"))
        delete_channel!(obj_tmp, ch=get_channel_bytype(obj_tmp, type="ref"))
        locs_tmp = obj_tmp.locs
        get_gtk_property(combo_flip, :active, Int64) == 0 && locs_flipx!(locs_tmp, polar=get_gtk_property(cb_polar, :active, Bool), cart=get_gtk_property(cb_cartesian, :active, Bool), spherical=get_gtk_property(cb_spherical, :active, Bool))
        get_gtk_property(combo_flip, :active, Int64) == 1 && locs_flipy!(locs_tmp, polar=get_gtk_property(cb_polar, :active, Bool), cart=get_gtk_property(cb_cartesian, :active, Bool), spherical=get_gtk_property(cb_spherical, :active, Bool))
        get_gtk_property(combo_flip, :active, Int64) == 2 && locs_flipz!(locs_tmp, polar=get_gtk_property(cb_polar, :active, Bool), cart=get_gtk_property(cb_cartesian, :active, Bool), spherical=get_gtk_property(cb_spherical, :active, Bool))
        locs[_find_bylabel(locs_tmp, locs_tmp[!, :labels]), :] = locs_tmp
        refresh = false
        _refresh_locs()
        refresh = true
        _refresh_plots()
    end

    signal_connect(bt_ax_rot, "clicked") do widget
        # do not modify "ref" and "eog" channels
        obj_tmp = deepcopy(obj_new)
        delete_channel!(obj_tmp, ch=get_channel_bytype(obj_tmp, type="ref"))
        delete_channel!(obj_tmp, ch=get_channel_bytype(obj_tmp, type="eog"))
        locs_tmp = obj_tmp.locs
        ax = get_gtk_property(combo_ax_rot, :active, Int64)
        if ax == 0
            locs_rotx!(locs_tmp, a=get_gtk_property(entry_ax_rot_degree, :value, Float64), polar=get_gtk_property(cb_polar, :active, Bool), cart=get_gtk_property(cb_cartesian, :active, Bool), spherical=get_gtk_property(cb_spherical, :active, Bool))
        elseif ax == 1
            locs_roty!(locs_tmp, a=get_gtk_property(entry_ax_rot_degree, :value, Float64), polar=get_gtk_property(cb_polar, :active, Bool), cart=get_gtk_property(cb_cartesian, :active, Bool), spherical=get_gtk_property(cb_spherical, :active, Bool))
        elseif ax == 2
            locs_rotz!(locs_tmp, a=get_gtk_property(entry_ax_rot_degree, :value, Float64), polar=get_gtk_property(cb_polar, :active, Bool), cart=get_gtk_property(cb_cartesian, :active, Bool), spherical=get_gtk_property(cb_spherical, :active, Bool))
        end
        locs[_find_bylabel(locs_tmp, locs_tmp[!, :labels]), :] = locs_tmp
        refresh = false
        _refresh_locs()
        refresh = true
        _refresh_plots()
    end

    signal_connect(bt_scale, "clicked") do widget
        # do not modify "ref" and "eog" channels
        obj_tmp = deepcopy(obj_new)
        delete_channel!(obj_tmp, ch=get_channel_bytype(obj_tmp, type="ref"))
        delete_channel!(obj_tmp, ch=get_channel_bytype(obj_tmp, type="eog"))
        locs_tmp = obj_tmp.locs
        locs_scale!(locs_tmp, r=get_gtk_property(entry_scale, :value, Float64), polar=get_gtk_property(cb_polar, :active, Bool), cart=get_gtk_property(cb_cartesian, :active, Bool), spherical=get_gtk_property(cb_spherical, :active, Bool))
        locs[_find_bylabel(locs_tmp, locs_tmp[!, :labels]), :] = locs_tmp
        refresh = false
        _refresh_locs()
        refresh = true
        _refresh_plots()
    end

    signal_connect(bt_normalize, "clicked") do widget
        # do not modify "ref" and "eog" channels
        obj_tmp = deepcopy(obj_new)
        delete_channel!(obj_tmp, ch=get_channel_bytype(obj_tmp, type="ref"))
        delete_channel!(obj_tmp, ch=get_channel_bytype(obj_tmp, type="eog"))
        locs_tmp = obj_tmp.locs
        locs_normalize!(locs_tmp, polar=get_gtk_property(cb_polar, :active, Bool), cart=get_gtk_property(cb_cartesian, :active, Bool), spherical=get_gtk_property(cb_spherical, :active, Bool))
        locs[_find_bylabel(locs_tmp, locs_tmp[!, :labels]), :] = locs_tmp
        refresh = false
        _refresh_locs()
        refresh = true
        _refresh_plots()
    end

    signal_connect(bt_transform, "clicked") do widget
        transform_type = get_gtk_property(combo_transform, :active, Int64)
        transform_type == 0 && locs_cart2pol!(locs)
        transform_type == 1 && locs_cart2sph!(locs)
        transform_type == 2 && locs_pol2cart!(locs)
        transform_type == 3 && locs_pol2sph!(locs)
        transform_type == 4 && locs_sph2cart!(locs)
        transform_type == 5 && locs_sph2pol!(locs)
        refresh = false
        _refresh_locs()
        refresh = true
        _refresh_plots()
    end

    signal_connect(bt_swapxy, "clicked") do widget
        # do not modify "ref" and "eog" channels
        obj_tmp = deepcopy(obj_new)
        delete_channel!(obj_tmp, ch=get_channel_bytype(obj_tmp, type="ref"))
        delete_channel!(obj_tmp, ch=get_channel_bytype(obj_tmp, type="eog"))
        locs_tmp = obj_tmp.locs
        locs_swapxy!(locs_tmp, polar=get_gtk_property(cb_polar, :active, Bool), cart=get_gtk_property(cb_cartesian, :active, Bool), spherical=get_gtk_property(cb_spherical, :active, Bool))
        locs[_find_bylabel(locs_tmp, locs_tmp[!, :labels]), :] = locs_tmp
        refresh = false
        _refresh_locs()
        refresh = true
        _refresh_plots()
    end

    ## TO DO: GENERATE
    signal_connect(bt_generate, "clicked") do widget
        info_dialog("This feature has not been implemented yet.")
    end

    signal_connect(bt_load, "clicked") do widget
        file_name = open_dialog("Pick locations file", GtkNullContainer(), (GtkFileFilter("*.ced, *.elc, *.locs, *.tsv, *.sfp, *.csd, *.geo, *.mat", name="All supported formats"), "*.ced, *.elc, *.locs, *.tsv, *.sfp, *.csd, *.geo, *.mat"))
        if file_name != ""
            if splitext(file_name)[2] in [".ced", ".elc", ".locs", ".tsv", ".sfp", ".csd", ".geo", ".mat"]
                if _has_locs(obj_new) && ask_dialog("Replace channel locations ?", "No", "Yes")
                    load_locs!(obj_new, file_name=file_name)
                    locs = obj_new.locs
                    locs_ch = _find_bylabel(locs, ch_labels[ch_signal])
                    refresh = false
                    _refresh_locs()
                    refresh = true
                    _refresh_plots()
                end
            else
                warn_dialog("Incorrect file name!")
            end
        end
    end

    signal_connect(bt_save, "clicked") do widget
        file_name = save_dialog("Pick locations file", GtkNullContainer(), (GtkFileFilter("*.ced, *.locs, *.tsv", name="All supported formats"), "*.ced, *.locs, *.tsv"))
        if file_name != ""
            if splitext(file_name)[2] in [".ced", ".locs", ".tsv"]
                export_locs(obj_new, file_name=file_name, overwrite=true)
            else
                warn_dialog("Incorrect file name!")
            end
        end
    end

    signal_connect(bt_apply, "clicked") do widget
        if ask_dialog("This operation will apply all changes made to channels and locs.\nPlease confirm.", "No", "Yes")
            obj.header = obj_new.header
            obj.data = obj_new.data
            obj.locs = obj_new.locs
            obj.components = obj_new.components
            obj.history = obj_new.history
            obj.markers = obj_new.markers
            Gtk.destroy(win)
            return nothing
        end
    end

    signal_connect(bt_cancel, "clicked") do widget
        Gtk.destroy(win)
        return nothing
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        s = event.state
        if s == 20
        end
    end

    return nothing

end