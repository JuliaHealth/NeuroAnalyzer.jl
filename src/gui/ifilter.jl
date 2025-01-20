export ifilter

"""
    ifilter(obj)

Interactive filter design.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function ifilter(obj::NeuroAnalyzer.NEURO)::Nothing

    fprototypes = [:butterworth, :chebyshev1, :chebyshev2, :elliptic, :fir, :firls, :iirnotch, :remez]
    fprototype = fprototypes[1]
    ftypes = [:lp, :hp, :bs, :bp]
    ftype = ftypes[1]
    cutoff = 0.1
    order = 8
    rp = 2
    rs = 20
    bw = 2
    w = nothing

    fs = sr(obj)

    p = plot_filter_response(fs=fs, n=epoch_len(obj), fprototype=fprototype, ftype=ftype, cutoff=cutoff, order=order, rp=rp, rs=rs, bw=bw, w=w)

    win = GtkWindow("NeuroAnalyzer: ifilter()", 1300, 800)
    set_gtk_property!(win, :border_width, 5)
    set_gtk_property!(win, :resizable, true)
    set_gtk_property!(win, :has_resize_grip, false)
    set_gtk_property!(win, :window_position, 3)
    set_gtk_property!(win, :startup_id, "org.neuroanalyzer")
    can = GtkCanvas(Int32(p.attr[:size][1]), Int32(p.attr[:size][2]))

    g = GtkGrid()
    g_opts = GtkGrid()
    set_gtk_property!(g, :column_homogeneous, false)
    set_gtk_property!(g_opts, :column_homogeneous, false)
    set_gtk_property!(g, :column_spacing, 5)
    set_gtk_property!(g, :row_spacing, 5)
    set_gtk_property!(g_opts, :row_spacing, 5)
    set_gtk_property!(g_opts, :column_spacing, 5)

    bt_refresh = GtkButton("Refresh")
    set_gtk_property!(bt_refresh, :tooltip_text, "Refresh the plot")

    bt_close = GtkButton("Close")
    set_gtk_property!(bt_close, :tooltip_text, "Close this window")

    cb_mono = GtkCheckButton()
    set_gtk_property!(cb_mono, :tooltip_text, "Use color or gray palette")

    combo_fprototype = GtkComboBoxText()
    for idx in uppercase.(String.(fprototypes))
        push!(combo_fprototype, idx)
    end
    set_gtk_property!(combo_fprototype, :active, 0)
    set_gtk_property!(combo_fprototype, :tooltip_text, "Filter prototype")

    combo_ftype = GtkComboBoxText()
    for idx in uppercase.(String.(ftypes))
        push!(combo_ftype, idx)
    end
    set_gtk_property!(combo_ftype, :active, 0)
    set_gtk_property!(combo_ftype, :tooltip_text, "Filter type")

    entry_cutoff1 = GtkSpinButton(0.1, fs / 2 - 0.1, 0.1)
    set_gtk_property!(entry_cutoff1, :value, cutoff)
    set_gtk_property!(entry_cutoff1, :tooltip_text, "Filter cutoff in Hz (lower)")

    entry_cutoff2 = GtkSpinButton(0.2, fs / 2, 0.1)
    set_gtk_property!(entry_cutoff2, :value, cutoff)
    set_gtk_property!(entry_cutoff2, :sensitive, false)
    set_gtk_property!(entry_cutoff2, :tooltip_text, "Filter cutoff in Hz (upper)")

    entry_order = GtkSpinButton(2, 256, 2)
    set_gtk_property!(entry_order, :value, order)
    set_gtk_property!(entry_order, :tooltip_text, "Filter order (6 dB/octave) for IIR filters\nNumber of taps for REMEZ filter\nAttenuation (× 4 dB) for FIR filter")

    entry_rp = GtkSpinButton(0.0025, 256.0, 0.0025)
    set_gtk_property!(entry_rp, :value, rp)
    set_gtk_property!(entry_rp, :tooltip_text, "Ripple amplitude in dB in the pass band; default: 0.0025 dB for ELLIPTIC, 2 dB for others")

    entry_rs = GtkSpinButton(0.0025, 256.0, 0.0025)
    set_gtk_property!(entry_rs, :value, rs)
    set_gtk_property!(entry_rs, :tooltip_text, "Ripple amplitude in dB in the stop band; default: 40 dB for ELLIPTIC, 20 dB for others")

    entry_bw = GtkSpinButton(1.0, 256.0, 1.0)
    set_gtk_property!(entry_bw, :value, bw)
    set_gtk_property!(entry_bw, :tooltip_text, "Bandwidth for IIRNOTCH and REMEZ filters")

    lab_fprototype = GtkLabel("Filter prototype:")
    set_gtk_property!(lab_fprototype, :halign, 2)
    lab_ftype = GtkLabel("Filter type:")
    set_gtk_property!(lab_ftype, :halign, 2)
    lab_mono = GtkLabel("Grayscale:")
    set_gtk_property!(lab_mono, :halign, 2)
    lab_cutoff1 = GtkLabel("Cutoff (L):")
    set_gtk_property!(lab_cutoff1, :halign, 2)
    lab_cutoff2 = GtkLabel("Cutoff (U):")
    set_gtk_property!(lab_cutoff2, :halign, 2)
    lab_order = GtkLabel("Order:")
    set_gtk_property!(lab_order, :halign, 2)
    lab_rp = GtkLabel("R_pass:")
    set_gtk_property!(lab_rp, :halign, 2)
    lab_rs = GtkLabel("R_stop:")
    set_gtk_property!(lab_rs, :halign, 2)
    lab_bw = GtkLabel("Bandwidth:")
    set_gtk_property!(lab_bw, :halign, 2)

    g_opts[1, 1] = lab_fprototype
    g_opts[1, 2] = lab_ftype
    g_opts[1, 3] = lab_cutoff1
    g_opts[1, 4] = lab_cutoff2
    g_opts[1, 5] = lab_order
    g_opts[1, 6] = lab_rp
    g_opts[1, 7] = lab_rs
    g_opts[1, 8] = lab_bw
    g_opts[1, 9] = lab_mono
    g_opts[2, 1] = combo_fprototype
    g_opts[2, 2] = combo_ftype
    g_opts[2, 3] = entry_cutoff1
    g_opts[2, 4] = entry_cutoff2
    g_opts[2, 5] = entry_order
    g_opts[2, 6] = entry_rp
    g_opts[2, 7] = entry_rs
    g_opts[2, 8] = entry_bw
    g_opts[2, 9] = cb_mono
    g_opts[1:2, 10] = bt_refresh
    g_opts[1:2, 11] = bt_close
    vbox = GtkBox(:v)
    push!(vbox, g_opts)

    g[1, 1] = vbox
    g[2, 1] = can
    push!(win, g)

    showall(win)

    @guarded draw(can) do widget
        if mod(get_gtk_property(entry_order, :value, Int64), 2) != 0
            warn_dialog("Order must be even!")
            Gtk.@sigatom begin
                set_gtk_property!(entry_order, :value, order)
            end
        else
            ctx = getgc(can)
            mono = get_gtk_property(cb_mono, :active, Bool)
            cutoff1 = get_gtk_property(entry_cutoff1, :value, Float64)
            cutoff2 = get_gtk_property(entry_cutoff2, :value, Float64)
            order = get_gtk_property(entry_order, :value, Int64)
            order = get_gtk_property(entry_order, :value, Int64)
            rp = get_gtk_property(entry_rp, :value, Float64)
            rs = get_gtk_property(entry_rs, :value, Float64)
            bw = get_gtk_property(entry_bw, :value, Float64)
            fprototype = fprototypes[get_gtk_property(combo_fprototype, :active, Int64) + 1]
            ftype = ftypes[get_gtk_property(combo_ftype, :active, Int64) + 1]
            fprototype === :iirnotch && (ftype = :bs)
            if ftype === :bp || ftype === :bs
                if fprototype !== :iirnotch
                    cutoff = (cutoff1, cutoff2)
                else
                    cutoff = cutoff1
                end
            else
                cutoff = cutoff1
            end
            p = plot_filter_response(fs=fs, n=epoch_len(obj), fprototype=fprototype, ftype=ftype, cutoff=cutoff, order=order, rp=rp, rs=rs, bw=bw, w=w, mono=mono)
            show(io, MIME("image/png"), p)
            img = read_from_png(io)
            set_source_surface(ctx, img, 0, 0)
            paint(ctx)
        end
    end

    signal_connect(combo_fprototype, "changed") do widget
        draw(can)
    end

    signal_connect(combo_ftype, "changed") do widget
        ftype = ftypes[get_gtk_property(combo_ftype, :active, Int64) + 1]
        if ftype === :bp || ftype === :bs
            Gtk.@sigatom begin
                set_gtk_property!(entry_cutoff2, :sensitive, true)
            end
        else
            Gtk.@sigatom begin
                set_gtk_property!(entry_cutoff2, :sensitive, false)
            end
        end
        draw(can)
    end

    signal_connect(entry_order, "changed") do widget
        draw(can)
    end

    signal_connect(entry_cutoff1, "changed") do widget
        Gtk.@sigatom begin
            GAccessor.range(entry_cutoff2, round(get_gtk_property(entry_cutoff1, :value, Float64) + 0.1, digits=1), fs / 2)
            set_gtk_property!(entry_cutoff2, :value, round(get_gtk_property(entry_cutoff1, :value, Float64) + 0.1, digits=1))
        end
        draw(can)
    end

    signal_connect(entry_cutoff2, "changed") do widget
        draw(can)
    end

    signal_connect(entry_rp, "changed") do widget
        draw(can)
    end

    signal_connect(entry_rs, "changed") do widget
        draw(can)
    end

    signal_connect(entry_bw, "changed") do widget
        draw(can)
    end

    signal_connect(cb_mono, "clicked") do widget
        draw(can)
    end

    signal_connect(bt_refresh, "clicked") do widget
        draw(can)
    end

    signal_connect(bt_close, "clicked") do widget
        Gtk.destroy(win)
        return nothing
    end

    signal_connect(win, "key-press-event") do widget, event
        k = event.keyval
        s = event.state
        if s == 0x00000004 || s == 0x00000014 # ctrl
            if k == 0x00000071 # q
                Gtk.destroy(win)
                return nothing
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