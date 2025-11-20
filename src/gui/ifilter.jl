export ifilter

"""
    ifilter(obj)

Interactive filter design.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `flt::Union{Nothing, Vector{Float64}, ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64}, Biquad{:z, Float64}}`

# Notes

The returned filter is based on sampling rate and epoch length of the OBJ used for designing the filter. Therefore, it should not be applied for objects of different sampling rate or epoch length.
"""
function ifilter(obj::NeuroAnalyzer.NEURO)::Union{Nothing, Vector{Float64}, ZeroPoleGain{:z, ComplexF64, ComplexF64, Float64}, Biquad{:z, Float64}}

    fprototypes = [:fir, :firls, :remez, :butterworth, :chebyshev1, :chebyshev2, :elliptic, :iirnotch]
    fprototype = fprototypes[1]
    ftypes = [:lp, :hp, :bs, :bp]
    ftype = ftypes[1]
    cutoff = 10
    order = 8
    rp = 2
    rs = 20
    bw = 1
    w = nothing

    fs = sr(obj)

    p = plot_filter_response(fs=fs, fprototype=fprototype, ftype=ftype, cutoff=cutoff, order=order, rp=rp, rs=rs, bw=bw, w=w)

    function activate(app)

        win = GtkApplicationWindow(app, "NeuroAnalyzer: ifilter()")
        set_gtk_property!(win, :startup_id, "org.neuroanalyzer")
        win.width_request = p.attr[:size][1] + 2
        win.height_request = p.attr[:size][2] + 2

        can = GtkCanvas(Int32(p.attr[:size][1]), Int32(p.attr[:size][2]))

        g = GtkGrid()
        g.column_homogeneous = false
        g.column_spacing = 5
        g.row_spacing = 5
        g_opts = GtkGrid()
        g_opts.column_homogeneous = false
        g_opts.column_spacing = 5
        g_opts.row_spacing = 5

        bt_refresh = GtkButton("Refresh")
        set_gtk_property!(bt_refresh, :tooltip_text, "Refresh the plot")

        bt_close = GtkButton("Close")
        set_gtk_property!(bt_close, :tooltip_text, "Close this window and return the filter")

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
        set_gtk_property!(entry_cutoff2, :value, cutoff + 0.1)
        set_gtk_property!(entry_cutoff2, :sensitive, false)
        set_gtk_property!(entry_cutoff2, :tooltip_text, "Filter cutoff in Hz (upper)")

        entry_order = GtkSpinButton(1, 10_000, 1)
        set_gtk_property!(entry_order, :value, order)
        set_gtk_property!(entry_order, :tooltip_text, "Filter order (number of taps)")

        entry_rp = GtkSpinButton(0.0025, 256.0, 0.0025)
        set_gtk_property!(entry_rp, :value, rp)
        set_gtk_property!(entry_rp, :tooltip_text, "Maximum ripple amplitude in dB; default: 0.0025 dB for ELLIPTIC, 2 dB for others")

        entry_rs = GtkSpinButton(0.0025, 256.0, 0.0025)
        set_gtk_property!(entry_rs, :value, rs)
        set_gtk_property!(entry_rs, :tooltip_text, "Minimum ripple attenuation in dB; default: 40 dB for ELLIPTIC, 20 dB for others")

        entry_bw = GtkSpinButton(0.1, cutoff, 0.1)
        set_gtk_property!(entry_bw, :value, bw)
        set_gtk_property!(entry_bw, :tooltip_text, "Transition band width in Hz")

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

        Gtk4.show(win)

        @guarded draw(can) do widget
            ctx = getgc(can)
            mono = cb_mono.active
            cutoff1 = entry_cutoff1.value
            cutoff2 = entry_cutoff2.value
            order = round(Int64, entry_order.value)
            rp = entry_rp.value
            rs = entry_rs.value
            bw = entry_bw.value
            fprototype = fprototypes[combo_fprototype.active + 1]
            ftype = ftypes[combo_ftype.active + 1]
            fprototype === :iirnotch && (ftype = :bs)
            fprototype = fprototypes[combo_fprototype.active + 1]
            ftype = ftypes[combo_ftype.active + 1]
            if ftype === :bp || ftype === :bs
                if fprototype !== :iirnotch
                    cutoff = (cutoff1, cutoff2)
                else
                    cutoff = cutoff1
                end
            else
                cutoff = cutoff1
            end
            if fprototype === :fir && ftype in [:hp, :bp, :bs] && mod(order, 2) == 0
                _warn("order must be odd. Filter was not generated.")
            else
                p = plot_filter_response(fs=fs, fprototype=fprototype, ftype=ftype, cutoff=cutoff, order=order, rp=rp, rs=rs, bw=bw, w=w, mono=mono)
                show(io, MIME("image/png"), p)
                img = read_from_png(io)
                set_source_surface(ctx, img, 0, 0)
                paint(ctx)
            end
        end

        signal_connect(combo_fprototype, "changed") do widget
            fprototype = fprototypes[combo_fprototype.active + 1]
            ftype = ftypes[combo_ftype.active + 1]
            if fprototype === :iirnotch
                combo_ftype.sensitive = false
                entry_order.sensitive = false
                entry_cutoff2.sensitive = false
            elseif fprototype !== :iirnotch
                combo_ftype.sensitive = true
                entry_order.sensitive = true
                if ftype === :bp || ftype === :bs
                    entry_cutoff2.sensitive = true
                end
            elseif fprototype === :fir
                entry_bw.sensitive = false
            elseif fprototype !== :fir
                entry_bw.sensitive = false
            elseif fprototype in [:fir, :firls, :remez, :iirnotch]
                entry_rp.sensitive = false
                entry_rs.sensitive = false
            elseif fprototype in [:butterworth, :chebyshev1, :chebyshev2, :elliptic]
                entry_rp.sensitive = true
                entry_rs.sensitive = true
            end
            draw(can)
        end

        signal_connect(combo_ftype, "changed") do widget
            ftype = ftypes[combo_ftype.active + 1]
            if ftype in [:bp, :bs]
                entry_cutoff2.sensitive = true
            else
                entry_cutoff2.sensitive = false
            end
            draw(can)
        end

        signal_connect(entry_order, "changed") do widget
            draw(can)
        end

        signal_connect(entry_cutoff1, "changed") do widget
            bw_adj = GtkAdjustment(entry_bw)
            bw_adj.lower = 0.1
            bw_adj.upper = entry_cutoff1.value - 0.1
            @idle_add Gtk4.adjustment(entry_bw, bw_adj)

            cutoff2_adj = GtkAdjustment(entry_cutoff2)
            cutoff2_adj.lower = round(entry_cutoff1.value + 0.1, digits=1)
            cutoff2_adj.upper = fs / 2
            @idle_add Gtk4.adjustment(entry_cutoff2, cutoff2_adj)

            @idle_add entry_cutoff2.value = round(entry_cutoff1.value + 0.1, digits=1)

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

        signal_connect(cb_mono, "toggled") do widget
            draw(can)
        end

        signal_connect(bt_refresh, "clicked") do widget
            draw(can)
        end

        signal_connect(bt_close, "clicked") do widget
            close(win)
        end

        help = "Keyboard shortcuts:\n\nCtrl + r\t\t\tRefresh\n\nCtrl + h\t\t\tThis info\nCtrl + q\t\t\tClose\n"

        win_key = Gtk4.GtkEventControllerKey(win)

        signal_connect(win_key, "key-pressed") do widget, keyval, keycode, state
            k = keyval
            if ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('q'))
                close(win)
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('h'))
                info_dialog(help, win) do
                    nothing
                end
            elseif ((ModifierType(state & Gtk4.MODIFIER_MASK) & mask_ctrl == mask_ctrl) && keyval == UInt('r'))
                draw(can)
            end
        end
    end

    app = GtkApplication()

    Gtk4.signal_connect(activate, app, :activate)

    Gtk4.GLib.stop_main_loop()

    run(app)

    if fprototype === :fir && ftype in [:hp, :bp, :bs] && mod(order, 2) == 0
        _warn("order must be odd. Filter was not generated.")
        return nothing
    else
        return filter_create(; fprototype=fprototype, ftype=ftype, cutoff=cutoff, fs=fs, order=order, rp=rp, rs=rs, bw=bw)
    end

end