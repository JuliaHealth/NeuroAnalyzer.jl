function _refresh_ica_can_set(obj_reconstructed::Vector{NeuroAnalyzer.NEURO}, ica_can_set::Vector{Gtk.GtkCanvas}, ic_idx::Vector{Int64}, time1::Float64, time2::Float64)::Nothing
    ica_set = Vector{Cairo.CairoSurfaceBase{UInt32}}()
    for idx in ic_idx
        p_tmp = plot_topo(obj_reconstructed[idx], ch=datatype(obj_reconstructed[1]), seg=(time1, time2), amethod=:mean, imethod=:sh, nmethod=:minmax, cb=false, large=false)
        cx_tmp = plot2canvas(p_tmp)
        push!(ica_set, cx_tmp)
    end
    for idx in ic_idx
        ica_can = ica_can_set[idx]
        @guarded draw(ica_can) do widget
            ctx_ica = getgc(ica_can)
            Cairo.set_source_surface(ctx_ica, ica_set[idx], 0, 0)
            Cairo.paint(ctx_ica)
            Cairo.move_to(ctx_ica, 10.0, 12.0)
            Cairo.set_source_rgb(ctx_ica, 0, 0, 0)
            Cairo.show_text(ctx_ica, "IC: $idx")
        end
    end
    return nothing
end

function _refresh_plots()::Nothing
    Gtk.@sigatom begin
        draw(can1)
        draw(can2)
        draw(can3)
        draw(can4)
    end
    return nothing
end

function _refresh_locs()::Nothing
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
            if isa(NeuroAnalyzer._find_bylabel(locs, ch_labels[current_channel]), Int64)
                set_gtk_property!(entry_loc_theta, :value, locs[_find_bylabel(locs, ch_labels[current_channel]), :loc_theta])
                set_gtk_property!(entry_loc_radius, :value, locs[_find_bylabel(locs, ch_labels[current_channel]), :loc_radius])
                set_gtk_property!(entry_loc_x, :value, locs[_find_bylabel(locs, ch_labels[current_channel]), :loc_x])
                set_gtk_property!(entry_loc_y, :value, locs[_find_bylabel(locs, ch_labels[current_channel]), :loc_y])
                set_gtk_property!(entry_loc_z, :value, locs[_find_bylabel(locs, ch_labels[current_channel]), :loc_z])
                set_gtk_property!(entry_loc_theta_sph, :value, locs[_find_bylabel(locs, ch_labels[current_channel]), :loc_theta_sph])
                set_gtk_property!(entry_loc_radius_sph, :value, locs[_find_bylabel(locs, ch_labels[current_channel]), :loc_radius_sph])
                set_gtk_property!(entry_loc_phi_sph, :value, locs[_find_bylabel(locs, ch_labels[current_channel]), :loc_phi_sph])
            end
        else
            set_gtk_property!(entry_loc_theta, :sensitive, false)
            set_gtk_property!(entry_loc_radius, :sensitive, false)
            set_gtk_property!(entry_loc_x, :sensitive, false)
            set_gtk_property!(entry_loc_y, :sensitive, false)
            set_gtk_property!(entry_loc_z, :sensitive, false)
            set_gtk_property!(entry_loc_theta_sph, :sensitive, false)
            set_gtk_property!(entry_loc_radius_sph, :sensitive, false)
            set_gtk_property!(entry_loc_phi_sph, :sensitive, false)
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
    return nothing
end
