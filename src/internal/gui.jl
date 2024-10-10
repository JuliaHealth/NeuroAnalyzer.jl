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
