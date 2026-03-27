_nill() = nothing

"""
    _refresh_ica_can_set(obj_reconstructed, ica_can_set, ic_idx, time1, time2)

Refresh ICA component plots in GUI canvas widgets.

# Arguments

- `obj_reconstructed::Vector{NeuroAnalyzer.NEURO}`: vector of reconstructed NEURO objects for each ICA component
- `ica_can_set::Vector{Gtk4.GtkCanvas}`: vector of GTK4 canvas widgets to display the ICA components
- `ic_idx::Vector{Int64}`: indices of ICA components to refresh
- `time1::Float64`: start time of the segment to plot (in seconds)
- `time2::Float64`: end time of the segment to plot (in seconds)

# Returns

- `Nothing`

# Notes

- This function updates the visualization of ICA components in a GUI interface.
- Each canvas is updated with:
    - The topographical plot of the ICA component
    - The component index label ("IC: X")
- The plots are generated using `plot_topo` with specific parameters for ICA visualization.
"""
function _refresh_ica_can_set(
    obj_reconstructed::Vector{NeuroAnalyzer.NEURO},
    ica_can_set::Vector{Gtk4.GtkCanvas},
    ic_idx::Vector{Int64},
    time1::Float64,
    time2::Float64,
)::Nothing

    # validate
    isempty(obj_reconstructed) &&
        throw(ArgumentError("Reconstructed objects vector cannot be empty"))
    isempty(ica_can_set) && throw(ArgumentError("Canvas vector cannot be empty"))
    isempty(ic_idx) && throw(ArgumentError("IC indices vector cannot be empty"))
    length(obj_reconstructed) == length(ica_can_set) == length(ic_idx) ||
        throw(ArgumentError("Input vectors must have equal lengths"))
    time1 < time2 || throw(ArgumentError("Start time must be less than end time"))
    all(1 .<= ic_idx .<= length(obj_reconstructed)) ||
        throw(ArgumentError("IC indices must be valid"))

    # create Cairo surfaces for each ICA component
    ica_set = Vector{Cairo.CairoSurfaceBase{UInt32}}()

    # denerate plots for each selected ICA component
    for idx in ic_idx
        # create topographical plot for the ICA component
        p_tmp = plot_topo(
            obj_reconstructed[idx];
            ch = datatype(obj_reconstructed[1]),  # use same channel type as first object
            seg = (time1, time2),                 # plot specified time segment
            amethod = :mean,                      # amplitude method
            imethod = :sh,                        # interpolation method
            nmethod = :minmax,                    # normalization method
            cb = false,                           # don't show color bar
            large = false,                        # don't use large plot size
        )
        # convert plot to Cairo surface for GTK display
        cx_tmp = plot2canvas(p_tmp)
        push!(ica_set, cx_tmp)
    end

    # update each canvas with the corresponding ICA component plot
    for (idx, can_idx) in enumerate(ic_idx)
        ica_can = ica_can_set[can_idx]

        # Safely draw on the canvas
        @guarded draw(ica_can) do widget
            ctx_ica = getgc(ica_can)
            Cairo.set_source_surface(ctx_ica, ica_set[idx], 0, 0)
            Cairo.paint(ctx_ica)

            # draw component index label
            Cairo.move_to(ctx_ica, 10.0, 12.0)
            Cairo.set_source_rgb(ctx_ica, 0, 0, 0) # black text
            return Cairo.show_text(ctx_ica, "IC: $can_idx")
        end
    end

    return nothing
end
