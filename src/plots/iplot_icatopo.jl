export iplot_icatopo

"""
    iplot_icatopo(obj; <keyword arguments>)

Interactive topographical plot of embedded ("ic" and "ic_mw") ICA components.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0`: component(s) to plot, default is all components
- `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
- `cb::Bool=false`: plot color bar
- `cb_label::String="[A.U.]"`: color bar label
- `amethod::Symbol=:mean`: averaging method:
    - `:mean`
    - `:median`
- `imethod::Symbol=:sh`: interpolation method:
    - `:sh`: Shepard
    - `:mq`: Multiquadratic
    - `:imq`: InverseMultiquadratic
    - `:tp`: ThinPlate
    - `:nn`: NearestNeighbour
    - `:ga`: Gaussian
- `nmethod::Symbol=:minmax`: method for normalization, see `normalize()`
- `plot_contours::Bools=true`: plot contours over topo plot
- `plot_electrodes::Bools=true`: plot electrodes over topo plot
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function iplot_icatopo(obj::NeuroAnalyzer.NEURO; ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0, seg::Tuple{Real, Real}=(0, 10), cb::Bool=false, cb_label::String="default", amethod::Symbol=:mean, imethod::Symbol=:sh, nmethod::Symbol=:minmax, plot_contours::Bool=true, plot_electrodes::Bool=true, kwargs...)

    _wip()

    @assert :ic in keys(obj.components) "OBJ does not contain :ic component. Perform ica_decompose() first."
    @assert :ic_mw in keys(obj.components) "OBJ does not contain :ic_mw component. Perform ica_decompose() first."
    
    ic = obj.components[:ic]
    ic_mw = obj.components[:ic_mw]

    # select component channels, default is all channels
    ic_idx == 0 && (ic_idx = _select_cidx(ic, ic_idx))
    _check_cidx(ic, ic_idx)

    p_topo = Vector{Plots.Plot{Plots.GRBackend}}()

    for idx in eachindex(ic_idx)
        obj_tmp = ica_reconstruct(obj, ic, ic_mw, ch=signal_channels(obj), ic_idx=ic_idx[idx], keep=true)
        p_tmp = plot_topo(obj_tmp, title="IC $(ic_idx[idx])", cb=cb, cb_label=cb_label, amethod=amethod, imethod=imethod, nmethod=nmethod, plot_contours=plot_contours, plot_electrodes=plot_electrodes, seg=seg, kwargs...)
        push!(p_topo, p_tmp)
    end
    
    if length(ic_idx) <= 4
        p = plot_compose(p_topo, layout=(1, 4))
    elseif length(ic_idx) <= 8
        p = plot_compose(p_topo, layout=(2, ceil(Int64, length(ic_idx) / 2)))
    elseif length(ic_idx) <= 12
        p = plot_compose(p_topo, layout=(3, ceil(Int64, length(ic_idx) / 3)))
    else
        p = plot_compose(p_topo, layout=(4, ceil(Int64, length(ic_idx) / 4)))
    end

    return p

end

"""
    iplot_icatopo(obj, ic, ic_mw; <keyword arguments>)

Interactive topographical plot of external ICA components.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ic::Matrix{Float64}`: components IC(1)..IC(n)
- `ic_mw::Matrix{Float64}`: weighting matrix IC(1)..IC(n)
- `ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0`: component(s) to plot, default is all components
- `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
- `cb::Bool=false`: plot color bar
- `cb_label::String="[A.U.]"`: color bar label
- `amethod::Symbol=:mean`: averaging method:
    - `:mean`
    - `:median`
- `imethod::Symbol=:sh`: interpolation method:
    - `:sh`: Shepard
    - `:mq`: Multiquadratic
    - `:imq`: InverseMultiquadratic
    - `:tp`: ThinPlate
    - `:nn`: NearestNeighbour
    - `:ga`: Gaussian
- `nmethod::Symbol=:minmax`: method for normalization, see `normalize()`
- `plot_contours::Bools=true`: plot contours over topo plot
- `plot_electrodes::Bools=true`: plot electrodes over topo plot
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function iplot_icatopo(obj::NeuroAnalyzer.NEURO, ic::Matrix{Float64}, ic_mw::Matrix{Float64}; ic_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0, seg::Tuple{Real, Real}=(0, 10), cb::Bool=false, cb_label::String="default", amethod::Symbol=:mean, imethod::Symbol=:sh, nmethod::Symbol=:minmax, plot_contours::Bool=true, plot_electrodes::Bool=true, kwargs...)

    _wip()

    # select component channels, default is all channels
    ic_idx == 0 && (ic_idx = _select_cidx(ic, ic_idx))
    _check_cidx(ic, ic_idx)

    p_topo = Vector{Plots.Plot{Plots.GRBackend}}()

    for idx in eachindex(ic_idx)
        obj_tmp = ica_reconstruct(obj, ic, ic_mw, ch=signal_channels(obj), ic_idx=ic_idx[idx], keep=true)
        p_tmp = plot_topo(obj_tmp, title="IC $(ic_idx[idx])", cb=cb, cb_label=cb_label, amethod=amethod, imethod=imethod, nmethod=nmethod, plot_contours=plot_contours, plot_electrodes=plot_electrodes, seg=seg, kwargs...)
        push!(p_topo, p_tmp)
    end

    if length(ic_idx) <= 4
        p = plot_compose(p_topo, layout=(1, 4))
    elseif length(ic_idx) <= 8
        p = plot_compose(p_topo, layout=(2, ceil(Int64, length(ic_idx) / 2)))
    elseif length(ic_idx) <= 12
        p = plot_compose(p_topo, layout=(3, ceil(Int64, length(ic_idx) / 3)))
    else
        p = plot_compose(p_topo, layout=(4, ceil(Int64, length(ic_idx) / 4)))
    end

    return p

end