export plinterpolate_channel
export plinterpolate_channel!
export plinterpolate

# sample range (within an epoch) of the time segment seg, in seconds of the epoch time
function _epoch_seg(obj::NeuroAnalyzer.NEURO, seg::Tuple{Real, Real})::UnitRange{Int64}
    seg[1] < seg[2] || throw(ArgumentError("seg[1] must be < seg[2]."))
    et = obj.epoch_time
    et[1] <= seg[1] && seg[2] <= et[end] ||
        throw(ArgumentError("seg must be within the epoch time [$(et[1]), $(et[end])] s."))
    return vsearch(seg[1], et):vsearch(seg[2], et)
end

"""
    plinterpolate_channel(obj; <keyword arguments>)

Interpolate channel using planar interpolation.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: channel to interpolate
- `ep::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}}=_c(nepochs(obj))`: epoch number(s) within to interpolate; default: all epochs
- `seg::Union{Nothing, Tuple{Real, Real}}=nothing`: time segment (from, to) in seconds of the epoch time to interpolate, e.g. a short artifact in a continuous recording; default: the whole epoch(s)
- `imethod::Symbol=:sh`: interpolation method:
    - `:sh`: Shepard
    - `:mq`: Multiquadratic
    - `:imq`: Inverse Multiquadratic
    - `:tp`: ThinPlate
    - `:nn`: Nearest Neighbour
    - `:ga`: Gaussian

The channel is interpolated from the other good signal channels (bad channels are not used), whose locations must be available; the interpolated value is computed exactly at the location of `ch`.
- `reset_bads::Bool=false`: if `true`, unmark the interpolated channel(s) as bad (`obj.header.recording[:bad_channel]`)

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function plinterpolate_channel(
    obj::NeuroAnalyzer.NEURO;
    ch::String,
    ep::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}} = _c(nepochs(obj)),
    seg::Union{Nothing, Tuple{Real, Real}} = nothing,
    imethod::Symbol = :sh,
    reset_bads::Bool = false,
)::NeuroAnalyzer.NEURO
    # validate
    _check_var(imethod, [:sh, :mq, :imq, :tp, :nn, :ga], "imethod")
    _has_locs(obj)
    _check_epochs(obj, ep)
    ep = _n2v(ep)
    # samples to interpolate within each epoch
    t_rng = isnothing(seg) ? (1:epoch_len(obj)) : _epoch_seg(obj, seg)

    # signal channels (names); the reference channels exclude ch and bad channels
    signal = get_channel(obj; type = datatype(obj))
    ch in signal ||
        throw(ArgumentError("ch must be a signal channel; cannot interpolate non-signal channels."))
    bad = labels(obj)[obj.header.recording[:bad_channel]]
    ref = setdiff(signal, [ch], bad)
    length(ref) >= 3 || throw(ArgumentError("At least 3 good signal channels are needed for interpolation."))

    # channel numbers and locations (matched by label)
    ch_num(l) = findfirst(==(l), labels(obj))
    loc_labels = lowercase.(String.(obj.locs[!, :label]))
    loc_idx = [findfirst(==(lowercase(l)), loc_labels) for l in vcat(ref, ch)]
    any(isnothing, loc_idx) &&
        throw(ArgumentError("Locations of all reference channels and of ch are required."))
    loc_x = Float64.(obj.locs[Vector{Int64}(loc_idx), :loc_x])
    loc_y = Float64.(obj.locs[Vector{Int64}(loc_idx), :loc_y])
    ref_num = ch_num.(ref)
    ch_idx = ch_num(ch)
    # target: the location of ch
    target = [loc_x[end], loc_y[end]]
    loc_x, loc_y = loc_x[1:(end - 1)], loc_y[1:(end - 1)]

    # create new dataset
    obj_tmp = deepcopy(obj)

    # initialize progress bar
    progbar = Progress(
        length(ep) * length(t_rng); dt = 1, barlen = 20, color = :white, enabled = progress_bar,
    )

    # value at the location of ch, interpolated from the reference channels at each sample
    @inbounds Threads.@threads :static for idx in CartesianIndices((length(t_rng), length(ep)))
        t_idx, ep_idx = t_rng[idx[1]], idx[2]
        itp = _interpolator(imethod, loc_x, loc_y, @view(obj.data[ref_num, t_idx, ep[ep_idx]]))
        obj_tmp.data[ch_idx, t_idx, ep[ep_idx]] = ScatteredInterpolation.evaluate(itp, target)[1]
        progress_bar && next!(progbar)
    end
    progress_bar && finish!(progbar)

    reset_bads && (obj_tmp.header.recording[:bad_channel][ch_idx] = false)

    push!(
        obj_tmp.history,
        "plinterpolate_channel(obj; ch=$ch, ep=$ep, seg=$seg, imethod=$imethod, reset_bads=$reset_bads)",
    )

    return obj_tmp
end

"""
    plinterpolate_channel!(obj; <keyword arguments>)

Interpolate channel using planar interpolation.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: channel to interpolate
- `ep::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}}=_c(nepochs(obj))`: epoch number(s) within to interpolate; default: all epochs
- `seg::Union{Nothing, Tuple{Real, Real}}=nothing`: time segment (from, to) in seconds of the epoch time to interpolate; default: the whole epoch(s)
- `imethod::Symbol=:sh`: interpolation method Shepard (`:sh`), Multiquadratic (`:mq`), InverseMultiquadratic (`:imq`), ThinPlate (`:tp`), NearestNeighbour (`:nn`), Gaussian (`:ga`)
- `reset_bads::Bool=false`: if `true`, unmark the interpolated channel(s) as bad (`obj.header.recording[:bad_channel]`)

# Returns

- `Nothing`
"""
function plinterpolate_channel!(
    obj::NeuroAnalyzer.NEURO;
    ch::String,
    ep::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}} = _c(nepochs(obj)),
    seg::Union{Nothing, Tuple{Real, Real}} = nothing,
    imethod::Symbol = :sh,
    reset_bads::Bool = false,
)::Nothing
    obj_tmp =
        plinterpolate_channel(
            obj; ch = ch, ep = ep, seg = seg, imethod = imethod, reset_bads = reset_bads,
        )
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    obj.history = obj_tmp.history
    obj_tmp = nothing

    return nothing
end

"""
    plinterpolate(s; <keyword arguments>)

Interpolate a 2-D map (planar interpolation) of channel values, leaving out channel `ch`: the map shows the values that the other channels predict at the location of `ch`.

# Arguments

- `s::AbstractVector`: one value per location (row) of `locs`
- `locs::DataFrame`: columns: channel, labels, loc_radius, loc_theta, loc_x, loc_y, loc_z, loc_radius_sph, loc_theta_sph, loc_phi_sph
- `ch::Int64`: location (row of `locs`) to leave out
- `imethod::Symbol=:sh`: interpolation method:
    - `:sh`: Shepard
    - `:mq`: Multiquadratic
    - `:imq`: Inverse Multiquadratic
    - `:tp`: ThinPlate
    - `:nn`: Nearest Neighbour
    - `:ga`: Gaussian
- `nmethod::Symbol=:minmax`: method for normalization, see `normalize()`
- `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use polar coordinates
- `ifactor::Int64=100`: interpolation quality (grid size)

# Returns

Named tuple:

- `int_s::Matrix{Float64}`: interpolated map
- `int_x::Vector{Float64}`: X-axis coordinates
- `int_y::Vector{Float64}`: Y-axis coordinates
"""
function plinterpolate(
    s::AbstractVector;
    locs::DataFrame,
    ch::Int64,
    imethod::Symbol = :sh,
    nmethod::Symbol = :minmax,
    cart::Bool = false,
    ifactor::Int64 = 100,
)::@NamedTuple{
    int_s::Matrix{Float64},
    int_x::Vector{Float64},
    int_y::Vector{Float64},
}
    # validate
    length(s) == DataFrames.nrow(locs) ||
        throw(ArgumentError("Length of s must equal the number of locations."))
    ch in eachindex(s) || throw(ArgumentError("ch must be in [1, $(length(s))]."))
    _check_var(imethod, [:sh, :mq, :imq, :tp, :nn, :ga], "imethod")

    keep = setdiff(eachindex(s), ch)
    if cart
        loc_x = Float64.(locs[keep, :loc_x])
        loc_y = Float64.(locs[keep, :loc_y])
    else
        xy = [pol2cart(locs[idx, :loc_radius], locs[idx, :loc_theta]) for idx in keep]
        loc_x = first.(xy)
        loc_y = last.(xy)
    end

    int_s, int_x, int_y = _interpolate2d(Float64.(s[keep]), loc_x, loc_y, ifactor, imethod, nmethod)

    return (; int_s, int_x, int_y)
end
