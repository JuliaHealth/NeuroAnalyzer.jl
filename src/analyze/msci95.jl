export msci95

"""
    msci95(obj; <keyword arguments>)

Calculate mean, standard deviation and 95% confidence interval.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `n::Int64=3`: number of bootstraps
- `method::Symbol=:normal`: use normal (`:normal`) method or `n`-times bootstrapping (`:boot`)

# Returns

Named tuple containing:
- `sm::Matrix{Float64}`: mean
- `ss::Matrix{Float64}`: standard deviation
- `su::Matrix{Float64}`: upper 95% CI
- `sl::Matrix{Float64}`: lower 95% CI
"""
function msci95(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, n::Int64=3, method::Symbol=:normal)::@NamedTuple{sm::Matrix{Float64}, ss::Matrix{Float64}, su::Matrix{Float64}, sl::Matrix{Float64}}

    ch = exclude_bads ? get_channel(obj, ch=ch, exclude="bad") : get_channel(obj, ch=ch, exclude="")
    sm, ss, su, sl = @views NeuroAnalyzer.msci95(obj.data[ch, :, :], n=n, method=method)

    return (sm=sm, ss=ss, su=su, sl=sl)

end

"""
    msci95(obj1, obj2; <keyword arguments>)

Calculate mean difference, standard deviation and 95% confidence interval.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2:NeuroAnalyzer.NEURO`
- `ch1::Union{String, Vector{String}}`: list of channels
- `ch2::Union{String, Vector{String}}`: list of channels
- `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1))`: default use all epochs
- `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2))`: default use all epochs

# Returns

Named tuple containing:
- `sm::Matrix{Float64}`: mean
- `ss::Matrix{Float64}`: standard deviation
- `su::Matrix{Float64}`: upper 95% CI bound
- `sl::Matrix{Float64}`: lower 95% CI bound
"""
function msci95(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch1::Union{String, Vector{String}}, ch2::Union{String, Vector{String}}, ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1)), ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2)))::@NamedTuple{sm::Matrix{Float64}, ss::Matrix{Float64}, su::Matrix{Float64}, sl::Matrix{Float64}}

    # check channels
    ch1 = exclude_bads ? get_channel(obj1, ch=ch1, exclude="bad") : get_channel(obj1, ch=ch1, exclude="")
    ch2 = exclude_bads ? get_channel(obj2, ch=ch2, exclude="bad") : get_channel(obj2, ch=ch2, exclude="")
    @assert length(ch1) == length(ch2) "ch1 and ch2 must have the same length."

    # check epochs
    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    @assert length(ep1) == length(ep2) "ep1 and ep2 must have the same length."
    @assert epoch_len(obj1) == epoch_len(obj2) "OBJ1 and OBJ2 must have the same epoch lengths."

    isa(ep1, Int64) && (ep1 = [ep1])
    isa(ep2, Int64) && (ep2 = [ep2])

    sm, ss, su, sl = @views NeuroAnalyzer.msci95(obj1.data[ch1, :, ep1], obj2.data[ch2, :, ep2])

    return (sm=sm, ss=ss, su=su, sl=sl)

end
