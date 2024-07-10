export vartest

"""
    vartest(obj; <keyword arguments>)

Calculate variance F-test.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

# Returns

Named tuple containing:
- `f::Array{Float64, 3}`
- `p::Array{Float64, 3}`
"""
function vartest(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj))

    _check_channels(obj, ch)
    isa(ch, Int64) && (ch = [ch])

    ch_n = length(ch)
    ep_n = nepochs(obj)

    f = zeros(ch_n, ch_n, ep_n)
    p = zeros(ch_n, ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
       Threads.@threads for ch_idx1 in 1:ch_n
            # create half of the matrix
            for ch_idx2 in 1:ch_idx1
                ftest = @views VarianceFTest(obj.data[ch[ch_idx1], :, ep_idx], obj.data[ch[ch_idx2], :, ep_idx])
                f[ch_idx1, ch_idx2, ep_idx] = ftest.F
                p[ch_idx1, ch_idx2, ep_idx] = pvalue(ftest)
            end
        end
    end

    # copy to the other half
    f = _copy_lt2ut(f)
    p = _copy_lt2ut(p)

    return (f=f, p=p)

end

"""
    vartest(obj1, obj2; <keyword arguments>)

Calculate variance F-test.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `ch1::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj1)`: index of channels, default is all signal channels
- `ch2::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj2)`: index of channels, default is all signal channels
- `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1))`: default use all epochs
- `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2))`: default use all epochs

# Returns

Named tuple containing:
- `f::Array{Float64, 3}`
- `p::Array{Float64, 3}`
"""
function vartest(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch1::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj1), ch2::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj2), ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1)), ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2)))

    _check_channels(obj1, ch1)
    _check_channels(obj2, ch2)
    @assert length(ch1) == length(ch2) "ch1 and ch2 must have the same length."

    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    @assert length(ep1) == length(ep2) "ep1 and ep2 must have the same length."
    @assert epoch_len(obj1) == epoch_len(obj2) "OBJ1 and OBJ2 must have the same epoch lengths."

    isa(ch1, Int64) && (ch1 = [ch1])
    isa(ch2, Int64) && (ch2 = [ch2])
    length(ep1) == 1 && (ep1 = [ep1])
    length(ep2) == 1 && (ep2 = [ep2])

    ep_n = length(ep1)
    ch_n = length(ch1)

    f = zeros(ch_n, ch_n, ep_n)
    p = zeros(ch_n, ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
       Threads.@threads for ch_idx1 in 1:ch_n
            for ch_idx2 in 1:ch_n
                ftest = @views VarianceFTest(obj1.data[ch1[ch_idx1], :, ep1[ep_idx]], obj2.data[ch2[ch_idx2], :, ep2[ep_idx]])
                f[ch_idx1, ch_idx2, ep_idx] = ftest.F
                p[ch_idx1, ch_idx2, ep_idx] = pvalue(ftest)
            end
        end
    end

    return (f=f, p=p)

end
