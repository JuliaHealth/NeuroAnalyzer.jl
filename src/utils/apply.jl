export apply

"""
    apply(obj; <keyword arguments>)

Apply a custom function to selected channels and epochs of a NEURO object.

The function string `f` must reference the signal data using the placeholder variable `obj`, which is bound to the actual channel data (a `Vector{Float64}` representing one channel × one epoch slice). The expression is parsed once and compiled into an anonymous function `obj -> f`.

Example: `f = "cumsum(obj)"` or `f = "obj .^ 2"`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `f::String`: Julia expression to evaluate, using `obj` as the signal placeholder

# Returns

- `Array{Float64, 3}`: result array, shape (channels, length of the result of `f`, epochs)
"""
function apply(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    f::String,
)::Array{Float64, 3}
    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))

    # number of channels
    ch_n = length(ch)
    # number of epochs
    ep_n = nepochs(obj)

    # compile the expression once into `obj -> f`; the signal is passed as an argument
    # instead of being pasted into the source text for every channel and epoch
    fn = try
        eval(Meta.parse("obj -> ($f)"))
    catch err
        throw(ArgumentError("Formula failed. Check expression `f`. Error: $err"))
    end

    # dry run to infer output shape
    # use the first selected channel (not hardcoded index 1) for a representative result
    # this also validates the formula before entering the main loop
    out_tmp = try
        Base.invokelatest(fn, obj.data[ch[1], :, 1])
    catch err
        throw(ArgumentError("Formula failed. Check expression `f`. Error: $err"))
    end

    # pre-allocate output
    result = zeros(eltype(out_tmp), ch_n, length(out_tmp), ep_n)

    # initialize progress bar
    progbar = Progress(
        ch_n * ep_n; dt = 1, barlen = 20, color = :white, enabled = progress_bar,
    )
    # compute over channels and epochs
    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        try
            result[ch_idx, :, ep_idx] = Base.invokelatest(fn, obj.data[ch[ch_idx], :, ep_idx])
        catch err
            throw(ArgumentError("Formula failed. Check expression `f`. Error: $err"))
        end
        # update progress bar
        progress_bar && next!(progbar)
    end
    progress_bar && finish!(progbar)

    return result
end
