export apply

"""
    apply(obj; <keyword arguments>)

Apply a custom function to selected channels and epochs of a NEURO object.

The function string `f` must reference the signal data using the placeholder variable `obj`, which is substituted at runtime with the actual channel data (a `Vector{Float64}` representing one channel × one epoch slice).

Example: `f = "cumsum(obj)"` or `f = "obj .^ 2"`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `f::String`: Julia expression to evaluate, using `obj` as the signal placeholder. 

# Returns

- `out::Array{Float64, 3}`: result array, shape `(channels, epoch length, epochs)`

# Throws

- `ArgumentError`: If the formula `f` produces an error on the dry-run evaluation.
- `ErrorException`: If the formula `f` fails for any channel/epoch combination.
"""
function apply(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    f::String
)::Array{Float64, 3}

    # resolve channel names to indices
    ch = get_channel(obj, ch = ch)

    # number of channels
    ch_n = length(ch)
    # number of epochs
    ep_n = nepochs(obj)

    # dry run to infer output shape
    # use the first selected channel (not hardcoded index 1) for a representative result
    # this also validates the formula before entering the main loop
    f_tmp = replace(f, "obj" => "$(obj.data[ch[1], :, 1])")
    local out_tmp
    try
        out_tmp = eval(Meta.parse(f_tmp))
    catch err
        throw(ArgumentError("Formula failed. Check expression `f`. Error: $err"))
    end

    # pre-allocate output
    out = zeros(eltype(out_tmp), ch_n, out_len, ep_n)

    # initialize progress bar
    progbar = Progress(
        ch_n * ep_n, dt = 1, barlen = 20, color = :white, enabled = progress_bar
    )
    # compute over channels and epochs
    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        f_tmp = replace(f, "obj" => "$(obj.data[ch[ch_idx], :, ep_idx])")
        try
            out[ch_idx, :, ep_idx] = eval(Meta.parse(f_tmp))
        catch err
            throw(ArgumentError("Formula failed. Check expression `f`. Error: $err"))
        end
        # update progress bar
        progress_bar && next!(progbar)
    end

    return out

end
