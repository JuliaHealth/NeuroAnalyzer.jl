export apply

"""
    apply(obj; channel, f)

Apply custom function.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `f::String`: function to be applied, e.g. `f="mean(obj, dims=3)"; OBJ signal is given using variable `obj` here.

# Returns

- `out::Array{Float64, 3}`
"""
function apply(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), f::String)

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    f_tmp = replace(f, "obj" => "$(obj.data[1, :, 1])")
    out_tmp = eval(Meta.parse(f_tmp))
    out = zeros(eltype(out_tmp), ch_n, length(out_tmp), ep_n)

    # initialize progress bar
    progress_bar == true && (p = Progress(ch_n * ep_n, 1))
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            f_tmp = replace(f, "obj" => "$(obj.data[channel[ch_idx], :, ep_idx])")
            try
                out[ch_idx, :, ep_idx] = eval(Meta.parse(f_tmp))

            catch
                @error "Formula is incorrect."
            end
            # update progress bar
            progress_bar == true && next!(p)
        end
    end
    return out
end
