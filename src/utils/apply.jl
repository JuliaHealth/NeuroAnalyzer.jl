export apply

"""
    apply(obj; ch, f)

Apply custom function.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `f::String`: function to be applied, e.g. `f="mean(obj, dims=3)"`; OBJ signal is given using variable `obj` here.

# Returns

- `out::Array{Float64, 3}`
"""
function apply(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), f::String)

    _check_channels(obj, ch)
    isa(ch, Int64) && (ch = [ch])
    ch_n = length(ch)
    ep_n = nepochs(obj)

    f_tmp = replace(f, "obj" => "$(obj.data[1, :, 1])")
    out_tmp = eval(Meta.parse(f_tmp))
    out = zeros(eltype(out_tmp), ch_n, length(out_tmp), ep_n)

    # initialize progress bar
    progress_bar && (progbar = Progress(ch_n * ep_n, dt=1, barlen=20, color=:white))
    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            f_tmp = replace(f, "obj" => "$(obj.data[ch[ch_idx], :, ep_idx])")
            try
                out[ch_idx, :, ep_idx] = eval(Meta.parse(f_tmp))
            catch
                @error "Formula is incorrect."
            end
            # update progress bar
            progress_bar && next!(progbar)
        end
    end
    
    return out

end
