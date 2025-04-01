export normalize

"""
    normalize(obj; <keyword arguments>)

Normalize channel(s).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `method::Symbol`:
    - `:zscore`: by z-score
    - `:minmax`: in [-1, +1]
    - `:log`: using log-transformation
    - `:log10`: using log10-transformation
    - `:neglog`: using -log-transformation
    - `:neglog10`: using -log10-transformation
    - `:neg`: in [-∞, 0]
    - `:pos`: in [0, +∞]
    - `:perc`: in percentages
    - `:gauss`: to Gaussian
    - `:invroot`: to inverse root: 1/sqrt(x)
    - `:n`: in [0, n], default is [0, 1]; <keyword arguments>) .+ n1`
    - `:softmax`: using softmax function: exp(x_i) / sum(exp(x))
    - `:sigmoid`: using sigmoid function: 1 /  1 + exp(-x_i)
    - `:mad`: by MAD
    - `:rank`: using tiedranks
    - `:none`
- `bych::Bool=false`: if true, normalize each channel separately

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function normalize(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, method::Symbol, bych::Bool=false)::NeuroAnalyzer.NEURO

    ch = get_channel(obj, ch=ch)
    ch_n = length(ch)
    ep_n = nepochs(obj)

    obj_new = deepcopy(obj)
    if bych
        @inbounds for ep_idx in 1:ep_n
            Threads.@threads :greedy for ch_idx in 1:ch_n
                @views obj_new.data[ch[ch_idx], :, ep_idx] = NeuroStats.normalize(obj_new.data[ch[ch_idx], :, ep_idx], method=method)
            end
        end
    else
        obj_new.data[ch, :, :] = NeuroStats.normalize(obj_new.data[ch, :, :], method=method, bych=false)
    end

    reset_components!(obj_new)
    push!(obj_new.history, "normalize(OBJ, ch=$ch, method=$method)")

    return obj_new

end

"""
    normalize!(obj; <keyword arguments>)

Normalize channel(s).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `method::Symbol`:
    - `:zscore`: by z-score
    - `:minmax`: in [-1, +1]
    - `:log`: using log-transformation
    - `:log10`: using log10-transformation
    - `:neglog`: using -log-transformation
    - `:neglog10`: using -log10-transformation
    - `:neg`: in [-∞, 0]
    - `:pos`: in [0, +∞]
    - `:perc`: in percentages
    - `:gauss`: to Gaussian
    - `:invroot`: to inverse root: 1/sqrt(x)
    - `:n`: in [0, n], default is [0, 1]
    - `:softmax`: using softmax function: exp(x_i) / sum(exp(x))
    - `:sigmoid`: using sigmoid function: 1 /  1 + exp(-x_i)
    - `:mad`: by MAD
    - `:rank`: using tiedranks
    - `:none`
- `bych::Bool=false`: if true, normalize each channel separately

# Returns

Nothing
"""
function normalize!(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, method::Symbol, bych::Bool=false)::Nothing

    obj_new = NeuroAnalyzer.normalize(obj, ch=ch, method=method, bych=bych)
    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history

    return nothing

end
