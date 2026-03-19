export vch

"""
    vch(obj; <keyword arguments>)

Calculate virtual channel using formula `f`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `f::String`: channel calculation formula, e.g. `"cz / mean(fp1 + fp2)"`; case of labels in the formula is ignored, all standard Julia math operators are available, channel labels must be the same as of the OBJ object

# Returns

- `vc::Array{Float64, 3}`: single channel × time × epochs
"""
function vch(obj::NeuroAnalyzer.NEURO; f::String)::Array{Float64, 3}

    ep_n = nepochs(obj)
    f = lowercase(f)
    clabels = lowercase.(labels(obj))
    vc = zeros(1, epoch_len(obj), ep_n)

#=
    # ── Step 1: find referenced channels once, outside the loop ──────────────────
    # occursin against the original formula; result never changes between epochs.
    active_idx = [i for i in eachindex(clabels) if occursin(clabels[i], f)]

    # ── Step 2: replace channel label strings with numbered variable names ────────
    # e.g. "fp1 + fp2" → "_x1 + _x2"
    # This avoids embedding array data as a string entirely.
    f_var = f
    for (k, ch_idx) in enumerate(active_idx)
        f_var = replace(f_var, clabels[ch_idx] => "_x$k")
    end

    # ── Step 3: parse and compile the formula into a real Julia function — once ───
    # eval is called exactly once here instead of ep_n times inside the loop.
    arg_list = join(["_x$k" for k in eachindex(active_idx)], ", ")
    try
        f_func = eval(Meta.parse("($arg_list) -> @. $f_var"))
    catch
        @error "Formula is incorrect, check channel labels and operators."
    end

    # ── Step 4: parallel loop — no eval, no parsing, no string alloc per epoch ───
    # @view avoids copying channel slices; ntuple builds a type-stable argument list.
    Threads.@threads :static for ep_idx in 1:ep_n
        args = ntuple(k -> @view(obj.data[active_idx[k], :, ep_idx]), length(active_idx))
        @inbounds vc[1, :, ep_idx] = f_func(args...)
    end
=#

    Threads.@threads :static for ep_idx in 1:ep_n
        f_tmp = f
        @inbounds for ch_idx in eachindex(clabels)
            occursin(clabels[ch_idx], f) && (
                f_tmp = replace(
                    f_tmp, clabels[ch_idx] => "$(obj.data[ch_idx, :, ep_idx])"
                )
            )
        end
        try
            @inbounds vc[1, :, ep_idx] = eval(Meta.parse("@. " * f_tmp))
        catch
            @error "Formula is incorrect, check channel labels and operators."
        end
    end

    return vc

end
