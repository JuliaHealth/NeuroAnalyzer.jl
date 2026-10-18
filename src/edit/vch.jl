export vch

"""
    vch(obj; <keyword arguments>)

Calculate a virtual channel using formula `f`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `f::String`: channel calculation formula, e.g. `"cz / mean(Fp1 + Fp2)"`; all standard Julia math operators are available, channel labels must be the same as of the OBJ object

# Returns

- `Array{Float64, 3}`: shape (1, time, epochs)
"""
function vch(obj::NeuroAnalyzer.NEURO; f::String)::Array{Float64, 3}
    # number of epochs
    ep_n = nepochs(obj)
    # channel labels
    clabels = labels(obj)
    # pre-allocate output
    vc = zeros(1, epoch_len(obj), ep_n)

    # a label matches only as a whole token, not inside another label or a function name
    # (e.g. "Fp1" does not match in "Fp10", a channel "E" does not match in "mean")
    token(l) = Regex("(?<![A-Za-z0-9_])\\Q" * l * "\\E(?![A-Za-z0-9_])")

    # ── step 1: find referenced channels once, outside the loop ──────────────────
    active_idx = [i for i in eachindex(clabels) if occursin(token(clabels[i]), f)]
    # longer labels first, so that a shorter label can never match part of a longer one
    active_idx = sort(active_idx; by = i -> length(clabels[i]), rev = true)

    # ── step 2: replace channel labels with numbered variable names ──────────────
    # e.g. "Fp1 + Fp2" → "_x1 + _x2"; this avoids embedding array data as a string
    f_var = f
    for (k, ch_idx) in enumerate(active_idx)
        f_var = replace(f_var, token(clabels[ch_idx]) => "_x$k")
    end

    # ── step 3: parse and compile the formula into a function - once ─────────────
    # (a local function: no global state shared between concurrent calls)
    arg_list = Base.join(["_x$k" for k in eachindex(active_idx)], ", ")
    f_func = try
        eval(Meta.parse("($arg_list) -> @. $f_var"))
    catch
        throw(ArgumentError("Formula is incorrect, check channel labels and operators."))
    end

    # ── step 4: parallel loop - no eval, no parsing, no string alloc per epoch ───
    # @view avoids copying channel slices; ntuple builds the argument list
    Threads.@threads :static for ep_idx = 1:ep_n
        args = ntuple(k -> @view(obj.data[active_idx[k], :, ep_idx]), length(active_idx))
        @inbounds vc[1, :, ep_idx] = Base.invokelatest(f_func, args...)
    end

    return vc
end
