export vch

"""
    vch(obj; f)

Calculate virtual channel using formula `f`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `f::String`: channel calculation formula, e.g. `"cz / mean(fp1 + fp2)"`; case of labels in the formula is ignored, all standard Julia math operators are available, channel labels must be the same as of the OBJ object

# Returns
 
- `vc::Array{Float64, 3}`: single channel × time × epochs
"""
function vch(obj::NeuroAnalyzer.NEURO; f::String)

    ep_n = nepochs(obj)
    f = lowercase(f)
    clabels = lowercase.(labels(obj))
    vc = zeros(1, epoch_len(obj), ep_n)
    Threads.@threads for ep_idx in 1:ep_n
        f_tmp = f
        @inbounds for ch_idx in eachindex(clabels)
            occursin(clabels[ch_idx], f) && (f_tmp = replace(f_tmp, clabels[ch_idx] => "$(obj.data[ch_idx, :, ep_idx])"))
        end
        try
            @inbounds vc[1, :, ep_idx] = eval(Meta.parse("@. " * f_tmp))
        catch
            @error "Formula is incorrect, check channel labels and operators."
        end
    end

    return vc

end

