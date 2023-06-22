export remove_powerline
export remove_powerline!

"""
    remove_powerline(obj; pl_frq, method)

Remove power line noise and harmonics.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `pl_frq::Real=50`: power line frequency
- `method::Symbol=:iir`: use IIR filter

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function remove_powerline(obj::NeuroAnalyzer.NEURO; pl_frq::Real=50, method::Symbol=:iir)
    
    _check_var(method, [:iir], "method")

    obj_new = deepcopy(obj)

    pl_frq < 0 && throw(ArgumentError("pl_freq must be > 0."))
    pl_frq > sr(obj) / 2 && throw(ArgumentError("pl_freq must be ≤ $(sr(obj) / 2)."))

    if method === :iir    
        p, f = psd(obj_new, ch=1)
        f_pl = vsearch(pl_frq + 10, f)
        p_tmp = p[1, :, 1][f_pl:end]
        p1 = maximum(p_tmp)
        f1 = vsearch(p1, p_tmp)
        pl_f1 = f[f_pl:end][f1]
        p_tmp[f1] = 0
        p2 = maximum(p_tmp)
        f2 = vsearch(p2, p_tmp)
        pl_f2 = f[f_pl:end][f2]

        _info("Removing power line noise at: $pl_frq Hz.")

        NeuroAnalyzer.filter!(obj_new, fprototype=:iirnotch, cutoff=pl_frq, bw=1.5)
        
        _info("Power line harmonics found at: $pl_f1 and $pl_f2 Hz.")

        NeuroAnalyzer.filter!(obj_new, fprototype=:iirnotch, cutoff=pl_f1, bw=.2)
        NeuroAnalyzer.filter!(obj_new, fprototype=:iirnotch, cutoff=pl_f2, bw=.23)

    end

    reset_components!(obj_new)
    push!(obj_new.history, "remove_powerline(OBJ, pl_frq=$pl_frq, method=:$method)")

    return obj_new

end

"""
    remove_powerline!(obj; pl_frq, method)

Remove power line noise and harmonics.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `pl_frq::Real=50`: power line frequency
- `method::Symbol=:iir`: use IIR filter
"""
function remove_powerline!(obj::NeuroAnalyzer.NEURO; pl_frq::Real=50, method::Symbol=:iir)

    obj_new = remove_powerline(obj, pl_frq=pl_frq, method=method)
    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history

    return nothing

end