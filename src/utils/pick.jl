export pick

"""
    pick(obj; p)

Return set of channel indices corresponding with `p` of electrodes

# Arguments

- `p::Vector{Symbol}`: pick of electrodes; picks may be combined, e.g. `[:left, :frontal]`
    - `:list`
    - `:central` (or `:c`)
    - `:left` (or `:l`)
    - `:right` (or `:r`)
    - `:frontal` (or `:f`)
    - `:temporal` (or `:t`)
    - `:parietal` (or `:p`)
    - `:occipital` (or `:o`)

# Returns

- `channels::Vector{Int64}`: channel numbers
"""
function pick(obj::NeuroAnalyzer.NEURO; p::Union{Symbol, Vector{Symbol}})

    _check_datatype(obj, "eeg")

    @assert length(labels(obj)) != 0 "OBJ does not contain channel labels."

    if p isa Vector{Symbol}
        for idx in p
            _check_var(idx, [:list, :central, :c, :left, :l, :right, :r, :frontal, :f, :temporal, :t, :parietal, :p, :occipital, :o], "p")
        end

        # convert picks to channel labels
        c = Vector{Char}()
        for idx in p
            (idx === :central || idx === :c) && push!(c, 'z')
            (idx === :frontal || idx === :f) && push!(c, 'F')
            (idx === :temporal || idx === :t) && push!(c, 'T')
            (idx === :parietal || idx === :p) && push!(c, 'P')
            (idx === :occipital || idx === :o) && push!(c, 'O')
        end
        
        # check which channels are in the picks list
        clabels = labels(obj)[get_channel_bytype(obj, type="eeg")]
        channels = Vector{Int64}()
        for idx1 in eachindex(clabels)
            for idx2 in eachindex(c)
                in(c[idx2], clabels[idx1]) && push!(channels, idx1)
            end
        end

        # check for both :l and :r
        for idx1 in eachindex(p)
            if (p[idx1] === :left || p[idx1] === :l)
                for idx2 in eachindex(p)
                    if (p[idx2] === :right || p[idx2] === :r)
                        return channels
                    end
                end
            end
            if (p[idx1] === :right || p[idx1] === :r)
                for idx2 in eachindex(p)
                    if (p[idx2] === :left || p[idx2] === :l)
                        return channels
                    end
                end
            end
        end

        clabels = labels(obj)[get_channel_bytype(obj, type="eeg")]
        clabels = clabels[channels]
        pat = nothing
        for idx in p
            # for :right remove lefts
            (idx === :right || idx === :r) && (pat = r"[z13579]$")
            # for :left remove rights
            (idx === :left || idx === :l) && (pat = r"[z02468]$")
        end
        if typeof(pat) == Regex
            for idx in length(clabels):-1:1
                typeof(match(pat, clabels[idx])) == RegexMatch && deleteat!(channels, idx)
            end
        end

        return channels
    else
        _check_var(p, [:central, :c, :left, :l, :right, :r, :frontal, :f, :temporal, :t, :parietal, :p, :occipital, :o], "p")

        c = Vector{Char}()
        (p === :central || p === :c) && (c = ['z'])
        (p === :left || p === :l) && (c = ['1', '3', '5', '7', '9'])
        (p === :right || p === :r) && (c = ['2', '4', '6', '8'])
        (p === :frontal || p === :f) && (c = ['F'])
        (p === :temporal || p === :t) && (c = ['T'])
        (p === :parietal || p === :p) && (c = ['P'])
        (p === :occipital || p === :o) && (c = ['O'])

        clabels = labels(obj)[get_channel_bytype(obj, type="eeg")]
        channels = Vector{Int64}()
        for idx1 in eachindex(c)
            for idx2 in eachindex(clabels)
                in(c[idx1], clabels[idx2]) && push!(channels, idx2)
            end
        end

        return channels
    end
end
