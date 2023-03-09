"""
    pick(obj; pick)

Return set of channel indices corresponding with `pick` of electrodes

# Arguments

- `pick::Vector{Symbol}`: pick of electrodes; picks may be combined, e.g. `[:left, :frontal]`
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
function pick(obj::NeuroAnalyzer.NEURO; pick::Union{Symbol, Vector{Symbol}})

    length(labels(obj)) == 0 && throw(ArgumentError("OBJ does not contain channel labels."))

    if typeof(pick) == Vector{Symbol}
        for idx in pick
            _check_var(idx, [:list, :central, :c, :left, :l, :right, :r, :frontal, :f, :temporal, :t, :parietal, :p, :occipital, :o], "pick")
        end

        # convert picks to channel labels
        c = Vector{Char}()
        for idx in pick
            (idx === :central || idx === :c) && push!(c, 'z')
            (idx === :frontal || idx === :f) && push!(c, 'F')
            (idx === :temporal || idx === :t) && push!(c, 'T')
            (idx === :parietal || idx === :p) && push!(c, 'P')
            (idx === :occipital || idx === :o) && push!(c, 'O')
        end
        
        # check which channels are in the picks list
        clabels = labels(obj)
        channels = Vector{Int64}()
        for idx1 in eachindex(clabels)
            for idx2 in eachindex(c)
                in(c[idx2], clabels[idx1]) && push!(channels, idx1)
            end
        end

        # check for both :l and :r
        for idx1 in eachindex(pick)
            if (pick[idx1] === :left || pick[idx1] === :l)
                for idx2 in eachindex(pick)
                    if (pick[idx2] === :right || pick[idx2] === :r)
                        return channels
                    end
                end
            end
            if (pick[idx1] === :right || pick[idx1] === :r)
                for idx2 in eachindex(pick)
                    if (pick[idx2] === :left || pick[idx2] === :l)
                        return channels
                    end
                end
            end
        end

        clabels = labels(obj)
        clabels = clabels[channels]
        pat = nothing
        for idx in pick
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
        _check_var(pick, [:central, :c, :left, :l, :right, :r, :frontal, :f, :temporal, :t, :parietal, :p, :occipital, :o], "pick")

        c = Vector{Char}()
        (pick === :central || pick === :c) && (c = ['z'])
        (pick === :left || pick === :l) && (c = ['1', '3', '5', '7', '9'])
        (pick === :right || pick === :r) && (c = ['2', '4', '6', '8'])
        (pick === :frontal || pick === :f) && (c = ['F'])
        (pick === :temporal || pick === :t) && (c = ['T'])
        (pick === :parietal || pick === :p) && (c = ['P'])
        (pick === :occipital || pick === :o) && (c = ['O'])

        clabels = labels(obj)
        channels = Vector{Int64}()
        for idx1 in eachindex(c)
            for idx2 in eachindex(clabels)
                in(c[idx1], clabels[idx2]) && push!(channels, idx2)
            end
        end

        return channels
    end
end
