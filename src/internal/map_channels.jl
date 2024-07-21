function _map_channels(ch::Union{String, Vector{String}}, chs=Vector{Int64})
    ch_orig = ch
    if ch isa Int64
        ch = vsearch(ch, chs)
    else
        for idx in eachindex(ch)
            ch[idx] = vsearch(ch[idx], chs)
        end
    end
    return ch, ch_orig
end
