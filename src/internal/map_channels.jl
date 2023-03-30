function _map_channels(ch::Union{Int64, Vector{Int64}, <:AbstractRange}, chs=Vector{Int64})
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
