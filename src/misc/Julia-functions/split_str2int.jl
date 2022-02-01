"""
    split_str2int(input::String, header::Bool=true, splitter::Char=' ')

Splits an `input` string into Int64 array removing `header` (if true) and splitting at `splitter` character.
"""

function split_str2int(input::String, header::Bool=true, splitter::Char=' ')
    s = split(input, splitter)
    if header
        s = s[2:end]
    end
    a = zeros(length(s))
    for n in 1:length(s)
        a[n] = parse(Int64, s[n])
    end
    a = Array{Int64}(a)
    return a
end
