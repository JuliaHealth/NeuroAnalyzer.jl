"""
    split_str2float(input::String, header::Bool=true, splitter::Char=' ')

Splits an `input` string into Float64 array removing `header` (if true) and splitting at `splitter` character.
"""

function split_str2float(input::String, header::Bool = true, splitter::Char = ' ')
    s = split(input, splitter)
    if header
        s = s[2:end]
    end
    a = zeros(length(s))
    for n in 1:length(s)
        a[n] = parse(Float64, s[n])
    end
    return a
end
