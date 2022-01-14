function result = split_str2float(input, header, splitter)
    % SPLIT_STR2INT(input, header, splitter)
    % Splits an input string into Float64 array removing header (if true) and splitting at splitter character.

    s = split(input, splitter);
    s = string(s);
    if header
        s = s(2:end);
    end
    result = zeros(length(s), 1);
    for n = 1:length(s)
        result(n) = str2double(s(n))
    end
end