function result = split_str2int(input, header, splitter)
    % SPLIT_STR2INT(input, header, splitter)
    % Splits an input string into Int64 array removing header (if true) and splitting at splitter character.

    if isempty(header)
        header = true;
    end
    if isempty(splitter)
        splitter = ' ';
    end
    s = split(input, splitter);
    s = string(s);
    if header
        s = s(2:end);
    end
    result = zeros(length(s), 1);
    for n = 1:length(s)
        result(n) = str2num(s(n))
    end
end