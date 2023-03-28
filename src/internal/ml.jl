function _split(df::DataFrame, ratio::Float64=0.8)
    n = nrow(df)
    idx = shuffle(1:n)
    train_idx = view(idx, 1:floor(Int, ratio * n))
    test_idx = view(idx, (floor(Int, ratio * n) + 1):n)
    return df[train_idx, :], df[test_idx, :]
end
