"""
    _split_df(df; <keyword arguments>)

Split a DataFrame into training and test subsets by random row sampling.

# Arguments

- `df::DataFrame`: input data
- `ratio::Float64=0.8`: fraction of rows assigned to the training set; must be in `(0.0, 1.0)` exclusive

# Returns

- `Tuple{DataFrame, DataFrame}`: `(train, test)` DataFrames
"""
function _split(df::DataFrame, ratio::Float64 = 0.8)::Tuple{DataFrame, DataFrame}
    _bin(ratio, (0.0, 1.0), "ratio")

    n         = DataFrames.nrow(df)
    n_train   = floor(Int, ratio * n)
    idx       = shuffle(1:n)
    train_idx = idx[1:n_train]
    test_idx  = idx[(n_train + 1):end]

    return df[train_idx, :], df[test_idx, :]
end
