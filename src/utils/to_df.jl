export to_df

"""
    to_df(obj)

Export object signal data as a `DataFrame`.

Each row corresponds to one time point across all epochs (epochs are concatenated). Columns are `:time` followed by one column per channel, named after the channel labels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object

# Returns

- `DataFrame`: Table with `signal_len(obj)` rows and `nchannels(obj) + 1` columns (`[:time, ch1, ch2, …]`).

# Notes

- For multi-epoch objects the time column contains `obj.time_pts`, which spans the full concatenated signal, not individual epoch times.
- Channel columns are named via `labels(obj)`; ensure the object has channel labels assigned before calling this function.

# See also
[`labels`](@ref), [`signal_len`](@ref)
"""
function to_df(obj::NeuroAnalyzer.NEURO)::DataFrame

    # reshape (channels × samples × epochs) → (channels × total_samples),
    # then transpose to (total_samples × channels) for row-per-timepoint layout
    data_matrix = reshape(obj.data, nchannels(obj), :)'
    df = DataFrame(hcat(obj.time_pts, data_matrix), :auto)
    DataFrames.rename!(df, vcat(:time, Symbol.(labels(obj))))

    return df

end
