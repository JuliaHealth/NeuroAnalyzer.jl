export to_df

"""
    to_df(obj)

Export OBJ data to DataFrame.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `df::DataFrame`: DataFrame containing time points and channels
"""
function to_df(obj::NeuroAnalyzer.NEURO)

    df = DataFrame(hcat(obj.time_pts, reshape(obj.data, channel_n(obj), :, 1)[:, :]'), :auto)
    DataFrames.rename!(df, vcat(:time, Symbol.(labels(obj))))

    return df

end
