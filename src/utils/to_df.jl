export to_df

"""
    to_df(obj)

Export object data as DataFrames.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `df::DataFrame`: DataFrame containing time points and channels
"""
function to_df(obj::NeuroAnalyzer.NEURO)::DataFrame

    df = DataFrame(hcat(obj.time_pts, reshape(obj.data, nchannels(obj), :, 1)[:, :]'), :auto)
    DataFrames.rename!(df, vcat(:time, Symbol.(labels(obj))))

    return df

end
