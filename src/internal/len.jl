function _len(obj::NeuroAnalyzer.NEURO, len::Int64, def_l::Int64)::Int64
    # return default length: one epoch (if len < def_l) or def_l seconds
    if len == 0
        if epoch_len(obj) > def_l * sr(obj)
            len = def_l * sr(obj)
        else
            len = epoch_len(obj)
        end
    end
    return len
end
