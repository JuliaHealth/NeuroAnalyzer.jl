function _len(obj::NeuroAnalyzer.NEURO, len::Int64, def_l::Int64)
    # return default length: one epoch (if epoch_len_seconds < def_l) or def_l seconds
    if len == 0
        if epoch_len(obj) > def_l * sr(obj)
            len = def_l * sr(obj)
        else
            len = epoch_len(obj)
        end
    end
    return len
end
