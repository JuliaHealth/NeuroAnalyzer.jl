function _free_gpumem(threshold::Real=0.95)::Nothing
    if na_gpu === :cuda
        m = CUDA.MemoryInfo()
        usedmem = m.total_bytes - m.free_bytes
        totalmem = m.total_bytes
        if usedmem / totalmem > threshold
            # CUDA.reclaim()
            GC.gc(true)
        end
    end
    return nothing
end
