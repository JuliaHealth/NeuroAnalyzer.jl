function _select_channels(obj::NeuroAnalyzer.NEURO, channel::Union{Int64, Vector{Int64}, <:AbstractRange}, def_chn::Int64=0)
    # select channels, default is all or def_chn
    def_chn > nchannels(obj) && (def_chn = nchannels(obj))
    def_chn == 0 && (def_chn = nchannels(obj))
    channel == 0 && (channel = 1:def_chn)
    typeof(channel) <: AbstractRange && (channel = collect(channel))
    length(channel) > 1 && sort!(channel)
    return channel
end

function _select_epochs(obj::NeuroAnalyzer.NEURO, epoch::Union{Int64, Vector{Int64}, <:AbstractRange}, def_ep::Int64=0)
    # select epochs, default is all or def_ep
    def_ep > nepochs(obj) && (def_ep = nepochs(obj))
    def_ep == 0 && (def_ep = nepochs(obj))
    epoch == 0 && (epoch = 1:def_ep)
    typeof(epoch) <: AbstractRange && (epoch = collect(epoch))
    length(epoch) > 1 && sort!(epoch)
    return epoch
end
