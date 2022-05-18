################################
#                              #
# Low-level internal functions #
#                              #
################################

################################

function _make_epochs(signal::Matrix{Float64}; epoch_n::Union{Int64, Nothing}=nothing, epoch_len::Union{Int64, Nothing}=nothing, average::Bool=false)

    (epoch_len === nothing && epoch_n === nothing) && throw(ArgumentError("Either epoch_n or epoch_len must be set."))
    (epoch_len !== nothing && epoch_n !== nothing) && throw(ArgumentError("Both epoch_n and epoch_len cannot be set."))
    (epoch_len !== nothing && epoch_len < 1) && throw(ArgumentError("epoch_len must be ≥ 1."))
    (epoch_n !== nothing && epoch_n < 1) && throw(ArgumentError("epoch_n must be ≥ 1."))

    channel_n, _ = size(signal)

    if epoch_n === nothing
        epoch_n = size(signal, 2) ÷ epoch_len
    else
        epoch_len = size(signal, 2) ÷ epoch_n
    end

    epochs = zeros(channel_n, epoch_len, epoch_n)

    idx1 = 1
    for idx2 in 1:epoch_len:(epoch_n * epoch_len - 1)
        epochs[:, :, idx1] = signal[:, idx2:(idx2 + epoch_len - 1), 1]
        idx1 += 1
    end

    if average == true
        epochs = mean(epochs, dims=3)[:, :]
    end

    return epochs
end

################################

"""
    eeg_add_component(eeg; c, v)

Add component name `c` of value `v` to `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`
- `c::Symbol`: component name
- `v::Any`: component value

# Returns

- `eeg::NeuroJ.EEG`
"""
function eeg_add_component(eeg::NeuroJ.EEG; c::Symbol, v::Any)

    eeg_new = deepcopy(eeg)
    c in eeg_new.eeg_header[:components] && throw(ArgumentError("Component $c already exists. Use eeg_delete_component() to remove it prior the operation."))
    push!(eeg_new.eeg_header[:components], c)
    push!(eeg_new.eeg_components, v)
    push!(eeg_new.eeg_header[:history], "eeg_add_component(EEG, c=$c, v=$v)")

    return eeg_new
end

"""
    eeg_add_component!(eeg; c, v)

Add component name `c` of value `v` to `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`
- `c::Symbol`: component name
- `v::Any`: component value
"""
function eeg_add_component!(eeg::NeuroJ.EEG; c::Symbol, v::Any)

    c in eeg.eeg_header[:components] && throw(ArgumentError("Component $c already exists. Use eeg_delete_component!() to remove it prior the operation."))
    push!(eeg.eeg_header[:components], c)
    push!(eeg.eeg_components, v)
    push!(eeg.eeg_header[:history], "eeg_add_component!(EEG, c=$c, v=$v)")

    nothing
end

"""
    eeg_list_components(eeg)

List `eeg` components.

# Arguments

- `eeg::NeuroJ.EEG`

# Returns

- `components::Vector{Symbol}`
"""
function eeg_list_components(eeg::NeuroJ.EEG)

    return eeg.eeg_header[:components]
end

"""
    eeg_extract_component(eeg, c)

Extract component `c` of `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`
- `c::Symbol`: component name

# Returns

- `component::Any`
"""
function eeg_extract_component(eeg::NeuroJ.EEG; c::Symbol)

    c in eeg.eeg_header[:components] || throw(ArgumentError("Component $c does not exist. Use eeg_list_component() to view existing components."))
    
    for idx in 1:length(eeg.eeg_header[:components])
        if c == eeg.eeg_header[:components][idx]
            return eeg.eeg_components[idx]
        end
    end

    nothing
end

"""
    eeg_delete_component(eeg; c)

Delete component `c` of `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`
- `c::Symbol`: component name

# Returns

- `eeg::NeuroJ.EEG`
"""
function eeg_delete_component(eeg::NeuroJ.EEG; c::Symbol)

    eeg_new = deepcopy(eeg)
    c in eeg_new.eeg_header[:components] || throw(ArgumentError("Component $c does not exist. Use eeg_list_component() to view existing components."))
    for idx in 1:length(eeg.eeg_header[:components])
        if c == eeg_new.eeg_header[:components][idx]
            deleteat!(eeg_new.eeg_components, idx)
            deleteat!(eeg_new.eeg_header[:components], idx)
            push!(eeg_new.eeg_header[:history], "eeg_delete_component(EEG, c=$c)")
            return eeg_new
        end
    end
end

"""
    eeg_delete_component!(eeg; c)

Delete component `c` of `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`
- `c::Symbol`: component name
"""
function eeg_delete_component!(eeg::NeuroJ.EEG; c::Symbol)

    c in eeg.eeg_header[:components] || throw(ArgumentError("Component $c does not exist. Use eeg_list_component() to view existing components."))
    
    for idx in length(eeg.eeg_header[:components]):-1:1
        if c == eeg.eeg_header[:components][idx]
            deleteat!(eeg.eeg_components, idx)
            deleteat!(eeg.eeg_header[:components], idx)
            push!(eeg.eeg_header[:history], "eeg_delete_component(EEG, c=$c)")
        end
    end

    nothing
end

"""
    eeg_reset_components(eeg)

Remove all `eeg` components.

# Arguments

- `eeg:EEG`

# Returns

- `eeg::NeuroJ.EEG`
"""
function eeg_reset_components(eeg::NeuroJ.EEG)

    eeg_new = deepcopy(eeg)
    eeg_new.eeg_header[:components] = []
    eeg_new.eeg_components = []

    return eeg_new
end

"""
    eeg_reset_components!(eeg)

Remove all `eeg` components.

# Arguments

- `eeg:EEG`
"""
function eeg_reset_components!(eeg::NeuroJ.EEG)

    eeg.eeg_header[:components] = []
    eeg.eeg_components = []

    nothing
end

"""
    eeg_component_idx(eeg, c)

Return index of `eeg` component.

# Arguments

- `eeg:EEG`
- `c::Symbol`: component name

# Return

- `c_idx::Int64`
"""
function eeg_component_idx(eeg::NeuroJ.EEG; c::Symbol)

    c in eeg.eeg_header[:components] || throw(ArgumentError("Component $c does not exist. Use eeg_list_component() to view existing components."))
    c_idx = findfirst(isequal(c), eeg.eeg_header[:components])

    return c_idx
end

"""
    eeg_component_type(eeg, c)

Return type of `eeg` components.

# Arguments

- `eeg:EEG`
- `c::Symbol`: component name

# Return

- `c_type::DataType`
"""
function eeg_component_type(eeg::NeuroJ.EEG; c::Symbol)

    c in eeg.eeg_header[:components] || throw(ArgumentError("Component $c does not exist. Use eeg_list_component() to view existing components."))
    c_idx = eeg_component_idx(eeg; c=c)
    c_type = typeof(eeg.eeg_components[c_idx])

    return c_type
end

"""
    eeg_rename_component(eeg, c_old, c_new)

Return type of `eeg` components.

# Arguments

- `eeg:EEG`
- `c_old::Symbol`: old component name
- `c_new::Symbol`: new component name

# Return

- `eeg_new:EEG`
"""
function eeg_rename_component(eeg::NeuroJ.EEG; c_old::Symbol, c_new::Symbol)

    c_old in eeg.eeg_header[:components] || throw(ArgumentError("Component $c_old does not exist. Use eeg_list_component() to view existing components."))
    c_new in eeg.eeg_header[:components] && throw(ArgumentError("Component $c_new already exists. Use eeg_list_component() to view existing components."))

    eeg_new = deepcopy(eeg)
    c_idx = eeg_component_idx(eeg, c=c_old)
    eeg_new.eeg_header[:components][c_idx] = c_new

    push!(eeg_new.eeg_header[:history], "eeg_rename_component(EEG, c_old=$c_old, c_new=$c_new)")

    return eeg_new
end

"""
    eeg_rename_component(eeg, c_old, c_new)

Return type of `eeg` components.

# Arguments

- `eeg:EEG`
- `c_old::Symbol`: old component name
- `c_new::Symbol`: new component name
"""
function eeg_rename_component!(eeg::NeuroJ.EEG; c_old::Symbol, c_new::Symbol)

    c_old in eeg.eeg_header[:components] || throw(ArgumentError("Component $c_old does not exist. Use eeg_list_component() to view existing components."))
    c_new in eeg.eeg_header[:components] && throw(ArgumentError("Component $c_new already exists. Use eeg_list_component() to view existing components."))

    c_idx = eeg_component_idx(eeg, c=c_old)
    eeg.eeg_header[:components][c_idx] = c_new

    push!(eeg.eeg_header[:history], "eeg_rename_component!(EEG, c_old=$c_old, c_new=$c_new)")

    nothing
end

"""
    eeg_delete_channel(eeg; channel)

Remove `channel` from the `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel index to be removed, vector of numbers or range

# Returns

- `eeg::NeuroJ.EEG`
"""
function eeg_delete_channel(eeg::NeuroJ.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange})

    typeof(channel) <: AbstractRange && (channel = collect(channel))
    channel_n = eeg_channel_n(eeg)
    length(channel) == channel_n && throw(ArgumentError("You cannot delete all channels."))

    length(channel) > 1 && (channel = sort!(channel, rev=true))

    if channel[end] < 1 || channel[1] > eeg_channel_n(eeg)
        throw(ArgumentError("channel does not match signal channels."))
    end

    eeg_new = deepcopy(eeg)

    # update headers
    eeg_new.eeg_header[:channel_n] = channel_n - length(channel)
    for idx1 in 1:length(channel)
        for idx2 in 1:channel_n
            if idx2 == channel[idx1]
                deleteat!(eeg_new.eeg_header[:labels], idx2)
                deleteat!(eeg_new.eeg_header[:channel_type], idx2)
                (eeg_new.eeg_header[:channel_locations] == true && length(eeg_new.eeg_header[:loc_theta]) > 0) && deleteat!(eeg_new.eeg_header[:loc_theta], idx2)
                (eeg_new.eeg_header[:channel_locations] == true && length(eeg_new.eeg_header[:loc_radius]) > 0) && deleteat!(eeg_new.eeg_header[:loc_radius], idx2)
                (eeg_new.eeg_header[:channel_locations] == true && length(eeg_new.eeg_header[:loc_x]) > 0) && deleteat!(eeg_new.eeg_header[:loc_x], idx2)
                (eeg_new.eeg_header[:channel_locations] == true && length(eeg_new.eeg_header[:loc_y]) > 0) && deleteat!(eeg_new.eeg_header[:loc_y], idx2)
                (eeg_new.eeg_header[:channel_locations] == true && length(eeg_new.eeg_header[:loc_z]) > 0) && deleteat!(eeg_new.eeg_header[:loc_z], idx2)
                (eeg_new.eeg_header[:channel_locations] == true && length(eeg_new.eeg_header[:loc_radius_sph]) > 0) && deleteat!(eeg_new.eeg_header[:loc_radius_sph], idx2)
                (eeg_new.eeg_header[:channel_locations] == true && length(eeg_new.eeg_header[:loc_theta_sph]) > 0) && deleteat!(eeg_new.eeg_header[:loc_theta_sph], idx2)
                (eeg_new.eeg_header[:channel_locations] == true && length(eeg_new.eeg_header[:loc_phi_sph]) > 0) && deleteat!(eeg_new.eeg_header[:loc_phi_sph], idx2)
                deleteat!(eeg_new.eeg_header[:transducers], idx2)
                deleteat!(eeg_new.eeg_header[:physical_dimension], idx2)
                deleteat!(eeg_new.eeg_header[:physical_minimum], idx2)
                deleteat!(eeg_new.eeg_header[:physical_maximum], idx2)
                deleteat!(eeg_new.eeg_header[:digital_minimum], idx2)
                deleteat!(eeg_new.eeg_header[:digital_maximum], idx2)
                deleteat!(eeg_new.eeg_header[:prefiltering], idx2)
                deleteat!(eeg_new.eeg_header[:samples_per_datarecord], idx2)
                deleteat!(eeg_new.eeg_header[:sampling_rate], idx2)
                deleteat!(eeg_new.eeg_header[:gain], idx2)
            end
        end 
    end

    # remove channel
    eeg_new.eeg_signals = eeg_new.eeg_signals[setdiff(1:end, (channel)), :, :]

    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_delete_channel(EEG, $channel)")

    return eeg_new
end

"""
    eeg_delete_channel!(eeg; channel)

Remove `channel` from the `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel index to be removed
"""
function eeg_delete_channel!(eeg::NeuroJ.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange})

    typeof(channel) <: AbstractRange && (channel = collect(channel))
    channel_n = eeg_channel_n(eeg)
    length(channel) == channel_n && throw(ArgumentError("You cannot delete all channels."))

    length(channel) > 1 && (channel = sort!(channel, rev=true))

    if channel[end] < 1 || channel[1] > eeg_channel_n(eeg)
        throw(ArgumentError("channel does not match signal channels."))
    end

    # update headers
    eeg.eeg_header[:channel_n] = channel_n - length(channel)
    for idx1 in 1:length(channel)
        for idx2 in 1:channel_n
            if idx2 == channel[idx1]
                deleteat!(eeg.eeg_header[:labels], idx2)
                deleteat!(eeg.eeg_header[:channel_type], idx2)
                (eeg.eeg_header[:channel_locations] == true && length(eeg.eeg_header[:loc_theta]) > 0) && deleteat!(eeg.eeg_header[:loc_theta], idx2)
                (eeg.eeg_header[:channel_locations] == true && length(eeg.eeg_header[:loc_radius]) > 0) && deleteat!(eeg.eeg_header[:loc_radius], idx2)
                (eeg.eeg_header[:channel_locations] == true && length(eeg.eeg_header[:loc_x]) > 0) && deleteat!(eeg.eeg_header[:loc_x], idx2)
                (eeg.eeg_header[:channel_locations] == true && length(eeg.eeg_header[:loc_y]) > 0) && deleteat!(eeg.eeg_header[:loc_y], idx2)
                (eeg.eeg_header[:channel_locations] == true && length(eeg.eeg_header[:loc_z]) > 0) && deleteat!(eeg.eeg_header[:loc_z], idx2)
                (eeg.eeg_header[:channel_locations] == true && length(eeg.eeg_header[:loc_radius_sph]) > 0) && deleteat!(eeg.eeg_header[:loc_radius_sph], idx2)
                (eeg.eeg_header[:channel_locations] == true && length(eeg.eeg_header[:loc_theta_sph]) > 0) && deleteat!(eeg.eeg_header[:loc_theta_sph], idx2)
                (eeg.eeg_header[:channel_locations] == true && length(eeg.eeg_header[:loc_phi_sph]) > 0) && deleteat!(eeg.eeg_header[:loc_phi_sph], idx2)
                deleteat!(eeg.eeg_header[:transducers], idx2)
                deleteat!(eeg.eeg_header[:physical_dimension], idx2)
                deleteat!(eeg.eeg_header[:physical_minimum], idx2)
                deleteat!(eeg.eeg_header[:physical_maximum], idx2)
                deleteat!(eeg.eeg_header[:digital_minimum], idx2)
                deleteat!(eeg.eeg_header[:digital_maximum], idx2)
                deleteat!(eeg.eeg_header[:prefiltering], idx2)
                deleteat!(eeg.eeg_header[:samples_per_datarecord], idx2)
                deleteat!(eeg.eeg_header[:sampling_rate], idx2)
                deleteat!(eeg.eeg_header[:gain], idx2)
            end
        end 
    end

    # remove channel
    eeg.eeg_signals = eeg.eeg_signals[setdiff(1:end, (channel)), :, :]

    eeg_reset_components!(eeg)
    push!(eeg.eeg_header[:history], "eeg_delete_channel!(EEG, $channel)")

    nothing
end

"""
    eeg_keep_channel(eeg; channel)

Keep `channels` in the `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel index to keep

# Returns

- `eeg::NeuroJ.EEG`
"""
function eeg_keep_channel(eeg::NeuroJ.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange})

    typeof(channel) <: AbstractRange && (channel = collect(channel))

    length(channel) > 1 && (channel = sort!(channel, rev=true))
    if channel[end] < 1 || channel[1] > eeg_channel_n(eeg)
        throw(ArgumentError("channel does not match signal channels."))
    end

    channel_list = collect(1:eeg_channel_n(eeg))
    channel_to_remove = setdiff(channel_list, channel)

    length(channel_to_remove) > 1 && (channel_to_remove = sort!(channel_to_remove, rev=true))

    eeg_new = deepcopy(eeg)
    channel_n = eeg_new.eeg_header[:channel_n]

    # update headers
    eeg_new.eeg_header[:channel_n] = channel_n - length(channel_to_remove)
    for idx1 in 1:length(channel_to_remove)
        for idx2 in channel_n:-1:1
            if idx2 == channel_to_remove[idx1]
                deleteat!(eeg_new.eeg_header[:labels], idx2)
                deleteat!(eeg_new.eeg_header[:channel_type], idx2)
                (eeg_new.eeg_header[:channel_locations] == true && length(eeg_new.eeg_header[:loc_theta]) > 0) && deleteat!(eeg_new.eeg_header[:loc_theta], idx2)
                (eeg_new.eeg_header[:channel_locations] == true && length(eeg_new.eeg_header[:loc_radius]) > 0) && deleteat!(eeg_new.eeg_header[:loc_radius], idx2)
                (eeg_new.eeg_header[:channel_locations] == true && length(eeg_new.eeg_header[:loc_x]) > 0) && deleteat!(eeg_new.eeg_header[:loc_x], idx2)
                (eeg_new.eeg_header[:channel_locations] == true && length(eeg_new.eeg_header[:loc_y]) > 0) && deleteat!(eeg_new.eeg_header[:loc_y], idx2)
                (eeg_new.eeg_header[:channel_locations] == true && length(eeg_new.eeg_header[:loc_z]) > 0) && deleteat!(eeg_new.eeg_header[:loc_z], idx2)
                (eeg_new.eeg_header[:channel_locations] == true && length(eeg_new.eeg_header[:loc_radius_sph]) > 0) && deleteat!(eeg_new.eeg_header[:loc_radius_sph], idx2)
                (eeg_new.eeg_header[:channel_locations] == true && length(eeg_new.eeg_header[:loc_theta_sph]) > 0) && deleteat!(eeg_new.eeg_header[:loc_theta_sph], idx2)
                (eeg_new.eeg_header[:channel_locations] == true && length(eeg_new.eeg_header[:loc_phi_sph]) > 0) && deleteat!(eeg_new.eeg_header[:loc_phi_sph], idx2)
                deleteat!(eeg_new.eeg_header[:transducers], idx2)
                deleteat!(eeg_new.eeg_header[:physical_dimension], idx2)
                deleteat!(eeg_new.eeg_header[:physical_minimum], idx2)
                deleteat!(eeg_new.eeg_header[:physical_maximum], idx2)
                deleteat!(eeg_new.eeg_header[:digital_minimum], idx2)
                deleteat!(eeg_new.eeg_header[:digital_maximum], idx2)
                deleteat!(eeg_new.eeg_header[:prefiltering], idx2)
                deleteat!(eeg_new.eeg_header[:samples_per_datarecord], idx2)
                deleteat!(eeg_new.eeg_header[:sampling_rate], idx2)
                deleteat!(eeg_new.eeg_header[:gain], idx2)
            end
        end
    end

    # remove channel
    eeg_new.eeg_signals = eeg_new.eeg_signals[setdiff(1:end, (channel_to_remove)), :, :]

    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_keep_channel(EEG, $channel)")

    return eeg_new
end

"""
    eeg_keep_channel!(eeg; channel)

Keep `channels` in the `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel index to keep
"""
function eeg_keep_channel!(eeg::NeuroJ.EEG; channel::Union{Int64, Vector{Int64}, AbstractRange})

    typeof(channel) <: AbstractRange && (channel = collect(channel))

    length(channel) > 1 && (channel = sort!(channel, rev=true))
    if channel[end] < 1 || channel[1] > eeg_channel_n(eeg)
        throw(ArgumentError("channel does not match signal channels."))
    end

    channel_list = collect(1:eeg_channel_n(eeg))
    channel_to_remove = setdiff(channel_list, channel)

    length(channel_to_remove) > 1 && (channel_to_remove = sort!(channel_to_remove, rev=true))

    channel_n = eeg_channel_n(eeg)

    # update headers
    eeg.eeg_header[:channel_n] = channel_n - length(channel_to_remove)
    for idx1 in 1:length(channel_to_remove)
        for idx2 in 1:channel_n
            if idx2 == channel_to_remove[idx1]
                deleteat!(eeg.eeg_header[:labels], idx2)
                deleteat!(eeg.eeg_header[:channel_type], idx2)
                (eeg.eeg_header[:channel_locations] == true && length(eeg.eeg_header[:loc_theta]) > 0) && deleteat!(eeg.eeg_header[:loc_theta], idx2)
                (eeg.eeg_header[:channel_locations] == true && length(eeg.eeg_header[:loc_radius]) > 0) && deleteat!(eeg.eeg_header[:loc_radius], idx2)
                (eeg.eeg_header[:channel_locations] == true && length(eeg.eeg_header[:loc_x]) > 0) && deleteat!(eeg.eeg_header[:loc_x], idx2)
                (eeg.eeg_header[:channel_locations] == true && length(eeg.eeg_header[:loc_y]) > 0) && deleteat!(eeg.eeg_header[:loc_y], idx2)
                (eeg.eeg_header[:channel_locations] == true && length(eeg.eeg_header[:loc_z]) > 0) && deleteat!(eeg.eeg_header[:loc_z], idx2)
                (eeg.eeg_header[:channel_locations] == true && length(eeg.eeg_header[:loc_radius_sph]) > 0) && deleteat!(eeg.eeg_header[:loc_radius_sph], idx2)
                (eeg.eeg_header[:channel_locations] == true && length(eeg.eeg_header[:loc_theta_sph]) > 0) && deleteat!(eeg.eeg_header[:loc_theta_sph], idx2)
                (eeg.eeg_header[:channel_locations] == true && length(eeg.eeg_header[:loc_phi_sph]) > 0) && deleteat!(eeg.eeg_header[:loc_phi_sph], idx2)
                deleteat!(eeg.eeg_header[:transducers], idx2)
                deleteat!(eeg.eeg_header[:physical_dimension], idx2)
                deleteat!(eeg.eeg_header[:physical_minimum], idx2)
                deleteat!(eeg.eeg_header[:physical_maximum], idx2)
                deleteat!(eeg.eeg_header[:digital_minimum], idx2)
                deleteat!(eeg.eeg_header[:digital_maximum], idx2)
                deleteat!(eeg.eeg_header[:prefiltering], idx2)
                deleteat!(eeg.eeg_header[:samples_per_datarecord], idx2)
                deleteat!(eeg.eeg_header[:sampling_rate], idx2)
                deleteat!(eeg.eeg_header[:gain], idx2)
            end
        end
    end

    # remove channel
    eeg.eeg_signals = eeg.eeg_signals[setdiff(1:end, (channel_to_remove)), :, :]

    eeg_reset_components!(eeg)
    push!(eeg.eeg_header[:history], "eeg_keep_channel!(EEG, channel=$channel)")

    nothing
end

"""
    eeg_get_channel(eeg; channel)

Return the `channel` index / name.

# Arguments

- `eeg::NeuroJ.EEG`
- `channel::Union{Int64, String}`: channel name

# Returns

- `channel_idx::Int64`
"""
function eeg_get_channel(eeg::NeuroJ.EEG; channel::Union{Int64, String})

    labels = eeg_labels(eeg)
    if typeof(channel) == String
        channel_idx = nothing
        for idx in 1:length(labels)
            if channel == labels[idx]
                channel_idx = idx
            end
        end
        if channel_idx === nothing
            throw(ArgumentError("channel name does not match signal labels."))
        end
        return channel_idx
    else
        if channel < 1 || channel > length(labels)
            throw(ArgumentError("channel index does not match signal channels."))
        end
        return labels[channel]
    end

    nothing
end

"""
    eeg_rename_channel(eeg; channel, new_name)

Renames the `eeg` `channel`.

# Arguments

- `eeg::NeuroJ.EEG`
- `channel::Union{Int64, String}`
- `new_name::String`

# Returns

- `eeg::NeuroJ.EEG`
"""
function eeg_rename_channel(eeg::NeuroJ.EEG; channel::Union{Int64, String}, new_name::String)

    # create new dataset
    eeg_new = deepcopy(eeg)
    labels = eeg_labels(eeg_new)
    
    if typeof(channel) == String
        channel_found = nothing
        for idx in 1:length(labels)
            if channel == labels[idx]
                labels[idx] = new_name
                channel_found = idx
            end
        end
        if channel_found === nothing
            throw(ArgumentError("channel name does not match signal labels."))
        end
    else
        if channel < 1 || channel > length(labels)
            throw(ArgumentError("channel index does not match signal channels."))
        else
            labels[channel] = new_name
        end
    end
    eeg_new.eeg_header[:labels] = labels
    
    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_rename_channel(EEG, channel=$channel, new_name=$new_name)")

    return eeg_new
end

"""
    eeg_rename_channel!(eeg; channel, new_name)

Renames the `eeg` `channel`.

# Arguments

- `eeg::NeuroJ.EEG`
- `channel::Union{Int64, String}`
- `new_name::String`
"""
function eeg_rename_channel!(eeg::NeuroJ.EEG; channel::Union{Int64, String}, new_name::String)

    labels = eeg_labels(eeg)
    
    if typeof(channel) == String
        channel_found = nothing
        for idx in 1:length(labels)
            if channel == labels[idx]
                labels[idx] = new_name
                channel_found = idx
            end
        end
        if channel_found === nothing
            throw(ArgumentError("channel name does not match signal labels."))
        end
    else
        if channel < 1 || channel > length(labels)
            throw(ArgumentError("channel index does not match signal channels."))
        else
            labels[channel] = new_name
        end
    end
    eeg.eeg_header[:labels] = labels
    
    push!(eeg.eeg_header[:history], "eeg_rename_channel!(EEG, channel=$channel, new_name=$new_name)")

    nothing
end

"""
    eeg_extract_channel(eeg; channel)

Extract `channel` number or name.

# Arguments

- `eeg::NeuroJ.EEG`
- `channel::Union{Int64, String}`

# Returns

- `channel::Vector{Float64}`
"""
function eeg_extract_channel(eeg::NeuroJ.EEG; channel::Union{Int64, String})

    labels = eeg_labels(eeg)
    if typeof(channel) == String
        channel_idx = nothing
        for idx in 1:length(labels)
            if channel == labels[idx]
                channel_idx = idx
            end
        end
        if channel_idx === nothing
            throw(ArgumentError("channel name does not match signal labels."))
        end
    else
        if channel < 1 || channel > length(labels)
            throw(ArgumentError("channel index does not match signal channels."))
        end
        channel_idx = channel
    end    
    eeg_channel = reshape(eeg.eeg_signals[channel_idx, :, :], 1, eeg_epoch_len(eeg), eeg_epoch_n(eeg))

    return eeg_channel
end

"""
    eeg_history(eeg)

Show processing history.

# Arguments

- `eeg::NeuroJ.EEG`

# Returns

- `eeg::NeuroJ.EEG`
"""
function eeg_history(eeg::NeuroJ.EEG)

    return eeg.eeg_header[:history]
end

"""
    eeg_labels(eeg)

Return `eeg` labels.

# Arguments

- `eeg::NeuroJ.EEG`

# Returns

- `labels::Vector{String}`
"""
function eeg_labels(eeg::NeuroJ.EEG)

    return eeg.eeg_header[:labels]
end

"""
    eeg_sr(eeg)

Return `eeg` sampling rate.

# Arguments

- `eeg::NeuroJ.EEG`

# Returns

- `sr::Int64`
"""
function eeg_sr(eeg::NeuroJ.EEG)

    return eeg.eeg_header[:sampling_rate][1]
end

"""
    eeg_channel_n(eeg; type=:eeg)

Return number of `eeg` channels of `type`.

# Arguments

- `eeg::NeuroJ.EEG`
- `type::Vector{Symbol}[:all, :eeg, :ecg, :eog, :emg]`

# Returns

- `channel_n::Int64`
"""
function eeg_channel_n(eeg::NeuroJ.EEG; type::Symbol=:all)

    channel_n = 0
    for idx in 1:eeg.eeg_header[:channel_n]
        eeg.eeg_header[:channel_type][idx] == string(type) && (channel_n += 1)
    end
    type === :all && (channel_n = size(eeg.eeg_signals, 1))

    return channel_n
end

"""
    eeg_epoch_n(eeg)

Return number of `eeg` epochs.

# Arguments

- `eeg::NeuroJ.EEG`

# Returns

- `epoch_n::Int64`
"""
function eeg_epoch_n(eeg::NeuroJ.EEG)

    epoch_n = eeg.eeg_header[:epoch_n]

    return epoch_n
end

"""
    eeg_signal_len(eeg)

Return length of `eeg` signal.

# Arguments

- `eeg::NeuroJ.EEG`

# Returns

- `signal_len::Int64`
"""
function eeg_signal_len(eeg::NeuroJ.EEG)

    return eeg.eeg_header[:eeg_duration_samples]
end

"""
    eeg_epoch_len(eeg)

Return length of `eeg` signal.

# Arguments

- `eeg::NeuroJ.EEG`

# Returns

- `epoch_len::Int64`
"""
function eeg_epoch_len(eeg::NeuroJ.EEG)

    epoch_len = eeg.eeg_header[:epoch_duration_samples]

    return epoch_len
end

"""
    eeg_info(eeg)

Show info.

# Arguments

- `eeg::NeuroJ.EEG`
"""
function eeg_info(eeg::NeuroJ.EEG)

    println("          EEG file name: $(eeg.eeg_header[:eeg_filename])")
    println("          EEG size [Mb]: $(eeg.eeg_header[:eeg_filesize_mb])")
    println("     Sampling rate (Hz): $(eeg_sr(eeg))")
    println("Signal length (samples): $(eeg_signal_len(eeg))")
    println("Signal length (seconds): $(round(eeg.eeg_header[:eeg_duration_seconds], digits=2))")
    println("     Number of channels: $(eeg_channel_n(eeg))")
    println("       Number of epochs: $(eeg_epoch_n(eeg))")
    println(" Epoch length (samples): $(eeg_epoch_len(eeg))")
    println(" Epoch length (seconds): $(round(eeg.eeg_header[:epoch_duration_seconds], digits=2))")
    if eeg.eeg_header[:reference] == ""
        println("         Reference type: unknown")
    else
        println("         Reference type: $(eeg.eeg_header[:reference])")
    end
    if length(eeg_labels(eeg)) == 0
        println("                 Labels: no")
    else
        println("                 Labels: yes")
    end
    if eeg.eeg_header[:channel_locations] == false
        println("      Channel locations: no")
    else
        println("      Channel locations: yes")
    end
    if eeg.eeg_header[:components] != []
        print("             Components: ")
        c = eeg_list_components(eeg)
        if length(c) == 1
            println(c[1])
        else
            for idx in 1:(length(c) - 1)
                print(c[idx], ", ")
            end
            println(c[end])
        end
    else
        println("             Components: no")
    end
    println("               Channels:")
    for idx in 1:length(eeg.eeg_header[:labels])
        println("                channel: $idx\tlabel: $(eeg.eeg_header[:labels][idx])\ttype: $(uppercase(eeg.eeg_header[:channel_type][idx]))")
    end

    nothing
end

"""
    eeg_epochs(eeg; epoch_n=nothing, epoch_len=nothing, average=false)

Splits `eeg` into epochs.

# Arguments

- `eeg::NeuroJ.EEG`
- `epoch_n::Union{Int64, Nothing}`: number of epochs
- `epoch_len::Union{Int64, Nothing}`: epoch length in samples
- `average::Bool`: average all epochs, return one averaged epoch; if false than return array of epochs, each row is one epoch

# Returns

- `eeg::NeuroJ.EEG`
"""
function eeg_epochs(eeg::NeuroJ.EEG; epoch_n::Union{Int64, Nothing}=nothing, epoch_len::Union{Int64, Nothing}=nothing, average::Bool=false)

    # unsplit epochs
    s_merged = reshape(eeg.eeg_signals,
                       eeg_channel_n(eeg),
                       eeg_epoch_len(eeg) * eeg_epoch_n(eeg))
    
    # split into epochs
    s_split = _make_epochs(s_merged, epoch_n=epoch_n, epoch_len=epoch_len, average=average)

    # convert into Array{Float64, 3}
    s_split = reshape(s_split, size(s_split, 1), size(s_split, 2), size(s_split, 3))

    # create new dataset
    epoch_n = size(s_split, 3)
    epoch_duration_samples = size(s_split, 2)
    epoch_duration_seconds = size(s_split, 2) / eeg.eeg_header[:sampling_rate][1]
    eeg_duration_samples = size(s_split, 2) * size(s_split, 3)
    eeg_duration_seconds = eeg_duration_samples / eeg.eeg_header[:sampling_rate][1]
    eeg_time = collect(0:(1 / eeg.eeg_header[:sampling_rate][1]):epoch_duration_seconds)
    eeg_time = eeg_time[1:(end - 1)]
    eeg_new = deepcopy(eeg)
    eeg_new.eeg_signals = s_split
    eeg_new.eeg_time = eeg_time

    # update epochs time
    fs = eeg_sr(eeg_new)
    ts = eeg.eeg_epochs_time[1, 1]
    new_epochs_time = linspace(ts, ts + (epoch_duration_samples / fs), epoch_duration_samples)
    eeg_new.eeg_epochs_time = repeat(new_epochs_time, 1, epoch_n)

    eeg_new.eeg_header[:eeg_duration_samples] = eeg_duration_samples
    eeg_new.eeg_header[:eeg_duration_seconds] = eeg_duration_seconds
    eeg_new.eeg_header[:epoch_n] = epoch_n
    eeg_new.eeg_header[:epoch_duration_samples] = epoch_duration_samples
    eeg_new.eeg_header[:epoch_duration_seconds] = epoch_duration_seconds

    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_epochs(EEG, epoch_n=$epoch_n, epoch_len=$epoch_len, average=$average)")

    return eeg_new
end

"""
    eeg_epochs!(eeg; epoch_n=nothing, epoch_len=nothing, average=false)

Splits `eeg` into epochs.

# Arguments

- `eeg::NeuroJ.EEG`
- `epoch_n::Union{Int64, Nothing}`: number of epochs
- `epoch_len::Union{Int64, Nothing}`: epoch length in samples
- `average::Bool`: average all epochs, return one averaged epoch
"""
function eeg_epochs!(eeg::NeuroJ.EEG; epoch_n::Union{Int64, Nothing}=nothing, epoch_len::Union{Int64, Nothing}=nothing, average::Bool=false)

    # unsplit epochs
    s_merged = reshape(eeg.eeg_signals,
                       eeg_channel_n(eeg),
                       eeg_epoch_len(eeg) * eeg_epoch_n(eeg))
    
    # split into epochs
    s_split = _make_epochs(s_merged, epoch_n=epoch_n, epoch_len=epoch_len, average=average)

    # convert into Array{Float64, 3}
    s_split = reshape(s_split, size(s_split, 1), size(s_split, 2), size(s_split, 3))

    # create new dataset
    epoch_n = size(s_split, 3)
    epoch_duration_samples = size(s_split, 2)
    epoch_duration_seconds = size(s_split, 2) / eeg.eeg_header[:sampling_rate][1]
    eeg_duration_samples = size(s_split, 2) * size(s_split, 3)
    eeg_duration_seconds = eeg_duration_samples / eeg.eeg_header[:sampling_rate][1]
    eeg_time = collect(0:(1 / eeg.eeg_header[:sampling_rate][1]):epoch_duration_seconds)
    eeg_time = eeg_time[1:(end - 1)]
    eeg.eeg_signals = s_split
    eeg.eeg_time = eeg_time

    # update epochs time
    fs = eeg_sr(eeg)
    ts = eeg.eeg_epochs_time[1, 1]
    new_epochs_time = linspace(ts, ts + (epoch_duration_samples / fs), epoch_duration_samples)
    eeg.eeg_epochs_time = repeat(new_epochs_time, 1, epoch_n)

    eeg.eeg_header[:eeg_duration_samples] = eeg_duration_samples
    eeg.eeg_header[:eeg_duration_seconds] = eeg_duration_seconds
    eeg.eeg_header[:epoch_n] = epoch_n
    eeg.eeg_header[:epoch_duration_samples] = epoch_duration_samples
    eeg.eeg_header[:epoch_duration_seconds] = epoch_duration_seconds

    eeg_reset_components!(eeg)
    push!(eeg.eeg_header[:history], "eeg_epochs!(EEG, epoch_n=$epoch_n, epoch_len=$epoch_len, average=$average)")

    nothing
end

"""
    eeg_extract_epoch(eeg; epoch)

Extract the `epoch` epoch.

# Arguments

- `eeg::NeuroJ.EEG`
- `epoch::Int64`: epoch index

# Returns

- `eeg::NeuroJ.EEG`
"""
function eeg_extract_epoch(eeg::NeuroJ.EEG; epoch::Int64)

    if epoch < 1 || epoch > eeg_epoch_n(eeg)
        throw(ArgumentError("epoch must be ≥ 1 and ≤ $(eeg_epoch_n(eeg))."))
    end

    s_new = reshape(eeg.eeg_signals[:, :, epoch], eeg_channel_n(eeg), eeg_signal_len(eeg), 1)
    eeg_new = deepcopy(eeg)
    eeg_new.eeg_signals = s_new
    eeg_new.eeg_epochs_time = repeat(eeg.eeg_epochs_time[:, epoch], 1, 1)
    eeg_new.eeg_header[:epoch_n] = 1
    eeg_new.eeg_header[:eeg_duration_samples] = eeg_new.eeg_header[:epoch_duration_samples]
    eeg_new.eeg_header[:eeg_duration_seconds] = eeg_new.eeg_header[:epoch_duration_seconds]

    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_get_epoch(EEG, epoch=$epoch)")

    return eeg_new
end

"""
    eeg_trim(eeg:EEG; len, offset=0, from=:start, keep_epochs::Bool=true)

Remove `len` samples from the beginning + `offset` (`from` = :start, default) or end (`from` = :end) of the `eeg`.

# Arguments

- `eeg:EEG`
- `len::Int64`: number of samples to remove
- `offset::Int64`: offset from which trimming starts, only works for `from` = :start
- `from::Symbol[:start, :end]`: trims from the signal start (default) or end
- `keep_epochs::Bool`: remove epochs containing signal to trim (keep_epochs=true) or remove signal and remove epoching

# Returns

- `eeg:EEG`
"""
function eeg_trim(eeg::NeuroJ.EEG; len::Int64, offset::Int64=1, from::Symbol=:start, keep_epochs::Bool=true)

    eeg_epoch_n(eeg) == 1 && (keep_epochs = false)

    if keep_epochs == false
        @warn "This operation will remove epochs, to keep epochs use keep_epochs=true."

        eeg_new = deepcopy(eeg)
        eeg_epoch_n(eeg) > 1 && (eeg_epochs!(eeg_new, epoch_n=1))
        channel_n = eeg_channel_n(eeg_new)
        epoch_n = eeg_epoch_n(eeg_new)
        s_trimmed = zeros(channel_n, (size(eeg_new.eeg_signals, 2) - len), epoch_n)
        @inbounds @simd for epoch_idx in 1:epoch_n
            Threads.@threads for idx in 1:channel_n
                s = @view eeg_new.eeg_signals[idx, :, epoch_idx]
                s_trimmed[idx, :, epoch_idx] = s_trim(s, len=len, offset=offset, from=from)
            end
        end
        t_trimmed = collect(0:(1 / eeg_sr(eeg)):(size(s_trimmed, 2) / eeg_sr(eeg)))[1:(end - 1)]
        eeg_new.eeg_signals = s_trimmed
        eeg_new.eeg_time = t_trimmed
        eeg_new.eeg_epochs_time = hcat(t_trimmed)
        eeg_new.eeg_header[:eeg_duration_samples] -= len
        eeg_new.eeg_header[:eeg_duration_seconds] -= len * (1 / eeg_sr(eeg))
        eeg_new.eeg_header[:epoch_duration_samples] -= len
        eeg_new.eeg_header[:epoch_duration_seconds] -= len * (1 / eeg_sr(eeg))
    else
        if from === :start
            epoch_from = floor(Int64, (offset / eeg_epoch_len(eeg)) + 1)
            epoch_to = ceil(Int64, ((offset + len) / eeg_epoch_len(eeg)) + 1)
        else
            epoch_from = floor(Int64, ((eeg.eeg_header[:eeg_duration_samples] - len) / eeg_epoch_len(eeg)) + 1)
            epoch_to = eeg_epoch_n(eeg)
        end
        eeg_new = eeg_delete_epoch(eeg, epoch=epoch_from:epoch_to)
    end

    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_trim(EEG, len=$len, offset=$offset, from=$from, keep_epochs=$keep_epochs)")

    return eeg_new
end

"""
    eeg_trim!(eeg:EEG; len, offset=0, from=:start, keep_epochs::Bool=true)

Remove `len` samples from the beginning + `offset` (`from` = :start, default) or end (`from` = :end) of the `eeg`.

# Arguments

- `eeg:EEG`
- `len::Int64`: number of samples to remove
- `offset::Int64`: offset from which trimming starts, only works for `from` = :start
- `from::Symbol[:start, :end]`: trims from the signal start (default) or end
- `keep_epochs::Bool`: remove epochs containing signal to trim (keep_epochs=true) or remove signal and remove epoching
"""
function eeg_trim!(eeg::NeuroJ.EEG; len::Int64, offset::Int64=1, from::Symbol=:start, keep_epochs::Bool=true)

    eeg_len = eeg_signal_len(eeg)
    from in [:start, :end] || throw(ArgumentError("from must be :start or :end."))
    len < 0 && throw(ArgumentError("len must be ≥ 1."))
    len >= eeg_len && throw(ArgumentError("len must be < $(eeg_len)."))
    offset < 0 && throw(ArgumentError("offset must be ≥ 1."))
    offset >= eeg_len - 1 && throw(ArgumentError("offset must be < $(eeg_len)."))
    (from ===:start && 1 + offset + len > eeg_len) && throw(ArgumentError("offset + len must be < $(eeg_len)."))

    eeg_epoch_n(eeg) == 1 && (keep_epochs = false)
    
    if keep_epochs == false
        @warn "This operation will remove epochs, to keep epochs use keep_epochs=true."
        eeg_epoch_n(eeg) > 1 && (eeg_epochs!(eeg, epoch_n=1))
        channel_n = eeg_channel_n(eeg)
        epoch_n = eeg_epoch_n(eeg)
        s_trimmed = zeros(channel_n, (size(eeg.eeg_signals, 2) - len), epoch_n)
        @inbounds @simd for epoch_idx in 1:epoch_n
            Threads.@threads for idx in 1:channel_n
                s = @view eeg.eeg_signals[idx, :, epoch_idx]
                s_trimmed[idx, :, epoch_idx] = s_trim(s, len=len, offset=offset, from=from)
            end
        end
        t_trimmed = collect(0:(1 / eeg_sr(eeg)):(size(s_trimmed, 2) / eeg_sr(eeg)))[1:(end - 1)]
        eeg.eeg_signals = s_trimmed
        eeg.eeg_time = t_trimmed
        eeg.eeg_epochs_time = hcat(t_trimmed)
        eeg.eeg_header[:eeg_duration_samples] -= len
        eeg.eeg_header[:eeg_duration_seconds] -= len * (1 / eeg_sr(eeg))
        eeg.eeg_header[:epoch_duration_samples] -= len
        eeg.eeg_header[:epoch_duration_seconds] -= len * (1 / eeg_sr(eeg))
    else
        if from === :start
            epoch_from = floor(Int64, (offset / eeg_epoch_len(eeg))) + 1
            epoch_to = floor(Int64, ((offset + len) / eeg_epoch_len(eeg))) + 1
        else
            epoch_from = floor(Int64, ((eeg.eeg_header[:eeg_duration_samples] - len) / eeg_epoch_len(eeg)) + 1)
            epoch_to = eeg_epoch_n(eeg)
        end
        eeg_delete_epoch!(eeg, epoch=epoch_from:epoch_to)
    end

    eeg_reset_components!(eeg)
    push!(eeg.eeg_header[:history], "eeg_trim!(EEG, len=$len, offset=$offset, from=$from, keep_epochs=$keep_epochs)")

    nothing
end

"""
    eeg_edit_header(eeg; field, value)

Change value of `eeg` `field` to `value`.

# Arguments

- `eeg::NeuroJ.EEG`
- `field::Symbol`
- `value::Any`

# Returns

- `eeg:EEG`
"""
function eeg_edit_header(eeg::NeuroJ.EEG; field::Symbol, value::Any)

    field === nothing && throw(ArgumentError("field cannot be empty."))
    value === nothing && throw(ArgumentError("value cannot be empty."))

    eeg_new = deepcopy(eeg)
    fields = keys(eeg_new.eeg_header)
    field in fields || throw(ArgumentError("$field does not exist."))
    typeof(eeg_new.eeg_header[field]) == typeof(value) || throw(ArgumentError("field type ($(typeof(eeg_new.eeg_header[field]))) does not mach value type ($(typeof(value)))."))
    eeg_new.eeg_header[field] = value
    push!(eeg_new.eeg_header[:history], "eeg_edit(EEG, field=$field, value=$value)")    

    return eeg_new
end

"""
    eeg_edit_header!(eeg; field, value)

Change value of `eeg` `field` to `value`.

# Arguments

- `eeg::NeuroJ.EEG`
- `field::Symbol`
- `value::Any`

# Returns

- `eeg:EEG`
"""
function eeg_edit_header!(eeg::NeuroJ.EEG; field::Symbol, value::Any)

    value === nothing && throw(ArgumentError("value cannot be empty."))

    fields = keys(eeg.eeg_header)
    field in fields || throw(ArgumentError("field does not exist."))
    typeof(eeg.eeg_header[field]) == typeof(value) || throw(ArgumentError("field type ($(typeof(eeg_new.eeg_header[field]))) does not mach value type ($(typeof(value)))."))
    eeg.eeg_header[field] = value
    push!(eeg.eeg_header[:history], "eeg_edit!(EEG, field=$field, value=$value)")    

    nothing
end

"""
    eeg_show_header(eeg)

Show keys and values of `eeg` header.

# Arguments

- `eeg::NeuroJ.EEG`
"""
function eeg_show_header(eeg::NeuroJ.EEG)

    for (key, value) in eeg.eeg_header
        println("$key: $value")
    end
end

"""
    eeg_delete_epoch(eeg; epoch)

Remove `epoch` from the `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}`: epoch index to be removed, vector of numbers or range

# Returns

- `eeg::NeuroJ.EEG`
"""
function eeg_delete_epoch(eeg::NeuroJ.EEG; epoch::Union{Int64, Vector{Int64}, AbstractRange})

    eeg_epoch_n(eeg) == 1 && throw(ArgumentError("You cannot delete the last epoch."))

    if typeof(epoch) <: AbstractRange
        epoch = collect(epoch)
    end

    length(epoch) == eeg_epoch_n(eeg) && throw(ArgumentError("You cannot delete all epochs."))

    length(epoch) > 1 && (epoch = sort!(epoch, rev=true))

    if epoch[end] < 1 || epoch[1] > eeg_epoch_n(eeg)
        throw(ArgumentError("epoch does not match signal epochs."))
    end

    eeg_new = deepcopy(eeg)

    # remove epoch
    eeg_new.eeg_signals = eeg_new.eeg_signals[:, :, setdiff(1:end, (epoch))]
    eeg_new.eeg_epochs_time = eeg_new.eeg_epochs_time[:, setdiff(1:end, (epoch))]
    # update headers
    eeg_new.eeg_header[:epoch_n] -= length(epoch)
    epoch_n = eeg_new.eeg_header[:epoch_n]
    eeg_new.eeg_header[:eeg_duration_samples] = epoch_n * size(eeg.eeg_signals, 2)
    eeg_new.eeg_header[:eeg_duration_seconds] = round((epoch_n * size(eeg.eeg_signals, 2)) / eeg_sr(eeg), digits=2)
    eeg_new.eeg_header[:epoch_duration_samples] = size(eeg.eeg_signals, 2)
    eeg_new.eeg_header[:epoch_duration_seconds] = round(size(eeg.eeg_signals, 2) / eeg_sr(eeg), digits=2)

    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_delete_epoch(EEG, $epoch)")
    
    return eeg_new
end

"""
    eeg_delete_epoch!(eeg; epoch)

Remove `epoch` from the `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}`: epoch index to be removed, vector of numbers or range
"""
function eeg_delete_epoch!(eeg::NeuroJ.EEG; epoch::Union{Int64, Vector{Int64}, AbstractRange})

    eeg_epoch_n(eeg) == 1 && throw(ArgumentError("You cannot delete the last epoch."))

    if typeof(epoch) <: AbstractRange
        epoch = collect(epoch)
    end

    length(epoch) == eeg_epoch_n(eeg) && throw(ArgumentError("You cannot delete all epochs."))

    length(epoch) > 1 && (epoch = sort!(epoch, rev=true))

    if epoch[end] < 1 || epoch[1] > eeg_epoch_n(eeg)
        throw(ArgumentError("epoch does not match signal epochs."))
    end

    # remove epoch
    eeg.eeg_signals = eeg.eeg_signals[:, :, setdiff(1:end, (epoch))]
    eeg.eeg_epochs_time = eeg.eeg_epochs_time[:, setdiff(1:end, (epoch))]

    # update headers
    eeg.eeg_header[:epoch_n] -= length(epoch)
    epoch_n = eeg.eeg_header[:epoch_n]
    eeg.eeg_header[:eeg_duration_samples] = epoch_n * size(eeg.eeg_signals, 2)
    eeg.eeg_header[:eeg_duration_seconds] = round((epoch_n * size(eeg.eeg_signals, 2)) / eeg_sr(eeg), digits=2)
    eeg.eeg_header[:epoch_duration_samples] = size(eeg.eeg_signals, 2)
    eeg.eeg_header[:epoch_duration_seconds] = round(size(eeg.eeg_signals, 2) / eeg_sr(eeg), digits=2)

    eeg_reset_components!(eeg)
    push!(eeg.eeg_header[:history], "eeg_delete_epoch!(EEG, $epoch)")

    nothing
end

"""
    eeg_keep_epoch(eeg; epoch)

Keep `epoch` in the `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}`: epoch index to keep, vector of numbers or range

# Returns

- `eeg::NeuroJ.EEG`
"""
function eeg_keep_epoch(eeg::NeuroJ.EEG; epoch::Union{Int64, Vector{Int64}, AbstractRange})

    eeg_epoch_n(eeg) == 1 && throw(ArgumentError("EEG contains only one epoch."))

    if typeof(epoch) <: AbstractRange
        epoch = collect(epoch)
    end

    length(epoch) > 1 && (epoch = sort!(epoch, rev=true))
    if epoch[end] < 1 || epoch[1] > eeg_epoch_n(eeg)
        throw(ArgumentError("epoch does not match signal epochs."))
    end

    epoch_list = collect(1:eeg_epoch_n(eeg))
    epoch_to_remove = setdiff(epoch_list, epoch)

    length(epoch_to_remove) > 1 && (epoch_to_remove = sort!(epoch_to_remove, rev=true))

    eeg_new = deepcopy(eeg)

    # remove epoch
    eeg_new.eeg_signals = eeg_new.eeg_signals[:, :, setdiff(1:end, (epoch_to_remove))]
    eeg_new.eeg_epochs_time = eeg_new.eeg_epochs_time[:, setdiff(1:end, (epoch))]

    # update headers
    eeg_new.eeg_header[:epoch_n] = eeg_new.eeg_header[:epoch_n] - length(epoch_to_remove)
    epoch_n = eeg_new.eeg_header[:epoch_n]
    eeg_new.eeg_header[:eeg_duration_samples] = epoch_n * eeg_epoch_len(eeg)
    eeg_new.eeg_header[:eeg_duration_seconds] = round(epoch_n * (eeg_epoch_len(eeg) / eeg_sr(eeg)), digits=2)
    eeg_new.eeg_header[:epoch_duration_samples] = eeg_epoch_len(eeg)
    eeg_new.eeg_header[:epoch_duration_seconds] = round(eeg_epoch_len(eeg) / eeg_sr(eeg), digits=2)

    eeg_reset_components!(eeg_new)
    push!(eeg_new.eeg_header[:history], "eeg_keep_epoch(EEG, $epoch)")    

    return eeg_new
end

"""
    eeg_keep_epoch!(eeg; epoch)

Keep `epoch` in the `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}`: epoch index to keep, vector of numbers or range
"""
function eeg_keep_epoch!(eeg::NeuroJ.EEG; epoch::Union{Int64, Vector{Int64}, AbstractRange})

    eeg_epoch_n(eeg) == 1 && throw(ArgumentError("EEG contains only one epoch."))

    if typeof(epoch) <: AbstractRange
        epoch = collect(epoch)
    end

    length(epoch) > 1 && (epoch = sort!(epoch, rev=true))
    if epoch[end] < 1 || epoch[1] > eeg_epoch_n(eeg)
        throw(ArgumentError("epoch does not match signal epochs."))
    end

    epoch_list = collect(1:eeg_epoch_n(eeg))
    epoch_to_remove = setdiff(epoch_list, epoch)

    length(epoch_to_remove) > 1 && (epoch_to_remove = sort!(epoch_to_remove, rev=true))

    # remove epoch
    eeg.eeg_signals = eeg.eeg_signals[:, :, setdiff(1:end, (epoch_to_remove))]
    eeg.eeg_epochs_time = eeg.eeg_epochs_time[:, setdiff(1:end, (epoch))]

    # update headers
    eeg.eeg_header[:epoch_n] = eeg_epoch_n(eeg) - length(epoch_to_remove)
    epoch_n = eeg_epoch_n(eeg)
    eeg.eeg_header[:eeg_duration_samples] = epoch_n * eeg_epoch_len(eeg)
    eeg.eeg_header[:eeg_duration_seconds] = round(epoch_n * (eeg_epoch_len(eeg) / eeg_sr(eeg)), digits=2)
    eeg.eeg_header[:epoch_duration_samples] = eeg_epoch_len(eeg)
    eeg.eeg_header[:epoch_duration_seconds] = round(eeg_epoch_len(eeg) / eeg_sr(eeg), digits=2)

    eeg_reset_components!(eeg)
    push!(eeg.eeg_header[:history], "eeg_keep_epoch(EEG, $epoch)")

    nothing
end

"""
    eeg_detect_bad_epochs(eeg; method=[:flat, :rmse, :rmsd, :euclid, :p2p], ch_t)

Detect bad `eeg` epochs based on:
- flat channel(s)
- RMSE
- RMSD
- Euclidean distance
- peak-to-peak amplitude

# Arguments

- `eeg::NeuroJ.EEG`
- `method::Vector{Symbol}=[:flat, :rmse, :rmsd, :euclid, :p2p]`
- `ch_t::Float64`: percentage of bad channels to mark the epoch as bad

# Returns

- `bad_epochs_idx::Vector{Int64}`
"""
function eeg_detect_bad_epochs(eeg::NeuroJ.EEG; method::Vector{Symbol}=[:flat, :rmse, :rmsd, :euclid, :p2p], ch_t::Float64=0.1)

    eeg_channel_n(eeg, type=:eeg) < eeg_channel_n(eeg, type=:all) && throw(ArgumentError("EEG contains non-eeg channels (e.g. ECG or EMG), remove them before processing."))

    for idx in method
        idx in [:flat, :rmse, :rmsd, :euclid, :p2p] || throw(ArgumentError("method must be :flat, :rmse, :rmsd, :euclid, :p2p"))
    end

    bad_epochs_idx = zeros(Int64, eeg_epoch_n(eeg))

    if :flat in method
        bad_epochs = s_detect_epoch_flat(eeg.eeg_signals)
        bad_epochs_idx[bad_epochs .> ch_t] .= 1
    end

    if :rmse in method
        bad_epochs = s_detect_epoch_rmse(eeg.eeg_signals)
        bad_epochs_idx[bad_epochs .> ch_t] .= 1
    end

    if :rmsd in method
        bad_epochs = s_detect_epoch_rmsd(eeg.eeg_signals)
        bad_epochs_idx[bad_epochs .> ch_t] .= 1
    end

    if :euclid in method
        bad_epochs = s_detect_epoch_euclid(eeg.eeg_signals)
        bad_epochs_idx[bad_epochs .> ch_t] .= 1
    end

    if :p2p in method
        bad_epochs = s_detect_epoch_p2p(eeg.eeg_signals)
        bad_epochs_idx[bad_epochs .> ch_t] .= 1
    end

    return bad_epochs_idx
end

"""
    eeg_add_labels(eeg::NeuroJ.EEG, labels::Vector{String})

Add `labels` to `eeg` channels.

# Arguments

- `eeg::NeuroJ.EEG`
- `labels::Vector{String}`

# Returns

- `eeg::NeuroJ.EEG`
"""
function eeg_add_labels(eeg::NeuroJ.EEG, labels::Vector{String})

    length(labels) == eeg_channel_n(eeg) || throw(ArgumentError("labels length must be $(eeg_channel_n(eeg))."))
    
    eeg_new = deepcopy(eeg)
    eeg_new.eeg_header[:labels] = labels

    push!(eeg_new.eeg_header[:history], "eeg_add_labels(EEG, labels=$labels")
 
    return eeg_new
end

"""
    eeg_add_labels!(eeg::NeuroJ.EEG, labels::Vector{String})

Add `labels` to `eeg` channels.

# Arguments

- `eeg::NeuroJ.EEG`
- `labels::Vector{String}`
"""
function eeg_add_labels!(eeg::NeuroJ.EEG, labels::Vector{String})

    length(labels) == eeg_channel_n(eeg) || throw(ArgumentError("labels length must be $(eeg_channel_n(eeg))."))
    eeg.eeg_header[:labels] = labels
    push!(eeg.eeg_header[:history], "eeg_add_labels(EEG, labels=$labels")
    
    nothing
end

"""
    eeg_edit_channel(eeg; channel, field, value)
Edits `eeg` `channel` properties.

# Arguments

- `eeg:EEG`
- `channel::Int64`
- `field::Any`
- `value::Any`

# Returns

- `eeg_new::NeuroJ.EEG`
"""
function eeg_edit_channel(eeg::NeuroJ.EEG; channel::Int64, field::Any, value::Any)
    
    field === nothing && throw(ArgumentError("field cannot be empty."))
    value === nothing && throw(ArgumentError("value cannot be empty."))
    (channel < 0 || channel > eeg_channel_n(eeg, type=:all)) && throw(ArgumentError("channel must be > 0 and ≤ $(eeg_channel_n(eeg, type=:all))."))
    
    field in [:channel_type, :loc_theta, :loc_radius, :loc_x, :loc_y, :loc_z, :loc_radius_sph, :loc_theta_sph, :loc_phi_sph, :labels] || throw(ArgumentError("field must be: :channel_type, :loc_theta, :loc_radius, :loc_x, :loc_y, :loc_z, :loc_radius_sph, :loc_theta_sph, :loc_phi_sph, :labels."))

    eeg_new = deepcopy(eeg)
    typeof(eeg_new.eeg_header[field][channel]) == typeof(value) || throw(ArgumentError("field type ($(eltype(eeg_new.eeg_header[field]))) does not mach value type ($(typeof(value)))."))
    eeg_new.eeg_header[field][channel] = value

    # add entry to :history field
    push!(eeg_new.eeg_header[:history], "eeg_edit_channel(EEG, channel=$channel, field=$field, value=$value)")   

    return eeg_new
end

"""
    eeg_edit_channel!(eeg; channel, field, value)

Edit `eeg` `channel` properties.

# Arguments

- `eeg:EEG`
- `channel::Int64`
- `field::Any`
- `value::Any`
"""
function eeg_edit_channel!(eeg::NeuroJ.EEG; channel::Int64, field::Any, value::Any)
    
    field === nothing && throw(ArgumentError("field cannot be empty."))
    value === nothing && throw(ArgumentError("value cannot be empty."))
    (channel < 0 || channel > eeg_channel_n(eeg, type=:all)) && throw(ArgumentError("channel must be > 0 and ≤ $(eeg_channel_n(eeg, type=:all))."))
    
    field in [:channel_type, :loc_theta, :loc_radius, :loc_x, :loc_y, :loc_z, :loc_radius_sph, :loc_theta_sph, :loc_phi_sph, :labels] || throw(ArgumentError("field must be: :channel_type, :loc_theta, :loc_radius, :loc_x, :loc_y, :loc_z, :loc_radius_sph, :loc_theta_sph, :loc_phi_sph, :labels."))

    typeof(eeg.eeg_header[field][channel]) == typeof(value) || throw(ArgumentError("field type ($(eltype(eeg.eeg_header[field]))) does not mach value type ($(typeof(value)))."))
    eeg.eeg_header[field][channel] = value

    push!(eeg.eeg_header[:history], "eeg_edit_channel(EEG, channel=$channel, field=$field, value=$value)")

    nothing
end

"""
    eeg_keep_eeg_channels(eeg::NeuroJ.EEG)

Keep only EEG channels of `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`

# Returns

- `eeg::NeuroJ.EEG`
"""
function eeg_keep_eeg_channels(eeg::NeuroJ.EEG)

    eeg_channel_n(eeg, type=:eeg) == eeg_channel_n(eeg, type=:all) && (return eeg)

    eeg_channels_idx = Vector{Int64}()
    for idx in 1:eeg_channel_n(eeg, type=:all)
        eeg.eeg_header[:channel_type][idx] == "eeg" && push!(eeg_channels_idx, idx)
    end
    eeg_new = eeg_keep_channel(eeg, channel=eeg_channels_idx)

    return eeg_new
end

"""
    eeg_keep_eeg_channels!(eeg::NeuroJ.EEG)

Keep only EEG channels of `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`
"""
function eeg_keep_eeg_channels!(eeg::NeuroJ.EEG)

    eeg_channel_n(eeg, type=:eeg) == eeg_channel_n(eeg, type=:all) && return

    eeg_channels_idx = Vector{Int64}()
    for idx in 1:eeg_channel_n(eeg, type=:all)
        eeg.eeg_header[:channel_type][idx] == "eeg" && push!(eeg_channels_idx, idx)
    end
    eeg_keep_channel!(eeg, channel=eeg_channels_idx)

    nothing
end

"""
    eeg_comment(eeg)

Return `eeg` comment.

# Arguments

- `eeg::NeuroJ.EEG`
"""
function eeg_comment(eeg::NeuroJ.EEG)

    return eeg.eeg_header[:comment]
end

"""
    eeg_copy(eeg::NeuroJ.EEG)

Make copy of `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`

# Returns

- `eeg::NeuroJ.EEG`

"""
function eeg_copy(eeg::NeuroJ.EEG)
    
    return deepcopy(eeg)
end

"""
    eeg_epochs_time(eeg; ts)

Edit `eeg` epochs time start.

# Arguments

- `eeg::NeuroJ.EEG`
- `ts::Real`: time start in seconds

# Returns

- `eeg::NeuroJ.EEG`
"""
function eeg_epochs_time(eeg::NeuroJ.EEG; ts::Real)

    epoch_n = eeg_epoch_n(eeg)
    epoch_len = eeg_epoch_len(eeg)
    fs = eeg_sr(eeg)
    new_epochs_time = linspace(ts, ts + (epoch_len / fs), epoch_len)
    eeg_new = deepcopy(eeg)
    eeg_new.eeg_epochs_time = repeat(new_epochs_time, 1, epoch_n)
    push!(eeg_new.eeg_header[:history], "eeg_epochs_time(EEG, ts=$ts)")

    return eeg_new
end

"""
    eeg_epochs_time!(eeg; ts)

Edit `eeg` epochs time start.

# Arguments

- `eeg::NeuroJ.EEG`
- `ts::Real`: time start in seconds

# Returns

- `eeg::NeuroJ.EEG`
"""
function eeg_epochs_time!(eeg::NeuroJ.EEG; ts::Real)

    epoch_n = eeg_epoch_n(eeg)
    epoch_len = eeg_epoch_len(eeg)
    fs = eeg_sr(eeg)
    new_epochs_time = linspace(ts, ts + (epoch_len / fs), epoch_len)
    eeg.eeg_epochs_time = repeat(new_epochs_time, 1, epoch_n)
    push!(eeg.eeg_header[:history], "eeg_epochs_time!(EEG, ts=$ts)")

    nothing
end