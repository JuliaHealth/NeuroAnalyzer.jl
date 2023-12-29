function _read_fiff_tag(fid::IOStream)
    tag = reinterpret(Int32, read(fid, sizeof(Int32) * 4))
    tag .= ntoh.(tag)
    tag_kind = tag[1]
    tag_type = tag[2]
    tag_size = tag[3]
    tag_next = tag[4]
    current_position = position(fid)

    if tag_next == 0
        seek(fid, current_position + tag_size)
    elseif tag_next > 0
        seek(fid, tag_next)
    end

    return tag_kind, tag_type, tag_size, tag_next
end

function _get_fiff_block_type(fid::IOStream, tags::Vector{NTuple{5, Int64}}, id)
    tag = tags[id]
    seek(fid, tag[1] + 16)
    buf = zeros(UInt8, tag[4])
    readbytes!(fid, buf, tag[4])
    return reinterpret(Int32, reverse(buf))
end

function _read_fiff_tag(fid::IOStream, fiff_blocks::Matrix{Int64}, tag_id::Int64)
    id = _find_fiff_tag(fiff_blocks, tag_id)
    if length(id) == 0
        return nothing
    elseif length(id) == 1
        out = _read_fiff_data(fid, fiff_blocks, id)
        if out isa String
            return out
        else
            return out[]
        end
    else
        out = Vector{Any}()
        for idx in eachindex(id)
            push!(out, _read_fiff_data(fid, fiff_blocks, id[idx]))
        end
        return out
    end
end

function _read_fiff_data(fid::IOStream, fiff_blocks::Matrix{Int64}, id::Union{Int64, Nothing})
    @assert id !== nothing && "id does not contain valid tag id."
    tag = fiff_blocks[id, :]
    seek(fid, tag[1] + 16)
    buf = zeros(UInt8, tag[4])
    readbytes!(fid, buf, tag[4])
    
    if tag[3] == 0
        return nothing
    elseif tag[3] == 1
        return reinterpret(Int8, buf)
    elseif tag[3] == 2
        return ntoh.(reinterpret(Int16, buf))
    elseif tag[3] == 3
        return ntoh.(reinterpret(Int32, buf))
    elseif tag[3] == 4
        return ntoh.(reinterpret(Float32, buf))
    elseif tag[3] == 5
        return ntoh.(reinterpret(Float64, buf))
    elseif tag[3] == 6
        return ntoh.(reinterpret(Float32, buf))
        # julian
    elseif tag[3] == 7
        return ntoh.(reinterpret(UInt16, buf))
    elseif tag[3] == 8
        return ntoh.(reinterpret(UInt32, buf))
    elseif tag[3] == 9
        return ntoh.(reinterpret(UInt64, buf))
    elseif tag[3] == 10
        return String(Char.(buf))
    elseif tag[3] == 11
        return ntoh.(reinterpret(Int64, buf))
    elseif tag[3] == 16
        return ntoh.(reinterpret(Int64, buf))
    elseif tag[3] == 20
        return ntoh.(reinterpret(ComplexF32, buf))
    elseif tag[3] == 21
        return ntoh.(reinterpret(ComplexF64, buf))
    elseif tag[3] == 21
        # old_pack
    elseif tag[3] == 30
        # ch_info_struct
        # recording order
        scan_no = ntoh.(reinterpret(Int32, buf[1:4]))[]
        # logical channels order
        log_no = ntoh.(reinterpret(Int32, buf[5:8]))[]
        # channel type
        kind = ntoh.(reinterpret(Int32, buf[9:12]))[]
        r = reinterpret(Float32, reverse(buf[13:16]))[]
        cal = reinterpret(Float32, reverse(buf[17:20]))[]
        coil_type = reinterpret(Int32, reverse(buf[21:24]))[]
        # coil coordinate system origin
        r01 = reinterpret(Float32, reverse(buf[25:28]))[]
        r02 = reinterpret(Float32, reverse(buf[29:32]))[]
        r03 = reinterpret(Float32, reverse(buf[33:36]))[]
        # coil coordinate system x-axis unit vector
        ex1 = reinterpret(Float32, reverse(buf[37:40]))[]
        ex2 = reinterpret(Float32, reverse(buf[41:44]))[]
        ex3 = reinterpret(Float32, reverse(buf[45:48]))[]
        # coil coordinate system y-axis unit vector
        ey1 = reinterpret(Float32, reverse(buf[49:52]))[]
        ey2 = reinterpret(Float32, reverse(buf[53:56]))[]
        ey3 = reinterpret(Float32, reverse(buf[57:60]))[]
        # coil coordinate system z-axis unit vector
        ez1 = reinterpret(Float32, reverse(buf[61:64]))[]
        ez2 = reinterpret(Float32, reverse(buf[65:68]))[]
        ez3 = reinterpret(Float32, reverse(buf[69:72]))[]
        unit = reinterpret(Int32, reverse(buf[73:76]))[]
        unit_mul = reinterpret(Int32, reverse(buf[77:80]))[]
        return(scan_no, log_no, kind, r, cal, coil_type, r01, r02, r03, ex1, ex2, ex2, ey1, ey2, ey2, ez1, ez2, ez2, unit, unit_mul)
    elseif tag[3] == 31
        # id_struct
        reverse(buf)
    elseif tag[3] == 32
        # dir_entry_struct
        reverse(buf)
    elseif tag[3] == 33
        # dig_point_struct
        reverse(buf)
    elseif tag[3] == 34
        # ch_pos_struct
        reverse(buf)
    elseif tag[3] == 35
        # coord_trans_struct
        reverse(buf)
    elseif tag[3] == 36
        # dig_string_struct
        reverse(buf)
    elseif tag[3] == 37
        # stream_segment_struct
        reverse(buf)
    else
        @error "Unknown tag type $(tag[3])."
    end
end

function _find_fiff_tag(fiff_blocks::Matrix{Int64}, id::Int64)
    ids = findall(x -> x == id, fiff_blocks[:, 2])
    if length(ids) == 1
        return ids[]
    else
        return ids
    end
end

function _extract_struct(s, id::Int64)
    @assert id >= 1 "id must be ≥ 1."
    @assert id <= length(s[1]) "id must be ≤ $(length(s[1]))."
    out = Vector{typeof(s[1][id])}()
    for channels_idx in eachindex(s)
        push!(out, s[channels_idx][id])
    end
    return out
end

function _fiff_channel_type(channel_types::Vector{Int32})
    channel_type = Vector{String}()
    for idx in eachindex(channel_types)
        if channel_types[idx] == 1
            push!(channel_type, "meg")
        elseif channel_types[idx] == 2
            push!(channel_type, "eeg")
        elseif channel_types[idx] == 3
            push!(channel_type, "stim")
        elseif channel_types[idx] == 102
            push!(channel_type, "bio")
        elseif channel_types[idx] == 201
            push!(channel_type, "mcg")
        elseif channel_types[idx] == 202
            push!(channel_type, "eog")
        elseif channel_types[idx] == 301
            push!(channel_type, "meg_ref")
        elseif channel_types[idx] == 302
            push!(channel_type, "emg")
        elseif channel_types[idx] == 402
            push!(channel_type, "ecg")
        elseif channel_types[idx] == 900
            push!(channel_type, "system")
        elseif channel_types[idx] == 910
            push!(channel_type, "ias")
        else
            push!(channel_type, "misc")
        end
    end
    return channel_type
end

function _create_fiff_block(file_name::String)

    fid = nothing
    try
        fid = open(file_name, "r")
    catch
        error("File $file_name cannot be loaded.")
    end

    # read file_id tag
    tag_kind = nothing
    tag_type = nothing
    tag_size = nothing
    tag_next = nothing
    try
        tag_kind, tag_type, tag_size, tag_next = _read_fiff_tag(fid)
    catch
        error("File $file_name first tag cannot be read.")
    end

    # check file_id tag
    @assert tag_kind == fiff_file_id "File $file_name is not a FIFF file."
    @assert tag_type == fiff_id_struct "File $file_name is not a FIFF file."
    @assert tag_size == 20 "File $file_name is not a FIFF file."

    # read dir_pointer tag
    try
        tag_kind, tag_type, tag_size, tag_next = _read_fiff_tag(fid)
    catch
        error("File $file_name dir_pointer tag cannot be read.")
    end

    # check dir_pointer tag
    @assert tag_kind == fiff_dir_pointer "File $file_name has no dir_pointer tag."

    # read tags
    seek(fid, 0)
    # tags: position in file, tag_id, data_type, data_size, next
    tags = Vector{Tuple{Int64, Int64, Int64, Int64, Int64}}()
    while tag_next != -1
        current_position = position(fid)
        tag_kind, tag_type, tag_size, tag_next = _read_fiff_tag(fid)
        push!(tags, (current_position, tag_kind, tag_type, tag_size, tag_next))
    end
    seek(fid, 0)

    # create list of tag IDs
    tag_pos = Vector{Int64}()
    for tag_idx in eachindex(tags)
        push!(tag_pos, tags[tag_idx][1])
    end
    tag_ids = Vector{Int64}()
    for tag_idx in eachindex(tags)
        push!(tag_ids, tags[tag_idx][2])
    end
    tag_type = Vector{Int64}()
    for tag_idx in eachindex(tags)
        push!(tag_type, tags[tag_idx][3])
    end
    tag_size = Vector{Int64}()
    for tag_idx in eachindex(tags)
        push!(tag_size, tags[tag_idx][4])
    end

    # create block structure
    block_level = ones(Int64, length(tags))
    block_number = ones(Int64, length(tags))
    block_type = Vector{Int64}()
    block_type_current = 999
    for tag_idx in eachindex(tags)
        if tags[tag_idx][2] == fiff_block_start
            block_level[tag_idx:end] .+= 1
            block_type_current = _get_fiff_block_type(fid, tags, tag_idx)[]
            push!(block_type, block_type_current)
        elseif tags[tag_idx][2] == fiff_block_end
            push!(block_type, block_type_current)
            block_level[tag_idx:end] .-= 1
        else
            push!(block_type, block_type_current)
        end
        tags[tag_idx][2] == fiff_block_start && (block_number[tag_idx:end] .+= 1)
    end

    return fid, hcat(tag_pos, tag_ids, tag_type, tag_size, block_number, block_level, block_type)
end

function _view_fiff_block(fiff_block::Matrix{Int64})
    for tag_idx in 1:size(fiff_block, 1)
        indent = repeat(" ", (fiff_block[tag_idx, 6] - 1))
        println(indent * "$(fiff_block[tag_idx, 6]) [$(fiff_block[tag_idx, 5])] $(fiff_block[tag_idx, 7])")
    end
end
