function _mesh_normalize_xyz(msh::GeometryBasics.AbstractMesh{3, Float32})::Float64
    msh_m = zeros(length(msh.position), 3)
    for idx in eachindex(msh.position)
        msh_m[idx, :] = msh.position[idx]
    end
    return maximum(abs.(msh_m))
end

function _mesh_normalize_xy(msh::GeometryBasics.AbstractMesh{3, Float32})::Float64
    msh_m = zeros(length(msh.position), 2)
    for idx in eachindex(msh.position)
        msh_m[idx, :] = msh.position[idx][1:2]
    end
    return maximum(abs.(msh_m))
end
