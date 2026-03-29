"""
    _mesh_normalize_xyz(msh)

Return the maximum absolute coordinate value across all x, y, z positions in `msh`.

Used as a uniform scaling factor to normalize the mesh to a unit sphere.
"""
function _mesh_normalize_xyz(msh::GeometryBasics.AbstractMesh{3, Float32})::Float64
    # Iterators.flatten streams over all coordinate values without allocation
    return maximum(abs, Iterators.flatten(msh.position))
end

"""
    _mesh_normalize_xy(msh)

Return the maximum absolute coordinate value across x and y positions in `msh`, ignoring z.

Used to normalize head meshes that should scale to the x-y plane.
"""
function _mesh_normalize_xy(msh::GeometryBasics.AbstractMesh{3, Float32})::Float32
    return maximum(p -> max(abs(p[1]), abs(p[2])), msh.position)
end