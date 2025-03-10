"""
    slice_mesh_with_plane!(mesh::BlockMesh, plane::String, position::Float64; export_file::Union{String, Nothing}=nothing)

Slice a mesh with a specified plane, removing all nodes beyond the cutting plane.

Arguments:
- `mesh`: The BlockMesh object to be sliced
- `plane`: Cutting plane orientation, one of "xy", "yz", or "xz"
- `position`: Relative position of the cut from 0.0 (beginning) to 1.0 (end)
- `export_file`: Optional filename to export the sliced mesh

Returns:
- The modified BlockMesh object with nodes beyond the cutting plane removed
"""
function slice_mesh_with_plane!(mesh::BlockMesh, plane::String, position::Float64; export_file::Union{String, Nothing}=nothing)
    # Validate input parameters
    if !(plane in ["z", "x", "y"])
        error("Invalid plane. Choose from: \"xy\", \"yz\", \"xz\"")
    end
    
    if position < 0.0 || position > 1.0
        error("Position must be between 0 and 1")
    end
    
    # Determine the axis perpendicular to the cutting plane
    cut_axis = if plane == "z" 
        3  # Z-axis is perpendicular to XY plane
    elseif plane == "x" 
        1  # X-axis is perpendicular to YZ plane
    else  # plane == "y"
        2  # Y-axis is perpendicular to XZ plane
    end
    
    # Calculate min and max coordinates along the cutting axis
    min_val = minimum(p[cut_axis] for p in mesh.X)
    max_val = maximum(p[cut_axis] for p in mesh.X)
    
    # Calculate the actual position of the cutting plane
    cut_pos = min_val + position * (max_val - min_val)
    
    @info "Cutting mesh with $(plane) plane at position $(position) (coordinate: $(cut_pos))"
    
    # Count original elements and nodes for reporting
    orig_nodes = length(mesh.X)
    orig_elements = length(mesh.IEN)
    
    # Define a function to check if a point is beyond the cutting plane
    is_beyond_plane(p) = p[cut_axis] > cut_pos
    
    # Filter elements to keep only those with all nodes not beyond the cutting plane
    # This removes all tetrahedra that have at least one node beyond the plane
    mesh.IEN = [tet for tet in mesh.IEN if all(node_idx -> !is_beyond_plane(mesh.X[node_idx]), tet)]
    
    # Update mesh data structures to remove unused nodes and rebuild connectivity
    cleanup_unused_nodes!(mesh)  # Removes nodes that are no longer part of any element
    create_INE!(mesh)            # Rebuilds the inverse node-to-element connectivity
    
    # Count remaining elements and nodes
    remaining_nodes = length(mesh.X)
    remaining_elements = length(mesh.IEN)
    
    @info "Slice complete. Removed $(orig_nodes - remaining_nodes) nodes and $(orig_elements - remaining_elements) elements."
    
    # Export mesh if requested
    if export_file !== nothing
        export_mesh_vtk(mesh, export_file)
        @info "Exported sliced mesh to $(export_file)"
    end
    
    # return mesh
end

# slice_mesh_with_plane!(mesh, "x", 0.6, export_file="sliced_mesh.vtu")
