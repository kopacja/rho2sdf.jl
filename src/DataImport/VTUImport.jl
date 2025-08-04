# src/DataImport/VTUImport.jl

"""
    import_vtu_mesh(vtu_file::String) -> Tuple{Vector{Vector{Float64}}, Vector{Vector{Int64}}, Vector{Float64}}

Import mesh data from VTU file format and return data structures compatible with Rho2sdf.

# Arguments
- `vtu_file::String`: Path to the VTU file

# Returns
- `Tuple`: (X, IEN, rho) where:
  - `X::Vector{Vector{Float64}}`: Node coordinates, each element is [x, y, z]
  - `IEN::Vector{Vector{Int64}}`: Element connectivity, each element contains node IDs
  - `rho::Vector{Float64}`: Element densities/volume fractions

# Supported Elements
- Hexahedron (VTK_HEXAHEDRON = 12) - 8 nodes
- Tetrahedron (VTK_TETRA = 10) - 4 nodes

"""
function import_vtu_mesh(vtu_file::String)
    if !isfile(vtu_file)
        error("VTU file not found: $vtu_file")
    end
    
    print_info("Importing VTU mesh from: $vtu_file")
    
    # Read VTU file using ReadVTK
    vtk_file = ReadVTK.VTKFile(vtu_file)
    
    # Extract node coordinates
    points = ReadVTK.get_points(vtk_file)
    n_nodes = size(points, 2)
    
    # Convert to Vector{Vector{Float64}} format expected by Rho2sdf
    X = [vec(points[:, i]) for i in 1:n_nodes]
    
    print_data("Imported $n_nodes nodes")
    
    # Extract cells (elements)
    vtk_cells = ReadVTK.get_cells(vtk_file)
    connectivity = vtk_cells.connectivity
    offsets = vtk_cells.offsets
    types = vtk_cells.types
    
    n_elements = length(types)
    print_data("Processing $n_elements elements")
    
    # Process elements
    IEN = Vector{Vector{Int64}}()
    start_indices = vcat(1, offsets[1:end-1] .+ 1)
    
    # Count supported element types
    supported_elements = 0
    unsupported_elements = 0
    
    for i in 1:n_elements
        vtk_type = types[i]
        conn_indices = start_indices[i]:offsets[i]
        cell_connectivity = connectivity[conn_indices]
        
        # Convert to 1-based indexing (VTK uses 0-based)
        cell_nodes = cell_connectivity # .+ 1
        
        # Check element type and node count
        if vtk_type == 12  # VTK_HEXAHEDRON
            if length(cell_nodes) != 8
                @warn "Hexahedron element $i has $(length(cell_nodes)) nodes, expected 8"
                unsupported_elements += 1
                continue
            end
            # VTK hexahedron node ordering matches our expectations
            push!(IEN, collect(cell_nodes))
            supported_elements += 1
            
        elseif vtk_type == 10  # VTK_TETRA
            if length(cell_nodes) != 4
                @warn "Tetrahedron element $i has $(length(cell_nodes)) nodes, expected 4"
                unsupported_elements += 1
                continue
            end
            # Add tetrahedron connectivity
            push!(IEN, collect(cell_nodes))
            supported_elements += 1
            
        else
            @warn "Unsupported VTK element type: $vtk_type (element $i), skipping"
            unsupported_elements += 1
        end
    end
    
    if supported_elements == 0
        error("No supported elements found in VTU file. Supported types: Hexahedron (12), Tetrahedron (10)")
    end
    
    if unsupported_elements > 0
        print_warning("Skipped $unsupported_elements unsupported elements")
    end
    
    print_success("Imported $supported_elements supported elements")
    
    # Extract density data from cell data
    rho = extract_density_from_vtu(vtk_file, supported_elements)
    
    return (X, IEN, rho)
end

"""
    extract_density_from_vtu(vtk_file, n_elements::Int) -> Vector{Float64}

Extract density/volume fraction data from VTU cell data.
"""
function extract_density_from_vtu(vtk_file, n_elements::Int)
    # Try to extract cell data
    cell_data = nothing
    try
        cell_data = ReadVTK.get_cell_data(vtk_file)
    catch e
        @warn "Could not read cell data: $e"
        # Return default densities if no cell data available
        print_warning("No cell data found, using default density = 1.0 for all elements")
        return ones(Float64, n_elements)
    end
    
    if cell_data === nothing
        print_warning("No cell data available, using default density = 1.0")
        return ones(Float64, n_elements)
    end
    
    # Look for density data with common field names
    density_field_names = [
        "density", "Density", "DENSITY",
        "rho", "Rho", "RHO", 
        "volfrac", "VolFrac", "vol_frac", "VOLFRAC",
        "material_density", "element_density",
        "topology", "design_variable"
    ]
    
    density_data = nothing
    found_field = nothing
    
    for field_name in density_field_names
        if field_name in keys(cell_data)
            try
                density_array = ReadVTK.get_data(cell_data[field_name])
                density_data = collect(Float64, density_array)
                found_field = field_name
                break
            catch e
                @warn "Could not extract data from field '$field_name': $e"
                continue
            end
        end
    end
    
    if density_data === nothing
        # List available fields for debugging
        available_fields = collect(keys(cell_data))
        print_warning("No density data found in common fields.")
        print_data("Available cell data fields: $(join(available_fields, ", "))")
        
        # Try to use the first numeric field as density
        for field_name in available_fields
            try
                test_data = ReadVTK.get_data(cell_data[field_name])
                if eltype(test_data) <: Number
                    density_data = collect(Float64, test_data)
                    found_field = field_name
                    print_warning("Using field '$field_name' as density data")
                    break
                end
            catch
                continue
            end
        end
    end
    
    if density_data === nothing
        print_warning("No numeric cell data found, using default density = 1.0")
        return ones(Float64, n_elements)
    end
    
    # Validate density data
    if length(density_data) != n_elements
        print_warning("Density data length ($(length(density_data))) doesn't match number of elements ($n_elements)")
        
        if length(density_data) > n_elements
            # Truncate if too long
            density_data = density_data[1:n_elements]
            print_warning("Truncated density data to match element count")
        else
            # Pad with ones if too short
            append!(density_data, ones(Float64, n_elements - length(density_data)))
            print_warning("Padded density data with 1.0 values")
        end
    end
    
    # Check density range
    min_density = minimum(density_data)
    max_density = maximum(density_data)
    
    if min_density < 0.0 || max_density > 1.0
        print_warning("Density values outside [0,1] range: [$min_density, $max_density]")
        print_warning("Consider normalizing density values for optimal results")
    end
    
    print_success("Extracted density data from field '$found_field'")
    print_data("Density range: [$min_density, $max_density]")
    
    return density_data
end

"""
    validate_vtu_mesh(X, IEN, rho) -> Bool

Validate imported mesh data for consistency.
"""
function validate_vtu_mesh(X, IEN, rho)
    print_info("Validating imported mesh data...")
    
    n_nodes = length(X)
    n_elements = length(IEN)
    n_densities = length(rho)
    
    # Check basic counts
    if n_elements != n_densities
        print_error("Mismatch: $n_elements elements but $n_densities density values")
        return false
    end
    
    # Check node dimension consistency
    if !isempty(X)
        node_dims = [length(node) for node in X]
        if !all(d -> d == 3, node_dims)
            print_error("All nodes must have 3D coordinates")
            return false
        end
    end
    
    # Check element connectivity
    max_node_id = 0
    min_node_id = typemax(Int)
    
    for (i, element) in enumerate(IEN)
        if isempty(element)
            print_error("Element $i has no nodes")
            return false
        end
        
        for node_id in element
            if node_id < 1 || node_id > n_nodes
                print_error("Element $i references invalid node ID: $node_id (valid range: 1-$n_nodes)")
                return false
            end
            max_node_id = max(max_node_id, node_id)
            min_node_id = min(min_node_id, node_id)
        end
    end
    
    # Check if all nodes are referenced
    if max_node_id != n_nodes
        print_warning("Not all nodes are referenced by elements (max ID: $max_node_id, total nodes: $n_nodes)")
    end
    
    if min_node_id != 1
        print_warning("Node numbering doesn't start from 1 (min ID: $min_node_id)")
    end
    
    print_success("Mesh validation completed successfully")
    print_data("Summary: $n_nodes nodes, $n_elements elements")
    
    return true
end
