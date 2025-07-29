"""
    create_test_cube_with_linear_density() -> Tuple{Mesh, Vector{Float64}}

Creates test cube geometry with linearly varying density.

# Geometry parameters:
- Cube centered at origin (0,0,0)
- Cube edge length: 10
- Discretization: 10×10×10 elements (each element has edge length 1)
- Total: 11×11×11 = 1331 nodes, 10×10×10 = 1000 elements

# Node density distribution:
- Cube center (0,0,0): density = 1.0
- Cube corners: density = 0.0  
- Other nodes: linear decrease with distance from center

# Returns:
- `Mesh`: Complete mesh structure
- `Vector{Float64}`: Nodal densities ρₙ
"""

function create_test_cube_with_linear_density()
        
    # Cube parameters
    side_length = 10.0
    n_elements_per_side = 10
    element_size = side_length / n_elements_per_side
    half_side = side_length / 2.0
    
    # Node and element counts
    n_nodes_per_side = n_elements_per_side + 1  # 11
    total_nodes = n_nodes_per_side^3  # 1331
    total_elements = n_elements_per_side^3  # 1000
    
    println("Creating test cube:")
    println("  - Edge length: $side_length")
    println("  - Elements per edge: $n_elements_per_side") 
    println("  - Element size: $element_size")
    println("  - Total nodes: $total_nodes")
    println("  - Total elements: $total_elements")
    
    # === NODE CREATION ===
    X = Vector{Vector{Float64}}()
    node_map = Dict{Tuple{Int,Int,Int}, Int}()
    
    node_id = 1
    for k in 0:n_elements_per_side  # z direction
        for j in 0:n_elements_per_side  # y direction
            for i in 0:n_elements_per_side  # x direction
                # Node coordinates (cube center is at origin)
                x = -half_side + i * element_size
                y = -half_side + j * element_size  
                z = -half_side + k * element_size
                
                push!(X, [x, y, z])
                node_map[(i, j, k)] = node_id
                node_id += 1
            end
        end
    end
    
    # === ELEMENT CREATION ===
    IEN = Vector{Vector{Int64}}()
    
    for k in 0:n_elements_per_side-1  # z direction
        for j in 0:n_elements_per_side-1  # y direction  
            for i in 0:n_elements_per_side-1  # x direction
                
                # Hexahedron nodes in correct order (compatible with hex8_shape)
                # Bottom face (z = k):
                n1 = node_map[(i,   j,   k)]      # (0,0,0) local
                n2 = node_map[(i+1, j,   k)]      # (1,0,0) local
                n3 = node_map[(i+1, j+1, k)]      # (1,1,0) local  
                n4 = node_map[(i,   j+1, k)]      # (0,1,0) local
                
                # Top face (z = k+1):
                n5 = node_map[(i,   j,   k+1)]    # (0,0,1) local
                n6 = node_map[(i+1, j,   k+1)]    # (1,0,1) local
                n7 = node_map[(i+1, j+1, k+1)]    # (1,1,1) local
                n8 = node_map[(i,   j+1, k+1)]    # (0,1,1) local
                
                # Element connectivity according to hex8 convention
                element_nodes = [n1, n2, n3, n4, n5, n6, n7, n8]
                push!(IEN, element_nodes)
            end
        end
    end
    
    # === NODAL DENSITY CALCULATION ===
    ρₙ = Vector{Float64}(undef, total_nodes)
    
    # Maximum distance from center to cube corner
    max_distance = sqrt(3) * half_side  # sqrt(5² + 5² + 5²) = 5√3
    
    for (idx, node_coords) in enumerate(X)
        # Distance of node from cube center (0,0,0)
        distance_from_center = norm(node_coords)
        
        # Linear density decrease with distance
        # density = 1 - (distance / max_distance)
        density = 1.0 - (distance_from_center / max_distance)
        
        # Ensure density is in range [0, 1]
        density = max(0.0, min(1.0, density))
        
        ρₙ[idx] = density
    end
    
    # === MESH STRUCTURE CREATION ===
    # Create element densities (required for Mesh constructor)
    # Use average of nodal densities for each element
    element_densities = Vector{Float64}(undef, total_elements)
    
    for (elem_idx, element_nodes) in enumerate(IEN)
        # Average density of element nodes
        avg_density = mean(ρₙ[node_id] for node_id in element_nodes)
        element_densities[elem_idx] = avg_density
    end
    
    # Create Mesh structure
    shape_func = coords -> shape_functions(HEX8, coords)
    mesh = Mesh(X, IEN, element_densities, shape_func; element_type=HEX8)
    
    # === STATISTICS OUTPUT ===
    println("\nDensity statistics:")
    println("  - Minimum density: $(round(minimum(ρₙ), digits=4))")
    println("  - Maximum density: $(round(maximum(ρₙ), digits=4))")  
    println("  - Average density: $(round(mean(ρₙ), digits=4))")
    
    println("\nMesh structure info:")
    println("  - Element type: $(mesh.element_type)")
    println("  - Total nodes: $(mesh.nnp)")
    println("  - Total elements: $(mesh.nel)")
    println("  - Nodes per element: $(mesh.nen)")
    println("  - Faces per element: $(mesh.nes)")

    return mesh, ρₙ
end
