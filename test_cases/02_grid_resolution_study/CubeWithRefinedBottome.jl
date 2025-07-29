"""
    create_test_cube_with_refined_bottom() -> Tuple{Mesh, Vector{Float64}}

Creates test cube geometry with refined bottom.
The bottom half has 2x refined elements, while the top half uses standard size.

# Geometry parameters:
- Cube centered at origin (0,0,0)
- Cube edge length: 10
- X,Z directions: 10 elements each (uniform)
- Y direction: 10 refined elements (bottom) + 5 coarse elements (top)

# Node density distribution:
- Cube center (0,0,0): density = 1.0
- Cube corners: density = 0.0  
- Other nodes: linear decrease with distance from center

# Returns:
- `Mesh`: Complete mesh structure
- `Vector{Float64}`: Nodal densities ρₙ
"""

function create_test_cube_with_refined_bottom()
    # Cube parameters
    side_length = 10.0
    half_side = side_length / 2.0
    
    # X and Z directions - unchanged
    n_elements_x = 10
    n_elements_z = 10
    element_size_x = 1.0
    element_size_z = 1.0
    
    # Y direction - variable refinement
    n_elements_y_bottom = 10  # Bottom half: 2x refined
    n_elements_y_top = 5      # Top half: original
    element_size_y_bottom = 0.5
    element_size_y_top = 1.0
    
    # Total counts
    n_nodes_x = n_elements_x + 1  # 11
    n_nodes_y = n_elements_y_bottom + n_elements_y_top + 1  # 16
    n_nodes_z = n_elements_z + 1  # 11
    
    total_nodes = n_nodes_x * n_nodes_y * n_nodes_z  # 1936
    total_elements = n_elements_x * (n_elements_y_bottom + n_elements_y_top) * n_elements_z  # 1500
    
    println("Creating refined test cube:")
    println("  - Total nodes: $total_nodes")
    println("  - Total elements: $total_elements")
    println("  - Y direction: $(n_elements_y_bottom) fine + $(n_elements_y_top) coarse elements")
    println("  - Bottom element size Y: $element_size_y_bottom")
    println("  - Top element size Y: $element_size_y_top")
    
    # === NODE CREATION ===
    X = Vector{Vector{Float64}}()
    node_map = Dict{Tuple{Int,Int,Int}, Int}()
    
    node_id = 1
    for k in 0:n_elements_z  # z direction
        for j in 0:(n_nodes_y-1)  # y direction
            for i in 0:n_elements_x  # x direction
                
                # X coordinate - uniform
                x = -half_side + i * element_size_x
                
                # Y coordinate - variable refinement
                if j <= n_elements_y_bottom  # Bottom part: j = 0 to 10
                    y = -half_side + j * element_size_y_bottom
                else  # Top part: j = 11 to 15
                    y = 0.0 + (j - n_elements_y_bottom) * element_size_y_top
                end
                
                # Z coordinate - uniform
                z = -half_side + k * element_size_z
                
                push!(X, [x, y, z])
                node_map[(i, j, k)] = node_id
                node_id += 1
            end
        end
    end
    
    # === ELEMENT CREATION ===
    IEN = Vector{Vector{Int64}}()
    
    for k in 0:(n_elements_z-1)  # z direction
        for j in 0:(n_elements_y_bottom + n_elements_y_top - 1)  # y direction: 0 to 14
            for i in 0:(n_elements_x-1)  # x direction
                
                # Hexahedron nodes in correct order
                n1 = node_map[(i,   j,   k)]
                n2 = node_map[(i+1, j,   k)]
                n3 = node_map[(i+1, j+1, k)]
                n4 = node_map[(i,   j+1, k)]
                n5 = node_map[(i,   j,   k+1)]
                n6 = node_map[(i+1, j,   k+1)]
                n7 = node_map[(i+1, j+1, k+1)]
                n8 = node_map[(i,   j+1, k+1)]
                
                element_nodes = [n1, n2, n3, n4, n5, n6, n7, n8]
                push!(IEN, element_nodes)
            end
        end
    end
    
    # === NODAL DENSITY CALCULATION ===
    ρₙ = Vector{Float64}(undef, total_nodes)
    max_distance = sqrt(3) * half_side
    
    for (idx, node_coords) in enumerate(X)
        distance_from_center = norm(node_coords)
        density = 1.0 - (distance_from_center / max_distance)
        density = max(0.0, min(1.0, density))
        ρₙ[idx] = density
    end
    
    # === ELEMENT DENSITIES ===
    element_densities = Vector{Float64}(undef, total_elements)
    
    for (elem_idx, element_nodes) in enumerate(IEN)
        avg_density = mean(ρₙ[node_id] for node_id in element_nodes)
        element_densities[elem_idx] = avg_density
    end
    
    # Create Mesh structure
    shape_func = coords -> shape_functions(HEX8, coords)
    mesh = Mesh(X, IEN, element_densities, shape_func; element_type=HEX8)
   
    # === VERIFICATION ===
    println("\nMesh verification:")
    println("  - Bottom elements (j=0-9): refined mesh with element size $(element_size_y_bottom)")
    println("  - Top elements (j=10-14): coarse mesh with element size $(element_size_y_top)")
    
    # Verify Y coordinates at transition
    transition_node_bottom = node_map[(0, n_elements_y_bottom, 0)]
    transition_node_top = node_map[(0, n_elements_y_bottom+1, 0)]
    
    println("  - Y coordinate at transition (bottom): $(X[transition_node_bottom][2])")
    println("  - Y coordinate at transition (top): $(X[transition_node_top][2])")
    
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
