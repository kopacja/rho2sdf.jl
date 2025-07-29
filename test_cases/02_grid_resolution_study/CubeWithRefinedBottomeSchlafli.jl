"""
    create_test_cube_with_refined_bottom_schlafli() -> Tuple{Mesh{TET4}, Vector{Float64}}

Creates test cube geometry with refined bottom using Schlafli tetrahedral discretization.
The bottom half has 2x refined elements, while the top half uses standard size.
Each hexahedral cell is converted to 6 tetrahedra using Schlafli orthoscheme.

# Geometry parameters:
- Cube centered at origin (0,0,0)
- Cube edge length: 10
- X,Z directions: 10 elements each (uniform)
- Y direction: 10 refined elements (bottom) + 5 coarse elements (top)
- Each hex cell → 6 tetrahedra via Schlafli scheme

# Returns:
- `Mesh{TET4}`: Complete tetrahedral mesh structure
- `Vector{Float64}`: Nodal densities ρₙ
"""

function create_test_cube_with_refined_bottom_schlafli()
    
    # Schlafli tetrahedral connectivity scheme
    schlafli_tet_connectivity = [
        [1, 2, 3, 7],  # Path 1: x, y, z
        [1, 6, 2, 7],  # Path 2: x, z, y
        [1, 3, 4, 7],  # Path 3: y, x, z
        [1, 4, 8, 7],  # Path 4: y, z, x
        [1, 5, 6, 7],  # Path 5: z, x, y
        [1, 8, 5, 7]   # Path 6: z, y, x
    ]
    
    # Cube parameters (same as refined bottom)
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
    total_hex_elements = n_elements_x * (n_elements_y_bottom + n_elements_y_top) * n_elements_z  # 1500
    total_tetrahedra = total_hex_elements * 6  # 9000
    
    println("Creating refined test cube with Schlafli tetrahedra:")
    println("  - Total nodes: $total_nodes")
    println("  - Total hexahedral cells: $total_hex_elements")
    println("  - Total tetrahedra: $total_tetrahedra")
    println("  - Y direction: $(n_elements_y_bottom) fine + $(n_elements_y_top) coarse elements")
    println("  - Bottom element size Y: $element_size_y_bottom")
    println("  - Top element size Y: $element_size_y_top")
    
    # === NODE CREATION (same as refined bottom) ===
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
    
    # === NODAL DENSITY CALCULATION ===
    ρₙ = Vector{Float64}(undef, total_nodes)
    max_distance = sqrt(3) * half_side
    
    for (idx, node_coords) in enumerate(X)
        distance_from_center = norm(node_coords)
        density = 1.0 - (distance_from_center / max_distance)
        density = max(0.0, min(1.0, density))
        ρₙ[idx] = density
    end
    
    # === TETRAHEDRAL ELEMENT CREATION ===
    IEN = Vector{Vector{Int64}}()
    element_densities = Vector{Float64}()
    
    for k in 0:(n_elements_z-1)  # z direction
        for j in 0:(n_elements_y_bottom + n_elements_y_top - 1)  # y direction: 0 to 14
            for i in 0:(n_elements_x-1)  # x direction
                
                # Get the 8 nodes of the hexahedral cell (same order as refined bottom)
                n1 = node_map[(i,   j,   k)]      # (0,0,0) local
                n2 = node_map[(i+1, j,   k)]      # (1,0,0) local
                n3 = node_map[(i+1, j+1, k)]      # (1,1,0) local  
                n4 = node_map[(i,   j+1, k)]      # (0,1,0) local
                n5 = node_map[(i,   j,   k+1)]    # (0,0,1) local
                n6 = node_map[(i+1, j,   k+1)]    # (1,0,1) local
                n7 = node_map[(i+1, j+1, k+1)]    # (1,1,1) local
                n8 = node_map[(i,   j+1, k+1)]    # (0,1,1) local
                
                # Hexahedral nodes in order expected by Schlafli scheme
                hex_nodes = [n1, n2, n3, n4, n5, n6, n7, n8]
                
                # Create 6 tetrahedra using Schlafli orthoscheme
                for tet_connectivity in schlafli_tet_connectivity
                    # Create tetrahedral element
                    tet_nodes = [hex_nodes[idx] for idx in tet_connectivity]
                    push!(IEN, tet_nodes)
                    
                    # Calculate tetrahedral element density as average of its nodes
                    tet_avg_density = mean(ρₙ[node_id] for node_id in tet_nodes)
                    push!(element_densities, tet_avg_density)
                end
            end
        end
    end
    
    # Create Mesh structure with TET4 elements
    shape_func = coords -> shape_functions(TET4, coords)
    mesh = Mesh(X, IEN, element_densities, shape_func; element_type=TET4)
    
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
    println("  - Minimum nodal density: $(round(minimum(ρₙ), digits=4))")
    println("  - Maximum nodal density: $(round(maximum(ρₙ), digits=4))")  
    println("  - Average nodal density: $(round(mean(ρₙ), digits=4))")
    
    println("\nMesh structure info:")
    println("  - Element type: $(mesh.element_type)")
    println("  - Total nodes: $(mesh.nnp)")
    println("  - Total elements: $(mesh.nel)")
    println("  - Nodes per element: $(mesh.nen)")
    println("  - Faces per element: $(mesh.nes)")
    
    return mesh, ρₙ
end
