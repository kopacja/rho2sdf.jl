"""
    create_test_cube_with_schlafli_tetrahedra() -> Tuple{Mesh{TET4}, Vector{Float64}}

Creates test cube geometry with linearly varying density using Schlafli scheme 
for tetrahedral discretization, returning a proper Mesh structure.

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

function create_test_cube_with_schlafli_tetrahedra()
    
    # Schlafli tetrahedral connectivity scheme
    schlafli_tet_connectivity = [
        [1, 2, 3, 7],  # Path 1: x, y, z
        [1, 6, 2, 7],  # Path 2: x, z, y
        [1, 3, 4, 7],  # Path 3: y, x, z
        [1, 4, 8, 7],  # Path 4: y, z, x
        [1, 5, 6, 7],  # Path 5: z, x, y
        [1, 8, 5, 7]   # Path 6: z, y, x
    ]
    
    # Cube parameters
    side_length = 10.0
    n_elements_per_side = 10
    element_size = side_length / n_elements_per_side
    half_side = side_length / 2.0
    
    # Node and element counts
    n_nodes_per_side = n_elements_per_side + 1  # 11
    total_nodes = n_nodes_per_side^3  # 1331
    total_hexes = n_elements_per_side^3  # 1000
    total_tetrahedra = total_hexes * 6  # 6000
    
    println("Creating test cube with Schlafli tetrahedra:")
    println("  - Edge length: $side_length")
    println("  - Hexahedral cells per edge: $n_elements_per_side") 
    println("  - Element size: $element_size")
    println("  - Total nodes: $total_nodes")
    println("  - Total hexahedral cells: $total_hexes")
    println("  - Total tetrahedra: $total_tetrahedra")
    
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
    
    # === NODAL DENSITY CALCULATION ===
    ρₙ = Vector{Float64}(undef, total_nodes)
    
    # Maximum distance from center to cube corner
    max_distance = sqrt(3) * half_side  # sqrt(5² + 5² + 5²) = 5√3
    
    for (idx, node_coords) in enumerate(X)
        # Distance of node from cube center (0,0,0)
        distance_from_center = norm(node_coords)
        
        # Linear density decrease with distance
        density = 1.0 - (distance_from_center / max_distance)
        
        # Ensure density is in range [0, 1]
        density = max(0.0, min(1.0, density))
        
        ρₙ[idx] = density
    end
    
    # === TETRAHEDRAL ELEMENT CREATION ===
    IEN = Vector{Vector{Int64}}()
    element_densities = Vector{Float64}()
    
    for k in 0:n_elements_per_side-1  # z direction
        for j in 0:n_elements_per_side-1  # y direction  
            for i in 0:n_elements_per_side-1  # x direction
                
                # Get the 8 nodes of the hexahedral cell
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
    
    # Create Mesh structure
    shape_func = coords -> shape_functions(TET4, coords)
    mesh = Mesh(X, IEN, element_densities, shape_func; element_type=TET4)
    
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
