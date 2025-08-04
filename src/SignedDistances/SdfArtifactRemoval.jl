# src/SignedDistances/SdfArtifactRemoval.jl
# SDF Artifact Removal using Union-Find Connected Components Analysis

"""
    UnionFind

Efficient Union-Find (Disjoint Set Union) data structure with path compression 
and union by rank for connected components analysis.
"""
mutable struct UnionFind
    parent::Vector{Int}
    rank::Vector{Int}
    size::Vector{Int}
    
    function UnionFind(n::Int)
        parent = collect(1:n)
        rank = zeros(Int, n)
        size = ones(Int, n)
        new(parent, rank, size)
    end
end

"""
    find_root!(uf::UnionFind, x::Int) -> Int

Find the root of element x with path compression for efficiency.
"""
function find_root!(uf::UnionFind, x::Int)
    if uf.parent[x] != x
        uf.parent[x] = find_root!(uf, uf.parent[x])  # Path compression
    end
    return uf.parent[x]
end

"""
    union!(uf::UnionFind, x::Int, y::Int) -> Bool

Union two components with union by rank optimization.
Returns true if union was performed, false if already in same component.
"""
function union!(uf::UnionFind, x::Int, y::Int)
    root_x = find_root!(uf, x)
    root_y = find_root!(uf, y)
    
    if root_x == root_y
        return false  # Already in same component
    end
    
    # Union by rank
    if uf.rank[root_x] < uf.rank[root_y]
        root_x, root_y = root_y, root_x
    end
    
    uf.parent[root_y] = root_x
    uf.size[root_x] += uf.size[root_y]
    
    if uf.rank[root_x] == uf.rank[root_y]
        uf.rank[root_x] += 1
    end
    
    return true
end

"""
    linear_to_3d(idx::Int, grid::Grid) -> Tuple{Int,Int,Int}

Convert linear index to 3D coordinates (i,j,k) for the grid.
"""
function linear_to_3d(idx::Int, grid::Grid)
    nx, ny, nz = grid.N .+ 1
    idx_0 = idx - 1  # Convert to 0-based indexing
    
    k = div(idx_0, nx * ny)
    remainder = idx_0 % (nx * ny)
    j = div(remainder, nx)
    i = remainder % nx
    
    return (i, j, k)
end

"""
    coord_3d_to_linear(i::Int, j::Int, k::Int, grid::Grid) -> Int

Convert 3D coordinates (i,j,k) to linear index for the grid.
"""
function coord_3d_to_linear(i::Int, j::Int, k::Int, grid::Grid)
    nx, ny, _ = grid.N .+ 1
    return k * nx * ny + j * nx + i + 1  # Convert to 1-based indexing
end

"""
    get_neighbors_3d(i::Int, j::Int, k::Int, grid::Grid) -> Vector{Tuple{Int,Int,Int}}

Get valid 3D neighbors (6-connectivity) for given coordinates.
"""
function get_neighbors_3d(i::Int, j::Int, k::Int, grid::Grid)
    nx, ny, nz = grid.N .+ 1
    neighbors = Tuple{Int,Int,Int}[]
    
    # 6-connectivity (face neighbors only)
    directions = [(-1,0,0), (1,0,0), (0,-1,0), (0,1,0), (0,0,-1), (0,0,1)]
    
    for (di, dj, dk) in directions
        ni, nj, nk = i + di, j + dj, k + dk
        if 0 <= ni < nx && 0 <= nj < ny && 0 <= nk < nz
            push!(neighbors, (ni, nj, nk))
        end
    end
    
    return neighbors
end

"""
    remove_sdf_artifacts!(sdf_values::Vector{Float64}, grid::Grid; 
                         threshold::Float64=0.0, min_component_ratio::Float64=0.01) -> Int

Remove small disconnected components from SDF field while preserving the main geometry.

# Arguments
- `sdf_values::Vector{Float64}`: SDF values (modified in-place)
- `grid::Grid`: Grid structure defining the spatial layout
- `threshold::Float64=0.0`: SDF threshold for interior detection (>= threshold is interior)
- `min_component_ratio::Float64=0.01`: Minimum size ratio for keeping components

# Returns
- `Int`: Number of nodes flipped from positive to negative

# Algorithm
1. Create binary mask for interior nodes (sdf >= threshold)
2. Run Union-Find connected components analysis  
3. Identify largest component
4. Flip sign of all nodes in smaller components (interior -> exterior)
"""
function remove_sdf_artifacts!(sdf_values::Vector{Float64}, grid::Grid; 
                              threshold::Float64=0.0, 
                              min_component_ratio::Float64=0.01)
    
    print_info("Removing SDF artifacts using Union-Find...")
    
    ngp = grid.ngp
    if length(sdf_values) != ngp
        error("SDF values length ($(length(sdf_values))) doesn't match grid points ($ngp)")
    end
    
    # Step 1: Create binary mask for interior nodes
    interior_mask = sdf_values .>= threshold
    interior_count = count(interior_mask)
    
    if interior_count == 0
        print_warning("No interior nodes found, nothing to process")
        return 0
    end
    
    println("Interior nodes detected: $interior_count / $ngp")
    
    # Step 2: Initialize Union-Find for all nodes
    uf = UnionFind(ngp)
    unions_performed = Atomic{Int}(0)
    
    # Step 3: Connect neighboring interior nodes
    @threads for linear_idx in 1:ngp
        if !interior_mask[linear_idx]
            continue  # Skip exterior nodes
        end
        
        i, j, k = linear_to_3d(linear_idx, grid)
        neighbors = get_neighbors_3d(i, j, k, grid)
        
        for (ni, nj, nk) in neighbors
            neighbor_idx = coord_3d_to_linear(ni, nj, nk, grid)
            
            # Connect only if neighbor is also interior
            if interior_mask[neighbor_idx]
                if union!(uf, linear_idx, neighbor_idx)
                    atomic_add!(unions_performed, 1)
                end
            end
        end
    end
    
    println("Union operations performed: $(unions_performed[])")
    
    # Step 4: Find component sizes and identify largest component
    component_sizes = Dict{Int, Int}()
    
    for linear_idx in 1:ngp
        if interior_mask[linear_idx]
            root = find_root!(uf, linear_idx)
            component_sizes[root] = get(component_sizes, root, 0) + 1
        end
    end
    
    if isempty(component_sizes)
        print_warning("No connected components found")
        return 0
    end
    
    largest_component_root = argmax(component_sizes)
    largest_component_size = component_sizes[largest_component_root]
    total_components = length(component_sizes)
    
    println("Connected components found: $total_components")
    println("Largest component size: $largest_component_size nodes")
    
    # Step 5: Determine minimum component size to keep
    min_component_size = max(1, round(Int, min_component_ratio * largest_component_size))
    println("Minimum component size to keep: $min_component_size nodes")
    
    # Step 6: Collect nodes to flip (small components)
    nodes_to_flip = Int[]
    small_components_count = 0
    small_components_total_size = 0
    
    for linear_idx in 1:ngp
        if interior_mask[linear_idx]
            root = find_root!(uf, linear_idx)
            component_size = component_sizes[root]
            
            # Keep largest component and components above minimum size
            if root != largest_component_root && component_size < min_component_size
                push!(nodes_to_flip, linear_idx)
                if component_sizes[root] > 0  # Count each component only once
                    small_components_count += 1
                    small_components_total_size += component_size
                    component_sizes[root] = 0  # Mark as counted
                end
            end
        end
    end
    
    # Step 7: Flip signs of nodes in small components
    if !isempty(nodes_to_flip)
        for idx in nodes_to_flip
            sdf_values[idx] = -abs(sdf_values[idx])  # Ensure negative (exterior)
        end
        
        println("Removed $small_components_count small components")
        println("Total nodes flipped: $(length(nodes_to_flip))")
        print_success("SDF artifacts successfully removed")
    else
        print_info("No small components found to remove")
    end
    
    return length(nodes_to_flip)
end

"""
    analyze_sdf_components(sdf_values::Vector{Float64}, grid::Grid; threshold::Float64=0.0)

Analyze connected components in SDF field without modifying the data.
Useful for understanding the structure before cleanup.

# Returns
- `Dict{Int, Int}`: Dictionary mapping component root to component size
"""
function analyze_sdf_components(sdf_values::Vector{Float64}, grid::Grid; threshold::Float64=0.0)
    print_info("Analyzing SDF components...")
    
    ngp = grid.ngp
    interior_mask = sdf_values .>= threshold
    interior_count = count(interior_mask)
    
    if interior_count == 0
        println("No interior nodes found")
        return Dict{Int, Int}()
    end
    
    uf = UnionFind(ngp)
    
    # Connect neighboring interior nodes
    for linear_idx in 1:ngp
        if !interior_mask[linear_idx]
            continue
        end
        
        i, j, k = linear_to_3d(linear_idx, grid)
        neighbors = get_neighbors_3d(i, j, k, grid)
        
        for (ni, nj, nk) in neighbors
            neighbor_idx = coord_3d_to_linear(ni, nj, nk, grid)
            if interior_mask[neighbor_idx]
                union!(uf, linear_idx, neighbor_idx)
            end
        end
    end
    
    # Calculate component sizes
    component_sizes = Dict{Int, Int}()
    for linear_idx in 1:ngp
        if interior_mask[linear_idx]
            root = find_root!(uf, linear_idx)
            component_sizes[root] = get(component_sizes, root, 0) + 1
        end
    end
    
    # Display analysis results
    if !isempty(component_sizes)
        sizes = collect(values(component_sizes))
        sort!(sizes, rev=true)
        
        println("Component analysis results:")
        println("  Total components: $(length(sizes))")
        println("  Largest component: $(sizes[1]) nodes")
        if length(sizes) > 1
            println("  Second largest: $(sizes[2]) nodes")
            println("  Smallest component: $(sizes[end]) nodes")
        end
    end
    
    return component_sizes
end
