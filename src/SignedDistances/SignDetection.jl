"""
    Sign_Detection_HEX8(mesh::Mesh{HEX8}, grid::Grid, points::Matrix, ρₙ::Vector, ρₜ::Float64) -> Vector{Float64}

Original Sign_Detection implementation specifically for HEX8 elements.
"""
function Sign_Detection_HEX8(mesh::Mesh{HEX8}, grid::Grid, points::Matrix, ρₙ::Vector{Float64}, ρₜ::Float64)
  X = mesh.X
  IEN = mesh.IEN
  sfce = mesh.sfce
  nel = mesh.nel

  ngp = grid.ngp
  signs = -1 * ones(ngp)

  # Pre-compute AABBs for all elements
  element_aabbs = Vector{NTuple{2,Vector{Float64}}}(undef, nel)
  @threads for el in 1:nel
    aabb = compute_aabb(@view X[:, IEN[:, el]])
    # Převedeme výstup z compute_aabb na požadovaný typ
    element_aabbs[el] = (vec(aabb[1]), vec(aabb[2]))
  end

  p_nodes = Progress(ngp, 1, "Processing grid nodes: ", 30)
  counter_nodes = Atomic{Int}(0)
  update_interval_nodes = max(1, div(ngp, 100))

  @threads for i in 1:ngp
    x = @view points[:, i]
    # Find potential elements that contain the point using AABB (Axis-Aligned Bounding Box) check
    candidate_elements = [el for el in 1:nel if is_point_inside_aabb(x, element_aabbs[el]...)]
    max_local = 10.0
    none = length(candidate_elements)
    ρₙₑ = ρₙ[IEN[:, candidate_elements]] # Get nodal densities for candidate elements

    # Skip if no elements or maximum density is below threshold
    if isempty(ρₙₑ) || maximum(ρₙₑ) < ρₜ
      continue
    end

    # Modified element loop with early exit condition
    for j in 1:none
      el = candidate_elements[j]
      Xₑ = @view X[:, IEN[:, el]]
      (_, local_coords) = find_local_coordinates(Xₑ, x)
      max_local_new = maximum(abs.(local_coords))

      # Check if point is inside element (with small tolerance)
      if max_local_new < 1.01 && max_local > max_local_new
        # If local coordinates indicate point is well inside element (< 0.95)
        # we can process this element and exit early
        if max_local_new < 0.95
          H = sfce(local_coords)  # Only need shape functions here
          ρₑ = ρₙ[IEN[:, el]]
          ρ = H ⋅ ρₑ
          if ρ >= ρₜ
            signs[i] = 1.0
          end
          break  # Exit the element loop early
        end

        # Otherwise process element normally
        H = sfce(local_coords)
        ρₑ = ρₙ[IEN[:, el]]
        ρ = H ⋅ ρₑ
        if ρ >= ρₜ
          signs[i] = 1.0
        end
        max_local = max_local_new
      end
    end

    # Update progress bar
    count = atomic_add!(counter_nodes, 1)
    if count % update_interval_nodes == 0 && Threads.threadid() == 1
      update!(p_nodes, count)
    end
  end

  finish!(p_nodes)
  return signs
end

"""
    Sign_Detection_TET4(mesh::Mesh{TET4}, grid::Grid, points::Matrix, ρₙ::Vector, ρₜ::Float64) -> Vector{Float64}

Determine sign (+1/-1) for each grid point based on tetrahedral mesh density field.
"""
function Sign_Detection_TET4(mesh::Mesh{TET4}, grid::Grid, points::Matrix, ρₙ::Vector{Float64}, ρₜ::Float64)
    X = mesh.X
    IEN = mesh.IEN
    nel = mesh.nel
    ngp = grid.ngp
    signs = -1.0 * ones(ngp)
    
    # Extract grid dimensions
    grid_dims = grid.N .+ 1
    grid_min = grid.AABB_min
    cell_size = grid.cell_size
    
    # Create spatial acceleration structure
    println("Building spatial acceleration structure for TET4 elements...")
    grid_tetrahedra = create_grid_tetrahedra_mapping_TET4(mesh, grid, grid_dims)
    
    # Process points in parallel
    p_nodes = Progress(ngp, 1, "Computing signs for TET4: ", 30)
    counter_nodes = Atomic{Int}(0)
    update_interval_nodes = max(1, div(ngp, 100))
    
    batch_size = 64
    num_batches = ceil(Int, ngp / batch_size)
    
    @threads for batch in 1:num_batches
        start_idx = (batch - 1) * batch_size + 1
        end_idx = min(batch * batch_size, ngp)
        
        for i in start_idx:end_idx
            x = @view points[:, i]
            grid_idx = point_to_grid_index(x, grid_min, cell_size, grid_dims)
            
            if any(grid_idx .< 1) || any(grid_idx .> grid_dims)
                continue
            end
            
            cell_tetrahedra = grid_tetrahedra[grid_idx...]
            
            # Check each candidate tetrahedron
            for el in cell_tetrahedra
                # Get tetrahedron vertices
                tetrahedron = @view X[:, IEN[:, el]]
                
                # Check if point is inside tetrahedron
                if is_point_in_tetrahedron(tetrahedron, x)
                    # Calculate density at this point using TET4 shape functions
                    (found, local_coords) = find_local_coordinates(tetrahedron, x, TET4)
                    
                    if found
                        # Get shape functions for TET4
                        N = shape_functions(TET4, local_coords)
                        
                        # Interpolate density
                        ρₑ = ρₙ[IEN[:, el]]
                        ρ = dot(N, ρₑ)
                        
                        # Set sign based on density threshold
                        if ρ >= ρₜ
                            signs[i] = 1.0
                            break
                        end
                    end
                end
            end
            
            # Update progress
            if i == end_idx
                count = atomic_add!(counter_nodes, end_idx - start_idx + 1)
                if Threads.threadid() == 1 && (count % update_interval_nodes <= batch_size)
                    update!(p_nodes, min(count, ngp))
                end
            end
        end
    end
    
    finish!(p_nodes)
    return signs
end

# Helper function specifically for TET4 mapping
function create_grid_tetrahedra_mapping_TET4(mesh::Mesh{TET4}, grid::Grid, grid_dims)
    # Same implementation as in stl2sdf but adapted for Mesh{TET4}
    X = mesh.X
    IEN = mesh.IEN
    nel = mesh.nel
    grid_min = grid.AABB_min
    cell_size = grid.cell_size
    
    num_threads = Threads.nthreads()
    
    local_grids = [
        [Vector{Int}() for _ in 1:grid_dims[1], _ in 1:grid_dims[2], _ in 1:grid_dims[3]]
        for _ in 1:num_threads
    ]
    
    @threads for el in 1:nel
        tid = Threads.threadid()
        
        tet_vertices = @view X[:, IEN[:, el]]
        
        min_bounds = vec(minimum(tet_vertices, dims=2))
        max_bounds = vec(maximum(tet_vertices, dims=2))
        
        min_idx = max.(1, floor.(Int, (min_bounds .- grid_min) ./ cell_size) .- 1)
        max_idx = min.(grid_dims, ceil.(Int, (max_bounds .- grid_min) ./ cell_size) .+ 1)
        
        for i in min_idx[1]:max_idx[1]
            for j in min_idx[2]:max_idx[2]
                for k in min_idx[3]:max_idx[3]
                    push!(local_grids[tid][i, j, k], el)
                end
            end
        end
    end
    
    # Merge thread-local results
    grid_tetrahedra = [Vector{Int}() for _ in 1:grid_dims[1], _ in 1:grid_dims[2], _ in 1:grid_dims[3]]
    
    for i in 1:grid_dims[1]
        for j in 1:grid_dims[2]
            for k in 1:grid_dims[3]
                for tid in 1:num_threads
                    append!(grid_tetrahedra[i, j, k], local_grids[tid][i, j, k])
                end
            end
        end
    end
    
    return grid_tetrahedra
end

# Reuse from stl2sdf
function is_point_in_tetrahedron(tetrahedron::AbstractMatrix{Float64}, point::AbstractVector{Float64}, tolerance::Float64=1e-10)
    # Quick AABB test
    min_bounds = minimum(tetrahedron, dims=2)
    max_bounds = maximum(tetrahedron, dims=2)
    
    if any(point .< min_bounds .- tolerance) || any(point .> max_bounds .+ tolerance)
        return false
    end
    
    # Get vertices
    v1 = tetrahedron[:, 1]
    v2 = tetrahedron[:, 2]
    v3 = tetrahedron[:, 3] 
    v4 = tetrahedron[:, 4]
    
    # Compute barycentric coordinates
    T = [v1 v2 v3 v4; 1 1 1 1]
    b = [point; 1.0]
    
    λ = T \ b
    
    return all(λ .>= -tolerance) && all(λ .<= 1.0 + tolerance)
end

"""
    point_to_grid_index(point, grid_min, cell_size, grid_dims) -> Vector{Int}

Convert a 3D point to grid indices.

# Arguments
- `point`: The 3D point coordinates
- `grid_min`: Minimum coordinates of the grid
- `cell_size`: Size of grid cells (scalar or vector)
- `grid_dims`: Dimensions of the grid

# Returns
- `Vector{Int}`: 3D index of the grid cell containing the point
"""
function point_to_grid_index(point, grid_min, cell_size, grid_dims)
    # Calculate cell indices
    if isa(cell_size, Number)
        idx = floor.(Int, (point .- grid_min) ./ cell_size) .+ 1
    else
        idx = floor.(Int, (point .- grid_min) ./ cell_size) .+ 1
    end
    
    # Ensure indices are within bounds
    return max.(1, min.(grid_dims, idx))
end

"""
    Sign_Detection(mesh::Mesh, grid::Grid, points::Matrix, ρₙ::Vector, ρₜ::Float64) -> Vector{Float64}

Main dispatch function that calls appropriate implementation based on element type.
"""
function Sign_Detection(mesh::Mesh{T}, grid::Grid, points::Matrix, ρₙ::Vector{Float64}, ρₜ::Float64) where {T<:AbstractElement}
    if T == HEX8
        return Sign_Detection_HEX8(mesh, grid, points, ρₙ, ρₜ)
    elseif T == TET4
        return Sign_Detection_TET4(mesh, grid, points, ρₙ, ρₜ)
    else
        error("Sign_Detection not implemented for element type: $T")
    end
end
