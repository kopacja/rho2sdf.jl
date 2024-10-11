function Sign_Detection(mesh::Mesh, grid::Grid, points::Matrix, ρₙ::Vector, ρₜ::Float64)
  X = mesh.X
  IEN = mesh.IEN
  INE = mesh.INE
  sfce = mesh.sfce

  ngp = grid.ngp
  signs = -1 * ones(ngp)

  p_nodes = Progress(ngp, 1, "Processing grid nodes: ", 30)
  counter_nodes = Atomic{Int}(0)
  update_interval_nodes = max(1, div(ngp, 100))

  @threads for i in 1:ngp
    x = @view points[:, i]
    max_local = 10.0
    IDmin = argmin(mapslices(norm, (X .- x), dims=1))[2] # Find ID of the node (in mesh) closest to grid point x
    none = length(INE[IDmin]) # Number of neighboring elements (connected to the node)
    ρₙₑ = ρₙ[IEN[:, INE[IDmin]]] # Nodal densities of none elements
    if maximum(ρₙₑ) < ρₜ # Empty elements -> skip to next i
      continue
    end
    # Loop through all elements that are part of the node X[:, IDmin]
    for j in 1:none
      el = INE[IDmin][j]
      Xₑ = @view X[:, IEN[:, el]] # Coordinates of the element nodes
      # Compute the Axis-Aligned Bounding Box (AABB)
      min_bounds, max_bounds = compute_aabb(Xₑ)
      # Check if the point x is inside the AABB
      if is_point_inside_aabb(x, min_bounds, max_bounds)
        (Conv, local_coords) = find_local_coordinates(sfce, Xₑ, x)
        max_local_new = maximum(abs.(local_coords))

        if max_local_new < 1.2 && max_local > max_local_new
          H, _, _, _ = sfce(local_coords) # Shape functions and their derivatives
          ρₑ = ρₙ[IEN[:, el]] # Nodal densities for one element
          ρ = H ⋅ ρₑ
          if ρ >= ρₜ
            signs[i] = 1.0
          end
          max_local = max_local_new
        end
      end
    end

    # Progress bar:
    count = atomic_add!(counter_nodes, 1)
    if count % update_interval_nodes == 0
      if Threads.threadid() == 1
        update!(p_nodes, count)
      end
    end
  end

  finish!(p_nodes)
  return signs
end
