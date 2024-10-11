function compute_sign_original(points, X, INE, IEN, ρₙ, ρₜ, ngp, sfce, nel)
  signs = -1 * ones(ngp)
  @threads for i in 1:ngp
    x = @view points[:, i]
    max_local = 10.0
    IDmin = argmin(mapslices(norm, (X .- x), dims=1))[2] # Find ID of the node (in mesh) closest to grid point x
    none = length(INE[IDmin]) # Number of neighboring elements (connected to the node)
    ρₙₑ = ρₙ[IEN[:, INE[IDmin]]] # Nodal densities of none elements
    if maximum(ρₙₑ) < ρₜ # Empty elements -> skip to next i
      continue
    end
    # Loop through all elements that are part of the node X[IDmin]
    for j in 1:none
      el = INE[IDmin][j]
      Xₑ = @view X[:, IEN[:, el]] # Coordinates of the element nodes
      # Compute the Axis-Aligned Bounding Box (AABB)
      min_bounds, max_bounds = compute_aabb(Xₑ)
      # Check if the point x is inside the AABB
      if is_point_inside_aabb(x, min_bounds, max_bounds)
        (Conv, local_coords) = find_local_coordinates(sfce, Xₑ, x) # Conv -> true/false (converged)
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
  end
  return signs
end
