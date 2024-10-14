function compute_sign_original(points, X, INE, IEN, ρₙ, ρₜ, ngp, sfce, nel)
  signs = -1 * ones(ngp)
  # Pre-compute AABBs for all elements
  element_aabbs = [compute_aabb(X[:, IEN[:, el]]) for el in 1:nel]
  @threads for i in 1:ngp
    x = @view points[:, i]
    candidate_elements = [el for el in 1:nel if is_point_inside_aabb(x, element_aabbs[el]...)]
    max_local = 10.0
    none = length(candidate_elements) # Number of neighboring elements (connected to the node)
    ρₙₑ = ρₙ[IEN[:, candidate_elements]] # Nodal densities of none elements
    if isempty(ρₙₑ) || maximum(ρₙₑ) < ρₜ # Empty elements -> skip to next i
      continue
    end
    # Loop through all elements for which their AABB contains x
    for j in 1:none
      el = candidate_elements[j]
      Xₑ = @view X[:, IEN[:, el]] # Coordinates of the element nodes
      (Conv, local_coords) = find_local_coordinates(sfce, Xₑ, x) # Conv -> true/false (converged)
      max_local_new = maximum(abs.(local_coords))

      if max_local_new < 1.01 && max_local > max_local_new
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
  return signs
end
