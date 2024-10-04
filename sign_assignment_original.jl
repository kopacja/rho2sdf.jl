function compute_sign_original(points, X, INE, IEN, ρₙ, ρₜ, ngp, sfce, nel)
  println("sign assignment by me")
  signs = -1 * ones(ngp)
  @threads for i in 1:ngp
    x = @view points[:, i]
    max_local = 10.0

    IDmin = argmin(mapslices(norm, (X .- x), dims=1))[2] # find ID of node (of mesh) that is closest to grid point x
    none = length(INE[IDmin]) # number of neighbour elements (that are connected to the node)
    ρₙₑ = ρₙ[IEN[:, INE[IDmin]]] # # nodal densities of none elements
    if maximum(ρₙₑ) < ρₜ # empty elements -> jump to another i
      continue
    end

    # cyklus přes všechny elementy, které jsou součástí uzlu X[IDmin]
    for j in 1:none
      el = INE[IDmin][j]

      Xₑ = @view X[:, IEN[:, el]] # coords of the element nodes
      # Compute the AABB
      min_bounds, max_bounds = compute_aabb(Xₑ)

      # Check if the point x is inside the AABB
      if is_point_inside_aabb(x, min_bounds, max_bounds)

        (Conv, local_coords) = find_local_coordinates(sfce, Xₑ, x)
        max_local_new = maximum(abs.(local_coords))
        
        if max_local_new < 1.2 && max_local > max_local_new

          H, _, _, _ = sfce(local_coords) # tvarové funkce a jejich derivace
          ρₑ = ρₙ[IEN[:, el]] # nodal densities for one element
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



