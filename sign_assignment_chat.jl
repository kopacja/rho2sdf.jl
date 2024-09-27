function compute_sign_chat(points, X, INE, IEN, ρₙ, ρₜ, ngp, sfce, nel)
  println("sign assignment by chat")
  ngp = size(points, 2)
  signs = zeros(Float64, ngp)

  # Pre-compute AABBs for all elements
  element_aabbs = [compute_aabb(X[:, IEN[:, el]]) for el in 1:nel]

  # Use multi-threading for parallelization
  @threads for i in 1:ngp
    x = @view points[:, i]

    # Find candidate elements
    candidate_elements = [el for el in 1:nel if is_point_inside_aabb(x, element_aabbs[el]...)]

    for el in candidate_elements
      Xₑ = @view X[:, IEN[:, el]]
      Conv, local_coords = find_local_coordinates(sfce, Xₑ, x)
      if Conv && all(abs.(local_coords) .<= 1.0001)
        H, _, _, _ = sfce(local_coords)
        ρₑ = @view ρₙ[IEN[:, el]]
        ρ = dot(H, ρₑ)
        if ρ >= ρₜ
          signs[i] = 1.0
          break  # Přidáno: ukončení cyklu po nalezení vnitřního bodu
        end
      end
    end

  end

  return signs
end
