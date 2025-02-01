function Sign_Detection(mesh::Mesh, grid::Grid, points::Matrix, ρₙ::Vector, ρₜ::Float64)
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
