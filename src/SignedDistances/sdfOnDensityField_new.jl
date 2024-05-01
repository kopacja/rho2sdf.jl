
# Reduction of the matrix K (eigenvectors) if the eigenvalues are close to zero:
function ReduceEigenvals(
  K::Matrix{Float64},
  r::Vector{Float64},
  Sign::Int,
  th::Float64=1.0e-6)

  Λ = real.(eigvals(K))
  Λ_min = minimum(abs.(Λ))

  if Λ_min < th
    Φ = real.(eigvecs(K))
    # Adjust the calculation if necessary
    idx_below_th = findall(x -> abs(x) < th, Λ)
    idx = setdiff(1:length(Λ), (idx_below_th))
    Φ_reduced = Φ[:, idx]
    ΔΞ̃_and_Δλ̃ = 1.0 ./ Λ[idx] .* (Φ_reduced' * r)
    ΔΞ_and_Δλ = (Φ_reduced .* Sign) * ΔΞ̃_and_Δλ̃

  else
    ΔΞ_and_Δλ = K \ (r .* Sign)
  end

  return ΔΞ_and_Δλ, Λ_min                         #??# Λ_min is not necessary - only for debug
end

# Selection of regular grid points that have been projected:
function SelectProjectedNodes(
  mesh::Mesh,
  grid::Grid,
  xp::Matrix{Float64},
  points::Matrix{Float64})

  ngp = grid.ngp # number of nodes in grid
  nsd = mesh.nsd # number of spacial dimensions

  # Assuming ngp is defined somewhere in your code
  # Preallocate arrays with maximum possible size
  max_size = ngp * 2  # Adjust this based on your knowledge of the data
  X = [zeros(Float64, nsd) for _ in 1:max_size]
  Xp = [zeros(Float64, nsd) for _ in 1:max_size]

  count = 0
  for i = 1:ngp
    if sum(abs.(xp[:, i])) > 1.0e-10
      count += 1
      X[count] = points[:, i]
      Xp[count] = xp[:, i]
    end
  end

  # If count is 0, indicating no points were added, handle gracefully
  if count == 0
    println("WARNING: no projected points!")
    return [], [], NaN, NaN
  end

  # Trim the unused preallocated space
  X = resize!(X, count)
  Xp = resize!(Xp, count)

  # Mean and max projected distance:
  mean_PD = mean(norm.(X - Xp))
  max_PD = maximum(norm.(X - Xp))

  return X, Xp, mean_PD, max_PD
end

# Correction of the sign of the distance function
# Prevention of dual contour. Out of the material -> negative sign
function SignCorrection4SDF(dist::Vector{Float64},
  grid::Grid,
  big::Float64)
  ngp = grid.ngp # number of nodes in grid

  Sign = -1
  for i in 1:ngp # For each grid point
    if dist[i] != big
      Sign = sign(dist[i])
    end
    if dist[i] == big && Sign == 1
      dist[i] = dist[i] * -1
    end
  end
  return dist
end

# Compute rho normal for node inside the element based on ρₑ
function RhoNorm(
  ρₑ::Vector{Float64},
  Ξ::Vector{Float64},
)

  (dρ_dΞ, d²ρ_dΞ², d³ρ_dΞ³) = ρ_derivatives(ρₑ, Ξ)
  #TODO: Can be optimized (d²ρ_dΞ², d³ρ_dΞ³ are not used)
  norm_dρ_dΞ = norm(dρ_dΞ)
  n = dρ_dΞ / norm_dρ_dΞ

  return n
end

# Update distance filed if the distance is smaller then the previous one
function WriteValue(
  dist_tmp::Float64,
  dist::Vector{Float64},
  xp::Matrix{Float64},
  xₚ::Vector{Float64},
  v::Int)

  if (abs(dist_tmp) < abs(dist[v]))
    dist[v] = dist_tmp
    xp[:, v] = xₚ
  end
  return dist, xp
end

function compute_coords(
  x::Vector,
  ρₜ::Float64,
  Xₑ::Matrix,
  ρₑ::Vector,
)
  model = Model(Ipopt.Optimizer)
  set_silent(model)
  # unset_silent(model)

  set_optimizer_attribute(model, "tol", 1e-6)
  set_optimizer_attribute(model, "max_iter", 50)
  set_optimizer_attribute(model, "acceptable_tol", 1e-6)

  @variable(model, ξ₁, lower_bound = -1.0, upper_bound = 1.0)
  @variable(model, ξ₂, lower_bound = -1.0, upper_bound = 1.0)
  @variable(model, ξ₃, lower_bound = -1.0, upper_bound = 1.0)
  # set_start_value(ξ₁, 0.1)
  # set_start_value(ξ₂, 0.1)
  # set_start_value(ξ₃, 0.1)

  N8 = [
    -1 / 8 * (ξ₁ - 1) * (ξ₂ - 1) * (ξ₃ - 1),
    1 / 8 * (ξ₁ + 1) * (ξ₂ - 1) * (ξ₃ - 1),
    -1 / 8 * (ξ₁ + 1) * (ξ₂ + 1) * (ξ₃ - 1),
    1 / 8 * (ξ₁ - 1) * (ξ₂ + 1) * (ξ₃ - 1),
    1 / 8 * (ξ₁ - 1) * (ξ₂ - 1) * (ξ₃ + 1),
    -1 / 8 * (ξ₁ + 1) * (ξ₂ - 1) * (ξ₃ + 1),
    1 / 8 * (ξ₁ + 1) * (ξ₂ + 1) * (ξ₃ + 1),
    -1 / 8 * (ξ₁ - 1) * (ξ₂ + 1) * (ξ₃ + 1)
  ]
  # @NLobjective(model, Min, sqrt(sum((x[i] - sum(Xₑ[i,k] * N8[k] for k in 1:length(N8)))^2 for i in 1:length(x))))
  @NLobjective(model, Min, sum((x[i] - sum(Xₑ[i, k] * N8[k] for k in 1:length(N8)))^2 for i in 1:length(x)))
  @NLconstraint(model, sum(ρₑ[k] * N8[k] for k in 1:length(N8)) == ρₜ)
  optimize!(model)

  return value.([ξ₁, ξ₂, ξ₃])
end


# MAIN FUNCTION for eval SDF
function evalSignedDistances(
  mesh::Mesh,
  grid::Grid,
  ρₙ::Vector{Float64},
  ρₜ::Float64,
)

  points = MeshGrid.generateGridPoints(grid) # uzly pravidelné mřížky
  linkedList = MeshGrid.LinkedList(grid, points) # pro rychlé vyhledávání

  head = linkedList.head # ID pravidelné bunky (pozice), index bodu z points
  next = linkedList.next # vec délky points, další uzly pro danou bunku, když -1 tak už další není
  N = linkedList.grid.N # Number of divisions along each axis of the grid
  AABB_min = linkedList.grid.AABB_min # Minimum coordinates of the Axis-Aligned Bounding Box (AABB)
  AABB_max = linkedList.grid.AABB_max # Maximum coordinates of the AABB
  δ = 2.5 * grid.cell_size # offset for mini AABB

  X = mesh.X   # vector of nodes positions
  IEN = mesh.IEN # ID element -> ID nodes
  INE = mesh.INE # ID node -> ID elements
  ISN = mesh.ISN # connectivity face - edges
  sfce = mesh.sfce # shape function handler
  nsd = mesh.nsd # number of spacial dimensions
  nel = mesh.nel # number of all elements
  nes = mesh.nes # number of element segments (faces)
  nsn = mesh.nsn # number of face nodes
  println("number of all elements: ", nel)

  ngp = grid.ngp # number of nodes in grid
  # big = -1.0e10
  big = 1.0e10
  dist = big * ones(ngp) # distance field initialization
  xp = zeros(nsd, ngp) # souřadnice bodů vrcholů (3xngp)

  # Tri mesh for pseudonormals:
  tri_mesh = Rho2sdf.extractSurfaceTriangularMesh(mesh)
  EN = NodePosition3D(tri_mesh)
  VPN, EPN = computePseudoNormals(tri_mesh)
  println("VPN", size(VPN))
  println("EPN", size(EPN))

  for el = 1:nel

    println("element ID: ", el)
    ρₑ = ρₙ[IEN[:, el]] # nodal densities for one element

    ρₑ_min = minimum(ρₑ)
    ρₑ_max = maximum(ρₑ)
    if (ρₑ_min >= ρₜ) # the boundary does not cross through the element
      commonEls = []

      # cycle through element faces (6)
      for sg = 1:nes
        commonEls = INE[IEN[mesh.ISN[sg][1], el]]
        for a = 2:nsn
          idx = findall(in(INE[IEN[ISN[sg][a], el]]), commonEls) # for how many elements does this face belong ?
          commonEls = commonEls[idx]
        end

        if (length(commonEls) == 1) # = is a part of the outer boundary of the body
          Xs = X[:, IEN[ISN[sg], el]]
          Xc = vec(mean(Xs, dims=2))

          for a = 1:nsn # cycle through number of all nodals belong to face

            # coordinates of nodes of the triangle
            x₁ = Xs[:, a]
            x₂ = Xs[:, (a%nsn)+1]
            x₃ = Xc
            # println(typeof(x₁))
            # exit()

            # coordinates of the vertices of the triangle
            Xt = [x₁ x₂ x₃]

            # finding coresponding triangle:
            ID_tri = find_triangle_position(EN, [x₁ x₂ x₃])

            #NOTE: From this part it is same as in sdfOnTriangularMesh ->

            # Triangle edges
            Et = calculate_triangle_edges(Xt)

            n = cross(Et[1], Et[2]) # norm of triangle
            n = n / norm(n) # unit norm

            # Nodes of mini AABB grid:
            Is = MeshGrid.calculateMiniAABB_grid(Xt, δ, N, AABB_min, AABB_max, nsd)

            for I ∈ Is # cycle through the nodes of the mini AABB grid
              ii = Int( # node ID
                I[3] * (N[1] + 1) * (N[2] + 1) + I[2] * (N[1] + 1) + I[1] + 1,
              )
              v = head[ii]
              while v != -1
                x = points[:, v]
                λ = barycentricCoordinates(x₁, x₂, x₃, n, x)

                xₚ = zeros(nsd) # projection

                isFaceOrEdge = false # projection check

                if (minimum(λ) >= 0.0) # xₚ is in the triangle, projection node x inside triangle 
                  xₚ = λ[1] * x₁ + λ[2] * x₂ + λ[3] * x₃
                  dist_tmp = dot(x - xₚ, n)

                  isFaceOrEdge = update_distance!(dist, dist_tmp, v, xp, xₚ, isFaceOrEdge)
                else

                  # Edges of the triangle:
                  for j = 1:3
                    L = norm(Et[j]) # length of j triangle edge
                    xᵥ = Xt[:, j]
                    P = dot(x - xᵥ, Et[j] / L) # skalar product of vector (vertex&node) and norm edge
                    if (P >= 0 && P <= L) # is the perpendicular projection of a node onto an edge in the edge interval?
                      xₚ = xᵥ + (Et[j] / L) * P
                      n_edge = n
                      n_edge = EPN[ID_tri][j]
                      dist_tmp = sign(dot(x - xₚ, n_edge)) * norm(x - xₚ)

                      isFaceOrEdge = update_distance!(dist, dist_tmp, v, xp, xₚ, isFaceOrEdge)
                    end
                  end
                end
                # Remaining cases:
                if (isFaceOrEdge == false)
                  dist_tmp, idx =
                    findmin([norm(x - x₁), norm(x - x₂), norm(x - x₃)]) # which node of the triangle is closer?
                  xₚ = Xt[:, idx] # the node of triangle
                  # n_vertex = n
                  n_vertex = VPN[tri_mesh.IEN[idx, ID_tri]]
                  dist_tmp = dist_tmp * sign(dot(x - xₚ, n_vertex))

                  isFaceOrEdge = update_distance!(dist, dist_tmp, v, xp, xₚ, isFaceOrEdge)
                end
                v = next[v]
              end
            end
          end
        end
      end
    else
      #WARNING:
      # continue
      #TODO: else -> elseif, delete if
      if (ρₑ_max > ρₜ) # The boundary (isocontour) goes through the element

        Xₑ = X[:, IEN[:, el]]

        Is = MeshGrid.calculateMiniAABB_grid(Xₑ, δ, N, AABB_min, AABB_max, nsd)

        for I ∈ Is

          ii = Int(
            I[3] * (N[1] + 1) * (N[2] + 1) + I[2] * (N[1] + 1) + I[1] + 1,
          )

          v = head[ii]
          while v != -1
            x = points[:, v]

            Ξ = zeros(Float64, 3)   # local coordinates

            Ξ = compute_coords(x, ρₜ, Xₑ, ρₑ)

            H, d¹N_dξ¹, d²N_dξ², d³N_dξ³ = sfce(Ξ)
            xₚ = Xₑ * H
            n = RhoNorm(ρₑ, Ξ)

            # Sign can be zero of norm is zero:
            if sign(dot(x - xₚ, n)) == 0
              println("Normal zero sign!")
              IDmin = argmin(mapslices(norm, (Xₑ .- x), dims=1))[2] # trying to find closest node for point x and looking for its density
              Sign = ρₑ[IDmin] < ρₜ ? -1 : 1
              dist_tmp = Sign * norm(x - xₚ)
            else
              dist_tmp = sign(dot(x - xₚ, n)) * norm(x - xₚ)
            end

            (dist, xp) = WriteValue(dist_tmp, dist, xp, xₚ, v)

            v = next[v]

          end
        end
      end
    end
  end
  println("typeof xp: ", typeof(xp))

  Xg, Xp, mean_PD, max_PD = SelectProjectedNodes(mesh, grid, xp, points)
  println("mean of projected distance: ", mean_PD)
  println("maximum projected distance: ", max_PD)

  nnp = size(Xg, 1)

  IEN = [[i; i + nnp] for i = 1:nnp]
  X = vec([Xg Xp])

  Rho2sdf.exportToVTU("lines.vtu", X, IEN, 3)

  IEN = [[i] for i = 1:nnp]
  Rho2sdf.exportToVTU("Xg.vtu", Xg, IEN, 1)
  Rho2sdf.exportToVTU("Xp.vtu", Xp, IEN, 1)

  # dist = SignCorrection4SDF(dist, grid, big)

  return dist, xp

end
