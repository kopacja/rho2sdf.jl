
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
    # println("WARNING: no projected points!")
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


# Update distance filed if the distance is smaller then the previous one
function WriteValue(
  dist_tmp::Float64,
  dist_local::Vector{Vector{Float64}},
  xp_local::Vector{Matrix{Float64}},
  xₚ::Vector{Float64},
  v::Int,
  tid::Int)

  if abs(dist_tmp) < abs(dist_local[tid][v])
    dist_local[tid][v] = dist_tmp
    isFaceOrEdge = true
    xp_local[tid][:, v] = xₚ  # Update the matrix column for vertex v, ???
  end

  if (abs(dist_tmp) < abs(dist_local[tid][v]))
    dist_local[tid][v] = dist_tmp
    xp_local[tid][:, v] = xₚ
  end
  return dist_local, xp_local
end


function compute_coords(
  x::Vector,
  ρₜ::Float64,
  Xₑ::Matrix,
  ρₑ::Vector,
  starting_points::Vector{Tuple{Float64,Float64,Float64}}=[
    (0.0, 0.0, 0.0),
    (-0.5, -0.5, -0.5),
    (0.5, -0.5, -0.5),
    (0.5, 0.5, -0.5),
    (-0.5, 0.5, -0.5),
    (-0.5, -0.5, 0.5),
    (0.5, -0.5, 0.5),
    (0.5, 0.5, 0.5),
    (-0.5, 0.5, 0.5),
  ]
)
  function objective(ξ::Vector, grad::Vector)
    ξ₁, ξ₂, ξ₃ = ξ
    N = [
      -1 / 8 * (ξ₁ - 1) * (ξ₂ - 1) * (ξ₃ - 1),
      1 / 8 * (ξ₁ + 1) * (ξ₂ - 1) * (ξ₃ - 1),
      -1 / 8 * (ξ₁ + 1) * (ξ₂ + 1) * (ξ₃ - 1),
      1 / 8 * (ξ₁ - 1) * (ξ₂ + 1) * (ξ₃ - 1),
      1 / 8 * (ξ₁ - 1) * (ξ₂ - 1) * (ξ₃ + 1),
      -1 / 8 * (ξ₁ + 1) * (ξ₂ - 1) * (ξ₃ + 1),
      1 / 8 * (ξ₁ + 1) * (ξ₂ + 1) * (ξ₃ + 1),
      -1 / 8 * (ξ₁ - 1) * (ξ₂ + 1) * (ξ₃ + 1)
    ]
    x_pred = Xₑ * N
    # return sum((x - x_pred) .^ 2)
    return sum(abs.(x - x_pred))
  end

  function constraint(ξ::Vector, grad::Vector)
    ξ₁, ξ₂, ξ₃ = ξ
    N = [
      -1 / 8 * (ξ₁ - 1) * (ξ₂ - 1) * (ξ₃ - 1),
      1 / 8 * (ξ₁ + 1) * (ξ₂ - 1) * (ξ₃ - 1),
      -1 / 8 * (ξ₁ + 1) * (ξ₂ + 1) * (ξ₃ - 1),
      1 / 8 * (ξ₁ - 1) * (ξ₂ + 1) * (ξ₃ - 1),
      1 / 8 * (ξ₁ - 1) * (ξ₂ - 1) * (ξ₃ + 1),
      -1 / 8 * (ξ₁ + 1) * (ξ₂ - 1) * (ξ₃ + 1),
      1 / 8 * (ξ₁ + 1) * (ξ₂ + 1) * (ξ₃ + 1),
      -1 / 8 * (ξ₁ - 1) * (ξ₂ + 1) * (ξ₃ + 1)
    ]
    return sum(ρₑ .* N) - ρₜ
  end

  best_solution = nothing
  best_objective = Inf

  for start_point in starting_points
    opt = Opt(:LN_COBYLA, 3)
    # opt = Opt(:LN_AUGLAG_EQ, 3)
    # local_opt = Opt(:LN_COBYLA, 3)
    opt.lower_bounds = [-1.0, -1.0, -1.0]
    opt.upper_bounds = [1.0, 1.0, 1.0]
    opt.xtol_rel = 1e-6
    opt.maxeval = 500
    opt.min_objective = objective
    opt.maxtime = 1.0
    # opt.local_optimizer = local_opt
    equality_constraint!(opt, constraint, 1e-8)

    (minf, minx, ret) = optimize(opt, collect(start_point))

    if minf < best_objective
      best_objective = minf
      best_solution = minx
    end
  end

  return best_solution
end


function find_local_coordinates(
  sfce::Function,
  Xₑ::Matrix,
  xₙ::Vector
)
  starting_points::Vector{Tuple{Float64,Float64,Float64}} = [
    (0.0, 0.0, 0.0),
    (-0.5, -0.5, -0.5),
    (0.5, -0.5, -0.5),
    (0.5, 0.5, -0.5),
    (-0.5, 0.5, -0.5),
    (-0.5, -0.5, 0.5),
    (0.5, -0.5, 0.5),
    (0.5, 0.5, 0.5),
    (-0.5, 0.5, 0.5)
  ]

  function objective(ξ::Vector, grad::Vector)
    ξ₁, ξ₂, ξ₃ = ξ
    N = [
      -1 / 8 * (ξ₁ - 1) * (ξ₂ - 1) * (ξ₃ - 1),
      1 / 8 * (ξ₁ + 1) * (ξ₂ - 1) * (ξ₃ - 1),
      -1 / 8 * (ξ₁ + 1) * (ξ₂ + 1) * (ξ₃ - 1),
      1 / 8 * (ξ₁ - 1) * (ξ₂ + 1) * (ξ₃ - 1),
      1 / 8 * (ξ₁ - 1) * (ξ₂ - 1) * (ξ₃ + 1),
      -1 / 8 * (ξ₁ + 1) * (ξ₂ - 1) * (ξ₃ + 1),
      1 / 8 * (ξ₁ + 1) * (ξ₂ + 1) * (ξ₃ + 1),
      -1 / 8 * (ξ₁ - 1) * (ξ₂ + 1) * (ξ₃ + 1)
    ]
    x = Xₑ * N
    R = x - xₙ
    return sum(abs.(R))
    # return sum(R[i]^2 for i in 1:3)
  end

  best_solution = nothing
  best_objective = Inf

  for start_point in starting_points
    opt = Opt(:LN_COBYLA, 3)
    # opt = Opt(:LD_MMA, 3) # -> dává blbé výsledky
    opt.lower_bounds = [-5.0, -5.0, -5.0]
    opt.upper_bounds = [5.0, 5.0, 5.0]
    opt.xtol_rel = 1e-6
    opt.maxeval = 500
    opt.min_objective = objective
    opt.maxtime = 1.0

    (minf, minx, ret) = optimize(opt, collect(start_point))

    if minf < best_objective
      best_objective = minf
      best_solution = minx
    end
  end

  if best_solution === nothing
    return (false, [10.0, 10.0, 10.0])
  else
    return (true, best_solution)
  end
end

###_________

# Function to create AABB from a set of points
function compute_aabb(points::Matrix{Float64})
  min_bounds = minimum(points, dims=2)
  max_bounds = maximum(points, dims=2)
  return min_bounds, max_bounds
end

# Function to check if a point is inside the AABB
function is_point_inside_aabb(x::Vector{Float64}, min_bounds, max_bounds)
  return all(min_bounds .<= x) && all(x .<= max_bounds)
end

function InOut(Xₑ, xₙ)
  V = collect(eachcol(Xₑ))
  # V = [Xₑ[i] for i in 1:length(Xₑ)]
  polytope = VPolytope(V)
  return inside = xₙ ∈ polytope
end

function IsProjectedOnFullSegment(
  sfce::Function,
  Xₑ::Matrix,
  xₚ::Vector,
  el::Int,
  IEN::Matrix,
  ρₙ::Vector,
  ρₜ::Float64,
  dist_local::Vector{Vector{Float64}},
  xp_local::Vector{Matrix{Float64}},
  v::Int,
  tid::Int,
  x::Vector,
)
  (_, local_coords) = find_local_coordinates(sfce, Xₑ, xₚ)
  max_local = maximum(local_coords)
  min_local = minimum(local_coords)

  if max_local < 1.001 && min_local > -1.001
    H, _, _, _ = sfce(local_coords) # shape functions and their derivatives
    ρₑ = view(ρₙ, IEN[:, el]) # view for better performance
    ρ = H ⋅ ρₑ

    if ρ >= ρₜ
      dist_tmp = norm(x - xₚ)
      (dist_local, xp_local) = WriteValue(dist_tmp, dist_local, xp_local, xₚ, v, tid)
      return true
    end
  else
    # println("Error! Projection outside the element face")
  end

  return false
end


function update_distance_parallel!(dist_local::Vector{Vector{Float64}},
  dist_tmp::Float64,
  v::Int,
  tid::Int,
  xp_local::Vector{Matrix{Float64}},
  xₚ::Vector{Float64},
  isFaceOrEdge::Bool)

  if abs(dist_tmp) < abs(dist_local[tid][v])
    dist_local[tid][v] = dist_tmp
    isFaceOrEdge = true
    xp_local[tid][:, v] = xₚ  # Update the matrix column for vertex v, ???
  end
  return isFaceOrEdge  # Optionally return whether the vertex was updated
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
  #TODO: δ nemůžu násobit grid cell size skalárem -> ale nejdelší hrana nepravidelného elementu * skalár!
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
  big = -1.0e10
  dist = big * ones(ngp) # distance field initialization
  xp = zeros(nsd, ngp) # souřadnice bodů vrcholů (3xngp)

  # Tri mesh for pseudonormals:
  tri_mesh = Rho2sdf.extractSurfaceTriangularMesh(mesh)
  EN = NodePosition3D(tri_mesh)
  VPN, EPN = computePseudoNormals(tri_mesh)

  nthreads = Threads.nthreads()

  dist_local = [fill(big, ngp) for _ in 1:nthreads]
  xp_local = [zeros(nsd, ngp) for _ in 1:nthreads]

  p_elements = Progress(nel, 1, "Processing elements: ", 30)
  p_nodes = Progress(ngp, 1, "Processing grid nodes: ", 30)

  # Atomic countery pro oba cykly
  counter_elements = Atomic{Int}(0)
  counter_nodes = Atomic{Int}(0)

  update_interval_elements = max(1, div(nel, 100))
  update_interval_nodes = max(1, div(ngp, 100))

  Threads.@threads for el in 1:nel
    println("Element id: ", el)
    tid = Threads.threadid()

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
                  dist_tmp = norm(x - xₚ)

                  isFaceOrEdge = update_distance_parallel!(dist_local, dist_tmp, v, tid, xp_local, xₚ, isFaceOrEdge)
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
                      dist_tmp = norm(x - xₚ)

                      isFaceOrEdge = update_distance_parallel!(dist_local, dist_tmp, v, tid, xp_local, xₚ, isFaceOrEdge)
                    end
                  end
                end
                # Remaining cases:
                if (isFaceOrEdge == false)
                  dist_tmp, idx =
                    findmin([norm(x - x₁), norm(x - x₂), norm(x - x₃)]) # which node of the triangle is closer?
                  xₚ = Xt[:, idx] # the node of triangle
                  dist_tmp = norm(x - xₚ)

                  isFaceOrEdge = update_distance_parallel!(dist_local, dist_tmp, v, tid, xp_local, xₚ, isFaceOrEdge)
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

        # NOTE: two isocountour in one element:

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

                  # NOTE: Projection is inside triangle:
                  if (minimum(λ) >= 0.0) # xₚ is in the triangle, projection node x inside triangle 
                    xₚ = λ[1] * x₁ + λ[2] * x₂ + λ[3] * x₃

                    isFaceOrEdge = IsProjectedOnFullSegment(sfce, Xₑ, xₚ, el, IEN, ρₙ, ρₜ, dist_local, xp_local, v, tid, x)
                  else

                    # NOTE: Projection is on the triangle edges:
                    for j = 1:3
                      L = norm(Et[j]) # length of j triangle edge
                      xᵥ = Xt[:, j]
                      P = dot(x - xᵥ, Et[j] / L) # skalar product of vector (vertex&node) and norm edge
                      if (P >= 0 && P <= L) # is the perpendicular projection of a node onto an edge in the edge interval?
                        xₚ = xᵥ + (Et[j] / L) * P

                        isFaceOrEdge = IsProjectedOnFullSegment(sfce, Xₑ, xₚ, el, IEN, ρₙ, ρₜ, dist_local, xp_local, v, tid, x)
                      end
                    end
                  end
                  # Remaining cases:
                  # NOTE: Projection is on the triangle vertices:
                  if (isFaceOrEdge == false)
                    idx = argmin([norm(x - x₁), norm(x - x₂), norm(x - x₃)]) # find the index of the closest node

                    xₚ = Xt[:, idx] # the node of the triangle

                    isFaceOrEdge = IsProjectedOnFullSegment(sfce, Xₑ, xₚ, el, IEN, ρₙ, ρₜ, dist_local, xp_local, v, tid, x)
                  end
                  v = next[v]
                end
              end
            end
          end
        end

        # NOTE: isocountour going trought element:

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
            H, _, _, _ = sfce(Ξ)
            xₚ = Xₑ * H
            dist_tmp = norm(x - xₚ)

            (dist_local, xp_local) = WriteValue(dist_tmp, dist_local, xp_local, xₚ, v, tid)

            v = next[v]

          end
        end
      end
    end
    count = atomic_add!(counter_elements, 1)
    if count % update_interval_elements == 0
      if Threads.threadid() == 1
        update!(p_elements, count)
      end
    end
  end
  finish!(p_elements)

  # Merging the results after parallel calculation:
  for i in 1:ngp
    min_dist, min_idx = findmin(abs.(getindex.(dist_local, i)))
    dist[i] = min_dist #* sign(dist_local[min_idx][i])
    xp[:, i] = xp_local[min_idx][:, i]
  end

  signs = -1 * ones(ngp)

  @threads for i in 1:ngp
    found = false  # Flag to indicate if we need to skip to the next i
    x = points[:, i]
    # NotConv = false
    max_local = 10.0

    IDmin = argmin(mapslices(norm, (X .- x), dims=1))[2] # find ID of node (of mesh) that is closest to grid point x
    none = length(INE[IDmin]) # number of neighbour elements (that are connected to the node)
    # ρₙₑ = ρₙ[IEN[:, none]] # # nodal densities of none elements
    ρₙₑ = ρₙ[IEN[:, INE[IDmin]]] # # nodal densities of none elements
    if maximum(ρₙₑ) < ρₜ # empty element -> jump to another i
      continue
    end

    for j in 1:none
      el = INE[IDmin][j]

      Xₑ = X[:, IEN[:, el]]
      # Compute the AABB
      min_bounds, max_bounds = compute_aabb(Xₑ)

      # Check if the point x is inside the AABB
      inside = is_point_inside_aabb(x, min_bounds, max_bounds)

      if inside
        # if InOut(Xₑ, x)

        (NotConv, local_coords) = find_local_coordinates(sfce, Xₑ, x)
        max_local_new = maximum(abs.(local_coords))
        # max_local = maximum(local_coords)
        # min_local = minimum(local_coords)
        #
        # if max_local < 1.0001 && min_local > -1.0001
        if max_local_new < 1.2 && max_local > max_local_new

          H, _, _, _ = sfce(local_coords) # tvarové funkce a jejich derivace
          ρₑ = ρₙ[IEN[:, el]] # nodal densities for one element
          ρ = H ⋅ ρₑ

          if ρ >= ρₜ
            # println("inside")
            signs[i] = 1.0
          end
          #   found = true
          #   break
          # else
          # if !NotConv
          # println("point: ", i)
          # end
          max_local = max_local_new
        end
      end
      # if !found && NotConv && mean(ρₙₑ) > ρₜ && j == none
      #   println("point: ", i)
      # end
      # if found
      #   break
      # end
    end

    # if !found
    #   println("node position not detected: ")
    #   println("id of grid point: ", i)
    # end
    count = atomic_add!(counter_nodes, 1)
    if count % update_interval_nodes == 0
      if Threads.threadid() == 1
        update!(p_nodes, count)
      end
    end
  end
  finish!(p_nodes)

  @save "Z_Chapadlo_xp.jld2" xp
  @save "Z_Chapadlo_dist.jld2" dist
  @save "Z_Chapadlo_mesh.jld2" mesh
  @save "Z_Chapadlo_grid.jld2" grid
  @save "Z_Chapadlo_points.jld2" points

  @save "Z_Chapadlo_signs.jld2" signs

  # @load "Z_Chapadlo_xp.jld2" xp
  # @load "Z_Chapadlo_dist.jld2" dist
  # @load "Z_Chapadlo_mesh.jld2" mesh
  # @load "Z_Chapadlo_grid.jld2" grid
  # @load "Z_Chapadlo_points.jld2" points


  # println("typeof xp: ", typeof(xp))

  Xg, Xp, mean_PD, max_PD = SelectProjectedNodes(mesh, grid, xp, points)
  # println("mean of projected distance: ", mean_PD)
  # println("maximum projected distance: ", max_PD)

  nnp = size(Xg, 1)

  IEN = [[i; i + nnp] for i = 1:nnp]
  X = vec([Xg Xp])

  Rho2sdf.exportToVTU("lines.vtu", X, IEN, 3)

  IEN = [[i] for i = 1:nnp]
  Rho2sdf.exportToVTU("Xg.vtu", Xg, IEN, 1)
  Rho2sdf.exportToVTU("Xp.vtu", Xp, IEN, 1)


  dist = abs.(dist) .* signs

  # @save "Z_Leg_dist.jld2" dist
  @save "Z_Chapadlo_cele_dist.jld2" dist
  # @save "Z_Leg_grid.jld2" grid
  @save "Z_Chapadlo_cele_grid.jld2" grid

  return dist, xp

end
