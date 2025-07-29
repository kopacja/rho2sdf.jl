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

  if (abs(dist_tmp) < abs(dist_local[tid][v]))
    dist_local[tid][v] = dist_tmp
    xp_local[tid][:, v] = xₚ
  end
  return dist_local, xp_local
end



# Function to create AABB from a set of points
function compute_aabb(points::SubArray)
  min_bounds = minimum(points, dims=2)
  max_bounds = maximum(points, dims=2)
  return min_bounds, max_bounds
end

# Function to check if a point is inside the AABB
function is_point_inside_aabb(x::SubArray, min_bounds, max_bounds)
  return all(min_bounds .<= x) && all(x .<= max_bounds)
end

function InOut(Xₑ, xₙ)
  V = collect(eachcol(Xₑ))
  # V = [Xₑ[i] for i in 1:length(Xₑ)]
  polytope = VPolytope(V)
  return inside = xₙ ∈ polytope
end

function IsProjectedOnFullSegment(
  mesh::Mesh{T},
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
) where {T<:AbstractElement}
  
  (_, local_coords) = find_local_coordinates(Xₑ, xₚ, mesh.element_type)
  
  # Validate coordinates based on element type
  coords_valid = if T == HEX8
    maximum(abs.(local_coords)) < 1.001
  elseif T == TET4
    validate_local_coords(TET4, local_coords) && sum(local_coords) ≤ 1.001
  else
    false
  end

  if coords_valid
    # Use shape functions appropriate for element type
    H = shape_functions(mesh.element_type, local_coords)
    ρₑ = view(ρₙ, IEN[:, el])
    ρ = H ⋅ ρₑ

    if ρ >= ρₜ
      dist_tmp = norm(x - xₚ)
      (dist_local, xp_local) = WriteValue(dist_tmp, dist_local, xp_local, xₚ, v, tid)
      return true
    end
  else
    # println("Error! Projection outside the element")
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
function evalDistances(
  mesh::Mesh,
  grid::Grid,
  points::Matrix,
  ρₙ::Vector{Float64},
  ρₜ::Float64;
  taskName::String="block",
  plot_projection_points_and_lines::Bool=false
)# where {T<:AbstractElement}

  print_info("\nComputing distace function...")

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
  # println("number of all elements: ", nel)

  ngp = grid.ngp # number of nodes in grid
  big = -1.0e10
  dist = big * ones(ngp) # distance field initialization
  xp = zeros(nsd, ngp) # souřadnice bodů vrcholů (3xngp)

  # Tri mesh for pseudonormals:
  tri_mesh = Rho2sdf.extractSurfaceTriangularMesh(mesh)
  EN = NodePosition3D(tri_mesh)
  _, EPN = computePseudoNormals(tri_mesh)

  nthreads = Threads.nthreads()

  dist_local = [fill(big, ngp) for _ in 1:nthreads]
  xp_local = [zeros(nsd, ngp) for _ in 1:nthreads]

  p_elements = Progress(nel, 1, "Processing elements: ", 30)

  # Atomic countery pro oba cykly
  counter_elements = Atomic{Int}(0)

  update_interval_elements = max(1, div(nel, 100))

  Threads.@threads for el in 1:nel
    # println("Element id: ", el)
    tid = Threads.threadid()

    ρₑ = ρₙ[IEN[:, el]] # nodal densities for one element

    ρₑ_min = minimum(ρₑ)
    ρₑ_max = maximum(ρₑ)
    if (ρₑ_min >= ρₜ) # the boundary does not cross through the element
      # commonEls = []
      #
      # # cycle through element faces (6)
      # for sg = 1:nes
      #   commonEls = INE[IEN[mesh.ISN[sg][1], el]]
      #   for a = 2:nsn
      #     idx = findall(in(INE[IEN[ISN[sg][a], el]]), commonEls) # for how many elements does this face belong ?
      #     commonEls = commonEls[idx]
      #   end
      #
      #   if (length(commonEls) == 1) # = is a part of the outer boundary of the body
      #     Xs = X[:, IEN[ISN[sg], el]]
      #     Xc = vec(mean(Xs, dims=2))
      #
      #     for a = 1:nsn # cycle through number of all nodals belong to face
      #
      #       # coordinates of nodes of the triangle
      #       x₁ = Xs[:, a]
      #       x₂ = Xs[:, (a%nsn)+1]
      #       x₃ = Xc
      #
      #       # coordinates of the vertices of the triangle
      #       Xt = [x₁ x₂ x₃]
      #
      #       # finding coresponding triangle:
      #       ID_tri = find_triangle_position(EN, [x₁ x₂ x₃])
      #
      #       #NOTE: From this part it is same as in sdfOnTriangularMesh ->
      #
      #       # Triangle edges
      #       Et = calculate_triangle_edges(Xt)
      #
      #       n = cross(Et[1], Et[2]) # norm of triangle
      #       n = n / norm(n) # unit norm
      #
      #       # Nodes of mini AABB grid:
      #       Is = MeshGrid.calculateMiniAABB_grid(Xt, δ, N, AABB_min, AABB_max, nsd)
      #
      #       for I ∈ Is # cycle through the nodes of the mini AABB grid
      #         ii = Int( # node ID
      #           I[3] * (N[1] + 1) * (N[2] + 1) + I[2] * (N[1] + 1) + I[1] + 1,
      #         )
      #         v = head[ii]
      #         while v != -1
      #           x = points[:, v]
      #           λ = barycentricCoordinates(x₁, x₂, x₃, n, x)
      #
      #           xₚ = zeros(nsd) # projection
      #
      #           isFaceOrEdge = false # projection check
      #
      #           if (minimum(λ) >= 0.0) # xₚ is in the triangle, projection node x inside triangle 
      #             xₚ = λ[1] * x₁ + λ[2] * x₂ + λ[3] * x₃
      #             dist_tmp = norm(x - xₚ)
      #
      #             isFaceOrEdge = update_distance_parallel!(dist_local, dist_tmp, v, tid, xp_local, xₚ, isFaceOrEdge)
      #           else
      #
      #             # Edges of the triangle:
      #             for j = 1:3
      #               L = norm(Et[j]) # length of j triangle edge
      #               xᵥ = Xt[:, j]
      #               P = dot(x - xᵥ, Et[j] / L) # skalar product of vector (vertex&node) and norm edge
      #               if (P >= 0 && P <= L) # is the perpendicular projection of a node onto an edge in the edge interval?
      #                 xₚ = xᵥ + (Et[j] / L) * P
      #                 n_edge = n
      #                 n_edge = EPN[ID_tri][j]
      #                 dist_tmp = norm(x - xₚ)
      #
      #                 isFaceOrEdge = update_distance_parallel!(dist_local, dist_tmp, v, tid, xp_local, xₚ, isFaceOrEdge)
      #               end
      #             end
      #           end
      #           # Remaining cases:
      #           if (isFaceOrEdge == false)
      #             dist_tmp, idx =
      #               findmin([norm(x - x₁), norm(x - x₂), norm(x - x₃)]) # which node of the triangle is closer?
      #             xₚ = Xt[:, idx] # the node of triangle
      #             dist_tmp = norm(x - xₚ)
      #
      #             isFaceOrEdge = update_distance_parallel!(dist_local, dist_tmp, v, tid, xp_local, xₚ, isFaceOrEdge)
      #           end
      #           v = next[v]
      #         end
      #       end
      #     end
      #   end
      # end
      process_boundary_faces!(mesh, el, grid, points, ρₙ, ρₜ, dist_local, xp_local, 
                          head, next, N, AABB_min, AABB_max, δ, EN, EPN, tid,
                          true)  # PŘIDÁN PARAMETR - true pro solidní elementy
    else
      #TODO: else -> elseif, delete if
      if (ρₑ_max > ρₜ) # The boundary (isocontour) goes through the element

        # Xₑ = X[:, IEN[:, el]]
        #
        # # NOTE: two isocountour in one element:
        #
        # # cycle through element faces (6)
        # for sg = 1:nes
        #   commonEls = INE[IEN[mesh.ISN[sg][1], el]]
        #   for a = 2:nsn
        #     idx = findall(in(INE[IEN[ISN[sg][a], el]]), commonEls) # for how many elements does this face belong ?
        #     commonEls = commonEls[idx]
        #   end
        #
        #   if (length(commonEls) == 1) # = is a part of the outer boundary of the body
        #     Xs = X[:, IEN[ISN[sg], el]]
        #     Xc = vec(mean(Xs, dims=2))
        #
        #     for a = 1:nsn # cycle through number of all nodals belong to face
        #
        #       # coordinates of nodes of the triangle
        #       x₁ = Xs[:, a]
        #       x₂ = Xs[:, (a%nsn)+1]
        #       x₃ = Xc
        #
        #       # coordinates of the vertices of the triangle
        #       Xt = [x₁ x₂ x₃]
        #
        #       #NOTE: From this part it is same as in sdfOnTriangularMesh ->
        #
        #       # Triangle edges
        #       Et = calculate_triangle_edges(Xt)
        #
        #       n = cross(Et[1], Et[2]) # norm of triangle
        #       n = n / norm(n) # unit norm
        #
        #       # Nodes of mini AABB grid:
        #       Is = MeshGrid.calculateMiniAABB_grid(Xt, δ, N, AABB_min, AABB_max, nsd)
        #
        #       for I ∈ Is # cycle through the nodes of the mini AABB grid
        #         ii = Int( # node ID
        #           I[3] * (N[1] + 1) * (N[2] + 1) + I[2] * (N[1] + 1) + I[1] + 1,
        #         )
        #         v = head[ii]
        #         while v != -1
        #           x = points[:, v]
        #           λ = barycentricCoordinates(x₁, x₂, x₃, n, x)
        #
        #           xₚ = zeros(nsd) # projection
        #
        #           isFaceOrEdge = false # projection check
        #
        #           # NOTE: Projection is inside triangle:
        #           if (minimum(λ) >= 0.0) # xₚ is in the triangle, projection node x inside triangle 
        #             xₚ = λ[1] * x₁ + λ[2] * x₂ + λ[3] * x₃
        #
        #             isFaceOrEdge = IsProjectedOnFullSegment(sfce, Xₑ, xₚ, el, IEN, ρₙ, ρₜ, dist_local, xp_local, v, tid, x)
        #           else
        #
        #             # NOTE: Projection is on the triangle edges:
        #             for j = 1:3
        #               L = norm(Et[j]) # length of j triangle edge
        #               xᵥ = Xt[:, j]
        #               P = dot(x - xᵥ, Et[j] / L) # skalar product of vector (vertex&node) and norm edge
        #               if (P >= 0 && P <= L) # is the perpendicular projection of a node onto an edge in the edge interval?
        #                 xₚ = xᵥ + (Et[j] / L) * P
        #
        #                 isFaceOrEdge = IsProjectedOnFullSegment(sfce, Xₑ, xₚ, el, IEN, ρₙ, ρₜ, dist_local, xp_local, v, tid, x)
        #               end
        #             end
        #           end
        #           # Remaining cases:
        #           # NOTE: Projection is on the triangle vertices:
        #           if (isFaceOrEdge == false)
        #             idx = argmin([norm(x - x₁), norm(x - x₂), norm(x - x₃)]) # find the index of the closest node
        #
        #             xₚ = Xt[:, idx] # the node of the triangle
        #
        #             isFaceOrEdge = IsProjectedOnFullSegment(sfce, Xₑ, xₚ, el, IEN, ρₙ, ρₜ, dist_local, xp_local, v, tid, x)
        #           end
        #           v = next[v]
        #         end
        #       end
        #     end
        #   end
        # end
        #
        # # NOTE: isocountour going trought element:
        #
        # Is = MeshGrid.calculateMiniAABB_grid(Xₑ, δ, N, AABB_min, AABB_max, nsd)
        #
        # for I ∈ Is
        #
        #   ii = Int(
        #     I[3] * (N[1] + 1) * (N[2] + 1) + I[2] * (N[1] + 1) + I[1] + 1,
        #   )
        #
        #   v = head[ii]
        #   while v != -1
        #     x = points[:, v]
        #
        #     Ξ = zeros(Float64, 3)   # local coordinates
        #
        #     Ξ = compute_coords_on_iso(x, ρₜ, Xₑ, ρₑ)
        #     H = sfce(Ξ)
        #     xₚ = Xₑ * H
        #     dist_tmp = norm(x - xₚ)
        #
        #     (dist_local, xp_local) = WriteValue(dist_tmp, dist_local, xp_local, xₚ, v, tid)
        #
        #     v = next[v]
        #
        #   end
        # end
        process_isocontour_element!(mesh, el, grid, points, ρₙ, ρₜ, dist_local, xp_local,
                                  head, next, N, AABB_min, AABB_max, δ, EN, EPN, tid)
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
    dist[i] = min_dist
    xp[:, i] = xp_local[min_idx][:, i]
  end

  # This section is used to check the DF results
  if plot_projection_points_and_lines
    Xg, Xp, mean_PD, max_PD = SelectProjectedNodes(mesh, grid, xp, points)

    #NOTE: This part is used for debugging:
    # println("mean of projected distance: ", mean_PD)
    # println("maximum projected distance: ", max_PD)

    nnp = size(Xg, 1)

    IEN = [[i; i + nnp] for i = 1:nnp]
    X = vec([Xg Xp])

    Rho2sdf.exportToVTU("lines_$(taskName).vtu", X, IEN, 3)

    IEN = [[i] for i = 1:nnp]
    Rho2sdf.exportToVTU("Xg_$(taskName).vtu", Xg, IEN, 1)
    Rho2sdf.exportToVTU("Xp_$(taskName).vtu", Xp, IEN, 1)
  end

  dist = abs.(dist)

  return dist, xp

end

# Helper function for boundary face processing
function process_boundary_faces!(mesh::Mesh{T}, el::Int, grid::Grid, points::Matrix,
                                ρₙ::Vector{Float64}, ρₜ::Float64, 
                                dist_local::Vector{Vector{Float64}}, xp_local::Vector{Matrix{Float64}},
                                head, next, N, AABB_min, AABB_max, δ, EN, EPN, tid::Int,
                                is_solid_element::Bool) where {T<:AbstractElement}  # PŘIDÁN PARAMETR
  
  # Process faces based on element type
  for sg in 1:mesh.nes
    # Check if face is on boundary
    commonEls = mesh.INE[mesh.IEN[mesh.ISN[sg][1], el]]
    for a in 2:mesh.nsn
      idx = findall(in(mesh.INE[mesh.IEN[mesh.ISN[sg][a], el]]), commonEls)
      commonEls = commonEls[idx]
    end

    if length(commonEls) == 1  # Boundary face
      Xs = mesh.X[:, mesh.IEN[mesh.ISN[sg], el]]
      Xc = vec(mean(Xs, dims=2))

      # Process triangular subdivisions of face
      for a in 1:mesh.nsn
        x₁ = Xs[:, a]
        x₂ = Xs[:, (a % mesh.nsn) + 1]
        x₃ = Xc

        Xt = [x₁ x₂ x₃]
        ID_tri = find_triangle_position(EN, [x₁ x₂ x₃])

        # Process triangle projection with solid element info
        process_triangle_projection!(Xt, grid, points, mesh, el, ρₙ, ρₜ, 
                                   dist_local, xp_local, head, next, N, 
                                   AABB_min, AABB_max, δ, EN, EPN, ID_tri, tid,
                                   is_solid_element)  # PŘEDÁN PARAMETR
      end
    end
  end
end

# Helper function for isocontour element processing  
function process_isocontour_element!(mesh::Mesh{T}, el::Int, grid::Grid, points::Matrix,
                                   ρₙ::Vector{Float64}, ρₜ::Float64,
                                   dist_local::Vector{Vector{Float64}}, xp_local::Vector{Matrix{Float64}},
                                   head, next, N, AABB_min, AABB_max, δ, EN, EPN, tid::Int) where {T<:AbstractElement}
  
  Xₑ = mesh.X[:, mesh.IEN[:, el]]
  ρₑ = ρₙ[mesh.IEN[:, el]]
  
  # First process boundary faces like solid elements
  process_boundary_faces!(mesh, el, grid, points, ρₙ, ρₜ, dist_local, xp_local,
                        head, next, N, AABB_min, AABB_max, δ, EN, EPN, tid,
                        false)  # PŘIDÁN PARAMETR - false pro elementy s isokonturou
  
  # Then process interior isocontour
  Is = MeshGrid.calculateMiniAABB_grid(Xₑ, δ, N, AABB_min, AABB_max, mesh.nsd)

  for I ∈ Is
    ii = Int(I[3] * (N[1] + 1) * (N[2] + 1) + I[2] * (N[1] + 1) + I[1] + 1)
    v = head[ii]
    
    while v != -1
      x = points[:, v]
      
      # Compute coordinates on isocontour using element-specific method
      Ξ = compute_coords_on_iso(x, ρₜ, Xₑ, ρₑ, mesh.element_type)
      H = shape_functions(mesh.element_type, Ξ)
      xₚ = Xₑ * H
      dist_tmp = norm(x - xₚ)

      (dist_local, xp_local) = WriteValue(dist_tmp, dist_local, xp_local, xₚ, v, tid)
      v = next[v]
    end
  end
end

# Complete triangle projection processing algorithm
function process_triangle_projection!(Xt, grid, points, mesh, el, ρₙ, ρₜ,
                                     dist_local, xp_local, head, next, N,
                                     AABB_min, AABB_max, δ, EN, EPN, ID_tri, tid,
                                     is_solid_element::Bool)  # PŘIDÁN PARAMETR
    
    nsd = mesh.nsd
    
    # Calculate triangle edges
    Et = calculate_triangle_edges(Xt)
    
    # Compute triangle normal vector
    n = cross(Et[1], Et[2])
    n = n / norm(n)  # unit normal
    
    # Extract triangle vertices
    x₁ = Xt[:, 1]
    x₂ = Xt[:, 2] 
    x₃ = Xt[:, 3]
    
    # Calculate mini AABB grid nodes for this triangle
    Is = MeshGrid.calculateMiniAABB_grid(Xt, δ, N, AABB_min, AABB_max, nsd)
    
    # Process all grid points in the mini AABB
    for I ∈ Is
        # Calculate grid node ID
        ii = Int(I[3] * (N[1] + 1) * (N[2] + 1) + I[2] * (N[1] + 1) + I[1] + 1)
        
        # Process all points associated with this grid node
        v = head[ii]
        while v != -1
            x = points[:, v]
            
            # Calculate barycentric coordinates
            λ = barycentricCoordinates(x₁, x₂, x₃, n, x)
            
            xₚ = zeros(nsd)  # projection point
            isFaceOrEdge = false  # projection success flag
            
            # Check if projection is inside the triangle
            if minimum(λ) >= 0.0
                # Point projects inside triangle - use barycentric interpolation
                xₚ = λ[1] * x₁ + λ[2] * x₂ + λ[3] * x₃
                dist_tmp = norm(x - xₚ)
                
                # OPRAVENÁ LOGIKA
                if is_solid_element
                    # Pro solidní elementy přímo aktualizujeme vzdálenost
                    isFaceOrEdge = update_distance_parallel!(dist_local, dist_tmp, v, tid, xp_local, xₚ, isFaceOrEdge)
                else
                    # Pro elementy s isokonturou musíme validovat projekci
                    isFaceOrEdge = IsProjectedOnFullSegment(mesh, mesh.X[:, mesh.IEN[:, el]], xₚ, el, mesh.IEN, ρₙ, ρₜ, 
                                                          dist_local, xp_local, v, tid, x)
                end
            else
                # Projection is outside triangle - check edges
                for j = 1:3
                    L = norm(Et[j])  # edge length
                    xᵥ = Xt[:, j]    # edge start vertex
                    
                    # Project point onto edge line
                    P = dot(x - xᵥ, Et[j] / L)  # scalar projection
                    
                    # Check if projection falls within edge segment
                    if P >= 0 && P <= L
                        xₚ = xᵥ + (Et[j] / L) * P  # projected point on edge
                        
                        # Use pseudo-normal for edge if available
                        if ID_tri > 0 && ID_tri <= length(EPN)
                            n_edge = EPN[ID_tri][j]
                        else
                            n_edge = n  # fallback to face normal
                        end
                        
                        dist_tmp = norm(x - xₚ)
                        
                        # OPRAVENÁ LOGIKA
                        if is_solid_element
                            isFaceOrEdge = update_distance_parallel!(dist_local, dist_tmp, v, tid, xp_local, xₚ, isFaceOrEdge)
                        else
                            isFaceOrEdge = IsProjectedOnFullSegment(mesh, mesh.X[:, mesh.IEN[:, el]], xₚ, el, mesh.IEN, ρₙ, ρₜ,
                                                                  dist_local, xp_local, v, tid, x)
                        end
                        
                        # Break on first valid edge projection
                        if isFaceOrEdge
                            break
                        end
                    end
                end
            end
            
            # If no face/edge projection worked, project to nearest vertex
            if !isFaceOrEdge
                distances_to_vertices = [norm(x - x₁), norm(x - x₂), norm(x - x₃)]
                dist_tmp, idx = findmin(distances_to_vertices)
                xₚ = Xt[:, idx]  # closest vertex
                
                # OPRAVENÁ LOGIKA
                if is_solid_element
                    isFaceOrEdge = update_distance_parallel!(dist_local, dist_tmp, v, tid, xp_local, xₚ, isFaceOrEdge)
                else
                    isFaceOrEdge = IsProjectedOnFullSegment(mesh, mesh.X[:, mesh.IEN[:, el]], xₚ, el, mesh.IEN, ρₙ, ρₜ,
                                                          dist_local, xp_local, v, tid, x)
                end
            end
            
            # Move to next point in linked list
            v = next[v]
        end
    end
end
