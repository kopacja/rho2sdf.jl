function barycentricCoordinates(
  x₁::Vector{Float64}, # coordinates of vertex one of the triangle
  x₂::Vector{Float64}, # coordinates of vertex two of the triangle
  x₃::Vector{Float64}, # coordinates of vertex tree of the triangle
  n::Vector{Float64},  # unit normal to the face of the triangle
  x::Vector{Float64})  # one node of the grid

    A = [
        (x₁[2]*n[3]-x₁[3]*n[2]) (x₂[2]*n[3]-x₂[3]*n[2]) (x₃[2]*n[3]-x₃[3]*n[2])
        (x₁[3]*n[1]-x₁[1]*n[3]) (x₂[3]*n[1]-x₂[1]*n[3]) (x₃[3]*n[1]-x₃[1]*n[3])
        (x₁[1]*n[2]-x₁[2]*n[1]) (x₂[1]*n[2]-x₂[2]*n[1]) (x₃[1]*n[2]-x₃[2]*n[1])
    ]
    b = [
        x[2] * n[3] - x[3] * n[2],
        x[3] * n[1] - x[1] * n[3],
        x[1] * n[2] - x[2] * n[1],
    ]
    
    n_max, i_max = findmax(abs.(n)) ##???
    A[i_max, :] = [1.0 1.0 1.0]
    b[i_max] = 1.0

    return λ = A \ b # barycentric coordinates
end


function calculate_triangle_edges(
    Xt::Matrix{Float64}) # Coordinates of the vertices of the triangle

    Et = Vector{Vector{Float64}}(undef, 3) # Preallocate with undefined values
    Et[1] = Xt[:, 2] - Xt[:, 1]
    Et[2] = Xt[:, 3] - Xt[:, 2]
    Et[3] = Xt[:, 1] - Xt[:, 3]
    return Et
end


function update_distance!(dist::Vector{Float64},
  dist_tmp::Float64,
  v::Int,
  xp::Matrix{Float64},
  xₚ::Vector{Float64},
  isFaceOrEdge::Bool)

    if abs(dist_tmp) < abs(dist[v])
        dist[v] = dist_tmp
        isFaceOrEdge = true
        xp[:, v] = xₚ  # Update the matrix column for vertex v, ???
    end
    return isFaceOrEdge  # Optionally return whether the vertex was updated
end


  
function evalSignedDistancesOnTriangularMesh(mesh::TriangularMesh, grid::Grid)

    points = MeshGrid.generateGridPoints(grid)
    linkedList = MeshGrid.LinkedList(grid, points)

    println("Init pseudo normals...")
    VPN, EPN = computePseudoNormals(mesh)
    println("...done.")

    head = linkedList.head # ID pravidelné bunky (pozice), index bodu z points
    next = linkedList.next # vec délky points, další uzly pro danou bunku, když -1 tak už další není
    N = linkedList.grid.N  # Number of divisions along each axis of the grid
    AABB_min = linkedList.grid.AABB_min # Minimum coordinates of the Axis-Aligned Bounding Box (AABB)
    AABB_max = linkedList.grid.AABB_max # Maximum coordinates of the AABB
    δ = 2.5 * grid.cell_size # offset for mini AABB
    
    X   = mesh.X   # vector of nodes positions
    IEN = mesh.IEN # ID element -> ID nodes
    nsd = mesh.nsd # number of spacial dimensions
    nel = mesh.nel # number of all elements

    ngp = grid.ngp
    big = 1.0e10
    # dists = [big * ones(ngp) for _ in 1:nthreads()]
    # xps = [zeros(nsd, ngp) for _ in 1:nthreads()]
    dist = big * ones(ngp)
    xp = zeros(nsd, ngp)

    for el = 1:nel
    # @threads for el = 1:nel
        # Access thread-local variables
        # dist = dists[threadid()]
        # xp = xps[threadid()]

        Xt = X[:, IEN[:, el]] # coordinates of the vertices of the triangle
        
        x₁, x₂, x₃ = Xt[:, 1], Xt[:, 2], Xt[:, 3] # coordinates of nodes of the triangle

        Et = calculate_triangle_edges(Xt)

        n = cross(Et[1], Et[2]) # norm of triangle
        n = n / norm(n) # unit norm

        # Nodes of mini AABB grid:
        Is = MeshGrid.calculateMiniAABB_grid(Xt, δ, N, AABB_min, AABB_max, nsd)

        for I ∈ Is # cycle through the nodes of the mini AABB grid
        # I is vector - node coords
            i = Int( # node ID
                I[3] * (N[1] + 1) * (N[2] + 1) + I[2] * (N[1] + 1) + I[1] + 1,
            )
            v = head[i]
            while v != -1
                x = points[:, v]
                λ = barycentricCoordinates(x₁, x₂, x₃, n, x)
        
                xₚ = zeros(nsd) # projection
                
                isFaceOrEdge = false # projection check

                # Inside triangle:
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
                            n_edge = EPN[el][j]
                            dist_tmp = sign(dot(x - xₚ, n_edge)) * norm(x - xₚ) ## hustý, ale nechápu, asi ok

                            isFaceOrEdge = update_distance!(dist, dist_tmp, v, xp, xₚ, isFaceOrEdge)
                        end
                    end
                end
                # Remaining cases:
                if (isFaceOrEdge == false) 
                    dist_tmp, idx =
                        findmin([norm(x - x₁), norm(x - x₂), norm(x - x₃)]) # which node of the triangle is closer?
                    xₚ = Xt[:, idx] # the node of triangle
                    n_vertex = VPN[IEN[idx, el]]
                    dist_tmp = dist_tmp * sign(dot(x - xₚ, n_vertex))

                    isFaceOrEdge = update_distance!(dist, dist_tmp, v, xp, xₚ, isFaceOrEdge)
                end
                v = next[v]
            end
        end
    end
    # dist = marchingCubes(dist, N.+1, big) 

    # Combine results from each thread
    # final_dist = reduce(min, dists)

    for i in eachindex(dist)
        if (abs(dist[i]) > norm(grid.cell_size))
            dist[i] = sign(dist[i]) * norm(grid.cell_size)
        end
    end

    Xg, Xp, mean_PD, max_PD = SelectProjectedNodes(mesh, grid, xp, points)
    println("mean of projected distance: ", mean_PD)
    println("maximum projected distance: ", max_PD)

    nnp = size(Xg, 1)

    IEN = [[i; i + nnp] for i = 1:nnp]
    X = vec([Xg Xp])

    Rho2sdf.exportToVTU("lines_STL.vtu", X, IEN, 3)

    IEN = [[i] for i = 1:nnp]
    Rho2sdf.exportToVTU("Xg_STL.vtu", Xg, IEN, 1)
    Rho2sdf.exportToVTU("Xp_STL.vtu", Xp, IEN, 1)

    return dist
end


