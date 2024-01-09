function barycentricCoordinates(x₁::Vector{Float64},
  x₂::Vector{Float64},
  x₃::Vector{Float64},
  n::Vector{Float64},
  x::Vector{Float64})

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
    
    n_max, i_max = findmax(abs.(n))
    A[i_max, :] = [1.0 1.0 1.0]
    b[i_max] = 1.0

    return λ = A \ b # baricentrické souřadnice
end


function calculate_triangle_edges(Xt::Matrix{Float64})
    Et = Vector{Vector{Float64}}(undef, 3) # Preallocate with undefined values
    Et[1] = Xt[:, 2] - Xt[:, 1]
    Et[2] = Xt[:, 3] - Xt[:, 2]
    Et[3] = Xt[:, 1] - Xt[:, 3]
    return Et
end

function evalSignedDistancesOnTriangularMesh(mesh::Mesh, grid::Grid)

    points = MeshGrid.generateGridPoints(grid)
    linkedList = MeshGrid.LinkedList(grid, points)

    println("Init pseudo normals...")
    VPN, EPN = computePseudoNormals(mesh)
    println("...done.")

    head = linkedList.head
    next = linkedList.next
    N = linkedList.grid.N  # Number of divisions along each axis of the grid
    AABB_min = linkedList.grid.AABB_min
    AABB_max = linkedList.grid.AABB_max
    δ = 2.5 * grid.cell_size
    X = mesh.X
    IEN = mesh.IEN
    # INE = mesh.INE
    nsd = mesh.nsd
    nel = mesh.nel

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
        
        # Triangle inside AABB 
        Xt_min = minimum(Xt, dims = 2) .- δ # bottom left corner coord
        Xt_max = maximum(Xt, dims = 2) .+ δ # top righ corner coord
        
        # Mini AABB for triangle:
        I_min = floor.(N .* (Xt_min .- AABB_min) ./ (AABB_max .- AABB_min)) # Triangle location (index) within the grid
        I_max = floor.(N .* (Xt_max .- AABB_min) ./ (AABB_max .- AABB_min)) # Triangle location (index) within the grid

        for j = 1:nsd # am I inside AABB?
            if (I_min[j] < 0)
                I_min[j] = 0
            end
            if (I_max[j] >= N[j])
                I_max[j] = N[j]
            end
        end

        x₁, x₂, x₃ = Xt[:, 1], Xt[:, 2], Xt[:, 3] # coordinates of nodes of the triangle

        Et = calculate_triangle_edges(Xt)

        n = cross(Et[1], Et[2]) # norm of triangle
        n = n / norm(n) # unit norm

        Is = Iterators.product( # step range of mini AABB
            I_min[1]:I_max[1],
            I_min[2]:I_max[2],
            I_min[3]:I_max[3],
        )
        for I ∈ Is # cycle through the nodes of the mini AABB grid
        # I is vector - node coords
            i = Int( # node ID
                I[3] * (N[1] + 1) * (N[2] + 1) + I[2] * (N[1] + 1) + I[1] + 1,
            )
            v = head[i]
            while v != -1
                x = points[:, v]
                λ = barycentricCoordinates(x₁, x₂, x₃, n, x)
        
                xₚ = zeros(nsd)
                # kam jsem se promítl?
                isFace = false
                isEdge = false
                isVertex = false
                if (minimum(λ) >= 0.0) # xₚ is in the triangle, projection node x inside triangle 
                    xₚ = λ[1] * x₁ + λ[2] * x₂ + λ[3] * x₃
                    dist_tmp = dot(x - xₚ, n)
                    if (abs(dist_tmp) < abs(dist[v]))
                        dist[v] = dist_tmp
                        isFace = true
                        xp[:, v] = xₚ # k čemu to je?
                    end
                else

                    for j = 1:3
                        L = norm(Et[j])
                        xᵥ = Xt[:, j]
                        P = dot(x - xᵥ, Et[j] / L)
                        if (P >= 0 && P <= L) # pohybuji se na intervalu hrany
                            xₚ = xᵥ + (Et[j] / L) * P
                            n_edge = EPN[el][j]
                            dist_tmp = sign(dot(x - xₚ, n_edge)) * norm(x - xₚ)

                            if (abs(dist_tmp) < abs(dist[v]))
                                dist[v] = dist_tmp
                                # isEdge = true
                                isVertex = true
                                xp[:, v] = xₚ # k čemu to je?
                            end
                        end
                    end
                end
                if (isFace == false && isEdge == false)
                    dist_tmp, idx =
                        findmin([norm(x - x₁), norm(x - x₂), norm(x - x₃)])
                    xₚ = Xt[:, idx]
                    n_vertex = VPN[IEN[idx, el]]
                    dist_tmp = dist_tmp * sign(dot(x - xₚ, n_vertex))
                    if (abs(dist_tmp) < abs(dist[v]))
                        dist[v] = dist_tmp
                        isVertex = true
                        xp[:, v] = xₚ # k čemu to je?
                    end
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

    return dist
end


