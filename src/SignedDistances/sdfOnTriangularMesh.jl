







function evalSignedDistancesOnTriangularMesh(mesh::Mesh, grid::Grid)

    points = generateGridPoints(grid)
    linkedList = LinkedList(grid, points)

    println("Init pseudo normals...")
    VPN, EPN = computePseudoNormals(mesh)
    println("...done.")

    head = linkedList.head
    next = linkedList.next
    N = linkedList.grid.N
    AABB_min = linkedList.grid.AABB_min
    AABB_max = linkedList.grid.AABB_max
    δ = 2.5 * grid.cell_size
    X = mesh.X
    IEN = mesh.IEN
    nsd = mesh.nsd
    nel = mesh.nel

    ngp = grid.ngp
    big = 1.0e10
    dist = big * ones(ngp)
    xp = zeros(nsd, ngp)

    for el in 1:nel
        Xt = X[:, IEN[:, el]] # souřadnice trojúhelníku
        Xt_min = minimum(Xt, dims=2) .- δ # trojúhelník vložený do BB a přifouknu o deltu
        Xt_max = maximum(Xt, dims=2) .+ δ

        I_min = floor.(N .* (Xt_min .- AABB_min) ./ (AABB_max .- AABB_min)) # index oddílu
        I_max = floor.(N .* (Xt_max .- AABB_min) ./ (AABB_max .- AABB_min))

        for j in 1:nsd # pokud jsem nepřetekl do mínusu nebo za max
            if (I_min[j] < 0)
                I_min[j] = 0
            end
            if (I_max[j] >= N[j])
                I_max[j] = N[j]
            end
        end

        x₁, x₂, x₃ = Xt[:, 1], Xt[:, 2], Xt[:, 3]
        Et = [x₂ - x₁, x₃ - x₂, x₁ - x₃] # vektory hran trojúhelníku
        n = cross(Et[1], Et[2]) # normála
        n = n / norm(n) # jednotková normála
        
        Is = Iterators.product(I_min[1]:I_max[1], I_min[2]:I_max[2], I_min[3]:I_max[3])
        for I in Is
            i = Int(I[3] * (N[1] + 1) * (N[2] + 1) + I[2] * (N[1] + 1) + I[1] + 1)
            v = head[i]
            while v != -1
                x = points[:, v]
                λ = barycentricCoordinates(x₁, x₂, x₃, n, x)
                xₚ = λ[1] * x₁ + λ[2] * x₂ + λ[3] * x₃
                dist_tmp = dot(x - xₚ, n)
                if (abs(dist_tmp) < abs(dist[v]))
                    dist[v] = dist_tmp
                    xp[:, v] = xₚ
                end
                v = next[v]
            end
        end
    end

    for i in 1:length(dist)
        if (abs(dist[i]) > norm(grid.cell_size))
            dist[i] = sign(dist[i]) * norm(grid.cell_size)
        end
    end

    return dist
end


function barycentricCoordinates(x₁, x₂, x₃, n, x)
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
    return λ = A \ b # baricentrické souřadnice
end


############x
using Base.Threads

function evalSignedDistancesOnTriangularMesh(mesh::Mesh, grid::Grid)
    # ... [previous code remains unchanged] ...

    # Create thread-local storage for variables modified inside the loop
    dists = [big * ones(ngp) for _ in 1:nthreads()]
    xps = [zeros(nsd, ngp) for _ in 1:nthreads()]

    @threads for el = 1:nel
        # Access thread-local variables
        dist = dists[threadid()]
        xp = xps[threadid()]

        # ... [rest of the loop code, using thread-local dist and xp] ...

    end

    # Combine results from each thread
    final_dist = reduce(min, dists)

    for i = 1:length(final_dist)
        if abs(final_dist[i]) > norm(grid.cell_size)
            final_dist[i] = sign(final_dist[i]) * norm(grid.cell_size)
        end
    end

    return final_dist
end

