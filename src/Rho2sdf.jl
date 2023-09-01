module Rho2sdf

using LinearAlgebra
using Statistics
using DelimitedFiles
using Einsum
using BenchmarkTools
using ProgressMeter

struct Mesh
    X::Matrix{Float64} # vector of nodes positions
    IEN::Matrix{Int64} # ID element -> ID nodes
    INE::Vector{Vector{Int64}} # ID node -> ID elements
    ISN::Vector{Vector{Int64}} # connectivity face - edges
    nsd::Int64 # number of spacial dimensions
    nnp::Int64 # number of all points
    nen::Int64 # number of element nodes
    nel::Int64 # number of all elements
    nes::Int64 # number of element segments (faces)
    nsn::Int64 # number of face nodes

    function Mesh(X::Vector{Vector{Float64}}, IEN::Vector{Vector{Int64}})
        IEN = reduce(hcat, IEN)
        INE = nodeToElementConnectivity(X, IEN)
        ISN = [
            [1, 4, 3, 2], # nodes define face 1
            [1, 2, 6, 5],
            [2, 3, 7, 6],
            [3, 4, 8, 7],
            [4, 1, 5, 8],
            [5, 6, 7, 8],
        ]
        X = reduce(hcat, X)
        nsd = size(X, 1)
        nnp = size(X, 2)
        nen = size(IEN, 1)
        nel = size(IEN, 2)
        nes = length(ISN)
        nsn = length(ISN[1])
        return new(X, IEN, INE, ISN, nsd, nnp, nen, nel, nes, nsn)
    end
end


struct Grid
    AABB_min::Vector{Float64}
    AABB_max::Vector{Float64}
    N::Vector{Int64}
    cell_size::Vector{Float64}
    ngp::Int64 # Number of Grid Points

    function Grid(
        AABB_min::Vector{Float64},
        AABB_max::Vector{Float64},
        N::Vector{Int64},
        margineCells::Int64 = 3,
    )

        cell_size = (AABB_max .- AABB_min) ./ N

        AABB_min = AABB_min .- margineCells * cell_size
        AABB_max = AABB_max .+ margineCells * cell_size

        AABB_size = AABB_max .- AABB_min

        N = N .+ 2 * margineCells
        ngp = prod(N .+ 1)

        return new(AABB_min, AABB_max, N, cell_size, ngp)

    end
end

mutable struct LinkedList # rozdělení pravidelné sítě na regiony
    grid::Grid
    head::Vector{Int64}
    next::Vector{Int64}
    N::Vector{Float64}

    function LinkedList(grid::Grid, X::Matrix{Float64})

        np = size(X, 2) # Number of Points
        N = grid.N
        AABB_min = grid.AABB_min
        AABB_max = grid.AABB_max

        head = -1 * ones(Int64, prod(N .+ 1))
        next = -1 * ones(Int64, np)

        I = floor.(N .* (X .- AABB_min) ./ (AABB_max .- AABB_min))
        I =
            Int.(
                I[3, :] .* (N[1] + 1) * (N[2] + 1) .+ I[2, :] .* (N[1] + 1) .+
                I[1, :] .+ 1,
            )


        for i = 1:np
            next[i] = head[I[i]]
            head[I[i]] = i
        end

        return new(grid, head, next, N)
    end
end

function nodeToElementConnectivity(
    X::Vector{Vector{Float64}},
    IEN::Matrix{Int64},
)::Vector{Vector{Int64}}
    INE = [Vector{Int64}() for _ = 1:length(X)]
    for el = 1:size(IEN, 2)
        for i = 1:size(IEN, 1)
            push!(INE[IEN[i, el]], el)
        end
    end
    return INE
end

function elementToNodalValues(
    mesh::Mesh,
    elVals::Vector{Float64},
)::Vector{Float64}
    nodalVals = zeros(Float64, mesh.nnp)

    for A = 1:mesh.nnp
        count = 0
        for el in mesh.INE[A]
            nodalVals[A] += elVals[el]
            count += 1
        end
        nodalVals[A] /= count
    end

    return nodalVals
end

function getMesh_AABB(X::Matrix{Float64})
    X_min = vec(minimum(X, dims = 2))
    X_max = vec(maximum(X, dims = 2))
    return X_min, X_max
end

function generateGridPoints(grid::Grid)::Matrix{Float64}

    X = zeros(3, grid.ngp)
    a = 1
    for k = 0:grid.N[3]
        for j = 0:grid.N[2]
            for i = 0:grid.N[1]
                X[:, a] = grid.AABB_min .+ grid.cell_size .* [i, j, k]
                a += 1
            end
        end
    end
    return X
end

function extractSurfaceTriangularMesh(mesh::Mesh, ρₙ::Vector{Float64})::Mesh
    X = mesh.X
    IEN = mesh.IEN # ID element -> nodes
    INE = mesh.INE # ID node -> ID elements
    ISN = mesh.ISN # connectivity face - edges
    # nsd = mesh.nsd # number of special dimension
    nel = mesh.nel # number of elements
    nes = mesh.nes # number of element segments (faces) 6
    nsn = mesh.nsn # number of segment nodes (kolik má stěna uzlů) 4

    X_new = Vector{Vector{Float64}}()
    push!(X_new, vec(X[:, 1]))

    IEN_new = Vector{Vector{Int64}}()

    for el = 1:nel
        # ρₑ = ρₙ[IEN[:, el]]
        commonEls = []
        for sg = 1:nes # 1:6 je face součástí pouze jednoho elementu?
            commonEls = INE[IEN[ISN[sg][1], el]]
            for a = 2:nsn # 2:4
                idx = findall(in(INE[IEN[ISN[sg][a], el]]), commonEls) # mají uzly (jedné stěny) společný pouze jeden element
                commonEls = commonEls[idx]
            end

            if (length(commonEls) == 1) # is a part of the outer boundary of the body

                Xs = X[:, IEN[ISN[sg], el]]
                Xc = mean(Xs, dims = 2) # center of the outside face

                for a = 1:nsn

                    IEN_el = zeros(Int64, 3)

                    x₁ = vec(Xs[:, a])
                    x₂ = vec(Xs[:, (a%nsn)+1])
                    x₃ = vec(Xc)

                    Xt = [x₁, x₂, x₃]
                    # Xt = reduce(hcat, Xt)

                    for i = 1:3
                        a = findfirst(x -> norm(x - Xt[i]) < 1.0e-5, X_new)
                        if (a === nothing)
                            push!(X_new, Xt[i])
                            IEN_el[i] = length(X_new)
                        else
                            IEN_el[i] = a[1]
                        end
                    end
                    push!(IEN_new, IEN_el)
                end # a = 1:nsn
            end # if (length(commonEls) == 1)
        end # sg
    end # el
    return Mesh(X_new, IEN_new)
end

function computePseudoNormals(mesh::Mesh)
    verticesPseudoNormals = [zeros(mesh.nsd) for _ = 1:mesh.nnp]
    edgesPseudoNormals = [[zeros(mesh.nsd) for _ = 1:3] for _ = 1:mesh.nel]


    for el = 1:mesh.nel
        As = mesh.IEN[:, el]
        Xt = mesh.X[:, As]
        n = cross(Xt[:, 2] - Xt[:, 1], Xt[:, 3] - Xt[:, 1])
        n = n / norm(n)

        for i = 1:3
            j = mod(i, 3) + 1
            el_adj = setdiff(intersect(mesh.INE[As[i]], mesh.INE[As[j]]), el)
            As_adj = mesh.IEN[:, el_adj]
            Xt_adj = mesh.X[:, As_adj]
            n_adj =
                cross(Xt_adj[:, 2] - Xt_adj[:, 1], Xt_adj[:, 3] - Xt_adj[:, 1])
            n_adj = π * n_adj / norm(n_adj)
            edgesPseudoNormals[el][i] = n + n_adj
        end
    end

    for A = 1:mesh.nnp
        pseudo_n = zeros(3)
        for k = 1:length(mesh.INE[A])
            el = mesh.INE[A][k]
            As = mesh.IEN[:, el]
            Xt = mesh.X[:, As]
            i = findfirst(x -> x == A, As)
            B = setdiff([1, 2, 3], i)
            e₁ = Xt[:, B[1]] - Xt[:, i]
            e₂ = Xt[:, B[2]] - Xt[:, i]
            α = acos(dot(e₁ / norm(e₁), e₂ / norm(e₂)))
            n = cross(Xt[:, 2] - Xt[:, 1], Xt[:, 3] - Xt[:, 1])
            n = n / norm(n)
            pseudo_n += α * n
        end
        verticesPseudoNormals[A] = pseudo_n
    end

    return verticesPseudoNormals, edgesPseudoNormals
end

function evalSignedDiscancesOnTriangularMesh(mesh::Mesh, grid::Grid)

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
    # INE = mesh.INE
    nsd = mesh.nsd
    nel = mesh.nel

    ngp = grid.ngp
    big = 1.0e10
    dist = big * ones(ngp)
    xp = zeros(nsd, ngp)

    for el = 1:nel
        Xt = X[:, IEN[:, el]] # souřadnice trojúhelníku
        Xt_min = minimum(Xt, dims = 2) .- δ # trojúhelník vložený do BB a přifouknu o deltu
        Xt_max = maximum(Xt, dims = 2) .+ δ

        I_min = floor.(N .* (Xt_min .- AABB_min) ./ (AABB_max .- AABB_min)) # index oddílu
        I_max = floor.(N .* (Xt_max .- AABB_min) ./ (AABB_max .- AABB_min))

        for j = 1:nsd # pokud jsem nepřetekl do mínusu nebo za max
            if (I_min[j] < 0)
                I_min[j] = 0
            end
            if (I_max[j] >= N[j])
                I_max[j] = N[j]
            end
        end

        x₁ = Xt[:, 1]
        x₂ = Xt[:, 2]
        x₃ = Xt[:, 3]

        Et = Vector{Vector{Float64}}() # vektory hran trojúhelníku
        push!(Et, Xt[:, 2] - Xt[:, 1])
        push!(Et, Xt[:, 3] - Xt[:, 2])
        push!(Et, Xt[:, 1] - Xt[:, 3])

        n = cross(Et[1], Et[2]) # normála
        n = n / norm(n) # jednotková normála

        Is = Iterators.product(
            I_min[1]:I_max[1],
            I_min[2]:I_max[2],
            I_min[3]:I_max[3],
        )
        for I ∈ Is
            i = Int(
                I[3] * (N[1] + 1) * (N[2] + 1) + I[2] * (N[1] + 1) + I[1] + 1,
            )
            v = head[i]
            while v != -1
                x = points[:, v]
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
                λ = A \ b # baricentrické souřadnice

                xₚ = zeros(nsd)
                # kam jsem se promítl?
                isFace = false
                isEdge = false
                isVertex = false
                if (minimum(λ) >= 0.0) # xₚ is in the triangle el
                    xₚ = λ[1] * x₁ + λ[2] * x₂ + λ[3] * x₃
                    dist_tmp = dot(x - xₚ, n)
                    if (abs(dist_tmp) < abs(dist[v]))
                        dist[v] = dist_tmp
                        isFace = true
                        xp[:, v] = xₚ
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
                                xp[:, v] = xₚ
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
                        xp[:, v] = xₚ
                    end
                end
                v = next[v]
            end
        end
    end
    # dist = marchingCubes(dist, N.+1, big)

    for i = 1:length(dist)
        if (abs(dist[i]) > norm(grid.cell_size))
            dist[i] = sign(dist[i]) * norm(grid.cell_size)
        end
    end

    return dist
end

function evalSignedDiscances(
    mesh::Mesh,
    grid::Grid,
    ρₙ::Vector{Float64},
    ρₜ::Float64,
)

    points = generateGridPoints(grid)
    linkedList = LinkedList(grid, points)

    head = linkedList.head
    next = linkedList.next
    N = linkedList.grid.N
    AABB_min = linkedList.grid.AABB_min
    AABB_max = linkedList.grid.AABB_max
    δ = 5.1 * grid.cell_size
    X = mesh.X
    IEN = mesh.IEN
    INE = mesh.INE
    ISN = mesh.ISN
    nsd = mesh.nsd
    nel = mesh.nel
    nes = mesh.nes
    nsn = mesh.nsn

    ngp = grid.ngp
    big = 1.0e10
    dist = big * ones(ngp)
    xp = zeros(nsd, ngp)

    for el = 1:nel
        ρₑ = ρₙ[IEN[:, el]]

        ρₑ_min = minimum(ρₑ)
        ρₑ_max = maximum(ρₑ)
        if (ρₑ_min >= ρₜ)
            continue # PRO PŘESKAKUJE HRANIČNÍ ELEMENTY (pouze pro ladění kodu)
            commonEls = []
            for sg = 1:nes
                commonEls = INE[IEN[mesh.ISN[sg][1], el]]
                for a = 2:nsn
                    idx = findall(in(INE[IEN[ISN[sg][a], el]]), commonEls)
                    commonEls = commonEls[idx]
                end

                # Mám vnější segment a chci jeho pseudonormály,...

                if (length(commonEls) == 1) # is a part of the outer boundary of the body
                    Xs = X[:, IEN[ISN[sg], el]]
                    Xc = mean(Xs, dims = 2)

                    for a = 1:nsn



                        As = IEN[ISN[sg][a], el]
                        adj_els = INE[As]
                        for adj_el in adj_els
                            common_adj_els = []
                            for adj_sg = 1:nes
                                common_adj_els =
                                    INE[IEN[mesh.ISN[adj_sg][1], adj_el]]
                                for b = 2:nsn
                                    idx = findall(
                                        in(INE[IEN[ISN[adj_sg][b], adj_el]]),
                                        common_adj_els,
                                    )
                                    common_adj_els = common_adj_els[idx]
                                end

                                if (
                                    length(common_adj_els) == 1 &&
                                    in(As, IEN[mesh.ISN[adj_sg], adj_el])
                                )
                                    # println("Adjacent element")

                                    adj_Xs = X[:, IEN[ISN[adj_sg], adj_el]]
                                    adj_Xc = mean(adj_Xs, dims = 2)

                                    as = indexin(
                                        As,
                                        IEN[mesh.ISN[adj_sg], adj_el],
                                    )
                                    a_prev = ((as + nsn - 1 - 1) % nsn) + 1
                                    a_next = ((as + nsn + 1 - 1) % nsn) + 1

                                    x_prev = X[
                                        :,
                                        IEN[mesh.ISN[adj_sg][a_prev], adj_el],
                                    ]
                                    xs = X[:, IEN[mesh.ISN[adj_sg][as], adj_el]]
                                    x_next = X[
                                        :,
                                        IEN[mesh.ISN[adj_sg][a_next], adj_el],
                                    ]

                                    Et = Vector{Vector{Float64}}()
                                    push!(Et, Xt[:, 2] - Xt[:, 1])
                                    push!(Et, Xt[:, 3] - Xt[:, 2])
                                    push!(Et, Xt[:, 1] - Xt[:, 3])

                                    n = cross(Et[1], Et[2])
                                    n = n / norm(n)



                                end
                            end
                        end





                        x₁ = Xs[:, a]
                        x₂ = Xs[:, (a%nsn)+1]
                        x₃ = Xc

                        Xt = [x₁, x₂, x₃]
                        Xt = reduce(hcat, Xt)

                        Xt_min = minimum(Xt, dims = 2) .- δ
                        Xt_max = maximum(Xt, dims = 2) .+ δ

                        I_min =
                            floor.(
                                N .* (Xt_min .- AABB_min) ./
                                (AABB_max .- AABB_min),
                            )
                        I_max =
                            floor.(
                                N .* (Xt_max .- AABB_min) ./
                                (AABB_max .- AABB_min),
                            )

                        for j = 1:nsd
                            if (I_min[j] < 0)
                                I_min[j] = 0
                            end
                            if (I_max[j] >= N[j])
                                I_max[j] = N[j]
                            end
                        end

                        Et = Vector{Vector{Float64}}()
                        push!(Et, Xt[:, 2] - Xt[:, 1])
                        push!(Et, Xt[:, 3] - Xt[:, 2])
                        push!(Et, Xt[:, 1] - Xt[:, 3])

                        n = cross(Et[1], Et[2])
                        n = n / norm(n)

                        Is = Iterators.product(
                            I_min[1]:I_max[1],
                            I_min[2]:I_max[2],
                            I_min[3]:I_max[3],
                        )
                        for I ∈ Is
                            ii = Int(
                                I[3] * (N[1] + 1) * (N[2] + 1) +
                                I[2] * (N[1] + 1) +
                                I[1] +
                                1,
                            )
                            v = head[ii]
                            while v != -1
                                x = points[:, v]
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
                                λ = A \ b

                                xₚ = zeros(nsd)
                                isFace = false
                                isEdge = false
                                isVertex = false
                                if (minimum(λ) >= 0.0) # xₚ is in the triangle el
                                    xₚ = λ[1] * x₁ + λ[2] * x₂ + λ[3] * x₃
                                    dist_tmp = dot(x - xₚ, n)
                                    if (abs(dist_tmp) < abs(dist[v]))
                                        dist[v] = dist_tmp
                                        isFace = true
                                        xp[:, v] = xₚ
                                    end
                                    # else
                                    #     for j = 1:3
                                    #         L = norm(Et[j])
                                    #         xᵥ = Xt[:, j]
                                    #         P = dot(x - xᵥ, Et[j] / L)
                                    #         if (P >= 0 && P <= L)
                                    #             xₚ = xᵥ + (Et[j] / L) * P
                                    #             #n_edge = EPN[el][j]
                                    #             dist_tmp = sign(dot(x - xₚ, n)) * norm(x - xₚ)
                                    #
                                    #             if (abs(dist_tmp) < abs(dist[v]))
                                    #                 dist[v] = dist_tmp
                                    #                 isVertex = true
                                    #                 xp[:, v] = xₚ
                                    #             end
                                    #         end
                                    #     end
                                end
                                # if (isFace == false && isEdge == false)
                                #     dist_tmp, idx = findmin([
                                #         norm(x - x₁),
                                #         norm(x - x₂),
                                #         norm(x - x₃),
                                #     ])
                                #     xₚ = Xt[:, idx]
                                #     #n_vertex = VPN[IEN[idx, el]]
                                #     dist_tmp =
                                #         dist_tmp * sign(dot(x - xₚ, n))
                                #     if (abs(dist_tmp) < abs(dist[v]))
                                #         dist[v] = dist_tmp
                                #         isVertex = true
                                #         xp[:, v] = xₚ
                                #     end
                                # end
                                v = next[v]
                            end
                        end
                    end
                    # Subdivide quad into two (or four?) triangles
                end
            end
        else # (ρₑ_min < ρₜ)
            # continue
            if (ρₑ_max > ρₜ)
                # println("Hranice prochází elementem...")

                Xₑ = X[:, IEN[:, el]]
                Xₑ_min = minimum(Xₑ, dims = 2) .- δ
                Xₑ_max = maximum(Xₑ, dims = 2) .+ δ

                I_min =
                    floor.(N .* (Xₑ_min .- AABB_min) ./ (AABB_max .- AABB_min),)
                I_max =
                    floor.(N .* (Xₑ_max .- AABB_min) ./ (AABB_max .- AABB_min),)

                for j = 1:nsd
                    if (I_min[j] < 0)
                        I_min[j] = 0
                    end
                    if (I_max[j] >= N[j])
                        I_max[j] = N[j]
                    end
                end

                Is = Iterators.product(
                    I_min[1]:I_max[1],
                    I_min[2]:I_max[2],
                    I_min[3]:I_max[3],
                )
                for I ∈ Is
                    ii = Int(
                        I[3] * (N[1] + 1) * (N[2] + 1) +
                        I[2] * (N[1] + 1) +
                        I[1] +
                        1,
                    )
                    v = head[ii]
                    while v != -1
                        x = points[:, v]

                        xₚ = [0.0, 0.0, 0.0] # odhad?
                        n = [0.0, 0.0, 0.0]
                        Ξ = [0.0, 0.0, 0.0]
                        λ = 1.0
                        Ξ_tol = 1e-2
                        Ξ_norm = 2 * Ξ_tol
                        r_tol = 1e-2
                        r_norm = 2 * r_tol
                        niter = 100
                        iter = 1

                        while (Ξ_norm ≥ Ξ_tol && iter ≤ niter)# || r_norm ≥ r_tol)

                            H, d¹N_dξ¹, d²N_dξ², d³N_dξ³ = sfce(Ξ) # tvarové funkce a jejich derivace

                            xₚ = Xₑ * H # Xₑ je souřadnice uzlů, H tvarové funkce -> souřadnice bodu KONTROLA!!
                            dx_dΞ = Xₑ * d¹N_dξ¹
                            dρ_dΞ = d¹N_dξ¹' * ρₑ

                            d²ρ_dΞ² = zeros(Float64, 3, 3)
                            for k = 1:length(H)
                                d²ρ_dΞ² += ρₑ[k] * d²N_dξ²[k, :, :]
                            end
                            
                            d³ρ_dΞ³ = zeros(Float64, 3, 3, 3)
                            for k = 1:length(H)
                                d³ρ_dΞ³ += ρₑ[k] * d³N_dξ³[k, :, :, :]
                            end

                            norm_dρ_dΞ = norm(dρ_dΞ) # ok
                            n = dρ_dΞ / norm_dρ_dΞ # ok

                            dn_dΞ = zeros(Float64, 3, 3)
                            @einsum dn_dΞ[i, j] := # dve tečky když matice neni alokovaná
                                d²ρ_dΞ²[i, j] / norm_dρ_dΞ -
                                (dρ_dΞ[i] * d²ρ_dΞ²[j, k] * dρ_dΞ[k]) / 
                                norm_dρ_dΞ^3  # ok

                            dd_dΞ = zeros(Float64, 3)
                            @einsum dd_dΞ[i] :=
                                -dx_dΞ[i, k] * n[k] +
                                (x[k] - xₚ[k]) * dn_dΞ[k, i] # ok


                            dL_dΞ = zeros(Float64, 3)
                            @einsum dL_dΞ[i] := dd_dΞ[i] + λ * dρ_dΞ[i] # ok

                            ρ = H ⋅ ρₑ # hustota v bodě
                            dL_dλ = ρ - ρₜ # ok

                            d²x_dΞ² = zeros(Float64, 3, 3, 3)
                            @einsum d²x_dΞ²[i, j, k] :=
                                Xₑ[i, m] * d²N_dξ²[m, j, k] # ok


                            d²n_dΞ² = zeros(Float64, 3, 3, 3)
                            @einsum d²n_dΞ²[i, j, k] :=
                                d³ρ_dΞ³[i, j, k] / norm_dρ_dΞ -
                                (d²ρ_dΞ²[i, j] * d²ρ_dΞ²[k, m] * dρ_dΞ[m]) /
                                norm_dρ_dΞ^3 -
                                d²ρ_dΞ²[i, j] * d²ρ_dΞ²[k, m] * dρ_dΞ[m] /
                                norm_dρ_dΞ^3 +
                                dρ_dΞ[i] * (
                                    d²ρ_dΞ²[j, m] * d²ρ_dΞ²[m, k] +
                                    d³ρ_dΞ³[j, k, m] * dρ_dΞ[m]
                                ) / norm_dρ_dΞ^3 +
                                3 * (
                                    dρ_dΞ[i] *
                                    dρ_dΞ[m] *
                                    d²ρ_dΞ²[m, j] *
                                    dρ_dΞ[l] *
                                    d²ρ_dΞ²[l, k]
                                ) / norm_dρ_dΞ^5 # ok

                            d²d_dΞ² = zeros(Float64, 3, 3)
                            @einsum d²d_dΞ²[i, j] :=
                                -d²x_dΞ²[i, j, k] * n[k] -
                                2 * dx_dΞ[i, k] * dn_dΞ[k, j] +
                                (x[k] - xₚ[k]) * d²n_dΞ²[k, i, j]
                            d²L_dΞ² = d²d_dΞ² + d²ρ_dΞ² * λ

                            d²L_dΞdλ = dρ_dΞ
                            d²L_dλ² = 0.0

                            K = [
                                d²L_dΞ² d²L_dΞdλ
                                d²L_dΞdλ' d²L_dλ²
                            ]
                            r = [dL_dΞ; dL_dλ] # r1 ... r4 vyčíslím v bode xi_temp
                            # r(xi_tem + )
                            ########################################
                            println("el:",el)
                            Ξ_tmp = zeros(Float64,4) # bod pro linearizaci ξ₁,ξ₂,ξ₃,λ
                            K_diff = zeros(Float64,4,4)
                            h = 1e-6    
                            sign = [1 -1]
                            ϵ = [h -h]
                            for kk in 1:4
                                K_diff_col = zeros(Float64,4)
                                for ll in 1:2
                                    ΔΞ_tmp = zeros(Float64,4)
                                    ΔΞ_tmp[kk] = ϵ[ll]

                                    H, d¹N_dξ¹, d²N_dξ², d³N_dξ³ = sfce(Ξ_tmp + ΔΞ_tmp)

                                    xₚ = Xₑ * H
                                    dx_dΞ = Xₑ * d¹N_dξ¹

                                    dρ_dΞ = d¹N_dξ¹' * ρₑ

                                    d²ρ_dΞ² = zeros(Float64,3,3)
                                    for m in 1:length(H)
                                        d²ρ_dΞ² +=  ρₑ[m] * d²N_dξ²[m,:,:]
                                    end


                                    norm_dρ_dΞ = norm(dρ_dΞ)
                                    n = dρ_dΞ / norm_dρ_dΞ

                                    dn_dΞ = zeros(Float64,3,3)
                                    @einsum dn_dΞ[i,j] := d²ρ_dΞ²[i,j] / norm_dρ_dΞ - (dρ_dΞ[i] * d²ρ_dΞ²[j,k] * dρ_dΞ[k] ) / norm_dρ_dΞ^(3/2)
                                    
                                    dd_dΞ = zeros(Float64,3)
                                    @einsum dd_dΞ[i] := -dx_dΞ[i,k] * n[k] + (x[k] - xₚ[k]) * dn_dΞ[k,i]

                                    dL_dΞ = zeros(Float64,3)
                                    @einsum dL_dΞ[i] := dd_dΞ[i] + (Ξ_tmp[end]+ΔΞ_tmp[end])*dρ_dΞ[i]
                                    
                                    ρ = H⋅ρₑ
                                    dL_dλ = ρ - ρₜ

                                    r_tmp = [dL_dΞ; dL_dλ] 
                                    
                                    K_diff_col = K_diff_col + sign[ll] .* r_tmp
                                end
                                K_diff[kk,:] =  K_diff_col ./ (2*h)
                            end
                            # println("K:",K)
                            # println("K_diff:",K_diff)
                            K = K_diff
                            ########################################
                            println("el:",el)

                            r_norm = norm(r)
                            Λ = real.(eigvals(K))
                            (Λ_min, idx_min) = findmin(Λ)

                            if (Λ_min < 1.0e-10)
                                Φ = real.(eigvecs(K))
                                idx = [1, 2, 3, 4]
                                deleteat!(idx, idx_min)
                                Φ = Φ[:, idx]
                                ΔΞ̃_and_Δλ̃ = 1.0 ./ Λ[idx] .* (Φ' * r)
                                ΔΞ_and_Δλ = Φ * ΔΞ̃_and_Δλ̃
                            else
                                ΔΞ_and_Δλ = K \ -r
                            end

                            max_abs_Ξ = maximum(abs.(ΔΞ_and_Δλ[1:end-1]))
                            if (max_abs_Ξ > 1.0)
                                ΔΞ_and_Δλ[1:end-1] =
                                    ΔΞ_and_Δλ[1:end-1] / max_abs_Ξ
                            end


                            Ξ = Ξ + ΔΞ_and_Δλ[1:end-1]
                            λ = λ + ΔΞ_and_Δλ[end]

                            Ξ_norm = norm(ΔΞ_and_Δλ)

                            iter = iter + 1
                        end

                        if (maximum(abs.(Ξ)) <= 1.0) # xₚ is in the element
                            dist_tmp = dot(x - xₚ, n)
                            if (abs(dist_tmp) < abs(dist[v]))
                                dist[v] = dist_tmp
                                xp[:, v] = xₚ
                            end
                        end

                        v = next[v]

                    end
                end
            end
        end
    end
    # dist = marchingCubes(dist, N.+1, big)

    X = Vector{Float64}(undef, 3)
    Xp = Vector{Float64}(undef, 3)
    for i = 1:size(xp, 2)
        if (sum(abs.(xp[:, i])) > 1.0e-10)
            X = [X points[:, i]]
            Xp = [Xp xp[:, i]]
        end
    end
    X = [X Xp]
    X = [X[:, i] for i = 1:size(X, 2)]

    nnp = size(Xp, 2)
    IEN = [[i; i + nnp] for i = 1:nnp]

    exportToVTU("xp.vtu", X, IEN)

    open("xp.csv", "w") do io
        writedlm(io, ['x' 'y' 'z'], ',')
        writedlm(io, xp', ',')
    end

    return dist
end

function exportToVTU(
    fileName::String,
    X::Vector{Vector{Float64}},
    IEN::Vector{Vector{Int64}},
)

    nnp = length(X)
    nsd = length(X[1])

    nel = length(IEN)
    nen = length(IEN[1])

    io = open(fileName, "w")

    println(
        io,
        "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">",
    )
    println(io, "  <UnstructuredGrid>")
    println(
        io,
        "    <Piece NumberOfPoints=\"",
        nnp,
        "\" NumberOfCells=\"",
        nel,
        "\">",
    )

    println(io, "	  <Points>")
    println(
        io,
        "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">",
    )
    for A = 1:nnp
        print(io, "          ")
        for i = 1:nsd
            if (abs(X[A][i]) < 1.0e-20)
                X[A][i] = 0.0
            end
            print(io, " ", X[A][i])
        end
        print(io, "\n")
    end
    println(io, "        </DataArray>")
    println(io, "	  </Points>")

    println(io, "      <Cells>")
    println(
        io,
        "		  <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">",
    )

    VTK_CODE = 5
    for el = 1:nel
        print(io, "         ")
        for a = 1:nen
            print(io, " ", IEN[el][a] - 1)
        end
        print(io, "\n")
    end

    println(io, "        </DataArray>")
    println(
        io,
        "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">",
    )
    for i = nen:nen:nen*nel
        println(io, "          ", i)
    end
    println(io, "        </DataArray>")
    println(
        io,
        "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">",
    )
    for el = 1:nel
        println(io, "          ", VTK_CODE)
    end
    println(io, "        </DataArray>")
    println(io, "      </Cells>")

    println(io, "    </Piece>")
    println(io, "  </UnstructuredGrid>")
    println(io, "</VTKFile>")

    close(io)
end

function exportStructuredPointsToVTK(
    fileName::String,
    grid::Grid,
    vals::Vector{Float64},
    valLabel::String,
)

    dim = grid.N .+ 1
    org = grid.AABB_min
    spacing = grid.cell_size

    io = open(fileName, "w")

    write(io, "# vtk DataFile Version 1.0\n")
    write(
        io,
        "Texture map for thresholding data (use boolean textures for 2D map)\n",
    )

    write(io, "ASCII\n\n")
    write(io, "DATASET STRUCTURED_POINTS\n")
    dim_x = dim[1]
    dim_y = dim[2]
    dim_z = dim[3]
    write(io, "DIMENSIONS $dim_x $dim_y $dim_z\n") # dimenze pravidelné sítě

    spacing_x = spacing[1]
    spacing_y = spacing[2]
    spacing_z = spacing[3]
    write(io, "SPACING $spacing_x $spacing_y $spacing_z\n") # krok sítě ve 3 směrech

    org_x = org[1]
    org_y = org[2]
    org_z = org[3]
    write(io, "ORIGIN $org_x $org_y $org_z\n\n") # souřadnice počátku

    n = prod(dim)
    write(io, "POINT_DATA $n\n") # počet uzlů pravidelné sítě
    write(io, "SCALARS $valLabel float 1\n") # druh hodnoty v uzlech (vzdálenost)
    write(io, "LOOKUP_TABLE default\n") # ??

    for val ∈ vals
        write(io, "$val\n") # hodnota vzdálenostní funkce v jednotlivých bodech
    end
    close(io)
end

function marchingCubes(
    dists::Vector{Float64},
    Nv::Vector{Int64},
    big::Float64,
)::Vector{Float64}

    idx = findall(x -> x < 0.0, dists)
    while true
        if isempty(idx)
            break
        end
        idx_new = Vector{Int64}()
        sameSign = false
        for I in idx
            I₁ = mod(mod(I - 1, Nv[1] * Nv[2]) + 1 - 1, Nv[1]) + 1
            I₂ = Int(ceil((mod(I - 1, Nv[1] * Nv[2]) + 1) / Nv[1]))
            I₃ = Int(ceil(I / (Nv[1] * Nv[2])))
            for i = I₁-1:I₁+1
                for j = I₂-1:I₂+1
                    for k = I₃-1:I₃+1
                        if minimum([i, j, k]) <= 0 ||
                           minimum(Nv - [i, j, k]) < 0
                            continue
                        end
                        I_adj = Int(
                            (k - 1) .* Nv[1] .* Nv[2] .+ (j - 1) .* Nv[1] .+
                            (i - 1) .+ 1,
                        )

                        if (sign(dists[I_adj]) == sign(dists[I]))
                            sameSign = true
                        end

                        if (abs(dists[I_adj] - big) < 1.0e-5)
                            dists[I_adj] = -big
                            push!(idx_new, I_adj)
                        end
                    end
                end
            end
        end
        if sameSign == false
            #dists[I] = -dists[I]
        end
        idx = idx_new
    end




    idx = findall(x -> x > 0.0, dists)
    while true
        if isempty(idx)
            break
        end
        idx_new = Vector{Int64}()

        sameSign = false
        for I in idx
            I₁ = mod(mod(I - 1, Nv[1] * Nv[2]) + 1 - 1, Nv[1]) + 1
            I₂ = Int(ceil((mod(I - 1, Nv[1] * Nv[2]) + 1) / Nv[1]))
            I₃ = Int(ceil(I / (Nv[1] * Nv[2])))
            for i = I₁-1:I₁+1
                for j = I₂-1:I₂+1
                    for k = I₃-1:I₃+1
                        if minimum([i, j, k]) <= 0 ||
                           minimum(Nv - [i, j, k]) < 0
                            continue
                        end
                        I_adj = Int(
                            (k - 1) .* Nv[1] .* Nv[2] .+ (j - 1) .* Nv[1] .+
                            (i - 1) .+ 1,
                        )

                        if (sign(dists[I_adj]) == sign(dists[I]))
                            sameSign = true
                        end
                        if (abs(dists[I_adj] + big) < 1.0e-5)
                            dists[I_adj] = +big
                            push!(idx_new, I_adj)
                        end
                    end
                end
            end
        end
        if sameSign == false
            #dists[I] = -dists[I]
        end
        idx = idx_new
    end



    return dists
end

function marchingCubesssss()
    # Kazdému vrcholu nastavit 0/1 podle toho, zda je < nebo >
    edgeTable = [
        0,
        265,
        515,
        778,
        1030,
        1295,
        1541,
        1804,
        2060,
        2309,
        2575,
        2822,
        3082,
        3331,
        3593,
        3840,
        400,
        153,
        915,
        666,
        1430,
        1183,
        1941,
        1692,
        2460,
        2197,
        2975,
        2710,
        3482,
        3219,
        3993,
        3728,
        560,
        825,
        51,
        314,
        1590,
        1855,
        1077,
        1340,
        2620,
        2869,
        2111,
        2358,
        3642,
        3891,
        3129,
        3376,
        928,
        681,
        419,
        170,
        1958,
        1711,
        1445,
        1196,
        2988,
        2725,
        2479,
        2214,
        4010,
        3747,
        3497,
        3232,
        1120,
        1385,
        1635,
        1898,
        102,
        367,
        613,
        876,
        3180,
        3429,
        3695,
        3942,
        2154,
        2403,
        2665,
        2912,
        1520,
        1273,
        2035,
        1786,
        502,
        255,
        1013,
        764,
        3580,
        3317,
        4095,
        3830,
        2554,
        2291,
        3065,
        2800,
        1616,
        1881,
        1107,
        1370,
        598,
        863,
        85,
        348,
        3676,
        3925,
        3167,
        3414,
        2650,
        2899,
        2137,
        2384,
        1984,
        1737,
        1475,
        1226,
        966,
        719,
        453,
        204,
        4044,
        3781,
        3535,
        3270,
        3018,
        2755,
        2505,
        2240,
        2240,
        2505,
        2755,
        3018,
        3270,
        3535,
        3781,
        4044,
        204,
        453,
        719,
        966,
        1226,
        1475,
        1737,
        1984,
        2384,
        2137,
        2899,
        2650,
        3414,
        3167,
        3925,
        3676,
        348,
        85,
        863,
        598,
        1370,
        1107,
        1881,
        1616,
        2800,
        3065,
        2291,
        2554,
        3830,
        4095,
        3317,
        3580,
        764,
        1013,
        255,
        502,
        1786,
        2035,
        1273,
        1520,
        2912,
        2665,
        2403,
        2154,
        3942,
        3695,
        3429,
        3180,
        876,
        613,
        367,
        102,
        1898,
        1635,
        1385,
        1120,
        3232,
        3497,
        3747,
        4010,
        2214,
        2479,
        2725,
        2988,
        1196,
        1445,
        1711,
        1958,
        170,
        419,
        681,
        928,
        3376,
        3129,
        3891,
        3642,
        2358,
        2111,
        2869,
        2620,
        1340,
        1077,
        1855,
        1590,
        314,
        51,
        825,
        560,
        3728,
        3993,
        3219,
        3482,
        2710,
        2975,
        2197,
        2460,
        1692,
        1941,
        1183,
        1430,
        666,
        915,
        153,
        400,
        3840,
        3593,
        3331,
        3082,
        2822,
        2575,
        2309,
        2060,
        1804,
        1541,
        1295,
        1030,
        778,
        515,
        265,
        0,
    ]


    triTable =
        [
            -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            0 8 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            0 1 9 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            1 8 3 9 8 1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            1 2 10 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            0 8 3 1 2 10 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            9 2 10 0 2 9 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            2 8 3 2 10 8 10 9 8 -1 -1 -1 -1 -1 -1 -1
            3 11 2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            0 11 2 8 11 0 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            1 9 0 2 3 11 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            1 11 2 1 9 11 9 8 11 -1 -1 -1 -1 -1 -1 -1
            3 10 1 11 10 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            0 10 1 0 8 10 8 11 10 -1 -1 -1 -1 -1 -1 -1
            3 9 0 3 11 9 11 10 9 -1 -1 -1 -1 -1 -1 -1
            9 8 10 10 8 11 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            4 7 8 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            4 3 0 7 3 4 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            0 1 9 8 4 7 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            4 1 9 4 7 1 7 3 1 -1 -1 -1 -1 -1 -1 -1
            1 2 10 8 4 7 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            3 4 7 3 0 4 1 2 10 -1 -1 -1 -1 -1 -1 -1
            9 2 10 9 0 2 8 4 7 -1 -1 -1 -1 -1 -1 -1
            2 10 9 2 9 7 2 7 3 7 9 4 -1 -1 -1 -1
            8 4 7 3 11 2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            11 4 7 11 2 4 2 0 4 -1 -1 -1 -1 -1 -1 -1
            9 0 1 8 4 7 2 3 11 -1 -1 -1 -1 -1 -1 -1
            4 7 11 9 4 11 9 11 2 9 2 1 -1 -1 -1 -1
            3 10 1 3 11 10 7 8 4 -1 -1 -1 -1 -1 -1 -1
            1 11 10 1 4 11 1 0 4 7 11 4 -1 -1 -1 -1
            4 7 8 9 0 11 9 11 10 11 0 3 -1 -1 -1 -1
            4 7 11 4 11 9 9 11 10 -1 -1 -1 -1 -1 -1 -1
            9 5 4 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            9 5 4 0 8 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            0 5 4 1 5 0 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            8 5 4 8 3 5 3 1 5 -1 -1 -1 -1 -1 -1 -1
            1 2 10 9 5 4 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            3 0 8 1 2 10 4 9 5 -1 -1 -1 -1 -1 -1 -1
            5 2 10 5 4 2 4 0 2 -1 -1 -1 -1 -1 -1 -1
            2 10 5 3 2 5 3 5 4 3 4 8 -1 -1 -1 -1
            9 5 4 2 3 11 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            0 11 2 0 8 11 4 9 5 -1 -1 -1 -1 -1 -1 -1
            0 5 4 0 1 5 2 3 11 -1 -1 -1 -1 -1 -1 -1
            2 1 5 2 5 8 2 8 11 4 8 5 -1 -1 -1 -1
            10 3 11 10 1 3 9 5 4 -1 -1 -1 -1 -1 -1 -1
            4 9 5 0 8 1 8 10 1 8 11 10 -1 -1 -1 -1
            5 4 0 5 0 11 5 11 10 11 0 3 -1 -1 -1 -1
            5 4 8 5 8 10 10 8 11 -1 -1 -1 -1 -1 -1 -1
            9 7 8 5 7 9 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            9 3 0 9 5 3 5 7 3 -1 -1 -1 -1 -1 -1 -1
            0 7 8 0 1 7 1 5 7 -1 -1 -1 -1 -1 -1 -1
            1 5 3 3 5 7 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            9 7 8 9 5 7 10 1 2 -1 -1 -1 -1 -1 -1 -1
            10 1 2 9 5 0 5 3 0 5 7 3 -1 -1 -1 -1
            8 0 2 8 2 5 8 5 7 10 5 2 -1 -1 -1 -1
            2 10 5 2 5 3 3 5 7 -1 -1 -1 -1 -1 -1 -1
            7 9 5 7 8 9 3 11 2 -1 -1 -1 -1 -1 -1 -1
            9 5 7 9 7 2 9 2 0 2 7 11 -1 -1 -1 -1
            2 3 11 0 1 8 1 7 8 1 5 7 -1 -1 -1 -1
            11 2 1 11 1 7 7 1 5 -1 -1 -1 -1 -1 -1 -1
            9 5 8 8 5 7 10 1 3 10 3 11 -1 -1 -1 -1
            5 7 0 5 0 9 7 11 0 1 0 10 11 10 0 -1
            11 10 0 11 0 3 10 5 0 8 0 7 5 7 0 -1
            11 10 5 7 11 5 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            10 6 5 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            0 8 3 5 10 6 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            9 0 1 5 10 6 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            1 8 3 1 9 8 5 10 6 -1 -1 -1 -1 -1 -1 -1
            1 6 5 2 6 1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            1 6 5 1 2 6 3 0 8 -1 -1 -1 -1 -1 -1 -1
            9 6 5 9 0 6 0 2 6 -1 -1 -1 -1 -1 -1 -1
            5 9 8 5 8 2 5 2 6 3 2 8 -1 -1 -1 -1
            2 3 11 10 6 5 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            11 0 8 11 2 0 10 6 5 -1 -1 -1 -1 -1 -1 -1
            0 1 9 2 3 11 5 10 6 -1 -1 -1 -1 -1 -1 -1
            5 10 6 1 9 2 9 11 2 9 8 11 -1 -1 -1 -1
            6 3 11 6 5 3 5 1 3 -1 -1 -1 -1 -1 -1 -1
            0 8 11 0 11 5 0 5 1 5 11 6 -1 -1 -1 -1
            3 11 6 0 3 6 0 6 5 0 5 9 -1 -1 -1 -1
            6 5 9 6 9 11 11 9 8 -1 -1 -1 -1 -1 -1 -1
            5 10 6 4 7 8 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            4 3 0 4 7 3 6 5 10 -1 -1 -1 -1 -1 -1 -1
            1 9 0 5 10 6 8 4 7 -1 -1 -1 -1 -1 -1 -1
            10 6 5 1 9 7 1 7 3 7 9 4 -1 -1 -1 -1
            6 1 2 6 5 1 4 7 8 -1 -1 -1 -1 -1 -1 -1
            1 2 5 5 2 6 3 0 4 3 4 7 -1 -1 -1 -1
            8 4 7 9 0 5 0 6 5 0 2 6 -1 -1 -1 -1
            7 3 9 7 9 4 3 2 9 5 9 6 2 6 9 -1
            3 11 2 7 8 4 10 6 5 -1 -1 -1 -1 -1 -1 -1
            5 10 6 4 7 2 4 2 0 2 7 11 -1 -1 -1 -1
            0 1 9 4 7 8 2 3 11 5 10 6 -1 -1 -1 -1
            9 2 1 9 11 2 9 4 11 7 11 4 5 10 6 -1
            8 4 7 3 11 5 3 5 1 5 11 6 -1 -1 -1 -1
            5 1 11 5 11 6 1 0 11 7 11 4 0 4 11 -1
            0 5 9 0 6 5 0 3 6 11 6 3 8 4 7 -1
            6 5 9 6 9 11 4 7 9 7 11 9 -1 -1 -1 -1
            10 4 9 6 4 10 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            4 10 6 4 9 10 0 8 3 -1 -1 -1 -1 -1 -1 -1
            10 0 1 10 6 0 6 4 0 -1 -1 -1 -1 -1 -1 -1
            8 3 1 8 1 6 8 6 4 6 1 10 -1 -1 -1 -1
            1 4 9 1 2 4 2 6 4 -1 -1 -1 -1 -1 -1 -1
            3 0 8 1 2 9 2 4 9 2 6 4 -1 -1 -1 -1
            0 2 4 4 2 6 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            8 3 2 8 2 4 4 2 6 -1 -1 -1 -1 -1 -1 -1
            10 4 9 10 6 4 11 2 3 -1 -1 -1 -1 -1 -1 -1
            0 8 2 2 8 11 4 9 10 4 10 6 -1 -1 -1 -1
            3 11 2 0 1 6 0 6 4 6 1 10 -1 -1 -1 -1
            6 4 1 6 1 10 4 8 1 2 1 11 8 11 1 -1
            9 6 4 9 3 6 9 1 3 11 6 3 -1 -1 -1 -1
            8 11 1 8 1 0 11 6 1 9 1 4 6 4 1 -1
            3 11 6 3 6 0 0 6 4 -1 -1 -1 -1 -1 -1 -1
            6 4 8 11 6 8 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            7 10 6 7 8 10 8 9 10 -1 -1 -1 -1 -1 -1 -1
            0 7 3 0 10 7 0 9 10 6 7 10 -1 -1 -1 -1
            10 6 7 1 10 7 1 7 8 1 8 0 -1 -1 -1 -1
            10 6 7 10 7 1 1 7 3 -1 -1 -1 -1 -1 -1 -1
            1 2 6 1 6 8 1 8 9 8 6 7 -1 -1 -1 -1
            2 6 9 2 9 1 6 7 9 0 9 3 7 3 9 -1
            7 8 0 7 0 6 6 0 2 -1 -1 -1 -1 -1 -1 -1
            7 3 2 6 7 2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            2 3 11 10 6 8 10 8 9 8 6 7 -1 -1 -1 -1
            2 0 7 2 7 11 0 9 7 6 7 10 9 10 7 -1
            1 8 0 1 7 8 1 10 7 6 7 10 2 3 11 -1
            11 2 1 11 1 7 10 6 1 6 7 1 -1 -1 -1 -1
            8 9 6 8 6 7 9 1 6 11 6 3 1 3 6 -1
            0 9 1 11 6 7 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            7 8 0 7 0 6 3 11 0 11 6 0 -1 -1 -1 -1
            7 11 6 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            7 6 11 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            3 0 8 11 7 6 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            0 1 9 11 7 6 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            8 1 9 8 3 1 11 7 6 -1 -1 -1 -1 -1 -1 -1
            10 1 2 6 11 7 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            1 2 10 3 0 8 6 11 7 -1 -1 -1 -1 -1 -1 -1
            2 9 0 2 10 9 6 11 7 -1 -1 -1 -1 -1 -1 -1
            6 11 7 2 10 3 10 8 3 10 9 8 -1 -1 -1 -1
            7 2 3 6 2 7 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            7 0 8 7 6 0 6 2 0 -1 -1 -1 -1 -1 -1 -1
            2 7 6 2 3 7 0 1 9 -1 -1 -1 -1 -1 -1 -1
            1 6 2 1 8 6 1 9 8 8 7 6 -1 -1 -1 -1
            10 7 6 10 1 7 1 3 7 -1 -1 -1 -1 -1 -1 -1
            10 7 6 1 7 10 1 8 7 1 0 8 -1 -1 -1 -1
            0 3 7 0 7 10 0 10 9 6 10 7 -1 -1 -1 -1
            7 6 10 7 10 8 8 10 9 -1 -1 -1 -1 -1 -1 -1
            6 8 4 11 8 6 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            3 6 11 3 0 6 0 4 6 -1 -1 -1 -1 -1 -1 -1
            8 6 11 8 4 6 9 0 1 -1 -1 -1 -1 -1 -1 -1
            9 4 6 9 6 3 9 3 1 11 3 6 -1 -1 -1 -1
            6 8 4 6 11 8 2 10 1 -1 -1 -1 -1 -1 -1 -1
            1 2 10 3 0 11 0 6 11 0 4 6 -1 -1 -1 -1
            4 11 8 4 6 11 0 2 9 2 10 9 -1 -1 -1 -1
            10 9 3 10 3 2 9 4 3 11 3 6 4 6 3 -1
            8 2 3 8 4 2 4 6 2 -1 -1 -1 -1 -1 -1 -1
            0 4 2 4 6 2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            1 9 0 2 3 4 2 4 6 4 3 8 -1 -1 -1 -1
            1 9 4 1 4 2 2 4 6 -1 -1 -1 -1 -1 -1 -1
            8 1 3 8 6 1 8 4 6 6 10 1 -1 -1 -1 -1
            10 1 0 10 0 6 6 0 4 -1 -1 -1 -1 -1 -1 -1
            4 6 3 4 3 8 6 10 3 0 3 9 10 9 3 -1
            10 9 4 6 10 4 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            4 9 5 7 6 11 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            0 8 3 4 9 5 11 7 6 -1 -1 -1 -1 -1 -1 -1
            5 0 1 5 4 0 7 6 11 -1 -1 -1 -1 -1 -1 -1
            11 7 6 8 3 4 3 5 4 3 1 5 -1 -1 -1 -1
            9 5 4 10 1 2 7 6 11 -1 -1 -1 -1 -1 -1 -1
            6 11 7 1 2 10 0 8 3 4 9 5 -1 -1 -1 -1
            7 6 11 5 4 10 4 2 10 4 0 2 -1 -1 -1 -1
            3 4 8 3 5 4 3 2 5 10 5 2 11 7 6 -1
            7 2 3 7 6 2 5 4 9 -1 -1 -1 -1 -1 -1 -1
            9 5 4 0 8 6 0 6 2 6 8 7 -1 -1 -1 -1
            3 6 2 3 7 6 1 5 0 5 4 0 -1 -1 -1 -1
            6 2 8 6 8 7 2 1 8 4 8 5 1 5 8 -1
            9 5 4 10 1 6 1 7 6 1 3 7 -1 -1 -1 -1
            1 6 10 1 7 6 1 0 7 8 7 0 9 5 4 -1
            4 0 10 4 10 5 0 3 10 6 10 7 3 7 10 -1
            7 6 10 7 10 8 5 4 10 4 8 10 -1 -1 -1 -1
            6 9 5 6 11 9 11 8 9 -1 -1 -1 -1 -1 -1 -1
            3 6 11 0 6 3 0 5 6 0 9 5 -1 -1 -1 -1
            0 11 8 0 5 11 0 1 5 5 6 11 -1 -1 -1 -1
            6 11 3 6 3 5 5 3 1 -1 -1 -1 -1 -1 -1 -1
            1 2 10 9 5 11 9 11 8 11 5 6 -1 -1 -1 -1
            0 11 3 0 6 11 0 9 6 5 6 9 1 2 10 -1
            11 8 5 11 5 6 8 0 5 10 5 2 0 2 5 -1
            6 11 3 6 3 5 2 10 3 10 5 3 -1 -1 -1 -1
            5 8 9 5 2 8 5 6 2 3 8 2 -1 -1 -1 -1
            9 5 6 9 6 0 0 6 2 -1 -1 -1 -1 -1 -1 -1
            1 5 8 1 8 0 5 6 8 3 8 2 6 2 8 -1
            1 5 6 2 1 6 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            1 3 6 1 6 10 3 8 6 5 6 9 8 9 6 -1
            10 1 0 10 0 6 9 5 0 5 6 0 -1 -1 -1 -1
            0 3 8 5 6 10 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            10 5 6 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            11 5 10 7 5 11 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            11 5 10 11 7 5 8 3 0 -1 -1 -1 -1 -1 -1 -1
            5 11 7 5 10 11 1 9 0 -1 -1 -1 -1 -1 -1 -1
            10 7 5 10 11 7 9 8 1 8 3 1 -1 -1 -1 -1
            11 1 2 11 7 1 7 5 1 -1 -1 -1 -1 -1 -1 -1
            0 8 3 1 2 7 1 7 5 7 2 11 -1 -1 -1 -1
            9 7 5 9 2 7 9 0 2 2 11 7 -1 -1 -1 -1
            7 5 2 7 2 11 5 9 2 3 2 8 9 8 2 -1
            2 5 10 2 3 5 3 7 5 -1 -1 -1 -1 -1 -1 -1
            8 2 0 8 5 2 8 7 5 10 2 5 -1 -1 -1 -1
            9 0 1 5 10 3 5 3 7 3 10 2 -1 -1 -1 -1
            9 8 2 9 2 1 8 7 2 10 2 5 7 5 2 -1
            1 3 5 3 7 5 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            0 8 7 0 7 1 1 7 5 -1 -1 -1 -1 -1 -1 -1
            9 0 3 9 3 5 5 3 7 -1 -1 -1 -1 -1 -1 -1
            9 8 7 5 9 7 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            5 8 4 5 10 8 10 11 8 -1 -1 -1 -1 -1 -1 -1
            5 0 4 5 11 0 5 10 11 11 3 0 -1 -1 -1 -1
            0 1 9 8 4 10 8 10 11 10 4 5 -1 -1 -1 -1
            10 11 4 10 4 5 11 3 4 9 4 1 3 1 4 -1
            2 5 1 2 8 5 2 11 8 4 5 8 -1 -1 -1 -1
            0 4 11 0 11 3 4 5 11 2 11 1 5 1 11 -1
            0 2 5 0 5 9 2 11 5 4 5 8 11 8 5 -1
            9 4 5 2 11 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            2 5 10 3 5 2 3 4 5 3 8 4 -1 -1 -1 -1
            5 10 2 5 2 4 4 2 0 -1 -1 -1 -1 -1 -1 -1
            3 10 2 3 5 10 3 8 5 4 5 8 0 1 9 -1
            5 10 2 5 2 4 1 9 2 9 4 2 -1 -1 -1 -1
            8 4 5 8 5 3 3 5 1 -1 -1 -1 -1 -1 -1 -1
            0 4 5 1 0 5 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            8 4 5 8 5 3 9 0 5 0 3 5 -1 -1 -1 -1
            9 4 5 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            4 11 7 4 9 11 9 10 11 -1 -1 -1 -1 -1 -1 -1
            0 8 3 4 9 7 9 11 7 9 10 11 -1 -1 -1 -1
            1 10 11 1 11 4 1 4 0 7 4 11 -1 -1 -1 -1
            3 1 4 3 4 8 1 10 4 7 4 11 10 11 4 -1
            4 11 7 9 11 4 9 2 11 9 1 2 -1 -1 -1 -1
            9 7 4 9 11 7 9 1 11 2 11 1 0 8 3 -1
            11 7 4 11 4 2 2 4 0 -1 -1 -1 -1 -1 -1 -1
            11 7 4 11 4 2 8 3 4 3 2 4 -1 -1 -1 -1
            2 9 10 2 7 9 2 3 7 7 4 9 -1 -1 -1 -1
            9 10 7 9 7 4 10 2 7 8 7 0 2 0 7 -1
            3 7 10 3 10 2 7 4 10 1 10 0 4 0 10 -1
            1 10 2 8 7 4 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            4 9 1 4 1 7 7 1 3 -1 -1 -1 -1 -1 -1 -1
            4 9 1 4 1 7 0 8 1 8 7 1 -1 -1 -1 -1
            4 0 3 7 4 3 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            4 8 7 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            9 10 8 10 11 8 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            3 0 9 3 9 11 11 9 10 -1 -1 -1 -1 -1 -1 -1
            0 1 10 0 10 8 8 10 11 -1 -1 -1 -1 -1 -1 -1
            3 1 10 11 3 10 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            1 2 11 1 11 9 9 11 8 -1 -1 -1 -1 -1 -1 -1
            3 0 9 3 9 11 1 2 9 2 11 9 -1 -1 -1 -1
            0 2 11 8 0 11 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            3 2 11 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            2 3 8 2 8 10 10 8 9 -1 -1 -1 -1 -1 -1 -1
            9 10 2 0 9 2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            2 3 8 2 8 10 0 1 8 1 10 8 -1 -1 -1 -1
            1 10 2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            1 3 8 9 1 8 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            0 9 1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            0 3 8 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
            -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
        ] .+ 1
end

function sfce(Ξ::Vector{Float64})

    ξ₁ = Ξ[1]
    ξ₂ = Ξ[2]
    ξ₃ = Ξ[3]

    N = Array{Float64}(undef, 8)
    N[1] = -1 / 8 * (ξ₁ - 1) * (ξ₂ - 1) * (ξ₃ - 1)
    N[2] = 1 / 8 * (ξ₁ + 1) * (ξ₂ - 1) * (ξ₃ - 1)
    N[3] = -1 / 8 * (ξ₁ + 1) * (ξ₂ + 1) * (ξ₃ - 1)
    N[4] = 1 / 8 * (ξ₁ - 1) * (ξ₂ + 1) * (ξ₃ - 1)

    N[5] = 1 / 8 * (ξ₁ - 1) * (ξ₂ - 1) * (ξ₃ + 1)
    N[6] = -1 / 8 * (ξ₁ + 1) * (ξ₂ - 1) * (ξ₃ + 1)
    N[7] = 1 / 8 * (ξ₁ + 1) * (ξ₂ + 1) * (ξ₃ + 1)
    N[8] = -1 / 8 * (ξ₁ - 1) * (ξ₂ + 1) * (ξ₃ + 1)

    d¹N_dξ¹ = Array{Float64}(undef, 8, 3)
    d¹N_dξ¹[1, :] = [
        -0.125 * (ξ₂ - 1) * (ξ₃ - 1)
        (0.125 - 0.125 * ξ₁) * (ξ₃ - 1)
        (0.125 - 0.125 * ξ₁) * (ξ₂ - 1)
    ]

    d¹N_dξ¹[2, :] = [
        0.125 * (ξ₂ - 1) * (ξ₃ - 1)
        (0.125 * ξ₁ + 0.125) * (ξ₃ - 1)
        (0.125 * ξ₁ + 0.125) * (ξ₂ - 1)
    ]

    d¹N_dξ¹[3, :] = [
        -0.125 * (ξ₂ + 1) * (ξ₃ - 1)
        (-0.125 * ξ₁ - 0.125) * (ξ₃ - 1)
        (-0.125 * ξ₁ - 0.125) * (ξ₂ + 1)
    ]

    d¹N_dξ¹[4, :] = [
        0.125 * (ξ₂ + 1) * (ξ₃ - 1)
        (0.125 * ξ₁ - 0.125) * (ξ₃ - 1)
        (0.125 * ξ₁ - 0.125) * (ξ₂ + 1)
    ]

    d¹N_dξ¹[5, :] = [
        0.125 * (ξ₂ - 1) * (ξ₃ + 1)
        (0.125 * ξ₁ - 0.125) * (ξ₃ + 1)
        (0.125 * ξ₁ - 0.125) * (ξ₂ - 1)
    ]

    d¹N_dξ¹[6, :] = [
        -0.125 * (ξ₂ - 1) * (ξ₃ + 1)
        (-0.125 * ξ₁ - 0.125) * (ξ₃ + 1)
        (-0.125 * ξ₁ - 0.125) * (ξ₂ - 1)
    ]

    d¹N_dξ¹[7, :] = [
        0.125 * (ξ₂ + 1) * (ξ₃ + 1)
        (0.125 * ξ₁ + 0.125) * (ξ₃ + 1)
        (0.125 * ξ₁ + 0.125) * (ξ₂ + 1)
    ]

    d¹N_dξ¹[8, :] = [
        -0.125 * (ξ₂ + 1) * (ξ₃ + 1)
        (0.125 - 0.125 * ξ₁) * (ξ₃ + 1)
        (0.125 - 0.125 * ξ₁) * (ξ₂ + 1)
    ]



    d²N_dξ² = Array{Float64}(undef, 8, 3, 3)
    d²N_dξ²[1, :, :] = [
        0 0.125*(1-ξ₃) 0.125*(1-ξ₂)
        0.125*(1-ξ₃) 0 0.125*(1-ξ₁)
        0.125*(1-ξ₂) 0.125*(1-ξ₁) 0
    ]

    d²N_dξ²[2, :, :] = [
        0 0.125*(ξ₃-1) 0.125*(ξ₂-1)
        0.125*(ξ₃-1) 0 0.125*(ξ₁+1)
        0.125*(ξ₂-1) 0.125*(ξ₁+1) 0
    ]

    d²N_dξ²[3, :, :] = [
        0 0.125*(1-ξ₃) -0.125*(ξ₂+1)
        0.125*(1-ξ₃) 0 -0.125*(ξ₁+1)
        -0.125*(ξ₂+1) -0.125*(ξ₁+1) 0
    ]

    d²N_dξ²[4, :, :] = [
        0 0.125*(ξ₃-1) 0.125*(ξ₂+1)
        0.125*(ξ₃-1) 0 0.125*(ξ₁-1)
        0.125*(ξ₂+1) 0.125*(ξ₁-1) 0
    ]

    d²N_dξ²[5, :, :] = [
        0 0.125*(ξ₃+1) 0.125*(ξ₂-1)
        0.125*(ξ₃+1) 0 0.125*(ξ₁-1)
        0.125*(ξ₂-1) 0.125*(ξ₁-1) 0
    ]

    d²N_dξ²[6, :, :] = [
        0 -0.125*(ξ₃+1) 0.125*(1-ξ₂)
        -0.125*(ξ₃+1) 0 -0.125*(ξ₁+1)
        0.125*(1-ξ₂) -0.125*(ξ₁+1) 0
    ]

    d²N_dξ²[7, :, :] = [
        0 0.125*(ξ₃+1) 0.125*(ξ₂+1)
        0.125*(ξ₃+1) 0 0.125*(ξ₁+1)
        0.125*(ξ₂+1) 0.125*(ξ₁+1) 0
    ]

    d²N_dξ²[8, :, :] = [
        0 -0.125*(ξ₃+1) -0.125*(ξ₂+1)
        -0.125*(ξ₃+1) 0 0.125*(1-ξ₁)
        -0.125*(ξ₂+1) 0.125*(1-ξ₁) 0
    ]



    d³N_dξ³ = Array{Float64}(undef, 8, 3, 3, 3)
    d³N_dξ³[1, :, :, :] = reshape([
        0 0 0
        0 0 -0.125000000000000
        0 -0.125000000000000 0
        0 0 -0.125000000000000
        0 0 0
        -0.125000000000000 0 0
        0 -0.125000000000000 0
        -0.125000000000000 0 0
        0 0 0
    ],(3,3,3))

    d³N_dξ³[2, :, :, :] = reshape([
        0 0 0
        0 0 0.125000000000000
        0 0.125000000000000 0
        0 0 0.125000000000000
        0 0 0
        0.125000000000000 0 0
        0 0.125000000000000 0
        0.125000000000000 0 0
        0 0 0
    ],(3,3,3))

    d³N_dξ³[3, :, :, :] = reshape([
        0 0 0
        0 0 -0.125000000000000
        0 -0.125000000000000 0
        0 0 -0.125000000000000
        0 0 0
        -0.125000000000000 0 0
        0 -0.125000000000000 0
        -0.125000000000000 0 0
        0 0 0
    ],(3,3,3))

    d³N_dξ³[4, :, :, :] = reshape([
        0 0 0
        0 0 0.125000000000000
        0 0.125000000000000 0
        0 0 0.125000000000000
        0 0 0
        0.125000000000000 0 0
        0 0.125000000000000 0
        0.125000000000000 0 0
        0 0 0
    ],(3,3,3))

    d³N_dξ³[5, :, :, :] = reshape([
        0 0 0
        0 0 0.125000000000000
        0 0.125000000000000 0
        0 0 0.125000000000000
        0 0 0
        0.125000000000000 0 0
        0 0.125000000000000 0
        0.125000000000000 0 0
        0 0 0
    ],(3,3,3))

    d³N_dξ³[6, :, :, :] = reshape([
        0 0 0
        0 0 -0.125000000000000
        0 -0.125000000000000 0
        0 0 -0.125000000000000
        0 0 0
        -0.125000000000000 0 0
        0 -0.125000000000000 0
        -0.125000000000000 0 0
        0 0 0
    ],(3,3,3))

    d³N_dξ³[7, :, :, :] = reshape([
        0 0 0
        0 0 0.125000000000000
        0 0.125000000000000 0
        0 0 0.125000000000000
        0 0 0
        0.125000000000000 0 0
        0 0.125000000000000 0
        0.125000000000000 0 0
        0 0 0
    ],(3,3,3))

    d³N_dξ³[8, :, :, :] = reshape([
        0 0 0
        0 0 -0.125000000000000
        0 -0.125000000000000 0
        0 0 -0.125000000000000
        0 0 0
        -0.125000000000000 0 0
        0 -0.125000000000000 0
        -0.125000000000000 0 0
        0 0 0
    ],(3,3,3))

    return N, d¹N_dξ¹, d²N_dξ², d³N_dξ³
end



end # module Rho2sdf
