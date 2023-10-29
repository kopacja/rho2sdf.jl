module SignedDistances

export evalSignedDiscancesOnTriangularMesh, evalSignedDiscances

using Einsum
using Statistics
using LinearAlgebra
using DelimitedFiles
using Rho2sdf.ShapeFunctions
using Rho2sdf.MeshGrid
using Rho2sdf

include("Derivatives.jl")


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
    # for el = 415
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
                        println("xxx:",typeof(x))

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
                            println("el:",el)
                            ########################################
                            (K, r) = AnalyticalDerivations(Ξ, Xₑ, ρₑ, λ, ρₜ, x)
                            # (K_diff) = NumericalDerivations(Ξ, Xₑ, ρₑ, λ, ρₜ, x)

                            # K = K_diff
                            # if (round.(K, digits=4))=! (round.(K_diff, digits=4))
                                # println("K:",K)
                                # println("K_diff:",K_diff)
                                # sleep(50)
                            # end
                            ########################################

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

    Rho2sdf.exportToVTU("xp.vtu", X, IEN)

    open("xp.csv", "w") do io
        writedlm(io, ['x' 'y' 'z'], ',')
        writedlm(io, xp', ',')
    end

    return dist
end


end
