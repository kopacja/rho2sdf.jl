

function evalSignedDistances(
    mesh::Mesh,
    grid::Grid,
    ρₙ::Vector{Float64},
    ρₜ::Float64,
)

    points = MeshGrid.generateGridPoints(grid) # uzly pravidelné mřížky
    linkedList = MeshGrid.LinkedList(grid, points) # pro rychlé vyhledávání

    head = linkedList.head # ID pravidelné bunky (pozice), index bodu z points
    next = linkedList.next # vel délky points, další uzly pro danou bunku, když -1 tak už další není
    N = linkedList.grid.N # počet buněk ve směru xyz
    AABB_min = linkedList.grid.AABB_min # bod boxu (minimum v xyz)
    AABB_max = linkedList.grid.AABB_max # bod boxu (maximum v xyz)
    δ = 5.1 * grid.cell_size # offset (přifouknutí celého boxu)
    X = mesh.X
    IEN = mesh.IEN
    INE = mesh.INE
    ISN = mesh.ISN
    nsd = mesh.nsd
    nel = mesh.nel
    nes = mesh.nes
    nsn = mesh.nsn

    ngp = grid.ngp # počet vrcholu prav sítě
    big = 1.0e10
    dist = big * ones(ngp) # inicializace dist fieldu
    xp = zeros(nsd, ngp) # souřadnice bodů vrcholů (3xngp)

    for el = 1:nel
    # for el = 415
        ρₑ = ρₙ[IEN[:, el]]

        ρₑ_min = minimum(ρₑ)
        ρₑ_max = maximum(ρₑ)
        if (ρₑ_min >= ρₜ) # hranice elementem neprochází
            # continue # PRO PŘESKAKUJE HRANIČNÍ ELEMENTY (pouze pro ladění kodu)
            commonEls = []
            for sg = 1:nes # je to hraniční element?
                commonEls = INE[IEN[mesh.ISN[sg][1], el]]
                for a = 2:nsn
                    idx = findall(in(INE[IEN[ISN[sg][a], el]]), commonEls)
                    commonEls = commonEls[idx]
                end

                # Mám vnější segment a chci jeho pseudonormály,...

                if (length(commonEls) == 1) # is a part of the outer boundary of the body
                    Xs = X[:, IEN[ISN[sg], el]]
                    Xc = mean(Xs, dims = 2)

                    for a = 1:nsn # cyklus přes uzly segmentu (face)

                        As = IEN[ISN[sg][a], el] # globální číslo uzlu
                        adj_els = INE[As] # všechny elementy které jsou součástí tohoto uzlu
                        for adj_el in adj_els
                            common_adj_els = []
                            for adj_sg = 1:nes # cyklus přes sousední stěny elementu
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

                                    a_prev = ((as[1] + nsn - 1 - 1) % nsn) + 1
                                    a_next = ((as[1] + nsn + 1 - 1) % nsn) + 1

                                    x_prev = X[
                                        :,
                                        IEN[mesh.ISN[adj_sg][a_prev], adj_el],
                                    ]
                                    xs = X[:, IEN[mesh.ISN[adj_sg][as], adj_el]]
                                    x_next = X[
                                        :,
                                        IEN[mesh.ISN[adj_sg][a_next], adj_el],
                                    ]

                                    Xt = [x_prev, xs, x_next]
                                    Xt = reduce(hcat, Xt)

                                    Et = calculate_triangle_edges(Xt)

                                    n = cross(Et[1], Et[2])
                                    n = n / norm(n)

                                end
                            end
                        end
                        
                        x₁ = Xs[:, a] # trojúhelník
                        x₂ = Xs[:, (a%nsn)+1]
                        x₃ = Xc

                        Xt = [x₁, x₂, x₃]
                        Xt = reduce(hcat, Xt)

                        # Hrany trojúhelníku
                        Et = calculate_triangle_edges(Xt)

                        # normála trojúhelníka
                        n = cross(Et[1], Et[2])
                        n = n / norm(n)

                        Is = MeshGrid.calculateMiniAABB_grid(Xt, δ, N, AABB_min, AABB_max, nsd)

                        for I ∈ Is # bunce přiřadím jedno číslo
                            ii = Int(
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
                    # Subdivide quad into two (or four?) triangles
                end
            end
        else # (ρₑ_min < ρₜ)
            # continue
            if (ρₑ_max > ρₜ)
                # println("Hranice prochází elementem...")

                Xₑ = X[:, IEN[:, el]]
                
                Is = MeshGrid.calculateMiniAABB_grid(Xₑ, δ, N, AABB_min, AABB_max, nsd)

                for I ∈ Is
                    ii = Int(
                        I[3] * (N[1] + 1) * (N[2] + 1) + I[2] * (N[1] + 1) + I[1] + 1,
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
                            sleep(0.5)
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


