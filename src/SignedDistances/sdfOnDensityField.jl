

function ReduceEigenvals(K::Matrix{Float64}, r::Vector{Float64}, sign::Int)
    Λ = real.(eigvals(K))
    (Λ_min, idx_min) = findmin(abs.(Λ))
                           
    if (abs(Λ_min) < 1.0e-6)
        Φ = real.(eigvecs(K))
        idx = [1, 2, 3, 4]
        deleteat!(idx, idx_min)
        Φ = Φ[:, idx]
        ΔΞ̃_and_Δλ̃ = 1.0 ./ Λ[idx] .* (Φ' * r)
        ΔΞ_and_Δλ = (Φ .* sign) * ΔΞ̃_and_Δλ̃
    else
        ΔΞ_and_Δλ = K \ (r .* sign)
    end
    return ΔΞ_and_Δλ, Λ_min
end

function ReturnLocalCoordsIntoTheElement(Ξ::Vector{Float64})
    Ξ_OutOfElement = 0
    Ξₘₐₓcomp = maximum(abs.(Ξ)) 

    if Ξₘₐₓcomp > 1
        Ξ = Ξ ./ Ξₘₐₓcomp
        Ξ_OutOfElement = Ξ_OutOfElement + 1
        if Ξ_OutOfElement > 5
            return Ξ, true  # Return a tuple with a flag indicating to break
        end
    end
    return Ξ, false
end

function SelectProjectedNodes(mesh::Mesh,
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

    # Trim the unused preallocated space
    X = resize!(X, count)
    Xp = resize!(Xp, count)

    # Mean and max projected distance:
    mean_PD = mean(norm.(X-Xp))
    max_PD = maximum(norm.(X-Xp))

    return X, Xp, mean_PD, max_PD
end


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
    δ = 5.1 * grid.cell_size # offset for mini AABB

    X   = mesh.X   # vector of nodes positions
    IEN = mesh.IEN # ID element -> ID nodes
    INE = mesh.INE # ID node -> ID elements
    ISN = mesh.ISN # connectivity face - edges
    nsd = mesh.nsd # number of spacial dimensions
    nel = mesh.nel # number of all elements
    nes = mesh.nes # number of element segments (faces)
    nsn = mesh.nsn # number of face nodes
    println("number of all elements: ", nel)

    ngp = grid.ngp # number of nodes in grid
    big = 1.0e10
    dist = big * ones(ngp) # distance field initialization
    xp = zeros(nsd, ngp) # souřadnice bodů vrcholů (3xngp)

    for el = 1:nel
        println("element ID: ", el)
        ρₑ = ρₙ[IEN[:, el]] # nodal densities for one element

        ρₑ_min = minimum(ρₑ)
        ρₑ_max = maximum(ρₑ)
        if (ρₑ_min >= ρₜ) # the boundary does not cross through the element
            # continue # PRO PŘESKAKUJE HRANIČNÍ ELEMENTY (pouze pro ladění kodu)
            commonEls = []

            # cycle through element faces (6)
            for sg = 1:nes 
                commonEls = INE[IEN[mesh.ISN[sg][1], el]] # 
                for a = 2:nsn
                    idx = findall(in(INE[IEN[ISN[sg][a], el]]), commonEls) # for how many elements does this face belong ?
                    commonEls = commonEls[idx]
                end

                if (length(commonEls) == 1) # = is a part of the outer boundary of the body
                    Xs = X[:, IEN[ISN[sg], el]]
                    Xc = vec(mean(Xs, dims = 2))

                    for a = 1:nsn # cycle through number of all nodals belong to face

                        # coordinates of nodes of the triangle
                        x₁ = Xs[:, a] 
                        x₂ = Xs[:, (a%nsn)+1]
                        x₃ = Xc
                        
                        # coordinates of the vertices of the triangle
                        Xt = [x₁, x₂, x₃]
                        Xt = reduce(hcat, Xt)

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
                                            # n_edge = EPN[el][j]
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
                                    n_vertex = n
                                    # n_vertex = VPN[IEN[idx, el]]
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
                # println("ρₑ_min: ", ρₑ_min)
                # println("ρₑ_max: ", ρₑ_max)

                Xₑ = X[:, IEN[:, el]]
                # println("Xₑ: ", Xₑ)
                
                # println("el: ", el)
                
                Is = MeshGrid.calculateMiniAABB_grid(Xₑ, δ, N, AABB_min, AABB_max, nsd)
                # println("length Is: ", length(Is))
                Isi = 0

                for I ∈ Is
                    Isi = Isi + 1
                    # println("where am I (Isi): ", Isi)

                    ii = Int(
                        I[3] * (N[1] + 1) * (N[2] + 1) + I[2] * (N[1] + 1) + I[1] + 1,
                    )
                    v = head[ii]
                    while v != -1
                        x = points[:, v]
                        # println("xxx:",typeof(x))

                        xₚ = [0.0, 0.0, 0.0] # estimation of point projection position
                        n = [0.0, 0.0, 0.0]  # estimation of norm
                        Ξ = [0.0, 0.0, 0.0]  # local coordinates
                        λ = 1.0              # Lagrange multiplier
                        Ξ_tol = 1e-2
                        Ξ_norm = 2 * Ξ_tol
                        Ξ_norm_old = 1000.0
                        r_tol = 1e-2
                        r_norm = 2 * r_tol   # 
                        niter = 15          # maximum number of iterations
                        iter = 1             # iteration form one 

                        while (Ξ_norm ≥ Ξ_tol && iter ≤ niter)# || r_norm ≥ r_tol)
                            # println("el: ",el)
                            # println("Ξ_norm: ",Ξ_norm)
                            # sleep(0.2)
                            ########################################
                            (K, r, n) = AnalyticalDerivations(Ξ, Xₑ, ρₑ, λ, ρₜ, x)
                            # (K_diff, r_tmp) = NumericalDerivations(Ξ, Xₑ, ρₑ, λ, ρₜ, x)

                            # K = K_diff
                            # r = r_tmp
                            # if (round.(K, digits=4))=! (round.(K_diff, digits=4))
                                # println("K:",K)
                                # println("K_diff:",K_diff)
                                # sleep(50)
                            # end
                            ########################################

                            ΔΞ_and_Δλ, Λ_min = ReduceEigenvals(K, r, 1) # sign (1/ -1) -> ΔΞ_and_Δλ = K \ (r .* sign)

                            max_abs_Ξ = maximum(abs.(ΔΞ_and_Δλ[1:end-1]))
                            if (max_abs_Ξ > 1.0)
                                ΔΞ_and_Δλ[1:end-1] =
                                    ΔΞ_and_Δλ[1:end-1] / max_abs_Ξ
                            end

                            # Coordinates update:
                            Ξ = Ξ - ΔΞ_and_Δλ[1:end-1]
                            λ = λ - ΔΞ_and_Δλ[end]

                            Ξ_norm = norm(ΔΞ_and_Δλ)


                            Ξ, should_break = ReturnLocalCoordsIntoTheElement(Ξ)
                            should_break && break


                            iter = iter + 1
                            
                            ####################################
                            # if Ξ_norm_old < (Ξ_norm * 0.6)
                            #     println("____Diverguje!____")
                            #     println("el: ",el)
                            #     println("length Is: ", length(Is))
                            #     println("where am I (Isi): ", Isi)
                            #     println("local coord Ξ: ",Ξ)
                            #     println("Ξ_norm: ",Ξ_norm)
                            #     println("Λ_min: ",Λ_min)
                            #     H, d¹N_dξ¹, d²N_dξ², d³N_dξ³ = C3D8_SFaD(Ξ) # tvarové funkce a jejich derivace
                            #     xₚ = Xₑ * H
                            #     println("xₚ: ",xₚ)
                            #     println("x: ",x)
                            #     dist_tmp = dot(x - xₚ, n)
                            #     println("dist: ", dist_tmp)
                            #     sleep(2.)
                            # end 
                            # Ξ_norm_old = Ξ_norm
                            ####################################

                        end
                        
                        if (maximum(abs.(Ξ)) <= 1.0) # xₚ is in the element
                            H, d¹N_dξ¹, d²N_dξ², d³N_dξ³ = C3D8_SFaD(Ξ) # tvarové funkce a jejich derivace
                            xₚ = Xₑ * H
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

    Xg, Xp, mean_PD, max_PD = SelectProjectedNodes(mesh, grid, xp, points)
    println("mean of projected distance: ", mean_PD)
    println("maximum projected distance: ", max_PD)

    nnp = Int(length(Xg)/2)
    IEN = [[i; i + nnp] for i = 1:nnp]
    
    # nnp₂ = Int(length(Xg)/2)
    # IEN₂ = [[i; i + nnp₂] for i = 1:nnp₂]

    X_combined = [Xg; Xp] 
    # X_combined_couples = [X Xp]

    Rho2sdf.exportToVTU("xp.vtu", X_combined, IEN)
    # Rho2sdf.exportToVTU("Xg.vtu", Xg, IEN₂)
    # Rho2sdf.exportToVTU("Xp.vtu", Xp, IEN₂)

    open("xp.csv", "w") do io
        writedlm(io, ['x' 'y' 'z'], ',')
        writedlm(io, xp', ',')
    end

    return dist
end

