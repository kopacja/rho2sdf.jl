

function ReduceEigenvals(K::Matrix{Float64}, r::Vector{Float64}, Sign::Int)
    Λ = real.(eigvals(K))
    (Λ_min, idx_min) = findmin(abs.(Λ))
                           
    if (abs(Λ_min) < 1.0e-6)
        Φ = real.(eigvecs(K))
        idx = [1, 2, 3, 4, 5]
        resize!(idx, length(Λ))
        deleteat!(idx, idx_min)
        Φ = Φ[:, idx]
        ΔΞ̃_and_Δλ̃ = 1.0 ./ Λ[idx] .* (Φ' * r)
        ΔΞ_and_Δλ = (Φ .* Sign) * ΔΞ̃_and_Δλ̃
    else
        ΔΞ_and_Δλ = K \ (r .* Sign)
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

    # Trim the unused preallocated space
    X = resize!(X, count)
    Xp = resize!(Xp, count)

    # Mean and max projected distance:
    mean_PD = mean(norm.(X-Xp))
    max_PD = maximum(norm.(X-Xp))

    return X, Xp, mean_PD, max_PD
end

function SignCorrection4SDF(dist::Vector{Float64},
    grid::Grid, 
    big::Float64)
    
    ngp = grid.ngp # number of nodes in grid
    Sign = -1
    for i in 1:ngp
        if dist[i] != big
            Sign = sign(dist[i])
        end
        if dist[i] == big && Sign == 1
           dist[i] = dist[i] * -1
        end
    end
    return dist
end


function evalSignedDistances(
    mesh::Mesh,
    grid::Grid,
    ρₙ::Vector{Float64},
    ρₜ::Float64,
)

    points = MeshGrid.generateGridPoints(grid) # uzly pravidelné mřížky
    # points = reshape([-0.5 -0.5 -0.5], 3,1)
    linkedList = MeshGrid.LinkedList(grid, points) # pro rychlé vyhledávání

    head = linkedList.head # ID pravidelné bunky (pozice), index bodu z points
    next = linkedList.next # vec délky points, další uzly pro danou bunku, když -1 tak už další není
    N = linkedList.grid.N # Number of divisions along each axis of the grid
    AABB_min = linkedList.grid.AABB_min # Minimum coordinates of the Axis-Aligned Bounding Box (AABB)
    AABB_max = linkedList.grid.AABB_max # Maximum coordinates of the AABB
    δ = 1.2 * grid.cell_size # offset for mini AABB
    
    X   = mesh.X   # vector of nodes positions
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
            if (ρₑ_max > ρₜ) # Hranice prochází elementem...

                Xₑ = X[:, IEN[:, el]]

                Is = MeshGrid.calculateMiniAABB_grid(Xₑ, δ, N, AABB_min, AABB_max, nsd)

                for I ∈ Is

                    ii = Int(
                        I[3] * (N[1] + 1) * (N[2] + 1) + I[2] * (N[1] + 1) + I[1] + 1,
                    )

                    v = head[ii]
                    while v != -1
                        x = points[:, v]

                        Ξ = zeros(Float64, 3)  # local coordinates
                        λ = 1.0              # Lagrange multiplier
                        Ξ_tol = 1e-2
                        Ξ_norm = 2 * Ξ_tol
                        Ξ_norm_old = 1000.0
                        r_tol = 1e-2
                        r_norm = 2 * r_tol   # 
                        niter = 10           # maximum number of iterations
                        iter = 1             # iteration form one 
 
                        while ((Ξ_norm ≥ Ξ_tol || r_norm ≥ r_tol) && iter ≤ niter)

                            K = Hessian(sfce, Ξ, λ, x, Xₑ, ρₑ)
                            r = Gradient(sfce, Ξ, λ, x, Xₑ, ρₑ, ρₜ)

                            r_norm = norm(r)

                            ΔΞ_and_Δλ = K \ -r
                            # (ΔΞ_and_Δλ, Λ_min) = ReduceEigenvals(K, r, -1)

                            Ξ += ΔΞ_and_Δλ[1:3]
                            λ += ΔΞ_and_Δλ[4]

                            Ξ_norm = norm(ΔΞ_and_Δλ)

                            iter = iter + 1
                        end
#=
                        # If projection is not inside the element it is a good idea to try
                        # to project on the edges and corners of the isosurface              
                        if (maximum(abs.(Ξ)) > 1.0) # xₚ is NOT in the element

                            # Let's loop  check whether there is a projection on the edges of the density isocontour.

                            # Each segment (face) have not three components ξ₁, ξ₂, ξ₃ but only two and the
                            # third one is known constant and must be fixed by constraint equation. 
                            # Following two vectors represents index (1, 2 or 3) of the fixed component and its
                            # value (-1 or 1):
                            idx = [3, 2, 1, 2, 1, 3] # Index of the third constant component 
                            Ξ_ = [-1, -1, 1, 1, -1, 1] #

                            # Loop over segments (nes=6 for hex element)
                            for sg = 1:nes
                                ρₛ = ρₑ[mesh.ISN[sg]]

                                ρₛ_min = minimum(ρₛ)
                                ρₛ_max = maximum(ρₛ)

                                if (ρₛ_min <= ρₜ && ρₛ_max >= ρₜ) # the boundary cross through the segment

                                    Ξ = zeros(3)
                                    λ = ones(2)
                                    Ξ_norm = 2 * Ξ_tol
                                    r_norm = 2 * r_tol
                                    iter = 1
                                    niter = 10
                                    while ((Ξ_norm ≥ Ξ_tol || r_norm ≥ r_tol) && iter ≤ niter)

                                        r4 = Gradient(sfce, Ξ, λ[1], x, Xₑ, ρₑ, ρₜ)
                                        K4 = Hessian(sfce, Ξ, λ[1], x, Xₑ, ρₑ)

                                        # Fifth equation, e.g. (ξ₁ - 1) = 0 etc., representing constraint of the fixed segment component is added into the residual and tangent matrix
                                        r = zeros(5)
                                        r[1:4] = r4
                                        r[5] = (Ξ[idx[sg]] - Ξ_[sg])
                                        r[idx[sg]] += λ[2]

                                        K = zeros(5, 5)
                                        K[1:4, 1:4] = K4
                                        K[5, idx[sg]] = 1.0
                                        K[idx[sg], 5] = 1.0

                                        r_norm = norm(r)

                                        ΔΞ_and_Δλ = K \ -r

                                        ΔΞ = ΔΞ_and_Δλ[1:3]
                                        Δλ = ΔΞ_and_Δλ[4:5]
                                        
                                        if (maximum(abs.(ΔΞ)) > 1.0)
                                            ΔΞ *= 0.2/maximum(abs.(ΔΞ))
                                        end

                                        Ξ += ΔΞ
                                        λ += Δλ

                                        Ξ_norm = norm(ΔΞ_and_Δλ)
                                        #println("SG: ", sg , ", iter: ", iter, ", Ξ_norm: ", Ξ_norm, ", r_norm: ", r_norm)
                                        iter = iter + 1
                                    end
                                end

                                if (maximum(abs.(Ξ)) <= 1.0) # xₚ is in the segment
                                    H, d¹N_dξ¹, d²N_dξ², d³N_dξ³ = sfce(Ξ)
                                    xₚ = Xₑ * H

                                    (dρ_dΞ, d²ρ_dΞ², d³ρ_dΞ³) = ρ_derivatives(ρₑ, Ξ)
                                    norm_dρ_dΞ = norm(dρ_dΞ)
                                    n = dρ_dΞ / norm_dρ_dΞ

                                    dist_tmp = sign(dot(x - xₚ, n)) * norm(x - xₚ)
                                    if (abs(dist_tmp) < abs(dist[v]))
                                        dist[v] = dist_tmp
                                        xp[:, v] = xₚ
                                    end
                                end
                            end # for sg

                        else # maximum(abs.(Ξ)) <= 1.0) # xₚ is in the element

                            H, d¹N_dξ¹, d²N_dξ², d³N_dξ³ = sfce(Ξ)
                            xₚ = Xₑ * H

                            (dρ_dΞ, d²ρ_dΞ², d³ρ_dΞ³) = ρ_derivatives(ρₑ, Ξ)
                            norm_dρ_dΞ = norm(dρ_dΞ)
                            n = dρ_dΞ / norm_dρ_dΞ

                            dist_tmp = dot(x - xₚ, n)
                            if (abs(dist_tmp) < abs(dist[v]))
                                dist[v] = dist_tmp
                                xp[:, v] = xₚ
                            end
                        end
=#

                        # The closed point could be a corner of the isocontour inside the element. 

                        # Let's loop over edges and check whether rho of the end points is below and above the threshold density
                        for a = 1:length(ρₑ)-1
                            ρ_min, min_idx = findmin([ρₑ[a], ρₑ[a+1]])
                            ρ_max, max_idx = findmax([ρₑ[a], ρₑ[a+1]])

                            a_min = a + min_idx - 1
                            a_max = a + max_idx - 1

                            if (ρ_min <= ρₜ && ρ_max >= ρₜ)

                                ratio = (ρₜ - ρ_min) / (ρ_max - ρ_min)
                                xₚ = Xₑ[:, a_min] + ratio .* (Xₑ[:, a_max] - Xₑ[:, a_min])

                                # Use normal vector at the center of the element
                                (dρ_dΞ, d²ρ_dΞ², d³ρ_dΞ³) = ρ_derivatives(ρₑ, [0.0, 0.0, 0.0])
                                norm_dρ_dΞ = norm(dρ_dΞ)
                                n = dρ_dΞ / norm_dρ_dΞ

                                dist_tmp = sign(dot(x - xₚ, n)) * norm(x - xₚ)
                                if (abs(dist_tmp) < abs(dist[v]))
                                    dist[v] = dist_tmp
                                    xp[:, v] = xₚ
                                end
                            end
                        end

                        v = next[v]

                    end
                end
            end
        end
    end

    # dist = marchingCubes(dist, N.+1, big)

    # Xg, Xp, mean_PD, max_PD = SelectProjectedNodes(mesh, grid, xp, points)
    # println("mean of projected distance: ", mean_PD)
    # println("maximum projected distance: ", max_PD)


    # nnp = Int(length(Xg)/2)
    # IEN = [[i; i + nnp] for i = 1:nnp]
    # 
    # nnp₂ = Int(length(Xg)/2)
    # IEN₂ = [[i; i + nnp₂] for i = 1:nnp₂]

    # X_combined = [Xg; Xp] 
    # # X_combined_couples = [X Xp]

    # nnp = length(Xg)
    # IEN = [[i; i + nnp] for i = 1:nnp]

   # Rho2sdf.exportToVTU("xp.vtu", X, IEN, 5)

    # Rho2sdf.exportToVTU("xp.vtu", X_combined, IEN)
    # Rho2sdf.exportToVTU("Xg.vtu", Xg, IEN₂)
    # Rho2sdf.exportToVTU("Xp.vtu", Xp, IEN₂)

    # open("xp.csv", "w") do io
    #     writedlm(io, ['x' 'y' 'z'], ',')
    #     writedlm(io, xp', ',')
    # end
    dist = SignCorrection4SDF(dist, grid, big)

    return dist, xp
end

