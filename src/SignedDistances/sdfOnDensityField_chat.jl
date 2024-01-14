# 1. Use meaningful variable names. Avoid using single-letter variables unless it's absolutely necessary as it makes the code harder to understand for others.
# 2. Break down complex operations into smaller functions. This can make the code easier to read and maintain. 
# 3. Implement error handling where appropriate, especially when dealing with external dependencies or libraries.
# 4. Use type annotations liberally in Julia, which will help prevent potential runtime errors due to incorrect types of arguments being passed in.
# 5. Utilize built-in functions whenever possible, as they are usually optimized for performance. 
# 6. Make use of comprehensions and map/reduce operations where appropriate. They can make your code more concise and faster than using loops.
# 7. Finally, consider adding comments to explain complex parts of the code that may not be immediately clear.

# Here's an example of how you might refactor this function:

function compute_distances(mesh::MeshGrid, grid::Grid, ρₙ::Vector{Float64}, ρₜ::Float64)
    # generate the grid points
    points = generate_grid_points(grid)
    
    # create a linked list for efficient lookups
    linked_list = LinkedList(grid, points)
    
    # retrieve necessary properties from the linked list
    head = linked_list.head 
    next = linked_list.next 
    N = linked_list.grid.N 
    AABB_min = linked_list.grid.AABB_min 
    AABB_max = linked_list.grid.AABB_max 
    
    # calculate other necessary parameters
    δ = 5.1 * grid.cell_size 
    X = mesh.X
    IEN = mesh.IEN
    INE = mesh.INE
    ISN = mesh.ISN
    nsd = mesh.nsd
    nel = mesh.nel
    nes = mesh.nes
    nsn = mesh.nsn
    
    # initialize distance field with a large number
    dist = fill(Inf, size(points)) 
    
    # process each element in the mesh
    for el = 1:nel
        ρₑ = ρₙ[IEN[:, el]]
        
        if all(<(ρₜ), ρₑ) 
            continue 
        end
        
        # find elements with common boundaries and calculate their pseudonormals
        common_elements = find_common_elements(INE, IEN, ISN, nes, el)
        
        for sg = 1:nes 
            common_elements[sg] = INE[IEN[ISN[sg][1], el]]
            
            for a = 2:nsn
                idx = findall(in.(INE[IEN[ISN[sg][a], el]]), common_elements[sg])
                common_elements[sg] = common_elements[sg][idx]
            end
        end
        
        # process boundary elements and calculate distances to them
        process_boundary_elements(X, IEN, ISN, mesh, grid, ρₙ, δ, AABB_min, AABB_max, N, nsd, nel, nes, nsn, el, dist)
        
    end
    
    # export the distances to a .vtu file and a .csv file
    Rho2sdf.exportToVTU("xp.vtu", X, IEN)
    open("xp.csv", "w") do io
        writedlm(io, ['x' 'y' 'z'], ',')
        writedlm(io, xp', ',')
    end
    
    return dist
end

#______
function evalSignedDistances(mesh::Mesh, grid::Grid, ρₙ::Vector{Float64}, ρₜ::Float64)
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

    for el in 1:nel
        ρₑ = ρₙ[IEN[:, el]]
        ρₑ_min, ρₑ_max = extrema(ρₑ)
        
        if (ρₑ_min >= ρₜ) # Element does not intersect boundary
            continue
        elseif (ρₑ_max > ρₜ) # Boundary intersects element
            Xₑ = X[:, IEN[:, el]]
            Xₑ_min, Xₑ_max = extrema(Xₑ, dims=2)
            
            I_min = floor.(N .* (Xₑ_min .- AABB_min) ./ (AABB_max .- AABB_min))
            I_max = floor.(N .* (Xₑ_max .- AABB_min) ./ (AABB_max .- AABB_min))
            
            I_min = max.(I_min, 0)
            I_max = min.(I_max, N)
            
            Is = Iterators.product(I_min[1]:I_max[1], I_min[2]:I_max[2], I_min[3]:I_max[3])
            for I in Is
                ii = Int(I[3] * (N[1] + 1) * (N[2] + 1) + I[2] * (N[1] + 1) + I[1] + 1)
                
                v = head[ii]
                while v != -1
                    x = points[:, v]
                    
                    xₚ = zeros(nsd) 
                    Ξ = zeros(nsd)
                    λ = 1.0
                    Ξ_tol = 1e-2
                    r_tol = 1e-2
                    niter = 100
                    
                    for iter in 1:niter
                        (K, r) = AnalyticalDerivations(Ξ, Xₑ, ρₑ, λ, ρₜ, x)
                        
                        r_norm = norm(r)
                        Λ = real.(eigvals(K))
                        (Λ_min, idx_min) = findmin(Λ)
                        
                        if Λ_min < 1.0e-10
                            Φ = real.(eigvecs(K))
                            
                            idx = setdiff([1, 2, 3, 4], idx_min)
                            Φ = Φ[:, idx]
                            ΔΞ̃_and_Δλ̃ = (1.0 ./ Λ[idx]) .* (transpose(Φ) * r)
                            ΔΞ_and_Δλ = Φ * ΔΞ̃_and_Δλ
                        else
                            ΔΞ_and_Δλ = K \ -r
                        end
                        
                        max_abs_Ξ = maximum(abs.(ΔΞ_and_Δλ[1:end-1]))
                        if (max_abs_Ξ > 1.0)
                            ΔΞ_and_Δλ[1:end-1] /= max_abs_Ξ
                        end
                        
                        Ξ += ΔΞ_and_Δλ[1:end-1]
                        λ += ΔΞ_and_Δλ[end]
                        
                        Ξ_norm = norm(ΔΞ_and_Δλ)
                        
                        if (maximum(abs.(Ξ)) <= 1.0 && r_norm < r_tol) # xₚ is in the element
                            dist_tmp = dot(x - xₚ, n)
                            if abs(dist_tmp) < abs(dist[v])
                                dist[v] = dist_tmp
                                xp[:, v] = xₚ
                            end
                        end
                    end
                    
                    v = next[v]
                end
            end
        else # Element is inside boundary
            commonEls = []
            
            for sg in 1:nes
                commonEls = INE[IEN[mesh.ISN[sg][1], el]]
                
                for a in 2:nsn
                    idx = findall(in(INE[IEN[ISN[sg][a], el]]), commonEls)
                    commonEls = commonEls[idx]
                end
                
                if length(commonEls) == 1 # is a part of the outer boundary of the body
                    Xs = X[:, IEN[ISN[sg], el]]
                    Xc = mean(Xs, dims=2)
                    
                    for a in 1:nsn
                        As = IEN[ISN[sg][a], el]
                        adj_els = INE[As]
                        
                        for adj_el in adj_els
                            common_adj_els = []
                            
                            for adj_sg in 1:nes
                                common_adj_els = INE[IEN[mesh.ISN[adj_sg][1], adj_el]]
                                
                                for b in 2:nsn
                                    idx = findall(in(INE[IEN[ISN[adj_sg][b], adj_el]]), common_adj_els)
                                    common_adj_els = common_adj_els[idx]
                                end
                                
                                if length(common_adj_els) == 1 && in(As, IEN[mesh.ISN[adj_sg], adj_el])
                                    adj_Xs = X[:, IEN[ISN[adj_sg], adj_el]]
                                    adj_Xc = mean(adj_Xs, dims=2)
                                    
                                    as = findall(in(As), IEN[mesh.ISN[adj_sg], adj_el])
                                    a_prev = ((as[1] + nsn - 1 - 1) % nsn) + 1
                                    a_next = ((as[1] + nsn + 1 - 1) % nsn) + 1
                                    
                                    x_prev = X[:, IEN[mesh.ISN[adj_sg][a_prev], adj_el]]
                                    xs = X[:, IEN[ISN[adj_sg][as], adj_el]]
                                    x_next = X[:, IEN[ISN[adj_sg][a_next], adj_el]]
                                    
                                    Xt = [x_prev, xs, x_next]
                                    Xt = reduce(hcat, Xt)
                                    
                                    Et = [[Xt[:, 2] - Xt[:, 1], Xt[:, 3] - Xt[:, 2], Xt[:, 1] - Xt[:, 3]]...]
                                    
                                    n = cross(Et[1], Et[2])
                                    n /= norm(n)
                                end
                            end
                        end
                        
                        x₁, x₂, x₃ = Xs[:, a], Xs[:, (a % nsn) + 1], Xc
                        Xt = [x₁, x₂, x₃]
                        Xt_min, Xt_max = minimum(Xt, dims=2), maximum(Xt, dims=2)
                        
                        Xt_min .-= δ
                        Xt_max .+= δ
                        
                        I_min = floor.(N .* (Xt_min .- AABB_min) ./ (AABB_max .- AABB_min))
                        I_max = floor.(N .* (Xt_max .- AABB_min) ./ (AABB_max .- AABB_min))
                        
                        I_min = max.(I_min, 0)
                        I_max = min.(I_max, N)
                        
                        Is = Iterators.product(I_min[1]:I_max[1], I_min[2]:I_max[2], I_min[3]:I_max[3])
                        
                        for I in Is
                            ii = Int(I[3] * (N[1] + 1) * (N[2] + 1) + I[2] * (N[1] + 1) + I[1] + 1)
                            
                            v = head[ii]
                            while v != -1
                                x = points[:, v]
                                
                                A = [(x₁[2]*n[3] - x₁[3]*n[2]) (x₂[2]*n[3] - x₂[3]*n[2]) (x₃[2]*n[3] - x₃[3]*n[2]);
                                      (x₁[3]*n[1] - x₁[1]*n[3]) (x₂[3]*n[1] - x₂[1]*n[3]) (x₃[3]*n[1] - x₃[1]*n[3]);
                                      (x₁[1]*n[2] - x₁[2]*n[1]) (x₂[1]*n[2] - x₂[2]*n[1]) (x₃[1]*n[2] - x₃[2]*n[1])]
                                b = [x[2] * n[3] - x[3] * n[2], x[3] * n[1] - x[1] * n[3], x[1] * n[2] - x[2] * n[1]]
                                
                                A[findmax(abs.(n))[2], :] = [1.0, 1.0, 1.0]
                                b[findmax(abs.(n))[2]] = 1.0
                                λ = A \ b
                                
                                xₚ = zeros(nsd)
                                isFace, isEdge, isVertex = false, false, false
                                
                                if minimum(λ) >= 0.0 # xₚ is in the triangle el
                                    xₚ = λ[1] * x₁ + λ[2] * x₂ + λ[3] * x₃
                                    
                                    dist_tmp = dot(x - xₚ, n)
                                    if abs(dist_tmp) < abs(dist[v]) # Checking distance (looking for the shortest one)
                                        dist[v] = dist_tmp
                                        isFace = true
                                        xp[:, v] = xₚ
                                    end
                                else
                                    for j in 1:3
                                        L = norm(Et[j])
                                        xᵥ = Xt[:, j]
                                        P = dot(x - xᵥ, Et[j] / L)
                                        
                                        if P >= 0 && P <= L
                                            xₚ = xᵥ + (Et[j] / L) * P
                                            
                                            dist_tmp = sign(dot(x - xₚ, n)) * norm(x - xₚ)
                                            if abs(dist_tmp) < abs(dist[v])
                                                dist[v] = dist_tmp
                                                isVertex = true
                                                xp[:, v] = xₚ
                                            end
                                        end
                                    end
                                end
                                
                                if !isFace && !isEdge
                                    dist_tmp, idx = findmin([norm(x - x₁), norm(x - x₂), norm(x - x₃)])
                                    xₚ = Xt[:, idx]
                                    
                                    dist_tmp *= sign(dot(x - xₚ, n))
                                    if abs(dist_tmp) < abs(dist[v])
                                        dist[v] = dist_tmp
                                        isVertex = true
                                        xp[:, v] = xₚ
                                    end
                                end
                                
                                v = next[v]
                            end
                        end
                    end
                end
            end
        end
    end
    
    X = Vector{Float64}(undef, 3)
    Xp = Vector{Float64}(undef, 3)
    
    for i in 1:size(xp, 2)
        if sum(abs.(xp[:, i])) > 1.0e-10
            X = [X; points[:, i]]
            Xp = [Xp; xp[:, i]]
        end
    end
    
    X = [X; Xp]
    X = [X[i] for i in 1:size(X, 2)]
    
    nnp = size(Xp, 2)
    IEN = [[i, i + nnp] for i in 1:nnp]
    
    Rho2sdf.exportToVTU("xp.vtu", X, IEN)
    
    open("xp.csv", "w") do io
        writedlm(io, ['x', 'y', 'z'], ',')
        writedlm(io, xp', ',')
    end
    
    return dist
end
