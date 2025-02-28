function interpolate_zero(p1::SVector{3,Float64}, p2::SVector{3,Float64},
                          f1::Float64, f2::Float64, mesh::BlockMesh; tol=mesh.grid_tol, max_iter=20)::Int64
                          
    # Ensure p1 corresponds to positive SDF and p2 to negative SDF
    if f1 < 0 || f2 >= 0
        p1, p2 = p2, p1
        f1, f2 = f2, f1
    end
    
    # Simple linear interpolation as first approximation
    t = f1 / (f1 - f2)
    t = clamp(t, 0.0, 1.0)  # Ensure t is in [0,1]
    mid = p1 * (1.0 - t) + p2 * t
    
    # Newton-Raphson refinement for better accuracy
    for iter in 1:max_iter
        f_mid = eval_sdf(mesh, mid)
        if abs(f_mid) < tol
            break
        end
        
        # Compute numerical gradient for Newton's method
        h = max(tol, 1e-6 * norm(p2 - p1))
        grad = SVector{3, Float64}(
            (eval_sdf(mesh, mid + SVector{3, Float64}(h, 0, 0)) - f_mid) / h,
            (eval_sdf(mesh, mid + SVector{3, Float64}(0, h, 0)) - f_mid) / h,
            (eval_sdf(mesh, mid + SVector{3, Float64}(0, 0, h)) - f_mid) / h
        )
        
        # Avoid division by very small numbers
        if norm(grad) < tol
            break
        end
        
        # Newton step
        step = f_mid * grad / dot(grad, grad)
        mid = mid - step
        
        # Ensure we stay on the line segment
        proj_t = dot(mid - p1, p2 - p1) / dot(p2 - p1, p2 - p1)
        if proj_t < 0.0 || proj_t > 1.0
            mid = p1 * (1.0 - clamp(proj_t, 0.0, 1.0)) + p2 * clamp(proj_t, 0.0, 1.0)
        end
    end

    # Quantize the point
    p_key = quantize(mid, tol)
    
    if haskey(mesh.node_hash, p_key)
        return mesh.node_hash[p_key]
    else
        push!(mesh.X, mid)
        push!(mesh.node_sdf, 0.0)  # Set exactly to 0
        new_index = length(mesh.X)
        mesh.node_hash[p_key] = new_index
        return new_index
    end
end

function edge_color(mesh::BlockMesh, v1::Int64, v2::Int64)::Symbol
    # Determine whether an edge is inside, outside, or on the boundary
    
    # Get SDF values of endpoints
    f1 = mesh.node_sdf[v1]
    f2 = mesh.node_sdf[v2]
    
    tol = mesh.grid_tol
    
    # Both endpoints have positive SDF (inside the body)
    if f1 > tol && f2 > tol
        return :positive
    end
    
    # Both endpoints have negative SDF (outside the body)
    if f1 < -tol && f2 < -tol
        return :negative
    end
    
    # One endpoint inside, one outside
    if (f1 > tol && f2 < -tol) || (f1 < -tol && f2 > tol)
        return :mixed
    end
    
    # At least one endpoint on boundary
    return :boundary
end

function quantize(p::SVector{3,Float64}, tol::Float64)::NTuple{3,Int64}
    # Quantize coordinates to a discrete grid for hashing purposes
    return (
        Int64(round(p[1] / tol)),
        Int64(round(p[2] / tol)),
        Int64(round(p[3] / tol))
    )
end

function stencil_1p3(mesh::BlockMesh, tet::Vector{Int64})::Vector{Vector{Int64}}
    # Get SDF values for all vertices
    f = [mesh.node_sdf[i] for i in tet]
    
    # Identify vertices by their position relative to the isosurface
    plus = map(x -> x > 0, f)      # vertices inside the body
    minus = map(x -> x < 0, f)     # vertices outside the body
    
    # Find indices of positive and negative vertices
    pos_idx = findfirst(identity, plus)
    neg_indices = findall(identity, minus)
    
    # If we don't have exactly one positive vertex, return empty result
    if isnothing(pos_idx) || length(neg_indices) != 3
        return Vector{Vector{Int64}}()
    end
    
    # Create interpolated points for each edge between positive and negative vertices
    ips = Int64[]
    
    for neg_idx in neg_indices
        # Interpolate between current negative vertex and the positive vertex
        ip = interpolate_zero(
            mesh.X[tet[pos_idx]],  # coordinates of positive vertex
            mesh.X[tet[neg_idx]],  # coordinates of negative vertex
            f[pos_idx],           # SDF value of positive vertex
            f[neg_idx],           # SDF value of negative vertex
            mesh
        )
        push!(ips, ip)
    end
    
    # Create one new tetrahedron
    new_tet = [tet[pos_idx], ips[1], ips[2], ips[3]]
    
    # Ensure correct orientation
    if tet_volume(mesh, new_tet) < 0
        # Swap two vertices to change orientation
        new_tet[3], new_tet[4] = new_tet[4], new_tet[3]
    end
    
    return [new_tet]
end

function stencil_3p1(mesh::BlockMesh, tet::Vector{Int64})::Vector{Vector{Int64}}
    # Get SDF values for all vertices
    f = [mesh.node_sdf[i] for i in tet]
    
    # Identify vertices by their position relative to the isosurface
    plus = map(x -> x > 0, f)      # vertices inside the body
    minus = map(x -> x < 0, f)     # vertices outside the body
    
    # Find indices of positive and negative vertices
    pos_indices = findall(identity, plus)
    neg_idx = findfirst(identity, minus)
    
    # If we don't have exactly three positive vertices, return empty result
    if isnothing(neg_idx) || length(pos_indices) != 3
        return Vector{Vector{Int64}}()
    end
    
    # Create interpolated points for each edge between positive and negative vertices
    ips = Int64[]
    
    for pos_idx in pos_indices
        # Interpolate between current positive vertex and the negative vertex
        ip = interpolate_zero(
            mesh.X[tet[pos_idx]],  # coordinates of positive vertex
            mesh.X[tet[neg_idx]],  # coordinates of negative vertex
            f[pos_idx],           # SDF value of positive vertex
            f[neg_idx],           # SDF value of negative vertex
            mesh
        )
        push!(ips, ip)
    end
    
    # Create three new tetrahedra
    new_tets = Vector{Int64}[]
    
    # Get global indices for vertices
    pos1, pos2, pos3 = [tet[i] for i in pos_indices]
    
    # Create tetrahedra connecting each pair of positive vertices with their corresponding
    # interpolated points on the isosurface
    # push!(new_tets, [pos1, pos2, pos3, ips[1]])
    # push!(new_tets, [pos2, pos3, ips[1], ips[2]])
    # push!(new_tets, [pos3, ips[1], ips[2], ips[3]])

    push!(new_tets, [pos1, pos2, pos3, ips[1]])
    push!(new_tets, [ips[1], ips[2], ips[3], pos3])
    push!(new_tets, [pos3, ips[1], pos2, ips[2]])
    
    # Ensure correct orientation
    for i in 1:length(new_tets)
        if tet_volume(mesh, new_tets[i]) < 0
            # Swap two vertices to change orientation
            new_tets[i][3], new_tets[i][4] = new_tets[i][4], new_tets[i][3]
        end
    end
    
    return new_tets
end


function stencil_2p2_variantA(mesh::BlockMesh, tet::Vector{Int64})::Vector{Vector{Int64}}
    # Get SDF values for all vertices
    f = [mesh.node_sdf[i] for i in tet]
    
    # Identify positive and negative vertices
    pos_locals = [i for (i, s) in enumerate(f) if s > 0]
    neg_locals = [i for (i, s) in enumerate(f) if s < 0]
    
    # Verify we have 2 positive and 2 negative vertices
    if length(pos_locals) != 2 || length(neg_locals) != 2
        return Vector{Vector{Int64}}()
    end
    
    # Get global indices of vertices
    pos1_idx = tet[pos_locals[1]]
    pos2_idx = tet[pos_locals[2]]
    neg1_idx = tet[neg_locals[1]]
    neg2_idx = tet[neg_locals[2]]
    
    # Calculate intersections for all edges between positive and negative vertices
    # There are 4 such edges (2 positive × 2 negative)
    
    # Edge 1: pos1 - neg1
    ip1 = interpolate_zero(
        mesh.X[pos1_idx],
        mesh.X[neg1_idx],
        f[pos_locals[1]],
        f[neg_locals[1]],
        mesh
    )
    
    # Edge 2: pos1 - neg2
    ip2 = interpolate_zero(
        mesh.X[pos1_idx],
        mesh.X[neg2_idx],
        f[pos_locals[1]],
        f[neg_locals[2]],
        mesh
    )
    
    # Edge 3: pos2 - neg1
    ip3 = interpolate_zero(
        mesh.X[pos2_idx],
        mesh.X[neg1_idx],
        f[pos_locals[2]],
        f[neg_locals[1]],
        mesh
    )
    
    # Edge 4: pos2 - neg2
    ip4 = interpolate_zero(
        mesh.X[pos2_idx],
        mesh.X[neg2_idx],
        f[pos_locals[2]],
        f[neg_locals[2]],
        mesh
    )
    
    # Calculate cost for both possible diagonals
    diag1 = (ip1, ip3)  # Diagonal 1: ip1-ip3
    diag2 = (ip2, ip4)  # Diagonal 2: ip2-ip4
    
    cost1 = diagonal_cost(mesh, [pos1_idx, pos2_idx], diag1)
    cost2 = diagonal_cost(mesh, [pos1_idx, pos2_idx], diag2)
    
    # Create tetrahedra based on selected diagonal
    new_tets = Vector{Int64}[]
    
    if cost1 <= cost2
        # Use diagonal ip1-ip3
        push!(new_tets, [pos1_idx, pos2_idx, ip1, ip3])
        push!(new_tets, [pos1_idx, ip1, ip2, ip3])
        push!(new_tets, [pos2_idx, ip3, ip4, ip1])
    else
        # Use diagonal ip2-ip4
        push!(new_tets, [pos1_idx, pos2_idx, ip2, ip4])
        push!(new_tets, [pos1_idx, ip1, ip2, ip4])
        push!(new_tets, [pos2_idx, ip3, ip4, ip2])
    end
    
    # Check and correct orientation of tetrahedra
    for i in 1:length(new_tets)
        if tet_volume(mesh, new_tets[i]) < 0
            # Swap two vertices to change orientation
            new_tets[i][3], new_tets[i][4] = new_tets[i][4], new_tets[i][3]
        end
    end
    
    return new_tets
end

function stencil_2p2_variantB(mesh::BlockMesh, tet::Vector{Int64})::Vector{Vector{Int64}}
    # Get SDF values for all vertices
    f = [mesh.node_sdf[i] for i in tet]
    
    # Identify vertices by their position relative to the isosurface
    # Use tolerance for detecting points on the surface
    tol = mesh.grid_tol
    pos_locals = [i for (i, s) in enumerate(f) if s > tol]              # inside body
    zero_locals = [i for (i, s) in enumerate(f) if abs(s) <= tol]       # on surface
    neg_locals = [i for (i, s) in enumerate(f) if s < -tol]             # outside body
    
    # Verify case variant B:
    # Two positive vertices, one on surface and one negative
    if length(pos_locals) != 2 || length(zero_locals) != 1 || length(neg_locals) != 1
        return Vector{Vector{Int64}}()
    end
    
    # Get global indices of vertices
    pos1_idx = tet[pos_locals[1]]
    pos2_idx = tet[pos_locals[2]]
    zero_idx = tet[zero_locals[1]]
    neg_idx = tet[neg_locals[1]]
    
    # Calculate intersection between positive vertices and negative vertex
    ip1 = interpolate_zero(
        mesh.X[pos1_idx],
        mesh.X[neg_idx],
        f[pos_locals[1]],
        f[neg_locals[1]],
        mesh
    )
    
    ip2 = interpolate_zero(
        mesh.X[pos2_idx],
        mesh.X[neg_idx],
        f[pos_locals[2]],
        f[neg_locals[1]],
        mesh
    )
    
    # Create two new tetrahedra
    tet1 = [pos1_idx, zero_idx, ip1, ip2]
    tet2 = [pos2_idx, zero_idx, ip2, ip1]
    
    # Check and correct orientation
    if tet_volume(mesh, tet1) < 0
        tet1 = [pos1_idx, zero_idx, ip2, ip1]
    end
    
    if tet_volume(mesh, tet2) < 0
        tet2 = [pos2_idx, zero_idx, ip1, ip2]
    end
    
    return [tet1, tet2]
end

function diagonal_cost(mesh::BlockMesh, pos_indices::Vector{Int64}, diag_pair::Tuple{Int64,Int64})
    # Get coordinates of diagonal endpoints
    p_diag1 = mesh.X[diag_pair[1]]
    p_diag2 = mesh.X[diag_pair[2]]
    
    # Base cost
    cost = 0.0
    
    # Add penalty if diagonal is not entirely inside the body
    col = edge_color(mesh, diag_pair[1], diag_pair[2])
    penalty = (col == :positive ? 0.0 : 1000.0)
    cost += penalty
    
    # Add cost based on distance from positive vertices
    for pi in pos_indices
        p_pos = mesh.X[pi]
        cost += norm(p_pos - p_diag1) + norm(p_pos - p_diag2)
    end
    
    # Add cost based on diagonal length (prefer shorter diagonals)
    diag_length = norm(p_diag1 - p_diag2)
    cost += diag_length
    
    return cost
end


function stencil_1p3_with_zero(mesh::BlockMesh, tet::Vector{Int64})::Vector{Vector{Int64}}
    f = [mesh.node_sdf[i] for i in tet]
    tol = mesh.grid_tol
    
    # Identify vertices by their position
    pos_locals = [i for (i, s) in enumerate(f) if s > tol]
    zero_locals = [i for (i, s) in enumerate(f) if abs(s) <= tol]
    neg_locals = [i for (i, s) in enumerate(f) if s < -tol]
    
    # Verify we have the correct configuration
    if length(pos_locals) != 1 || length(zero_locals) != 1 || length(neg_locals) != 2
        return Vector{Vector{Int64}}()
    end
    
    # Get global indices
    pos_idx = tet[pos_locals[1]]
    zero_idx = tet[zero_locals[1]]
    neg1_idx = tet[neg_locals[1]]
    neg2_idx = tet[neg_locals[2]]
    
    # Calculate intersections
    ip1 = interpolate_zero(
        mesh.X[pos_idx],
        mesh.X[neg1_idx],
        f[pos_locals[1]],
        f[neg_locals[1]],
        mesh
    )
    
    ip2 = interpolate_zero(
        mesh.X[pos_idx],
        mesh.X[neg2_idx],
        f[pos_locals[1]],
        f[neg_locals[2]],
        mesh
    )
    
    # Create tetrahedron
    new_tet = [pos_idx, zero_idx, ip1, ip2]
    
    # Check orientation
    if tet_volume(mesh, new_tet) < 0
        new_tet = [pos_idx, zero_idx, ip2, ip1]
    end
    
    return [new_tet]
end

function stencil_1p1_with_two_zeros(mesh::BlockMesh, tet::Vector{Int64})::Vector{Vector{Int64}}
    f = [mesh.node_sdf[i] for i in tet]
    tol = mesh.grid_tol
    
    # Identify vertices by their position
    pos_locals = [i for (i, s) in enumerate(f) if s > tol]
    zero_locals = [i for (i, s) in enumerate(f) if abs(s) <= tol]
    neg_locals = [i for (i, s) in enumerate(f) if s < -tol]
    
    # Verify we have the correct configuration
    if length(pos_locals) != 1 || length(zero_locals) != 2 || length(neg_locals) != 1
        return Vector{Vector{Int64}}()
    end
    
    # Get global indices
    pos_idx = tet[pos_locals[1]]
    zero1_idx = tet[zero_locals[1]]
    zero2_idx = tet[zero_locals[2]]
    neg_idx = tet[neg_locals[1]]
    
    # Calculate intersection
    ip = interpolate_zero(
        mesh.X[pos_idx],
        mesh.X[neg_idx],
        f[pos_locals[1]],
        f[neg_locals[1]],
        mesh
    )
    
    # Create tetrahedron
    new_tet = [pos_idx, zero1_idx, zero2_idx, ip]
    
    # Check orientation
    if tet_volume(mesh, new_tet) < 0
        new_tet = [pos_idx, zero1_idx, ip, zero2_idx]
    end
    
    return [new_tet]
end

function handle_special_cases(mesh::BlockMesh, tet::Vector{Int64}, n_pos::Int, n_neg::Int, n_zero::Int)::Vector{Vector{Int64}}
    # This function handles any remaining cases that don't fit the standard stencils
    
    # Keep tetrahedron if it's fully on or inside the body
    if n_neg == 0
        return [tet]
    end
    
    # Handle other potential cases
    if n_pos > 0 && n_zero > 0
        # There's at least one positive vertex and one on the boundary,
        # We'll try to keep it with proper orientation
        result_tet = copy(tet)
        if tet_volume(mesh, result_tet) < 0
            # Swap two vertices to fix orientation
            result_tet[3], result_tet[4] = result_tet[4], result_tet[3]
        end
        return [result_tet]
    end
    
    # For cases we're not sure about, log them for future inspection
    @debug "Unhandled stencil case: $(n_pos) positive, $(n_neg) negative, $(n_zero) zero"
    
    # Default: discard the tetrahedron
    return Vector{Vector{Int64}}()
end

function tet_volume(mesh::BlockMesh, tet::Vector{Int64})::Float64
    # Get coordinates of tetrahedron vertices
    p1 = mesh.X[tet[1]]
    p2 = mesh.X[tet[2]]
    p3 = mesh.X[tet[3]]
    p4 = mesh.X[tet[4]]
    
    # Calculate vectors from first vertex to others
    v1 = p2 - p1
    v2 = p3 - p1
    v3 = p4 - p1
    
    # Calculate volume using triple product
    # Volume = (1/6) * |v1 · (v2 × v3)|
    volume = dot(v1, cross(v2, v3)) / 6.0
    
    return volume
end

function stencil_1p_3zero(mesh::BlockMesh, tet::Vector{Int64})::Vector{Vector{Int64}}
    f = [mesh.node_sdf[i] for i in tet]
    tol = mesh.grid_tol
    
    # Identify vertices by their position
    pos_locals = [i for (i, s) in enumerate(f) if s > tol]
    zero_locals = [i for (i, s) in enumerate(f) if abs(s) <= tol]
    
    # Verify we have the correct configuration
    if length(pos_locals) != 1 || length(zero_locals) != 3
        return Vector{Vector{Int64}}()
    end
    
    # Get global indices
    pos_idx = tet[pos_locals[1]]
    zero1_idx = tet[zero_locals[1]]
    zero2_idx = tet[zero_locals[2]]
    zero3_idx = tet[zero_locals[3]]
    
    # Create tetrahedron - keep it unchanged as it's fully within the body
    new_tet = [pos_idx, zero1_idx, zero2_idx, zero3_idx]
    
    # Check orientation
    if tet_volume(mesh, new_tet) < 0
        # Swap two vertices to correct orientation
        new_tet[3], new_tet[4] = new_tet[4], new_tet[3]
    end
    
    return [new_tet]
end

function apply_stencil(mesh::BlockMesh, tet::Vector{Int64})::Vector{Vector{Int64}}
    f = [mesh.node_sdf[i] for i in tet]
    tol = mesh.grid_tol
    
    # Count vertices with different SDF values
    n_pos = count(x -> x > tol, f)                 # Strictly inside
    n_neg = count(x -> x < -tol, f)                # Strictly outside
    n_zero = count(x -> abs(x) <= tol, f)          # On the boundary
    n_pos_zero = count(x -> x > -tol, f)           # Inside or on boundary
    
    # All vertices are on or inside the body
    if n_pos_zero == 4
        return [tet]  # Keep the tetrahedron unchanged
        # TODO: Chybí elementy se 3 uzly na hladině a s 2 uzly na hladině. -> je to ok
    end
    
    # All vertices are outside the body - discard the tetrahedron
    if n_pos == 0 #&& n_zero == 0
        return Vector{Vector{Int64}}()
    end

    if n_pos == 1 && n_zero == 3
        # println("zero!")
    end
    
    # Select appropriate stencil based on vertex counts
    if n_pos == 1 && n_neg == 3
        return stencil_1p3(mesh, tet)
    elseif n_pos == 3 && n_neg == 1
        return stencil_3p1(mesh, tet)
    elseif n_pos == 2 && n_neg == 2
        # return stencil_2p2_variantA(mesh, tet)
        return Vector{Vector{Int64}}()
    elseif n_pos == 2 && n_neg == 1 && n_zero == 1
        return Vector{Vector{Int64}}()
        # return stencil_2p2_variantB(mesh, tet) # -> asi ok?
    elseif n_pos == 1 && n_neg == 2 && n_zero == 1
        return stencil_1p3_with_zero(mesh, tet)
    elseif n_pos == 1 && n_neg == 1 && n_zero == 2
        return stencil_1p1_with_two_zeros(mesh, tet)
    elseif n_pos == 1 && n_zero == 3
        println("zero!")
        return [tet]
        # return stencil_1p_3zero(mesh, tet)
    else
        # Handle any remaining cases
        return handle_special_cases(mesh, tet, n_pos, n_neg, n_zero)
    end
end