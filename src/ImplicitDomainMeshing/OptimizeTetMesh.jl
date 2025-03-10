# ----------------------------
# Quality metric for tetrahedra
# ----------------------------
function compute_tet_quality(v1::SVector{3,Float64}, v2::SVector{3,Float64}, 
                            v3::SVector{3,Float64}, v4::SVector{3,Float64})
    # Calculate edge vectors
    edges = [
        v2 - v1, v3 - v1, v4 - v1,
        v3 - v2, v4 - v2,
        v4 - v3
    ]
    
    # Calculate edge lengths squared
    edge_lengths_squared = [dot(e, e) for e in edges]
    
    # Calculate volume
    volume = dot(cross(v2 - v1, v3 - v1), v4 - v1) / 6.0
    
    if volume <= 0
        return -1.0  # Inverted or degenerate tetrahedron
    end
    
    # Calculate aspect ratio (using the formula from the paper)
    max_edge_length_squared = maximum(edge_lengths_squared)
    max_edge_length = sqrt(max_edge_length_squared)
    
    # Scale factor for better numerical behavior
    volume_scaled = volume * 8.48528137423857  # 6 * sqrt(2)
    
    # Quality metric: normalized by dividing by the cube of the longest edge
    quality = volume_scaled / (max_edge_length * max_edge_length * max_edge_length)
    
    return quality
end

# Compute quality for a tetrahedron in the mesh given by its index
function compute_element_quality(mesh::BlockMesh, tet_idx::Int)
    tet = mesh.IEN[tet_idx]
    v1, v2, v3, v4 = mesh.X[tet[1]], mesh.X[tet[2]], mesh.X[tet[3]], mesh.X[tet[4]]
    return compute_tet_quality(v1, v2, v3, v4)
end

# ----------------------------
# Check if moving a node would cause element inversion
# ----------------------------
function would_cause_inversion(mesh::BlockMesh, node_idx::Int, new_position::SVector{3,Float64})
    # Store original position
    original_position = mesh.X[node_idx]
    
    # Temporarily move the node
    mesh.X[node_idx] = new_position
    
    # Check all elements containing this node
    inversion_found = false
    for elem_idx in mesh.INE[node_idx]
        tet = mesh.IEN[elem_idx]
        v1, v2, v3, v4 = mesh.X[tet[1]], mesh.X[tet[2]], mesh.X[tet[3]], mesh.X[tet[4]]
        
        # Calculate volume (proportional to determinant)
        volume = dot(cross(v2 - v1, v3 - v1), v4 - v1)
        
        if volume <= 0
            inversion_found = true
            break
        end
    end
    
    # Restore original position
    mesh.X[node_idx] = original_position
    
    return inversion_found
end

# ----------------------------
# Laplacian smoothing with quality control
# ----------------------------
function laplacian_smooth!(mesh::BlockMesh; max_iterations::Int=20, damping_factor::Float64=0.5, 
                          quality_threshold::Float64=0.0, tolerance::Float64=1e-6)
    @info "Starting Laplacian smoothing optimization..."
    
    # Identify boundary nodes (nodes on the isosurface)
    boundary_nodes = Set{Int}()
    for i in 1:length(mesh.X)
        if abs(mesh.node_sdf[i]) < tolerance
            push!(boundary_nodes, i)
        end
    end
    
    # Perform smoothing iterations
    for iter in 1:max_iterations
        max_displacement = 0.0
        nodes_moved = 0
        
        # For each node
        for node_idx in 1:length(mesh.X)
            # Skip boundary nodes - we want to preserve the isosurface
            if node_idx in boundary_nodes
                continue
            end
            
            # Get all adjacent nodes from the elements
            adjacent_nodes = Set{Int}()
            for elem_idx in mesh.INE[node_idx]
                for n in mesh.IEN[elem_idx]
                    if n != node_idx
                        push!(adjacent_nodes, n)
                    end
                end
            end
            
            # Skip isolated nodes
            if isempty(adjacent_nodes)
                continue
            end
            
            # Calculate barycenter of adjacent nodes
            barycenter = sum(mesh.X[n] for n in adjacent_nodes) / length(adjacent_nodes)
            
            # Calculate proposed new position with damping
            new_position = (1.0 - damping_factor) * mesh.X[node_idx] + damping_factor * barycenter
            
            # Skip if the displacement is negligible
            displacement = norm(new_position - mesh.X[node_idx])
            if displacement < tolerance
                continue
            end
            
            # Check if the move would cause element inversion
            if would_cause_inversion(mesh, node_idx, new_position)
                # Try reduced damping if normal damping causes inversion
                reduced_damping = damping_factor * 0.5
                new_position = (1.0 - reduced_damping) * mesh.X[node_idx] + reduced_damping * barycenter
                
                # Skip if still causing inversion
                if would_cause_inversion(mesh, node_idx, new_position)
                    continue
                end
                
                # Recalculate displacement with reduced damping
                displacement = norm(new_position - mesh.X[node_idx])
            end
            
            # Check if the move would degrade quality too much
            current_min_quality = minimum(compute_element_quality(mesh, elem_idx) for elem_idx in mesh.INE[node_idx])
            
            # Temporarily apply the move
            old_position = mesh.X[node_idx]
            mesh.X[node_idx] = new_position
            
            # Compute new quality
            new_min_quality = minimum(compute_element_quality(mesh, elem_idx) for elem_idx in mesh.INE[node_idx])
            
            # Revert if quality decreases below threshold
            if new_min_quality < quality_threshold || new_min_quality < current_min_quality * 0.9
                mesh.X[node_idx] = old_position
                continue
            end
            
            # Update maximum displacement
            max_displacement = max(max_displacement, displacement)
            nodes_moved += 1
        end
        
        @info "Iteration $iter: Max displacement = $max_displacement, Nodes moved = $nodes_moved"
        
        # Check for convergence
        if max_displacement < tolerance || nodes_moved == 0
            @info "Smoothing converged after $iter iterations"
            break
        end
    end
    
    @info "Laplacian smoothing completed"
end

# ----------------------------
# Improved mesh optimization with adaptive damping
# ----------------------------
function optimize_mesh_adaptive!(mesh::BlockMesh; max_iterations::Int=20, 
                               initial_damping::Float64=0.5, min_damping::Float64=0.1,
                               quality_threshold::Float64=0.01, tolerance::Float64=1e-6)
    @info "Starting adaptive mesh optimization..."
    
    # Identify boundary nodes (nodes on the isosurface)
    boundary_nodes = Set{Int}()
    for i in 1:length(mesh.X)
        if abs(mesh.node_sdf[i]) < tolerance
            push!(boundary_nodes, i)
        end
    end
    
    # Compute initial quality for all elements
    element_qualities = [compute_element_quality(mesh, i) for i in 1:length(mesh.IEN)]
    
    # Track quality improvement
    initial_min_quality = minimum(element_qualities)
    initial_avg_quality = sum(element_qualities) / length(element_qualities)
    
    @info "Initial mesh quality: Min = $initial_min_quality, Avg = $initial_avg_quality"
    
    # For each iteration
    for iter in 1:max_iterations
        max_displacement = 0.0
        nodes_moved = 0
        
        # Process nodes in random order to avoid bias
        node_indices = collect(1:length(mesh.X))
        shuffle!(node_indices)
        
        for node_idx in node_indices
            # Skip boundary nodes - preserve the isosurface
            if node_idx in boundary_nodes
                continue
            end
            
            # Get adjacent nodes
            adjacent_nodes = Set{Int}()
            for elem_idx in mesh.INE[node_idx]
                for n in mesh.IEN[elem_idx]
                    if n != node_idx
                        push!(adjacent_nodes, n)
                    end
                end
            end
            
            isempty(adjacent_nodes) && continue
            
            # Calculate barycenter of adjacent nodes
            barycenter = sum(mesh.X[n] for n in adjacent_nodes) / length(adjacent_nodes)
            
            # Get incident elements and their qualities
            incident_elements = mesh.INE[node_idx]
            element_min_quality = minimum(element_qualities[e] for e in incident_elements)
            
            # Adaptive damping based on minimum quality
            # Lower quality elements get less aggressive movement
            damping = initial_damping
            if element_min_quality < quality_threshold
                # Scale damping down for low quality elements
                quality_ratio = max(element_min_quality / quality_threshold, 0.1)
                damping = max(initial_damping * quality_ratio, min_damping)
            end
            
            # Calculate proposed new position with adaptive damping
            new_position = (1.0 - damping) * mesh.X[node_idx] + damping * barycenter
            
            # Skip if displacement is negligible
            displacement = norm(new_position - mesh.X[node_idx])
            displacement < tolerance && continue
            
            # Binary search for maximum safe damping if inversion would occur
            if would_cause_inversion(mesh, node_idx, new_position)
                low = 0.0
                high = damping
                safe_position = mesh.X[node_idx]
                
                # Binary search for maximum safe damping
                for _ in 1:10  # Limit binary search iterations
                    mid = (low + high) / 2
                    test_position = (1.0 - mid) * mesh.X[node_idx] + mid * barycenter
                    
                    if would_cause_inversion(mesh, node_idx, test_position)
                        high = mid
                    else
                        low = mid
                        safe_position = test_position
                    end
                end
                
                # If we couldn't find a safe damping, skip this node
                if low < min_damping
                    continue
                end
                
                new_position = safe_position
                displacement = norm(new_position - mesh.X[node_idx])
            end
            
            # Store old position and apply the move
            old_position = mesh.X[node_idx]
            mesh.X[node_idx] = new_position
            
            # Update element qualities for affected elements
            quality_improved = true
            for elem_idx in incident_elements
                new_quality = compute_element_quality(mesh, elem_idx)
                if new_quality < element_qualities[elem_idx] * 0.95  # Allow 5% quality reduction
                    quality_improved = false
                    break
                end
                # Update quality
                element_qualities[elem_idx] = new_quality
            end
            
            # Revert if quality significantly decreases
            if !quality_improved
                mesh.X[node_idx] = old_position
                # Restore original qualities
                for elem_idx in incident_elements
                    element_qualities[elem_idx] = compute_element_quality(mesh, elem_idx)
                end
                continue
            end
            
            # Track maximum displacement and count moved nodes
            max_displacement = max(max_displacement, displacement)
            nodes_moved += 1
        end
        
        # Check for convergence
        @info "Iteration $iter: Max displacement = $max_displacement, Nodes moved = $nodes_moved"
        
        if max_displacement < tolerance || nodes_moved == 0
            @info "Optimization converged after $iter iterations"
            break
        end
    end
    
    # Final quality assessment
    final_qualities = [compute_element_quality(mesh, i) for i in 1:length(mesh.IEN)]
    final_min_quality = minimum(final_qualities)
    final_avg_quality = sum(final_qualities) / length(final_qualities)
    quality_improvement = (final_avg_quality - initial_avg_quality) / initial_avg_quality * 100
    
    @info "Mesh optimization complete. Quality improvement: $quality_improvement%"
    @info "Final mesh quality: Min = $final_min_quality, Avg = $final_avg_quality"
    
    # Update connectivity 
    update_connectivity!(mesh)
    
    return mesh
end

# ----------------------------
# Surface-preserving optimization
# ----------------------------
function optimize_surface_nodes!(mesh::BlockMesh; max_iterations::Int=10, 
                               damping_factor::Float64=0.3, tolerance::Float64=1e-6)
    @info "Starting surface node optimization..."
    
    # Identify surface nodes (on the isosurface)
    surface_nodes = Set{Int}()
    for i in 1:length(mesh.X)
        if abs(mesh.node_sdf[i]) < tolerance
            push!(surface_nodes, i)
        end
    end
    
    @info "Found $(length(surface_nodes)) surface nodes to optimize"
    
    # For each iteration
    for iter in 1:max_iterations
        max_displacement = 0.0
        nodes_moved = 0
        
        # For each surface node
        for node_idx in surface_nodes
            # Get adjacent nodes that are also on the surface
            adjacent_surface_nodes = Set{Int}()
            
            # Get all elements containing this node
            for elem_idx in mesh.INE[node_idx]
                tet = mesh.IEN[elem_idx]
                
                # Find other nodes from this element that are on the surface
                for n in tet
                    if n != node_idx && n in surface_nodes
                        push!(adjacent_surface_nodes, n)
                    end
                end
            end
            
            # Need at least 2 adjacent surface nodes to do tangential smoothing
            if length(adjacent_surface_nodes) < 2
                continue
            end
            
            # Calculate the normal vector at this node
            normal = compute_gradient(mesh, node_idx)
            normal_length = norm(normal)
            
            if normal_length < tolerance
                continue
            end
            
            # Normalize the normal vector
            normal = normal / normal_length
            
            # Calculate barycenter of adjacent surface nodes
            barycenter = sum(mesh.X[n] for n in adjacent_surface_nodes) / length(adjacent_surface_nodes)
            
            # Calculate the vector from current node to barycenter
            move_vector = barycenter - mesh.X[node_idx]
            
            # Project the move vector onto the tangent plane
            projected_move = move_vector - dot(move_vector, normal) * normal
            
            # Apply damping to the tangential movement
            new_position = mesh.X[node_idx] + damping_factor * projected_move
            
            # Calculate displacement
            displacement = norm(projected_move) * damping_factor
            
            # Skip if displacement is negligible
            if displacement < tolerance
                continue
            end
            
            # Check if this move would invert any elements
            if would_cause_inversion(mesh, node_idx, new_position)
                continue
            end
            
            # Apply the move
            mesh.X[node_idx] = new_position
            
            # Project back onto the isosurface using the SDF
            warp_node_to_isocontour!(mesh, node_idx, 5)  # Use fewer iterations here
            
            # Count the node as moved
            max_displacement = max(max_displacement, displacement)
            nodes_moved += 1
        end
        
        @info "Surface optimization iteration $iter: Max displacement = $max_displacement, Nodes moved = $nodes_moved"
        
        # Check for convergence
        if max_displacement < tolerance || nodes_moved == 0
            @info "Surface optimization converged after $iter iterations"
            break
        end
    end
    
    @info "Surface node optimization complete"
    return mesh
end

# ----------------------------
# Complete mesh optimization pipeline
# ----------------------------
function optimize_mesh!(mesh::BlockMesh; 
                               interior_iterations::Int=20, 
                               surface_iterations::Int=10,
                               adaptive_iterations::Int=15,
                               damping_factor::Float64=0.5)
    @info "Starting complete mesh optimization..."
    
    # First optimize the interior with standard Laplacian smoothing
    laplacian_smooth!(mesh, max_iterations=interior_iterations, damping_factor=damping_factor)
    
    # Then optimize the surface nodes while preserving the isosurface
    optimize_surface_nodes!(mesh, max_iterations=surface_iterations, damping_factor=damping_factor*0.6)
    
    # Finally, apply adaptive optimization with quality-based damping
    optimize_mesh_adaptive!(mesh, max_iterations=adaptive_iterations, initial_damping=damping_factor)
    
    # Final connectivity update
    update_connectivity!(mesh)
    
    # Calculate and report final mesh quality
    quality_stats = [compute_element_quality(mesh, i) for i in 1:length(mesh.IEN)]
    min_quality = minimum(quality_stats)
    max_quality = maximum(quality_stats)
    avg_quality = sum(quality_stats) / length(quality_stats)
    med_quality = sort(quality_stats)[div(length(quality_stats),2)]
    
    @info "Complete mesh optimization finished."
    @info "Final mesh quality statistics:"
    @info "  Min quality: $min_quality"
    @info "  Median quality: $med_quality"
    @info "  Avg quality: $avg_quality"
    @info "  Max quality: $max_quality"
    
    return mesh
end
