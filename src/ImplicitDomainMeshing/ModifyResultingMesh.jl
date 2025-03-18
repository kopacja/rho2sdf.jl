# Function to check if a point lies on a bounded plane
function is_on_plane(plane::BoundedPlane, point::Vector{Float64}, tolerance::Float64=1e-10)
    # Check if the point lies in the plane
    dist_to_plane = abs(dot(plane.normal, point - plane.point))
    if dist_to_plane > tolerance
        return false
    end
    
    # Project the point to the plane's coordinate system
    vec_to_point = point - plane.point
    
    # Coordinates of the point in the plane
    u_coord = dot(vec_to_point, plane.u)
    v_coord = dot(vec_to_point, plane.v)
    
    # Check if the projected point lies within the specified shape
    return is_in_shape(plane.shape, u_coord, v_coord)
end

# Check if a point (in plane coordinates) lies within a rectangle
function is_in_shape(shape::Rectangle, u::Float64, v::Float64)
    return abs(u) <= shape.width/2 && abs(v) <= shape.height/2
end

# Check if a point (in plane coordinates) lies within a circle
function is_in_shape(shape::Circle, u::Float64, v::Float64)
    return u^2 + v^2 <= shape.radius^2
end

# Check if a point (in plane coordinates) lies within an ellipse
function is_in_shape(shape::Ellipse, u::Float64, v::Float64)
    return (u/shape.a)^2 + (v/shape.b)^2 <= 1
end

function distance_to_bounded_plane(plane::BoundedPlane, point::Vector{Float32})
    # Vector from the plane point to the current point
    vec_to_point = point - plane.point
    
    # Calculate dot product to determine which side of the plane
    dot_product = dot(plane.normal, vec_to_point)
    
    # Distance to infinite plane (absolute value for checks)
    dist_to_infinite_plane = abs(dot_product)
    
    # Project point onto the infinite plane
    projected_point = point - dot_product * plane.normal
    
    # Check if the projected point lies within the bounded plane
    if is_on_plane(plane, projected_point)
        # If the projection lies within the bounded plane, return distance to the infinite plane
        # with sign - negative in normal direction, positive in opposite direction
        return dot_product > 0 ? -dist_to_infinite_plane : dist_to_infinite_plane
    else
        # If the projection is outside the bounded plane, return a high value
        # while maintaining the correct sign
        return dot_product > 0 ? -1.0e10 : 1.0e10
    end
end

# Function to evaluate the distance of a point from planes (planes_sdf)
function eval_planes_sdf(mesh::BlockMesh, p::SVector{3,Float64}, plane_definitions::Vector{PlaneDefinition})
    # Create bounded planes from definitions
    planes = [BoundedPlane(def.normal, def.point, def.shape) for def in plane_definitions]
    
    # Initialize with a high positive value
    min_dist = 1.0e10
    
    # Check distance to each plane
    for plane in planes
        # Convert to Vector{Float32} for compatibility with the distance_to_bounded_plane function
        p_float32 = Vector{Float32}([p[1], p[2], p[3]])
        dist = distance_to_bounded_plane(plane, p_float32)
        
        # Update minimum distance (comparing absolute values)
        if abs(dist) < abs(min_dist)
            min_dist = dist
        end
    end
    
    return min_dist
end

# Approximate the gradient of planes_sdf to determine warp direction
function approximate_planes_gradient(mesh::BlockMesh, p::SVector{3,Float64}, plane_definitions::Vector{PlaneDefinition}; h::Float64=1e-3)
    dx = SVector{3,Float64}(h, 0.0, 0.0)
    dy = SVector{3,Float64}(0.0, h, 0.0)
    dz = SVector{3,Float64}(0.0, 0.0, h)
    
    # Calculate partial derivatives using central differences
    df_dx = (eval_planes_sdf(mesh, p + dx, plane_definitions) - eval_planes_sdf(mesh, p - dx, plane_definitions)) / (2 * h)
    df_dy = (eval_planes_sdf(mesh, p + dy, plane_definitions) - eval_planes_sdf(mesh, p - dy, plane_definitions)) / (2 * h)
    df_dz = (eval_planes_sdf(mesh, p + dz, plane_definitions) - eval_planes_sdf(mesh, p - dz, plane_definitions)) / (2 * h)
    
    return SVector{3,Float64}(df_dx, df_dy, df_dz)
end

# Function to move a node to the zero level of planes_sdf
function warp_node_to_planes_isocontour!(mesh::BlockMesh, node_index::Int, plane_definitions::Vector{PlaneDefinition}, max_iter::Int)
    tol = mesh.grid_tol
    current_position = mesh.X[node_index]
    
    for iter in 1:max_iter
        # Evaluate planes_sdf at current position
        f = eval_planes_sdf(mesh, current_position, plane_definitions)
        
        # If we're close enough to the isosurface, end
        abs2(f) < tol * tol && break
        
        # Calculate gradient for displacement direction
        grad = approximate_planes_gradient(mesh, current_position, plane_definitions)
        norm_grad_squared = sum(abs2, grad)
        
        # If gradient is too small, end
        norm_grad_squared < 1e-16 && break
        
        # Newton step
        dp = (f / norm_grad_squared) * grad
        current_position -= dp
    end
    
    # Calculate current planes_sdf value after warping
    current_sdf = eval_planes_sdf(mesh, current_position, plane_definitions)
    
    # Update node position and its SDF value
    mesh.X[node_index] = current_position
    if abs(current_sdf) < tol*2
        mesh.node_sdf[node_index] = 0.0
    else
        mesh.node_sdf[node_index] = current_sdf
    end
end

# Main function for modifying the mesh according to planes_sdf
function warp_mesh_by_planes_sdf!(mesh::BlockMesh, plane_definitions::Vector{PlaneDefinition}, warp_param::Float64; max_iter::Int=20)
    # Check if any plane definitions were provided
    if isempty(plane_definitions)
        @info "No plane definitions provided, skipping plane-based warping."
        return
    end
    
    # Calculate the longest edge and then the threshold for displacement - same logic as in the warp! function
    max_edge = longest_edge(mesh)
    threshold_sdf = warp_param * max_edge
    
    @info "Planes warping: max edge = $max_edge, threshold_sdf = $threshold_sdf"
    
    # Preparation for warping - find nodes near planes (in one direction)
    nodes_to_warp = Int[]
    for i in 1:length(mesh.X)
        plane_sdf = eval_planes_sdf(mesh, mesh.X[i], plane_definitions)
        # If the node is close enough to the plane (from both sides), add it to nodes to warp
        # if abs(plane_sdf) < threshold_sdf
        if plane_sdf < threshold_sdf    
            push!(nodes_to_warp, i)
        end
    end
    
    @info "Found $(length(nodes_to_warp)) nodes near planes for warping"
    
    # Warp nodes near planes to the zero level
    for node_idx in nodes_to_warp
        warp_node_to_planes_isocontour!(mesh, node_idx, plane_definitions, max_iter)
    end
    
    # Update mesh topology
    update_connectivity!(mesh)
    
    # Remove elements that have all nodes at the zero level or don't have at least one node with positive SDF
    new_IEN = Vector{Vector{Int64}}()
    for tet in mesh.IEN
        # Get SDF values for all nodes of the tetrahedron
        tet_sdf = [eval_planes_sdf(mesh, mesh.X[i], plane_definitions) for i in tet]
        
        # Count nodes at zero level
        zero_nodes = count(x -> abs(x) < mesh.grid_tol, tet_sdf)
        
        # Count nodes with positive SDF value
        pos_nodes = count(x -> x > 0, tet_sdf)
        
        # Keep the tetrahedron only if:
        # 1. It has at least one node with positive value (relative to planes_sdf)
        # 2. It doesn't have all nodes at the zero level
        if pos_nodes > 0 && zero_nodes < 4
            push!(new_IEN, tet)
        end
    end
    
    # Update connectivity
    mesh.IEN = new_IEN
    @info "After removing unsuitable elements: $(length(mesh.IEN)) tetrahedra"
    
    # Final topology update
    update_connectivity!(mesh)
    
    @info "Mesh modification according to planes_sdf completed"
end

#TODO: Create a mesh trimming that will support element cutting and subsequent mesh optimization.

# # Function to warp nodes onto plane surfaces based on element connectivity
# function adjust_nodes_to_planes!(mesh::BlockMesh, plane_definitions::Vector{PlaneDefinition})
#     @info "Starting node warping to plane surfaces..."
    
#     # Track nodes that have been processed to avoid duplicating work
#     processed_nodes = Set{Int}()
#     nodes_moved = 0
#     nodes_skipped = 0
    
#     # Helper function to check if moving a node would cause element inversion
#     function would_cause_inversion(node_idx, new_position)
#         # Store original position
#         original_position = mesh.X[node_idx]
        
#         # Temporarily move the node
#         mesh.X[node_idx] = new_position
        
#         # Check all elements containing this node
#         inversion_found = false
#         for elem_idx in mesh.INE[node_idx]
#             tet = mesh.IEN[elem_idx]
#             vertices = [mesh.X[v] for v in tet]
            
#             # Calculate Jacobian determinant (proportional to signed volume)
#             jacobian_det = dot(
#                 cross(vertices[2] - vertices[1], vertices[3] - vertices[1]),
#                 vertices[4] - vertices[1]
#             )
            
#             if jacobian_det < 0
#                 inversion_found = true
#                 break
#             end
#         end
        
#         # Restore original position
#         mesh.X[node_idx] = original_position
        
#         return inversion_found
#     end
    
#     # Helper function for interpolating to plane surface
#     function interpolate_to_plane(p1::SVector{3,Float64}, p2::SVector{3,Float64}, 
#                                  f1::Float64, f2::Float64)
#         # Ensure f1 and f2 have opposite signs (crossing the plane)
#         if sign(f1) == sign(f2)
#             error("Points must be on opposite sides of the plane")
#         end
        
#         # Linear interpolation to find zero crossing
#         t = f1 / (f1 - f2)
#         return p1 + t * (p2 - p1)
#     end
    
#     # First pass: Find nodes with negative plane_sdf and process them
#     for node_idx in 1:length(mesh.X)
#         # Skip already processed nodes
#         node_idx in processed_nodes && continue
        
#         # Evaluate plane_sdf for this node
#         node_plane_sdf = eval_planes_sdf(mesh, mesh.X[node_idx], plane_definitions)
        
#         # Only process nodes with negative plane_sdf values (inside cutting planes)
#         if node_plane_sdf < 0
#             # Get all elements containing this node
#             elements = mesh.INE[node_idx]
            
#             # Skip nodes that aren't part of any elements
#             isempty(elements) && continue
            
#             # Try to find an element with exactly 3 negative nodes
#             edge_to_warp = nothing
            
#             for elem_idx in elements
#                 tet = mesh.IEN[elem_idx]
                
#                 # Get plane_sdf values for all nodes in this tetrahedron
#                 tet_plane_sdf = [eval_planes_sdf(mesh, mesh.X[i], plane_definitions) for i in tet]
                
#                 # Count negative nodes in this tetrahedron
#                 neg_nodes = [i for (i, sdf) in zip(tet, tet_plane_sdf) if sdf < 0]
                
#                 if length(neg_nodes) == 3 && node_idx in neg_nodes
#                     # This element has exactly 3 negative nodes including our target node
                    
#                     # Find the positive node
#                     pos_node = first([i for (i, sdf) in zip(tet, tet_plane_sdf) if sdf >= 0])
                    
#                     # Calculate the plane_sdf values
#                     pos_sdf = eval_planes_sdf(mesh, mesh.X[pos_node], plane_definitions)
#                     neg_sdf = node_plane_sdf
                    
#                     # Check if the edge between positive and negative node crosses zero
#                     if sign(pos_sdf) != sign(neg_sdf)
#                         # Calculate approximate zero crossing using linear interpolation
#                         zero_point = interpolate_to_plane(
#                             mesh.X[pos_node], 
#                             mesh.X[node_idx], 
#                             pos_sdf, 
#                             neg_sdf
#                         )
                        
#                         # Check if moving to this point would cause inversion
#                         if !would_cause_inversion(node_idx, zero_point)
#                             # Found a suitable edge to warp along
#                             edge_to_warp = (pos_node, node_idx)
#                             break
#                         end
#                     end
#                 end
#             end
            
#             # If no suitable edge found yet, try all edges connecting to positive nodes
#             if edge_to_warp === nothing
#                 for elem_idx in elements
#                     tet = mesh.IEN[elem_idx]
        
#                     # Look for edges connecting our negative node to positive nodes
#                     for other_node in tet
#                         # Skip if it's the same node
#                         other_node == node_idx && continue
            
#                         # Evaluate plane_sdf for the other node
#                         other_plane_sdf = eval_planes_sdf(mesh, mesh.X[other_node], plane_definitions)
                        
#                         # Check if the other node has positive plane_sdf (crosses zero level)
#                         if other_plane_sdf > 0
#                             # Calculate approximate zero crossing
#                             pos_sdf = other_plane_sdf
#                             neg_sdf = node_plane_sdf
                            
#                             # Make sure signs are opposite (to cross the zero level)
#                             if sign(pos_sdf) == sign(neg_sdf)
#                                 continue
#                             end
                            
#                             # Calculate the plane intersection point
#                             zero_point = interpolate_to_plane(
#                                 mesh.X[other_node],
#                                 mesh.X[node_idx],
#                                 pos_sdf,
#                                 neg_sdf
#                             )
                            
#                             # Check if moving to this point would cause inversion
#                             if !would_cause_inversion(node_idx, zero_point)
#                                 # Found a suitable edge to warp along
#                                 edge_to_warp = (other_node, node_idx)
#                                 break
#                             end
#                         end
#                     end
        
#                     # If we found a suitable edge, no need to check more elements
#                     edge_to_warp !== nothing && break
#                 end
#             end
            
#             # If we found an edge to warp along, perform the warping
#             if edge_to_warp !== nothing
#                 pos_node, neg_node = edge_to_warp
                
#                 # Get plane_sdf values for interpolation
#                 pos_sdf = eval_planes_sdf(mesh, mesh.X[pos_node], plane_definitions)
#                 neg_sdf = eval_planes_sdf(mesh, mesh.X[neg_node], plane_definitions)
                
#                 # Calculate precise intersection point
#                 intersection_point = interpolate_to_plane(
#                     mesh.X[pos_node],
#                     mesh.X[neg_node],
#                     pos_sdf,
#                     neg_sdf
#                 )
                
#                 # Quantize the point to avoid duplicates
#                 p_key = quantize(intersection_point, mesh.grid_tol)
                
#                 # Check if this point already exists in the mesh
#                 if haskey(mesh.node_hash, p_key)
#                     # Get existing node index
#                     new_node_idx = mesh.node_hash[p_key]
                    
#                     # Update all elements that contain the old node
#                     for elem_idx in elements
#                         tet = mesh.IEN[elem_idx]
#                         for i in 1:length(tet)
#                             if tet[i] == neg_node
#                                 mesh.IEN[elem_idx][i] = new_node_idx
#                             end
#                         end
#                     end
#                     nodes_moved += 1
#                 else
#                     # Update the current node's position
#                     mesh.X[neg_node] = intersection_point
                    
#                     # Update the node's hash entry
#                     mesh.node_hash[p_key] = neg_node
                    
#                     nodes_moved += 1
#                 end
#             else
#                 # No suitable edge was found - all would cause inversion
#                 @warn "Node $node_idx couldn't be moved to plane surface - all potential movements would cause element inversion"
#                 nodes_skipped += 1
#             end
            
#             # Mark this node as processed
#             push!(processed_nodes, node_idx)
#         end
#     end
    
#     @info "Node warping to planes complete. $nodes_moved nodes were moved to plane surfaces. $nodes_skipped nodes were skipped to prevent element inversion."
    
#     # Odstranění elementů, které mají všechny uzly na nulové hladině nebo nemají alespoň jeden uzel s kladnou SDF
#     new_IEN = Vector{Vector{Int64}}()
#     for tet in mesh.IEN
#         # Získáme SDF hodnoty pro všechny uzly tetraedru
#         tet_sdf = [eval_planes_sdf(mesh, mesh.X[i], plane_definitions) for i in tet]
        
#         # Počet uzlů na nulové hladině
#         zero_nodes = count(x -> abs(x) < mesh.grid_tol, tet_sdf)
        
#         # Počet uzlů s kladnou SDF hodnotou
#         pos_nodes = count(x -> x > 0, tet_sdf)
        
#         # Zachováme tetraedron pouze pokud:
#         # 1. Má alespoň jeden uzel s kladnou hodnotou (vzhledem k planes_sdf)
#         # 2. Nemá všechny uzly na nulové hladině
#         if pos_nodes > 0 && zero_nodes < 4
#             push!(new_IEN, tet)
#         end
#     end
    
#     # Aktualizace konektivity
#     mesh.IEN = new_IEN
#     @info "After removing elements inside cutting planes: $(length(mesh.IEN)) tetrahedra remain"
    
#     # After warping, we need to update the mesh connectivity
#     update_connectivity!(mesh)
#   end
