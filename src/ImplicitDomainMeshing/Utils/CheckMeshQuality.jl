using LinearAlgebra
using Statistics
using StaticArrays
using CairoMakie

"""
    compute_quality_metrics!(mesh::BlockMesh)

Computes geometric quality metrics for each tetrahedron in the mesh. 
Returns statistics for aspect ratio, dihedral angles, and shape quality.
"""
function compute_quality_metrics!(mesh::BlockMesh)
    # Pre-allocate arrays for metrics
    aspect_ratios = Float64[]
    min_dihedral_angles = Float64[]
    max_dihedral_angles = Float64[]
    shape_qualities = Float64[]
    
    for tet in mesh.IEN
        # Get node coordinates for this tetrahedron
        vertices = [mesh.X[i] for i in tet]
        
        # Compute edges
        edges = [vertices[j] - vertices[i] for i in 1:3 for j in i+1:4]
        
        # Compute aspect ratio (longest edge / shortest edge)
        edge_lengths = [norm(edge) for edge in edges]
        aspect_ratio = maximum(edge_lengths) / minimum(edge_lengths)
        push!(aspect_ratios, aspect_ratio)
        
        # Compute face normals
        face_normals = [
            normalize(cross(vertices[2]-vertices[1], vertices[3]-vertices[1])),
            normalize(cross(vertices[2]-vertices[1], vertices[4]-vertices[1])),
            normalize(cross(vertices[3]-vertices[1], vertices[4]-vertices[1])),
            normalize(cross(vertices[3]-vertices[2], vertices[4]-vertices[2]))
        ]
        
        # Compute dihedral angles
        dihedral_angles = [acos(clamp(dot(face_normals[i], face_normals[j]), -1.0, 1.0)) 
                           for i in 1:3 for j in i+1:4]
        push!(min_dihedral_angles, minimum(dihedral_angles) * 180/π)
        push!(max_dihedral_angles, maximum(dihedral_angles) * 180/π)
        
        # Compute shape quality metric (ratio of inscribed to circumscribed sphere radii)
        volume = abs(dot(cross(vertices[2]-vertices[1], vertices[3]-vertices[1]), vertices[4]-vertices[1])) / 6.0
        surface_area = sum([0.5 * norm(cross(vertices[j]-vertices[i], vertices[k]-vertices[i])) 
                           for (i,j,k) in [(1,2,3), (1,2,4), (1,3,4), (2,3,4)]])
        
        # Mean ratio quality metric (normalized shape measure)
        mean_ratio = 12.0 * (3.0 * volume)^(2.0/3.0) / sum([dot(edge, edge) for edge in edges])
        push!(shape_qualities, mean_ratio)
    end
    
    # Return statistics on quality metrics
    return Dict(
        "aspect_ratio" => (minimum=minimum(aspect_ratios), 
                          maximum=maximum(aspect_ratios), 
                          mean=mean(aspect_ratios)),
        "min_dihedral_angle" => (minimum=minimum(min_dihedral_angles), 
                               maximum=maximum(min_dihedral_angles), 
                               mean=mean(min_dihedral_angles)),
        "max_dihedral_angle" => (minimum=minimum(max_dihedral_angles), 
                               maximum=maximum(max_dihedral_angles), 
                               mean=mean(max_dihedral_angles)),
        "shape_quality" => (minimum=minimum(shape_qualities), 
                          maximum=maximum(shape_qualities), 
                          mean=mean(shape_qualities))
    )
end

"""
    check_inverted_elements(mesh::BlockMesh)

Identifies inverted tetrahedral elements by checking for negative Jacobian determinants.
Returns a list of indices of inverted elements.
"""
function check_inverted_elements(mesh::BlockMesh)
    inverted_elements = Int[]
    
    for (i, tet) in enumerate(mesh.IEN)
        vertices = [mesh.X[j] for j in tet]
        
        # Compute Jacobian determinant (6 times the volume)
        jacobian_det = dot(cross(vertices[2]-vertices[1], vertices[3]-vertices[1]), vertices[4]-vertices[1])
        
        # If negative, the element is inverted
        if jacobian_det < 0
            push!(inverted_elements, i)
        end
    end
    
    return inverted_elements
end

"""
    check_boundary_integrity(mesh::BlockMesh)

Checks the topological integrity of the mesh boundary by identifying boundary faces
and verifying that the boundary is watertight (closed).
Returns a dictionary with boundary information.
"""
function check_boundary_integrity(mesh::BlockMesh)
    # Extract all triangular faces (each face appears once or twice)
    faces = Dict{Set{Int}, Vector{Int}}()
    
    for (elem_idx, tet) in enumerate(mesh.IEN)
        # For each face in the tetrahedron
        for face in [(tet[1], tet[2], tet[3]), 
                    (tet[1], tet[2], tet[4]), 
                    (tet[1], tet[3], tet[4]), 
                    (tet[2], tet[3], tet[4])]
            face_set = Set(face)
            if haskey(faces, face_set)
                push!(faces[face_set], elem_idx)
            else
                faces[face_set] = [elem_idx]
            end
        end
    end
    
    # Boundary faces appear exactly once
    boundary_faces = [face for (face, elems) in faces if length(elems) == 1]
    
    # Check if boundary is closed (watertight)
    is_closed = true
    boundary_edges = Dict{Set{Int}, Int}()
    
    for face in boundary_faces
        face_vec = collect(face)
        for i in 1:3
            edge = Set([face_vec[i], face_vec[mod1(i+1, 3)]])
            boundary_edges[edge] = get(boundary_edges, edge, 0) + 1
        end
    end
    
    # Every edge in a closed boundary must appear exactly twice
    non_manifold_edges = [edge for (edge, count) in boundary_edges if count != 2]
    
    return Dict(
        "boundary_faces" => length(boundary_faces),
        "is_closed" => isempty(non_manifold_edges),
        "non_manifold_edges" => non_manifold_edges
    )
end

"""
    check_element_orientation_consistency(mesh::BlockMesh)

Checks the orientation consistency of mesh elements by verifying that
for each interior face shared by two tetrahedra, the 4th vertices of both
tetrahedra lie on opposite sides of the face.

Returns a list of inconsistently oriented element pairs.
"""
function check_element_orientation_consistency(mesh::BlockMesh)
    # Create a dictionary mapping faces to their parent tetrahedra
    faces_to_tets = Dict{Set{Int}, Vector{Int}}()
    
    # Populate the dictionary
    for (tet_idx, tet) in enumerate(mesh.IEN)
        # For each face in the tetrahedron
        for face in [(tet[1], tet[2], tet[3]), 
                     (tet[1], tet[2], tet[4]), 
                     (tet[1], tet[3], tet[4]), 
                     (tet[2], tet[3], tet[4])]
            face_set = Set(face)
            if haskey(faces_to_tets, face_set)
                push!(faces_to_tets[face_set], tet_idx)
            else
                faces_to_tets[face_set] = [tet_idx]
            end
        end
    end
    
    # List to store inconsistent orientations
    inconsistent_orientations = Vector{Tuple{Int, Int, Set{Int}}}()
    
    # Check orientation consistency for each interior face
    for (face_set, tet_indices) in faces_to_tets
        # Skip boundary faces (belong to only one tetrahedron)
        if length(tet_indices) != 2
            continue
        end
        
        # Get the two tetrahedra sharing this face
        tet1_idx, tet2_idx = tet_indices
        tet1 = mesh.IEN[tet1_idx]
        tet2 = mesh.IEN[tet2_idx]
        
        # Identify vertices of the face
        face_verts = collect(face_set)
        
        # Find the vertices that are not part of the face
        tet1_opposite_vertex = setdiff(tet1, face_verts)[1]
        tet2_opposite_vertex = setdiff(tet2, face_verts)[1]
        
        # Get coordinates of vertices
        face_coords = [mesh.X[v] for v in face_verts]
        tet1_opposite_coords = mesh.X[tet1_opposite_vertex]
        tet2_opposite_coords = mesh.X[tet2_opposite_vertex]
        
        # Calculate face normal (using first three vertices)
        edge1 = face_coords[2] - face_coords[1]
        edge2 = face_coords[3] - face_coords[1]
        face_normal = normalize(cross(edge1, edge2))
        
        # Vectors from a face vertex to opposite vertices
        vec_to_tet1_opposite = tet1_opposite_coords - face_coords[1]
        vec_to_tet2_opposite = tet2_opposite_coords - face_coords[1]
        
        # Check if the opposite vertices are on opposite sides of the face
        dot_product1 = dot(face_normal, vec_to_tet1_opposite)
        dot_product2 = dot(face_normal, vec_to_tet2_opposite)
        
        # If dot products have the same sign, the orientation is inconsistent
        if dot_product1 * dot_product2 > 0
            push!(inconsistent_orientations, (tet1_idx, tet2_idx, face_set))
        end
    end
    
    return inconsistent_orientations
end

"""
    verify_mesh_intersections(mesh::BlockMesh)

Performs a comprehensive check for element intersections using a three-phase approach:
1. Topological phase: Build adjacency graphs
2. Spatial phase: Filter potential collisions using AABBs
3. Geometric phase: Precise intersection testing

Returns detailed statistics and lists of problematic elements.
"""
function verify_mesh_intersections(mesh::BlockMesh)
    @info "Verifying mesh integrity - checking for element intersections..."
    
    # Output statistics
    stats = Dict(
        "tested_pairs" => 0,
        "intersecting_pairs" => 0,
        "invalid_topology" => 0,
        "time_topology" => 0.0,
        "time_spatial" => 0.0,
        "time_geometric" => 0.0
    )
    
    # 1. Topological phase - build adjacency graphs
    time_start = time()
    
    # Create dictionary of edges - key: pair of vertices, value: list of elements
    edges_to_tets = Dict{Set{Int}, Vector{Int}}()
    # Create dictionary of faces - key: triplet of vertices, value: list of elements
    faces_to_tets = Dict{Set{Int}, Vector{Int}}()
    # Create set of nodes for each element
    tet_vertices = [Set(tet) for tet in mesh.IEN]
    
    for (tet_idx, tet) in enumerate(mesh.IEN)
        # Register all edges of the tetrahedron
        for i in 1:4
            for j in (i+1):4
                edge = Set([tet[i], tet[j]])
                if haskey(edges_to_tets, edge)
                    push!(edges_to_tets[edge], tet_idx)
                else
                    edges_to_tets[edge] = [tet_idx]
                end
            end
        end
        
        # Register all faces of the tetrahedron
        faces = [
            Set([tet[1], tet[2], tet[3]]),
            Set([tet[1], tet[2], tet[4]]),
            Set([tet[1], tet[3], tet[4]]),
            Set([tet[2], tet[3], tet[4]])
        ]
        
        for face in faces
            if haskey(faces_to_tets, face)
                push!(faces_to_tets[face], tet_idx)
            else
                faces_to_tets[face] = [tet_idx]
            end
        end
    end
    
    # Check if any face belongs to more than two tetrahedra (topological invalidity)
    invalid_faces = filter(pair -> length(pair.second) > 2, faces_to_tets)
    
    if !isempty(invalid_faces)
        stats["invalid_topology"] = length(invalid_faces)
        @warn "Found $(length(invalid_faces)) faces shared by more than two tetrahedra!"
    end
    
    # Create list of elements that share any vertex (potential intersections)
    # For each element, store a set of elements that share at least one vertex
    vertex_sharing_elements = [Set{Int}() for _ in 1:length(mesh.IEN)]
    
    # For each vertex, store a list of elements that contain it
    vertex_to_tets = [Vector{Int}() for _ in 1:length(mesh.X)]
    for (tet_idx, tet) in enumerate(mesh.IEN)
        for v in tet
            push!(vertex_to_tets[v], tet_idx)
        end
    end
    
    # For each element, find elements that share vertices with it
    for (tet_idx, tet) in enumerate(mesh.IEN)
        for v in tet
            for other_tet in vertex_to_tets[v]
                if other_tet != tet_idx
                    push!(vertex_sharing_elements[tet_idx], other_tet)
                end
            end
        end
    end
    
    # Create adjacency graph (elements sharing a face)
    neighbors = [Set{Int}() for _ in 1:length(mesh.IEN)]
    
    for (_, tets) in faces_to_tets
        if length(tets) == 2
            push!(neighbors[tets[1]], tets[2])
            push!(neighbors[tets[2]], tets[1])
        end
    end
    
    stats["time_topology"] = time() - time_start
    
    # 2. Spatial phase - filtering candidates for detailed tests
    time_start = time()
    
    # Calculate AABB (axis-aligned bounding box) for each tetrahedron
    tet_aabbs = Vector{Tuple{SVector{3,Float64}, SVector{3,Float64}}}(undef, length(mesh.IEN))
    
    for (tet_idx, tet) in enumerate(mesh.IEN)
        verts = [mesh.X[v] for v in tet]
        min_point = SVector{3,Float64}(
            minimum(v[1] for v in verts),
            minimum(v[2] for v in verts),
            minimum(v[3] for v in verts)
        )
        max_point = SVector{3,Float64}(
            maximum(v[1] for v in verts),
            maximum(v[2] for v in verts),
            maximum(v[3] for v in verts)
        )
        tet_aabbs[tet_idx] = (min_point, max_point)
    end
    
    # Function to check intersection of two AABBs
    function aabb_intersect(box1, box2)
        # Use smaller tolerance to eliminate false positives
        tol = mesh.grid_tol * 0.1
        return all(box1[1][i] - tol <= box2[2][i] && box1[2][i] + tol >= box2[1][i] for i in 1:3)
    end
    
    # Create list of potential collisions
    # Tetrahedra must: 
    # 1. Not be neighbors (not share a face)
    # 2. Their AABBs must intersect
    # 3. Must not share an edge or vertex (legitimate for adjacent tetrahedra)
    potential_collisions = Vector{Tuple{Int,Int}}()
    
    for i in 1:length(mesh.IEN)
        # Check only with larger indices (to avoid duplicates)
        for j in (i+1):length(mesh.IEN)
            # Skip neighbors (sharing a face)
            if j in neighbors[i]
                continue
            end
            
            # Skip elements sharing any vertex (which is acceptable)
            if j in vertex_sharing_elements[i]
                continue
            end
            
            # Check for AABB intersection
            if aabb_intersect(tet_aabbs[i], tet_aabbs[j])
                push!(potential_collisions, (i, j))
            end
            
            stats["tested_pairs"] += 1
        end
    end
    
    stats["time_spatial"] = time() - time_start
    
    # 3. Geometric phase - precise intersection checks
    time_start = time()
    
    # List of actual intersections
    intersections = Vector{Tuple{Int,Int}}()
    
    # Function to compute tetrahedron volume
    function tetrahedron_volume(vertices)
        return abs(dot(cross(vertices[2]-vertices[1], vertices[3]-vertices[1]), vertices[4]-vertices[1])) / 6.0
    end
    
    # Improved function to test if a point is inside a tetrahedron
    function robust_point_in_tetrahedron(point, vertices)
        # Using more robust method based on sign of subtetrahedra volumes
        # This approach is less susceptible to numerical errors
        
        # Volume of the entire tetrahedron
        vol_total = tetrahedron_volume(vertices)
        
        # If total volume is almost zero, tetrahedron is degenerate
        if vol_total < 1e-12
            return false
        end
        
        # Volume of tetrahedra created by the point and three vertices of the original tetrahedron
        sub_vols = [
            tetrahedron_volume([point, vertices[2], vertices[3], vertices[4]]),
            tetrahedron_volume([vertices[1], point, vertices[3], vertices[4]]),
            tetrahedron_volume([vertices[1], vertices[2], point, vertices[4]]),
            tetrahedron_volume([vertices[1], vertices[2], vertices[3], point])
        ]
        
        # Point is inside the tetrahedron if sum of subtetrahedra volumes equals volume of entire tetrahedron
        # Use relative test due to numerical errors
        sum_vol = sum(sub_vols)
        is_inside = abs(sum_vol - vol_total) < 1e-9 * vol_total
        
        return is_inside
    end
    
    # Check for actual geometric intersections between non-adjacent elements
    for (i, j) in potential_collisions
        tet_i = [mesh.X[v] for v in mesh.IEN[i]]
        tet_j = [mesh.X[v] for v in mesh.IEN[j]]
        
        # Check if any vertex of one tetrahedron lies strictly inside the other
        intersection_found = false
        
        # Check vertices of first tetrahedron against second
        for v in tet_i
            if robust_point_in_tetrahedron(v, tet_j)
                intersection_found = true
                break
            end
        end
        
        # Check vertices of second tetrahedron against first
        if !intersection_found
            for v in tet_j
                if robust_point_in_tetrahedron(v, tet_i)
                    intersection_found = true
                    break
                end
            end
        end
        
        # If intersection found, record it
        if intersection_found
            push!(intersections, (i, j))
        end
    end
    
    stats["intersecting_pairs"] = length(intersections)
    stats["time_geometric"] = time() - time_start
    
    # Results of check
    if isempty(intersections) && stats["invalid_topology"] == 0
        @info "Mesh integrity check completed: No element intersections found."
    else
        @warn "Found $(length(intersections)) pairs of intersecting tetrahedra!"
        @warn "Found $(stats["invalid_topology"]) topological invalidities!"
    end
    
    # Prepare output with intersection details for visualization
    detailed_result = Dict(
        "stats" => stats,
        "intersections" => intersections,
        "invalid_faces" => collect(keys(invalid_faces))
    )
    
    return detailed_result
end

"""
    compute_all_dihedral_angles(mesh::BlockMesh)

Computes all interior dihedral angles for every tetrahedron in the mesh.
Returns a vector containing all angles in degrees.
The interior dihedral angle is the angle between two adjacent faces measured 
inside the tetrahedron (ideal angle in a regular tetrahedron is ~70.53°).
"""
function compute_all_dihedral_angles(mesh::BlockMesh)
    all_angles = Float64[]
    
    for tet in mesh.IEN
        # Get node coordinates for this tetrahedron
        vertices = [mesh.X[i] for i in tet]
        
        # Compute face normals (outward facing)
        face_normals = [
            normalize(cross(vertices[2]-vertices[1], vertices[3]-vertices[1])),
            normalize(cross(vertices[2]-vertices[1], vertices[4]-vertices[1])),
            normalize(cross(vertices[3]-vertices[1], vertices[4]-vertices[1])),
            normalize(cross(vertices[3]-vertices[2], vertices[4]-vertices[2]))
        ]
        
        # For a tetrahedron, we need to ensure normals are consistently outward-facing
        # We'll check if the normal points away from the opposite vertex
        opposite_vertices = [
            vertices[4],  # opposite to face 1-2-3
            vertices[3],  # opposite to face 1-2-4
            vertices[2],  # opposite to face 1-3-4
            vertices[1]   # opposite to face 2-3-4
        ]
        
        # Ensure normals point outward (away from the tetrahedron)
        for i in 1:4
            face_center = (vertices[setdiff(1:4, [i])[1]] + 
                         vertices[setdiff(1:4, [i])[2]] + 
                         vertices[setdiff(1:4, [i])[3]]) / 3
            
            vector_to_opposite = opposite_vertices[i] - face_center
            
            # If normal points inward, flip it
            if dot(face_normals[i], vector_to_opposite) > 0
                face_normals[i] = -face_normals[i]
            end
        end
        
        # Compute interior dihedral angles (the angle between faces inside the tetrahedron)
        # For interior angles, we use (π - angle between outward normals)
        # Equivalently: acos(-dot(normal1, normal2))
        tet_angles = [acos(-clamp(dot(face_normals[i], face_normals[j]), -1.0, 1.0)) * 180/π 
                      for i in 1:3 for j in i+1:4]
        
        append!(all_angles, tet_angles)
    end
    
    return all_angles
end

"""
    plot_dihedral_angles_histogram(angles::Vector{Float64}, output_path::String)

Creates a histogram of dihedral angles using the Makie library.
Highlights critical angle thresholds and saves the plot to the specified output path.
"""
function plot_dihedral_angles_histogram(angles::Vector{Float64}, output_path::String)
    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1], 
              xlabel = "Interior Dihedral Angle (degrees)",
              ylabel = "Frequency",
              title = "Histogram of Interior Dihedral Angles")
    
    # Create histogram
    hist!(ax, angles, bins=50, color=:skyblue, strokewidth=1, strokecolor=:black)
    
    # Add vertical lines for critical and ideal angles
    vlines!(ax, [10], color=:red, linestyle=:dash, linewidth=2, label="10° (critical low)")
    vlines!(ax, [30], color=:orange, linestyle=:dash, linewidth=2, label="30° (warning low)")
    vlines!(ax, [70.53], color=:green, linestyle=:solid, linewidth=2, label="70.53° (ideal)")
    vlines!(ax, [150], color=:orange, linestyle=:dash, linewidth=2, label="150° (warning high)")
    vlines!(ax, [170], color=:red, linestyle=:dash, linewidth=2, label="170° (critical high)")
    
    # Add statistics
    angle_stats = [
        "Min: $(round(minimum(angles), digits=2))°",
        "Max: $(round(maximum(angles), digits=2))°",
        "Mean: $(round(mean(angles), digits=2))°",
        "Median: $(round(median(angles), digits=2))°"
    ]
    
    critical_low_count = count(θ -> θ < 10, angles)
    warning_low_count = count(θ -> θ >= 10 && θ < 30, angles)
    warning_high_count = count(θ -> θ > 150 && θ <= 170, angles)
    critical_high_count = count(θ -> θ > 170, angles)
    
    push!(angle_stats, "Critical low angles (<10°): $critical_low_count")
    push!(angle_stats, "Warning low angles (10°-30°): $warning_low_count")
    push!(angle_stats, "Warning high angles (150°-170°): $warning_high_count")
    push!(angle_stats, "Critical high angles (>170°): $critical_high_count")
    
    # Add statistics text box
    text!(ax, 90, 0.9*maximum(ax.finallimits.val.origin[2] + ax.finallimits.val.widths[2]), 
          text=join(angle_stats, "\n"),
          align=(:center, :top),
          fontsize=12,
          color=:black,
          space=:data)
    
    axislegend(ax, position=:rt)
    
    # Save the figure
    save(output_path, fig)
    
    return fig
end

"""
    assess_mesh_quality(mesh::BlockMesh, output_prefix::String)

Main function that assesses mesh quality using multiple checks:
- Geometric quality metrics
- Element validity (inverted elements)
- Boundary integrity
- Element orientation consistency
- Element intersections
- Dihedral angle distribution (with histogram visualization)

Returns a comprehensive report and detailed results for further analysis.
"""
function assess_mesh_quality(mesh::BlockMesh, output_prefix::String)
    # Calculate geometric quality metrics
    quality_metrics = compute_quality_metrics!(mesh)
    
    # Check for inverted elements
    inverted = check_inverted_elements(mesh)
    
    # Check topological integrity
    boundary_check = check_boundary_integrity(mesh)
    
    # Check element orientation consistency
    orientation_check = check_element_orientation_consistency(mesh)
    
    # Check for element intersections
    intersection_check = verify_mesh_intersections(mesh)
    
    # Compute all dihedral angles and plot histogram
    all_angles = compute_all_dihedral_angles(mesh)
    plot_dihedral_angles_histogram(all_angles, "$(output_prefix)_dihedral_angles.png")
    
    # Calculate angle statistics
    critical_low_count = count(θ -> θ < 10, all_angles)
    warning_low_count = count(θ -> θ >= 10 && θ < 30, all_angles)
    warning_high_count = count(θ -> θ > 150 && θ <= 170, all_angles)
    critical_high_count = count(θ -> θ > 170, all_angles)
    
    total_critical = critical_low_count + critical_high_count
    total_warning = warning_low_count + warning_high_count
    
    percent_critical = round(100 * total_critical / length(all_angles), digits=2)
    percent_warning = round(100 * total_warning / length(all_angles), digits=2)
    
    # Calculate deviation from ideal angle
    ideal_angle = 70.53 # degrees
    angles_near_ideal = count(θ -> abs(θ - ideal_angle) < 5, all_angles)
    percent_near_ideal = round(100 * angles_near_ideal / length(all_angles), digits=2)
    
    # Print summary report
    println("Mesh Quality Assessment Report")
    println("=============================")
    println("Total elements: $(length(mesh.IEN))")
    println("Total nodes: $(length(mesh.X))")
    println()
    
    println("Geometric Quality:")
    println("  Aspect ratio: range [$(quality_metrics["aspect_ratio"].minimum), $(quality_metrics["aspect_ratio"].maximum)], mean $(quality_metrics["aspect_ratio"].mean)")
    println("  Min dihedral angle: range [$(quality_metrics["min_dihedral_angle"].minimum), $(quality_metrics["min_dihedral_angle"].maximum)], mean $(quality_metrics["min_dihedral_angle"].mean)")
    println("  Shape quality: range [$(quality_metrics["shape_quality"].minimum), $(quality_metrics["shape_quality"].maximum)], mean $(quality_metrics["shape_quality"].mean)")
    println()
    
    println("Dihedral Angle Analysis:")
    println("  Total angles: $(length(all_angles))")
    println("  Critical low angles (<10°): $critical_low_count")
    println("  Warning low angles (10°-30°): $warning_low_count")
    println("  Warning high angles (150°-170°): $warning_high_count")
    println("  Critical high angles (>170°): $critical_high_count")
    println("  Total critical angles: $total_critical ($percent_critical%)")
    println("  Total warning angles: $total_warning ($percent_warning%)")
    println("  Angles near ideal (65.53°-75.53°): $angles_near_ideal ($percent_near_ideal%)")
    println("  Histogram saved to: $(output_prefix)_dihedral_angles.png")
    println()
    
    println("Element Validity:")
    println("  Inverted elements: $(length(inverted))")
    if !isempty(inverted)
        println("  First 5 inverted elements: $(inverted[1:min(5, length(inverted))])")
    end
    println()
    
    println("Topological Integrity:")
    println("  Boundary faces: $(boundary_check["boundary_faces"])")
    println("  Watertight boundary: $(boundary_check["is_closed"])")
    println("  Elements with inconsistent orientation: $(length(orientation_check))")
    println("  Intersecting element pairs: $(intersection_check["stats"]["intersecting_pairs"])")
    println()
    
    # Return comprehensive results
    return Dict(
        "quality_metrics" => quality_metrics,
        "inverted_elements" => inverted,
        "boundary_check" => boundary_check,
        "orientation_check" => orientation_check,
        "intersection_check" => intersection_check,
        "dihedral_angles" => all_angles,
        "angle_statistics" => Dict(
            "critical_low_count" => critical_low_count,
            "warning_low_count" => warning_low_count,
            "warning_high_count" => warning_high_count,
            "critical_high_count" => critical_high_count,
            "total_critical" => total_critical,
            "total_warning" => total_warning,
            "percent_critical" => percent_critical,
            "percent_warning" => percent_warning,
            "angles_near_ideal" => angles_near_ideal,
            "percent_near_ideal" => percent_near_ideal
        )
    )
end
