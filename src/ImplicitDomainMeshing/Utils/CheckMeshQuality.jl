using Statistics

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

function check_non_manifold_vertices(mesh::BlockMesh)
    # For each vertex, get all boundary faces containing it
    vertex_boundary_faces = [Set() for _ in 1:length(mesh.X)]
    
    # First, identify boundary faces
    faces = Dict{Set{Int}, Vector{Int}}()
    for (elem_idx, tet) in enumerate(mesh.IEN)
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
    
    # For each boundary face, add it to the corresponding vertex sets
    for face in boundary_faces
        for vertex in face
            push!(vertex_boundary_faces[vertex], face)
        end
    end
    
    # Check if the boundary faces at each vertex form a single, connected ring
    non_manifold_vertices = Int[]
    
    for (idx, faces) in enumerate(vertex_boundary_faces)
        isempty(faces) && continue  # Skip interior vertices
        
        # A vertex is non-manifold if its star is not a topological disk
        # This is a simplified check - a complete check needs more complex topology analysis
        edges = Set()
        for face in faces
            face_vec = collect(face)
            for i in 1:3
                if face_vec[i] == idx
                    edge = Set([face_vec[mod1(i+1, 3)], face_vec[mod1(i+2, 3)]])
                    if edge in edges
                        edges = setdiff(edges, [edge])
                    else
                        push!(edges, edge)
                    end
                end
            end
        end
        
        if !isempty(edges)
            push!(non_manifold_vertices, idx)
        end
    end
    
    return non_manifold_vertices
end

function export_quality_metrics_vtk(mesh::BlockMesh, filename::String)
    # Compute per-element quality metrics
    aspect_ratios = Float64[]
    dihedral_angles = Float64[]
    shape_qualities = Float64[]
    
    for tet in mesh.IEN
        # Implement quality metrics as in compute_quality_metrics function
        # ... (calculation code)
        
        push!(aspect_ratios, aspect_ratio)
        push!(dihedral_angles, minimum(dihedral_angles))
        push!(shape_qualities, mean_ratio)
    end
    
    # Export the mesh with quality metrics as cell data
    npoints = length(mesh.X)
    points = zeros(Float64, 3, npoints)
    for i in 1:npoints
        points[:, i] = mesh.X[i]
    end
    
    cells = [MeshCell(VTKCellTypes.VTK_TETRA, tet) for tet in mesh.IEN]
    vtkfile = vtk_grid(filename, points, cells)
    
    # Add quality metrics as cell data
    vtk_cell_data(vtkfile, aspect_ratios, "aspect_ratio")
    vtk_cell_data(vtkfile, dihedral_angles, "min_dihedral_angle")
    vtk_cell_data(vtkfile, shape_qualities, "shape_quality")
    
    # Add SDF values as point data
    vtk_point_data(vtkfile, mesh.node_sdf, "sdf")
    
    vtk_save(vtkfile)
end


function assess_mesh_quality(mesh::BlockMesh, output_prefix::String)
    # Calculate geometric quality metrics
    quality_metrics = compute_quality_metrics!(mesh)
    
    # Check for inverted elements
    inverted = check_inverted_elements(mesh)
    
    # Check topological integrity
    boundary_check = check_boundary_integrity(mesh)
    non_manifold_check = check_non_manifold_vertices(mesh)
    
    # Export mesh with quality metrics for visualization
    # export_quality_metrics_vtk(mesh, "$(output_prefix)_quality.vtu")
    
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
    
    println("Element Validity:")
    println("  Inverted elements: $(length(inverted))")
    if !isempty(inverted)
        println("  First 5 inverted elements: $(inverted[1:min(5, length(inverted))])")
    end
    println()
    
    println("Topological Integrity:")
    println("  Boundary faces: $(boundary_check["boundary_faces"])")
    println("  Watertight boundary: $(boundary_check["is_closed"])")
    println("  Non-manifold vertices: $(length(non_manifold_check))")
    println()
    
    # Return comprehensive results
    return Dict(
        "quality_metrics" => quality_metrics,
        "inverted_elements" => inverted,
        "boundary_check" => boundary_check,
        "non_manifold_vertices" => non_manifold_check
    )
end

assess_mesh_quality(mesh, "IS-simple")


#_____________________________________

function edge_triangle_intersection(edge_start, edge_end, tri_verts)
    # Směrový vektor hrany
    edge_dir = edge_end - edge_start
    
    # Normála trojúhelníku
    tri_normal = normalize(cross(tri_verts[2] - tri_verts[1], tri_verts[3] - tri_verts[1]))
    
    # Kontrola, zda hrana protíná rovinu trojúhelníku
    denom = dot(tri_normal, edge_dir)
    
    # Hrany rovnoběžné s rovinou ignorujeme
    if abs(denom) < 1e-10
        return false
    end
    
    # Parametr průsečíku hrany s rovinou
    t = dot(tri_normal, tri_verts[1] - edge_start) / denom
    
    # Průsečík musí ležet na hraně
    if t < 0 || t > 1
        return false
    end
    
    # Průsečík s rovinou
    intersection = edge_start + t * edge_dir
    
    # Kontrola, zda průsečík leží uvnitř trojúhelníku (barycentrické souřadnice)
    v0 = tri_verts[3] - tri_verts[1]
    v1 = tri_verts[2] - tri_verts[1]
    v2 = intersection - tri_verts[1]
    
    # Skalární součiny pro barycentrické souřadnice
    dot00 = dot(v0, v0)
    dot01 = dot(v0, v1)
    dot02 = dot(v0, v2)
    dot11 = dot(v1, v1)
    dot12 = dot(v1, v2)
    
    # Barycentrické souřadnice
    inv_denom = 1.0 / (dot00 * dot11 - dot01 * dot01)
    u = (dot11 * dot02 - dot01 * dot12) * inv_denom
    v = (dot00 * dot12 - dot01 * dot02) * inv_denom
    
    # Test, zda je průsečík uvnitř trojúhelníku
    return (u >= 0) && (v >= 0) && (u + v <= 1)
end

function verify_mesh_intersections(mesh::BlockMesh)
    @info "Verifikace integrity sítě - kontrola průniků elementů..."
    
    # Výstupní statistiky
    stats = Dict(
        "tested_pairs" => 0,
        "intersecting_pairs" => 0,
        "invalid_topology" => 0,
        "time_topology" => 0.0,
        "time_spatial" => 0.0,
        "time_geometric" => 0.0
    )
    
    # 1. Topologická fáze - sestavení grafu sousednosti
    time_start = time()
    
    # Vytvoření slovníku hran - klíč: dvojice vrcholů, hodnota: seznam elementů
    edges_to_tets = Dict{Set{Int}, Vector{Int}}()
    # Vytvoření slovníku stěn - klíč: trojice vrcholů, hodnota: seznam elementů
    faces_to_tets = Dict{Set{Int}, Vector{Int}}()
    # Vytvoření množiny uzlů pro každý element
    tet_vertices = [Set(tet) for tet in mesh.IEN]
    
    for (tet_idx, tet) in enumerate(mesh.IEN)
        # Registrace všech hran tetraedru
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
        
        # Registrace všech stěn tetraedru
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
    
    # Kontrola, zda nějaká stěna patří více než dvěma tetraedrům (topologická nevalidita)
    invalid_faces = filter(pair -> length(pair.second) > 2, faces_to_tets)
    
    if !isempty(invalid_faces)
        stats["invalid_topology"] = length(invalid_faces)
        @warn "Nalezeno $(length(invalid_faces)) stěn sdílených více než dvěma tetraedry!"
    end
    
    # Vytvoření seznamu elementů, které sdílí nějaký vrchol (potenciálně by se mohly protínat)
    # Pro každý element si uložíme množinu elementů, které s ním sdílí alespoň jeden vrchol
    vertex_sharing_elements = [Set{Int}() for _ in 1:length(mesh.IEN)]
    
    # Pro každý vrchol si uložíme seznam elementů, které ho obsahují
    vertex_to_tets = [Vector{Int}() for _ in 1:length(mesh.X)]
    for (tet_idx, tet) in enumerate(mesh.IEN)
        for v in tet
            push!(vertex_to_tets[v], tet_idx)
        end
    end
    
    # Pro každý element najdeme elementy, které s ním sdílí vrcholy
    for (tet_idx, tet) in enumerate(mesh.IEN)
        for v in tet
            for other_tet in vertex_to_tets[v]
                if other_tet != tet_idx
                    push!(vertex_sharing_elements[tet_idx], other_tet)
                end
            end
        end
    end
    
    # Vytvoření grafu sousednosti (elementy sdílející stěnu)
    neighbors = [Set{Int}() for _ in 1:length(mesh.IEN)]
    
    for (_, tets) in faces_to_tets
        if length(tets) == 2
            push!(neighbors[tets[1]], tets[2])
            push!(neighbors[tets[2]], tets[1])
        end
    end
    
    stats["time_topology"] = time() - time_start
    
    # 2. Prostorová fáze - filtrování kandidátů pro detailní testy
    time_start = time()
    
    # Výpočet AABB (axis-aligned bounding box) pro každý tetraedr
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
    
    # Funkce pro kontrolu průniku dvou AABB
    function aabb_intersect(box1, box2)
        # Použijeme menší toleranci, abychom eliminovali false positives
        tol = mesh.grid_tol * 0.1
        return all(box1[1][i] - tol <= box2[2][i] && box1[2][i] + tol >= box2[1][i] for i in 1:3)
    end
    
    # Vytvoření seznamu potenciálních kolizí
    # Tetraedry musí: 
    # 1. Nebýt sousedy (nesdílet stěnu)
    # 2. Jejich AABB se musí protínat
    # 3. Nesmí sdílet hranu nebo vrchol (to je legitimní pro sousedící tetraedry)
    potential_collisions = Vector{Tuple{Int,Int}}()
    
    for i in 1:length(mesh.IEN)
        # Kontrolujeme pouze s většími indexy (vyhnutí se duplicitám)
        for j in (i+1):length(mesh.IEN)
            # Přeskočím sousedy (sdílející stěnu)
            if j in neighbors[i]
                continue
            end
            
            # Přeskočím elementy sdílející jakýkoliv vrchol (to je v pořádku)
            if j in vertex_sharing_elements[i]
                continue
            end
            
            # Kontrola průniku AABB
            if aabb_intersect(tet_aabbs[i], tet_aabbs[j])
                push!(potential_collisions, (i, j))
            end
            
            stats["tested_pairs"] += 1
        end
    end
    
    stats["time_spatial"] = time() - time_start
    
    # 3. Geometrická fáze - přesná kontrola průniků
    time_start = time()
    
    # Seznam skutečných průniků
    intersections = Vector{Tuple{Int,Int}}()
    
    # Funkce pro výpočet objemu tetraedru
    function tetrahedron_volume(vertices)
        return abs(dot(cross(vertices[2]-vertices[1], vertices[3]-vertices[1]), vertices[4]-vertices[1])) / 6.0
    end
    
    # Vylepšená funkce pro test, zda bod leží uvnitř tetraedru
    function robust_point_in_tetrahedron(point, vertices)
        # Použití robustnější metody založené na znaménku objemů podtetraedrů
        # Tento přístup je méně náchylný na numerické chyby
        
        # Objem celého tetraedru
        vol_total = tetrahedron_volume(vertices)
        
        # Je-li celkový objem téměř nulový, tetraedr je degenerovaný
        if vol_total < 1e-12
            return false
        end
        
        # Objem tetraedrů vytvořených bodem a třemi vrcholy původního tetraedru
        sub_vols = [
            tetrahedron_volume([point, vertices[2], vertices[3], vertices[4]]),
            tetrahedron_volume([vertices[1], point, vertices[3], vertices[4]]),
            tetrahedron_volume([vertices[1], vertices[2], point, vertices[4]]),
            tetrahedron_volume([vertices[1], vertices[2], vertices[3], point])
        ]
        
        # Bod je uvnitř tetraedru, pokud je součet objemů podtetraedrů roven objemu celého tetraedru
        # Použijeme relativní test kvůli numerickým chybám
        sum_vol = sum(sub_vols)
        is_inside = abs(sum_vol - vol_total) < 1e-9 * vol_total
        
        return is_inside
    end
    
    # Kontrola skutečných geometrických průniků mezi nesousedícími elementy
    for (i, j) in potential_collisions
        tet_i = [mesh.X[v] for v in mesh.IEN[i]]
        tet_j = [mesh.X[v] for v in mesh.IEN[j]]
        
        # Kontrola, zda některý vrchol jednoho tetraedru leží striktně uvnitř druhého
        intersection_found = false
        
        # Kontrola vrcholů prvního tetraedru vůči druhému
        for v in tet_i
            if robust_point_in_tetrahedron(v, tet_j)
                intersection_found = true
                break
            end
        end
        
        # Kontrola vrcholů druhého tetraedru vůči prvnímu
        if !intersection_found
            for v in tet_j
                if robust_point_in_tetrahedron(v, tet_i)
                    intersection_found = true
                    break
                end
            end
        end
        
        # Pokud byl nalezen průnik, zaznamenáme ho
        if intersection_found
            push!(intersections, (i, j))
        end
    end
    
    stats["intersecting_pairs"] = length(intersections)
    stats["time_geometric"] = time() - time_start
    
    # Výsledky kontroly
    if isempty(intersections) && stats["invalid_topology"] == 0
        @info "Kontrola integrity sítě dokončena: Žádné průniky elementů nebyly nalezeny."
    else
        @warn "Nalezeno $(length(intersections)) párů protínajících se tetraedrů!"
        @warn "Nalezeno $(stats["invalid_topology"]) topologických nevalidit!"
    end
    
    # Příprava výstupu s detaily průniků pro vizualizaci
    detailed_result = Dict(
        "stats" => stats,
        "intersections" => intersections,
        "invalid_faces" => collect(keys(invalid_faces))
    )
    
    return detailed_result
end

intersection_data = verify_mesh_intersections(mesh)

function export_intersection_visualization(mesh::BlockMesh, intersection_data, filename::String)
    # Extrakce problematických elementů
    problem_elements = Set{Int}()
    
    # Přidání protínajících se tetraedrů
    for (i, j) in intersection_data["intersections"]
        push!(problem_elements, i)
        push!(problem_elements, j)
    end
    
    # Přidání tetraedrů s nevalidní topologií
    for face in intersection_data["invalid_faces"]
        face_tets = faces_to_tets[face]
        for tet_idx in face_tets
            push!(problem_elements, tet_idx)
        end
    end
    
    # Vytvoření elementového pole pro vizualizaci (1 = problém, 0 = OK)
    element_status = zeros(Int, length(mesh.IEN))
    for elem_idx in problem_elements
        element_status[elem_idx] = 1
    end
    
    # Export do VTK
    npoints = length(mesh.X)
    points = zeros(Float64, 3, npoints)
    for i in 1:npoints
        points[:, i] = mesh.X[i]
    end
    
    cells = [MeshCell(VTKCellTypes.VTK_TETRA, tet) for tet in mesh.IEN]
    vtkfile = vtk_grid(filename, points, cells)
    
    # Přidání statusu elementů
    vtk_cell_data(vtkfile, element_status, "problem_element")
    
    # Přidání SDF hodnot
    vtk_point_data(vtkfile, mesh.node_sdf, "sdf")
    
    vtk_save(vtkfile)
    
    @info "Vizualizace průniků exportována do: $filename"
end

export_intersection_visualization(mesh, intersection_data, "průniky")