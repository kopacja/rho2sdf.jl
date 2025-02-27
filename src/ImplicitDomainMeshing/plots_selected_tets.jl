function export_specific_tetrahedra(mesh::BlockMesh, tet_ids::Vector{Int}, filename::String)
    # Ověření vstupních hodnot
    if isempty(tet_ids)
        @warn "Nebyly zadány žádné tetrahedry k exportu."
        return
    end
    
    # Ověření rozsahu ID tetraedrů
    max_tet_id = length(mesh.IEN)
    valid_tet_ids = filter(id -> 1 <= id <= max_tet_id, tet_ids)
    
    if length(valid_tet_ids) < length(tet_ids)
        @warn "Některá ID tetraedrů byla mimo rozsah (1-$max_tet_id)."
    end
    
    if isempty(valid_tet_ids)
        @error "Žádné z požadovaných ID tetraedrů není platné."
        return
    end
    
    # Shromáždění uzlů použitých ve vybraných tetraedrech
    used_nodes = Set{Int}()
    for tet_id in valid_tet_ids
        union!(used_nodes, mesh.IEN[tet_id])
    end
    
    # Vytvoření mapování ze starých ID uzlů na nová
    node_map = Dict{Int,Int}()
    for (new_id, old_id) in enumerate(sort(collect(used_nodes)))
        node_map[old_id] = new_id
    end
    
    # Připravení dat pro VTK export
    n_nodes = length(used_nodes)
    points = zeros(Float64, 3, n_nodes)
    
    # Vytvoření pole bodů
    for (old_id, new_id) in node_map
        points[:, new_id] = mesh.X[old_id]
    end
    
    # Vytvoření buněk s přemapovanými indexy uzlů
    cells = MeshCell[]
    for tet_id in valid_tet_ids
        new_tet = [node_map[node_id] for node_id in mesh.IEN[tet_id]]
        push!(cells, MeshCell(VTKCellTypes.VTK_TETRA, new_tet))
    end
    
    # Extrakce SDF hodnot pro exportované uzly
    sdf_values = zeros(Float64, n_nodes)
    for (old_id, new_id) in node_map
        sdf_values[new_id] = mesh.node_sdf[old_id]
    end
    
    # Vytvoř VTK soubor
    vtkfile = vtk_grid(filename, points, cells)
    
    # Přidej SDF data
    vtk_point_data(vtkfile, sdf_values, "sdf")
    
    # Přidej ID jednotlivých tetraedrů jako celltdata
    original_tet_ids = collect(valid_tet_ids)
    vtk_cell_data(vtkfile, original_tet_ids, "original_tet_id")
    
    # Ulož VTK soubor
    vtk_save(vtkfile)
    
    @info "Exportováno $(length(valid_tet_ids)) tetraedrů a $n_nodes uzlů do souboru $filename"
    
    # Vrať užitečné informace pro další analýzu
    return Dict(
        "exported_tet_ids" => valid_tet_ids,
        "node_map" => node_map,
        "n_nodes" => n_nodes
    )
end

filename = "vybrane_elementy"
tet_ids = [2621,2644].+1

export_specific_tetrahedra(mesh, tet_ids, filename)


function detect_close_nodes(mesh::BlockMesh; proximity_factor::Float64=0.1)
    # Tolerance pro detekci blízkých uzlů (0.1 * velikost hrany mřížky)
    proximity_threshold = proximity_factor * mesh.grid_step
    close_pairs = Vector{Tuple{Int,Int,Float64}}()
    
    # Vytvořím prostorové hašování pro efektivní vyhledávání
    spatial_hash = Dict{NTuple{3,Int},Vector{Int}}()
    
    # Určím faktor kvantizace
    quant_factor = proximity_threshold / 3.0
    
    # Vložím uzly do prostorového haše
    for (i, p) in enumerate(mesh.X)
        # Kvantizace pozice pro prostorový haš
        cell = (
            floor(Int, p[1] / quant_factor),
            floor(Int, p[2] / quant_factor),
            floor(Int, p[3] / quant_factor)
        )
        
        if !haskey(spatial_hash, cell)
            spatial_hash[cell] = Int[]
        end
        push!(spatial_hash[cell], i)
    end
    
    # Zkontroluje každou buňku a jejích 26 sousedů
    for (cell, node_indices) in spatial_hash
        # Zkontroluje uzly v rámci stejné buňky
        for i in 1:length(node_indices)
            for j in i+1:length(node_indices)
                node_i, node_j = node_indices[i], node_indices[j]
                distance = norm(mesh.X[node_i] - mesh.X[node_j])
                if distance < proximity_threshold
                    push!(close_pairs, (node_i, node_j, distance))
                end
            end
        end
        
        # Zkontroluje sousední buňky
        for dx in -1:1
            for dy in -1:1
                for dz in -1:1
                    (dx == 0 && dy == 0 && dz == 0) && continue
                    
                    neighbor = (cell[1] + dx, cell[2] + dy, cell[3] + dz)
                    haskey(spatial_hash, neighbor) || continue
                    
                    # Porovnání uzlů mezi buňkami
                    for node_i in node_indices
                        for node_j in spatial_hash[neighbor]
                            distance = norm(mesh.X[node_i] - mesh.X[node_j])
                            if distance < proximity_threshold
                                push!(close_pairs, (node_i, node_j, distance))
                            end
                        end
                    end
                end
            end
        end
    end
    
    return close_pairs
end

function check_mesh_for_duplicates(mesh::BlockMesh; proximity_factor::Float64=0.1)
    close_pairs = detect_close_nodes(mesh; proximity_factor=proximity_factor)
    
    if isempty(close_pairs)
        @info "Mesh je v pořádku: Nebyly nalezeny žádné blízké uzly (faktor blízkosti: $proximity_factor)"
        return false
    else
        @info "PROBLÉM: Nalezeno $(length(close_pairs)) párů blízkých uzlů (faktor blízkosti: $proximity_factor)"
        # Zobrazím několik příkladů
        if length(close_pairs) > 5
            @info "Prvních 5 párů blízkých uzlů (uzel_i, uzel_j, vzdálenost):"
            for i in 1:5
                @info "  $(close_pairs[i])"
            end
        else
            @info "Všechny páry blízkých uzlů:"
            for pair in close_pairs
                @info "  $pair"
            end
        end
        return true
    end
end

detect_close_nodes(mesh)
check_mesh_for_duplicates(mesh)