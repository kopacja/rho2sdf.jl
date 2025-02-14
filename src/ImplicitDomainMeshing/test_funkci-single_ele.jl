function slice_ambiguous_tetrahedra!(mesh::BlockMesh)
    new_IEN = Vector{Vector{Int64}}()
    found_valid_result = false
    original_tet = nothing  # Pro uložení původního tetraedru
    
    for tet in mesh.IEN
        new_tets = apply_stencil(mesh, tet)
        
        if !isempty(new_tets)
            for nt in new_tets
                t_sdf = [mesh.node_sdf[i] for i in nt]
                if any(x -> x >= 0, t_sdf)
                    push!(new_IEN, nt)
                    if !found_valid_result
                        original_tet = tet  # Uložíme původní tetraedron
                        found_valid_result = true
                        
                        # Vyexportujeme samostatný VTK soubor s původním tetraedrem
                        original_mesh = deepcopy(mesh)  # Vytvoříme kopii meshe
                        original_mesh.IEN = [tet]  # Ponecháme pouze původní tetraedron
                        update_connectivity!(original_mesh)  # Aktualizujeme topologii
                        export_mesh_vtk(original_mesh, "original_tetrahedron.vtu")
                        
                        println("Původní tetraedron: ", tet)
                        println("Rozdělené tetraedry: ", new_tets)
                    end
                end
            end
            if found_valid_result
                break
            end
        end
    end
    
    mesh.IEN = new_IEN
    @info "Po kompletním řezu tetraedrů: $(length(mesh.IEN)) tetraedrů"
end

mesh = BlockMesh()
@time generate_mesh!(mesh, "A15")
@time warp!(mesh)
update_connectivity!(mesh)

# Tento krok vytvoří dva soubory:
# 1. original_tetrahedron.vtu - obsahuje pouze původní tetraedron
# 2. block-mesh.vtu - obsahuje rozdělené tetraedry
slice_ambiguous_tetrahedra!(mesh)
update_connectivity!(mesh)
export_mesh_vtk(mesh, "block-mesh.vtu")