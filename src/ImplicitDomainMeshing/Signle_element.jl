using LinearAlgebra
using WriteVTK
using StaticArrays

# Ponechávám původní struktury TetrahedronElement a SphereMesh

struct TetrahedronElement
    vertices::SMatrix{4,3,Float64}
    sdf_values::SVector{4,Float64}
    
    function TetrahedronElement(vertices::AbstractMatrix{Float64}, sdf_values::AbstractVector{Float64})
        size(vertices, 1) == 4 || throw(ArgumentError("Vertices must have 4 rows"))
        size(vertices, 2) == 3 || throw(ArgumentError("Vertices must have 3 columns"))
        length(sdf_values) == 4 || throw(ArgumentError("SDF values must have length 4"))
        return new(SMatrix{4,3}(vertices), SVector{4}(sdf_values))
    end
end

# Funkce pro vytvoření jednoho pravidelného elementu
function create_regular_element()
    # Vytvoření vrcholů krychle
    vertices = Matrix{Float64}([
        0.0 0.0 0.0;  # 1: front-bottom-left
        1.0 0.0 0.0;  # 2: front-bottom-right
        1.0 1.0 0.0;  # 3: front-top-right
        0.0 1.0 0.0;  # 4: front-top-left
        0.0 0.0 1.0;  # 5: back-bottom-left
        1.0 0.0 1.0;  # 6: back-bottom-right
        1.0 1.0 1.0;  # 7: back-top-right
        0.0 1.0 1.0   # 8: back-top-left
    ])
    
    # Definice tetraedrů (1-based indexing)
    tets = SMatrix{5,4,Int}([
        1 2 4 6;  # První čtyřstěn
        2 3 4 6;  # Druhý čtyřstěn
        1 4 6 8;  # Třetí čtyřstěn
        1 6 5 8;  # Čtvrtý čtyřstěn
        3 4 7 6;  # Pátý čtyřstěn
        4 6 7 8   # Šestý čtyřstěn
    ])
    
    # Vytvoříme fiktivní SDF hodnoty (všechny kladné pro jednoduchý test)
    sdf_values = ones(8)
    
    # Vytvoření tetraedrů
    tetrahedra = TetrahedronElement[]
    for i in 1:size(tets, 1)
        tet_vertices = vertices[tets[i,:], :]
        tet_sdf = sdf_values[tets[i,:]]
        push!(tetrahedra, TetrahedronElement(tet_vertices, tet_sdf))
    end
    
    return tetrahedra
end

# Export do VTK formátu
function export_single_element_vtk(tetrahedra::Vector{TetrahedronElement}, filename::String)
    if !endswith(filename, ".vtu")
        filename = filename * ".vtu"
    end
    
    points_dict = Dict{NTuple{3,Float64}, Int}()
    points = Matrix{Float64}(undef, 3, 0)
    sdf_values = Float64[]
    
    cells = Vector{MeshCell{VTKCellTypes.VTKCellType, Vector{Int}}}()
    
    # Zpracování všech tetraedrů
    for tet in tetrahedra
        tet_points = zeros(Int, 4)
        
        for i in 1:4
            vertex = Tuple(tet.vertices[i,:])
            
            if !haskey(points_dict, vertex)
                points = hcat(points, [vertex[1], vertex[2], vertex[3]])
                push!(sdf_values, tet.sdf_values[i])
                points_dict[vertex] = size(points, 2)
            end
            
            tet_points[i] = points_dict[vertex]
        end
        
        push!(cells, MeshCell(VTKCellTypes.VTK_TETRA, tet_points))
    end
    
    if size(points, 2) == 0
        @warn "No points to export!"
        return
    end
    
    # Export základní krychle (wireframe)
    cube_vertices = [
        0 0 0;
        1 0 0;
        1 1 0;
        0 1 0;
        0 0 1;
        1 0 1;
        1 1 1;
        0 1 1
    ]'  # Transpozice pro správný formát
    
    cube_cells = [
        MeshCell(VTKCellTypes.VTK_LINE, [1,2]),
        MeshCell(VTKCellTypes.VTK_LINE, [2,3]),
        MeshCell(VTKCellTypes.VTK_LINE, [3,4]),
        MeshCell(VTKCellTypes.VTK_LINE, [4,1]),
        MeshCell(VTKCellTypes.VTK_LINE, [5,6]),
        MeshCell(VTKCellTypes.VTK_LINE, [6,7]),
        MeshCell(VTKCellTypes.VTK_LINE, [7,8]),
        MeshCell(VTKCellTypes.VTK_LINE, [8,5]),
        MeshCell(VTKCellTypes.VTK_LINE, [1,5]),
        MeshCell(VTKCellTypes.VTK_LINE, [2,6]),
        MeshCell(VTKCellTypes.VTK_LINE, [3,7]),
        MeshCell(VTKCellTypes.VTK_LINE, [4,8])
    ]
    
    # Export drátěného modelu krychle
    vtk_cube = vtk_grid("cube_wireframe", cube_vertices, cube_cells)
    vtk_save(vtk_cube)
    
    # Export tetraedrů
    vtk = vtk_grid(filename, points, cells)
    vtk["element_index"] = collect(1:length(cells))  # Pro rozlišení jednotlivých tetraedrů
    vtk_save(vtk)
    
    @info "Exportováno $(size(points, 2)) bodů a $(length(cells)) tetraedrů"
end

# Hlavní testovací funkce
function test_single_element()
    @info "Vytváření tetraedrální diskretizace jednoho elementu..."
    tetrahedra = create_regular_element()
    
    @info "Export do VTK..."
    export_single_element_vtk(tetrahedra, "single_element_test")
    
    @info "Hotovo! Vytvořeno $(length(tetrahedra)) tetraedrů"
end

# Spuštění testu
test_single_element()
