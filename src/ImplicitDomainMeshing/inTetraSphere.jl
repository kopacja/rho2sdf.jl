using LinearAlgebra
using WriteVTK
using StaticArrays
using JLD2

mutable struct BlockMesh
    nx::Int
    ny::Int
    nz::Int
    grid::Array{Vector{Float64}, 3}   # 3D pole vektorů souřadnic uzlů
    SDF::Array{Float64,3}             # SDF hodnoty
    X::Vector{Vector{Float64}}        # Souřadnice použitých uzlů
    IEN::Vector{Vector{Int64}}        # Konektivita elementů
    node_sdf::Vector{Float64}         # SDF hodnoty pro použité uzly
    node_map::Dict{Int64, Int64}      # původní_index => nový_index

    function BlockMesh()
        # Načtení dat
        @load "data/Z_block_FineGrid.jld2" fine_grid
        # @load "src/ImplicitDomainMeshing/data/Z_block_FineGrid.jld2" fine_grid
        @load "data/Z_block_FineSDF.jld2" fine_sdf
        # @load "src/ImplicitDomainMeshing/data/Z_block_FineSDF.jld2" fine_sdf
        
        # Konverze Float32 na Float64 pro konzistenci
        grid = Array{Vector{Float64}, 3}(undef, size(fine_grid))
        for i in eachindex(fine_grid)
            grid[i] = Vector{Float64}(fine_grid[i])
        end
        
        sdf = Float64.(fine_sdf)
        
        nx, ny, nz = size(grid)
        
        mesh = new(nx, ny, nz)
        mesh.grid = grid
        mesh.SDF = sdf
        mesh.node_map = Dict{Int64, Int64}()
        mesh.X = Vector{Vector{Float64}}()
        mesh.IEN = Vector{Vector{Int64}}()
        mesh.node_sdf = Vector{Float64}()
        return mesh
    end
end

# Pomocná funkce pro získání původního globálního indexu uzlu z i,j,k indexů
function get_original_node_id(mesh::BlockMesh, i::Int, j::Int, k::Int)
    return i + (j-1)*mesh.nx + (k-1)*mesh.nx*mesh.ny
end

# Pomocná funkce pro získání nebo vytvoření nového indexu uzlu
function get_or_create_node!(mesh::BlockMesh, i::Int, j::Int, k::Int)
    orig_id = get_original_node_id(mesh, i, j, k)
    
    if !haskey(mesh.node_map, orig_id)
        # Vytvoření nového uzlu
        push!(mesh.X, mesh.grid[i,j,k])
        push!(mesh.node_sdf, mesh.SDF[i,j,k])
        mesh.node_map[orig_id] = length(mesh.X)
    end
    
    return mesh.node_map[orig_id]
end

# Získání SDF hodnot pro vrcholy buňky
function get_cell_sdf_values(mesh::BlockMesh, i::Int, j::Int, k::Int)
    1 <= i < mesh.nx || throw(BoundsError(mesh.SDF, i))
    1 <= j < mesh.ny || throw(BoundsError(mesh.SDF, j))
    1 <= k < mesh.nz || throw(BoundsError(mesh.SDF, k))
    
    return SVector{8}(
        mesh.SDF[i,j,k],        # front-bottom-left
        mesh.SDF[i+1,j,k],      # front-bottom-right
        mesh.SDF[i+1,j+1,k],    # front-top-right
        mesh.SDF[i,j+1,k],      # front-top-left
        mesh.SDF[i,j,k+1],      # back-bottom-left
        mesh.SDF[i+1,j,k+1],    # back-bottom-right
        mesh.SDF[i+1,j+1,k+1],  # back-top-right
        mesh.SDF[i,j+1,k+1]     # back-top-left
    )
end

# Zpracování jedné buňky a vytvoření tetraedrů
function process_cell!(mesh::BlockMesh, i::Int, j::Int, k::Int)
    # Získání SDF hodnot
    sdf_values = get_cell_sdf_values(mesh, i, j, k)
    
    # Kontrola, zda jsou všechny hodnoty SDF kladné
    if !any(>=(0), sdf_values)
        return
    end
    
    # Získání/vytvoření indexů vrcholů buňky
    node_ids = [
        get_or_create_node!(mesh, i, j, k),      # front-bottom-left
        get_or_create_node!(mesh, i+1, j, k),    # front-bottom-right
        get_or_create_node!(mesh, i+1, j+1, k),  # front-top-right
        get_or_create_node!(mesh, i, j+1, k),    # front-top-left
        get_or_create_node!(mesh, i, j, k+1),    # back-bottom-left
        get_or_create_node!(mesh, i+1, j, k+1),  # back-bottom-right
        get_or_create_node!(mesh, i+1, j+1, k+1),# back-top-right
        get_or_create_node!(mesh, i, j+1, k+1)   # back-top-left
    ]
    
    # Definice tetraedrů
    tet_connectivity = [
        [1, 2, 4, 6],  # První čtyřstěn
        [2, 3, 4, 6],  # Druhý čtyřstěn
        [1, 4, 6, 8],  # Třetí čtyřstěn
        [1, 6, 5, 8],  # Čtvrtý čtyřstěn
        [3, 4, 7, 6],  # Pátý čtyřstěn
        [4, 6, 7, 8]   # Šestý čtyřstěn
    ]
    
    # Vytvoření elementů s novými ID uzlů
    for tet in tet_connectivity
        push!(mesh.IEN, [node_ids[i] for i in tet])
    end
end

# Generování sítě
function generate_mesh!(mesh::BlockMesh)
    empty!(mesh.X)
    empty!(mesh.IEN)
    empty!(mesh.node_sdf)
    empty!(mesh.node_map)
    
    for i in 1:mesh.nx-1
        for j in 1:mesh.ny-1
            for k in 1:mesh.nz-1
                process_cell!(mesh, i, j, k)
            end
        end
    end
    
    @info "Vytvořeno $(length(mesh.X)) uzlů a $(length(mesh.IEN)) tetraedrů"
end

# Export do VTK formátu
function export_vtk(mesh::BlockMesh, filename::String)
    if !endswith(filename, ".vtu")
        filename = filename * ".vtu"
    end
    
    # Převod uzlů do formátu pro VTK (3×N matice)
    points = reduce(hcat, mesh.X)
    
    # Vytvoření VTK buněk
    cells = [MeshCell(VTKCellTypes.VTK_TETRA, tet) for tet in mesh.IEN]
    
    # Export
    vtk = vtk_grid(filename, points, cells)
    vtk["SDF"] = mesh.node_sdf
    vtk["radius"] = fill(0.5, length(mesh.node_sdf))
    vtk_save(vtk)
    
    @info "VTK soubor uložen jako: $(abspath(filename))"
    @info "Exportováno $(length(mesh.X)) uzlů a $(length(mesh.IEN)) tetraedrů"
end

# Hlavní funkce
function main()
    mesh = BlockMesh()
    @info "Generování sítě..."
    generate_mesh!(mesh)
    @info "Export do VTK..."
    export_vtk(mesh, "sphere_interior_mesh")
end

# Spuštění
main()
