using LinearAlgebra
using WriteVTK
using StaticArrays

mutable struct SphereMesh
    nx::Int
    ny::Int
    nz::Int
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    SDF::Array{Float64,3}
    X::Vector{Vector{Float64}}    # Souřadnice použitých uzlů
    IEN::Vector{Vector{Int64}}    # Konektivita elementů
    node_sdf::Vector{Float64}     # SDF hodnoty pro použité uzly
    
    # Pomocná mapa pro sledování použitých uzlů a jejich nových indexů
    node_map::Dict{Int64, Int64}  # původní_index => nový_index

    function SphereMesh(nx::Int=30, ny::Int=30, nz::Int=30)
        nx > 1 || throw(ArgumentError("nx must be > 1"))
        ny > 1 || throw(ArgumentError("ny must be > 1"))
        nz > 1 || throw(ArgumentError("nz must be > 1"))
        
        mesh = new(nx, ny, nz)
        setup_grid!(mesh)
        mesh.node_map = Dict{Int64, Int64}()
        mesh.X = Vector{Vector{Float64}}()
        mesh.IEN = Vector{Vector{Int64}}()
        mesh.node_sdf = Vector{Float64}()
        return mesh
    end
end

# Inicializace výpočetní mřížky
function setup_grid!(mesh::SphereMesh)
    mesh.x = collect(range(-1, 1, length=mesh.nx))
    mesh.y = collect(range(-1, 1, length=mesh.ny))
    mesh.z = collect(range(-1, 1, length=mesh.nz))
    
    mesh.SDF = Array{Float64}(undef, mesh.nx, mesh.ny, mesh.nz)
    radius = 0.5
    
    for i in 1:mesh.nx, j in 1:mesh.ny, k in 1:mesh.nz
        mesh.SDF[i,j,k] = radius - sqrt(mesh.x[i]^2 + mesh.y[j]^2 + mesh.z[k]^2)
    end
end

# Pomocná funkce pro získání původního globálního indexu uzlu z i,j,k indexů
function get_original_node_id(mesh::SphereMesh, i::Int, j::Int, k::Int)
    return i + (j-1)*mesh.nx + (k-1)*mesh.nx*mesh.ny
end

# Pomocná funkce pro získání nebo vytvoření nového indexu uzlu
function get_or_create_node!(mesh::SphereMesh, i::Int, j::Int, k::Int)
    orig_id = get_original_node_id(mesh, i, j, k)
    
    if !haskey(mesh.node_map, orig_id)
        # Vytvoření nového uzlu
        push!(mesh.X, [mesh.x[i], mesh.y[j], mesh.z[k]])
        push!(mesh.node_sdf, mesh.SDF[i,j,k])
        mesh.node_map[orig_id] = length(mesh.X)
    end
    
    return mesh.node_map[orig_id]
end

# Získání SDF hodnot pro vrcholy buňky
function get_cell_sdf_values(mesh::SphereMesh, i::Int, j::Int, k::Int)
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
function process_cell!(mesh::SphereMesh, i::Int, j::Int, k::Int)
    # Získání SDF hodnot
    sdf_values = get_cell_sdf_values(mesh, i, j, k)
    
    # Kontrola, zda jsou všechny hodnoty SDF kladné
    if !all(>=(0), sdf_values)
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
function generate_mesh!(mesh::SphereMesh)
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
function export_vtk(mesh::SphereMesh, filename::String)
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
    mesh = SphereMesh(30, 30, 30)
    @info "Generování sítě..."
    generate_mesh!(mesh)
    @info "Export do VTK..."
    export_vtk(mesh, "sphere_interior_mesh")
end

# Spuštění
main()
