using LinearAlgebra
using WriteVTK
using StaticArrays

# Reprezentace čtyřstěnu
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

# Hlavní třída pro generování sítě koule
mutable struct SphereMesh
    nx::Int
    ny::Int
    nz::Int
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    SDF::Array{Float64,3}
    tetrahedra::Vector{TetrahedronElement}

    function SphereMesh(nx::Int=30, ny::Int=30, nz::Int=30)
        nx > 1 || throw(ArgumentError("nx must be > 1"))
        ny > 1 || throw(ArgumentError("ny must be > 1"))
        nz > 1 || throw(ArgumentError("nz must be > 1"))
        
        mesh = new(nx, ny, nz)
        mesh.tetrahedra = TetrahedronElement[]
        setup_grid!(mesh)
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

# Vytvoření tetraedrů pro jeden element
function tetrahedralize_regular_cell(mesh::SphereMesh, i::Int, j::Int, k::Int)
    vertices = Matrix{Float64}(undef, 8, 3)
    
    vertices[1,:] .= [mesh.x[i],   mesh.y[j],   mesh.z[k]]    # front-bottom-left
    vertices[2,:] .= [mesh.x[i+1], mesh.y[j],   mesh.z[k]]    # front-bottom-right
    vertices[3,:] .= [mesh.x[i+1], mesh.y[j+1], mesh.z[k]]    # front-top-right
    vertices[4,:] .= [mesh.x[i],   mesh.y[j+1], mesh.z[k]]    # front-top-left
    vertices[5,:] .= [mesh.x[i],   mesh.y[j],   mesh.z[k+1]]  # back-bottom-left
    vertices[6,:] .= [mesh.x[i+1], mesh.y[j],   mesh.z[k+1]]  # back-bottom-right
    vertices[7,:] .= [mesh.x[i+1], mesh.y[j+1], mesh.z[k+1]]  # back-top-right
    vertices[8,:] .= [mesh.x[i],   mesh.y[j+1], mesh.z[k+1]]  # back-top-left
    
    tets = SMatrix{6,4,Int}([
        1 2 4 6;  # První čtyřstěn
        2 3 4 6;  # Druhý čtyřstěn
        1 4 6 8;  # Třetí čtyřstěn
        1 6 5 8;  # Čtvrtý čtyřstěn
        3 4 7 6;  # Pátý čtyřstěn
        4 6 7 8   # Šestý čtyřstěn
    ])
    
    return vertices, tets
end

# Zpracování jedné buňky
function process_cell(mesh::SphereMesh, i::Int, j::Int, k::Int)
    vertices, tets = tetrahedralize_regular_cell(mesh, i, j, k)
    sdf_values = get_cell_sdf_values(mesh, i, j, k)
    
    # Kontrola, zda jsou všechny hodnoty SDF kladné
    if !all(>=(0), sdf_values)
        return TetrahedronElement[]
    end
    
    # Pokud jsme zde, element je celý uvnitř koule
    result_tets = TetrahedronElement[]
    
    for tet_idx in 1:size(tets, 1)
        tet_indices = view(tets, tet_idx, :)
        tet_vertices = vertices[tet_indices, :]
        tet_sdf = Vector(sdf_values[tet_indices])
        push!(result_tets, TetrahedronElement(tet_vertices, tet_sdf))
    end
    
    return result_tets
end

# Generování sítě
function generate_mesh!(mesh::SphereMesh)
    empty!(mesh.tetrahedra)
    sizehint!(mesh.tetrahedra, mesh.nx * mesh.ny * mesh.nz * 6)
    
    for i in 1:mesh.nx-1
        for j in 1:mesh.ny-1
            for k in 1:mesh.nz-1
                append!(mesh.tetrahedra, process_cell(mesh, i, j, k))
            end
        end
    end
end

# Export do VTK formátu
function export_vtk(mesh::SphereMesh, filename::String)
    if !endswith(filename, ".vtu")
        filename = filename * ".vtu"
    end
    
    points_dict = Dict{NTuple{3,Float64}, Int}()
    points = Matrix{Float64}(undef, 3, 0)
    sdf_values = Float64[]
    
    cells = Vector{MeshCell{VTKCellTypes.VTKCellType, Vector{Int}}}()
    sizehint!(cells, length(mesh.tetrahedra))
    
    for tet in mesh.tetrahedra
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
    
    vtk = vtk_grid(filename, points, cells)
    vtk["SDF"] = sdf_values
    vtk["radius"] = fill(0.5, length(sdf_values))
    vtk_save(vtk)
    
    @info "VTK soubor uložen jako: $(abspath(filename))"
    @info "Exported $(size(points, 2)) points and $(length(cells)) tetrahedra"
end

# Hlavní funkce
function main()
    mesh = SphereMesh(30, 30, 30)
    @info "Generování sítě..."
    generate_mesh!(mesh)
    @info "Export do VTK..."
    export_vtk(mesh, "sphere_interior_mesh")
    @info "Hotovo! Počet tetraedrů: $(length(mesh.tetrahedra))"
end

# Spuštění
main()
