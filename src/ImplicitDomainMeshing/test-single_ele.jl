"""
    test_single_tetrahedron(tet_nodes::Vector{SVector{3,Float64}}, node_sdf::Vector{Float64}, 
                           grid_tol::Float64=1e-6)

Testuje aplikaci šablony na jeden tetrahedrální element.

Vstupy:
- tet_nodes: vektor 4 bodů definujících tetrahedron (každý bod je SVector{3,Float64})
- node_sdf: vektor 4 SDF hodnot odpovídajících vrcholům tetrahedronu
- grid_tol: tolerance pro geometrické operace (volitelné)

Výstupy:
- Vytvoří dva VTK soubory:
  1. "original_test_tet.vtu" - původní tetrahedron
  2. "subdivided_test_tet.vtu" - výsledek po aplikaci šablony
"""

using WriteVTK
using StaticArrays
using Logging
using JLD2
using LinearAlgebra

# ----------------------------
# Block mesh structure (extended with node_hash)
# ----------------------------
mutable struct BlockMesh
  nx::Int
  ny::Int
  nz::Int
  grid::Array{SVector{3,Float64},3}          # 3D array of node coordinates (basic grid) using static vectors
  grid_step::Float64
  grid_tol::Float64
  SDF::Array{Float64,3}                      # SDF values
  X::Vector{SVector{3,Float64}}    # List of physical node coordinates (nodes used in mesh)
  IEN::Vector{Vector{Int64}}                 # Tetrahedral connectivity (elements)
  INE::Vector{Vector{Int64}}                 # Inverse connectivity: for each node, list of adjacent elements
  node_sdf::Vector{Float64}                  # SDF values at nodes
  node_map::Dict{Int64,Int64}                # Mapping from original grid node index -> new node index
  cell_center_map::Dict{Tuple{Int,Int,Int},Int64}  # Mapping for cell centers (Steiner points)
  node_hash::Dict{NTuple{3,Float64},Int64}    # Global dictionary for merging nodes

  function BlockMesh()
    @load "src/ImplicitDomainMeshing/data/Z_block_FineGrid_B-0.2_smooth-1.jld2" fine_grid
    @load "src/ImplicitDomainMeshing/data/Z_block_FineSDF_B-0.2_smooth-1.jld2" fine_sdf

    # Convert fine_grid into a 3D array of SVectors
    grid = Array{SVector{3,Float64},3}(undef, size(fine_grid))
    for i in eachindex(fine_grid)
      # Assume fine_grid[i] is an array of Float64; convert to SVector
      grid[i] = SVector{3,Float64}(fine_grid[i]...)
    end
    step = maximum(abs.(grid[1, 1, 1] - grid[2, 2, 2]))
    sdf = Float64.(fine_sdf)
    nx, ny, nz = size(grid)

    mesh = new(nx, ny, nz)
    mesh.grid = grid
    mesh.grid_step = step
    mesh.grid_tol = 1e-6 * step
    mesh.SDF = sdf
    mesh.node_map = Dict{Int64,Int64}()
    mesh.cell_center_map = Dict{Tuple{Int,Int,Int},Int64}()
    mesh.X = Vector{SVector{3,Float64}}()
    mesh.IEN = Vector{Vector{Int64}}()
    mesh.INE = Vector{Vector{Int64}}()
    mesh.node_sdf = Vector{Float64}()
    mesh.node_hash = Dict{NTuple{3,Float64},Int64}()
    return mesh
  end
end

function export_mesh_vtk(mesh::BlockMesh, filename::String)
    npoints = length(mesh.X)
    points = zeros(Float64, 3, npoints)
    @inbounds for i in 1:npoints
      points[:, i] = mesh.X[i]
    end
    cells = [MeshCell(VTKCellTypes.VTK_TETRA, tet) for tet in mesh.IEN]
    vtkfile = vtk_grid(filename, points, cells)
    vtk_point_data(vtkfile, mesh.node_sdf, "sdf")
    vtk_save(vtkfile)
  end

  function eval_sdf(mesh::BlockMesh, p::SVector{3,Float64})
    # Získáme minimální a maximální souřadnice mřížky
    vmin = mesh.grid[1, 1, 1]
    vmax = mesh.grid[end, end, end]
    # Normalizace souřadnic bodu p do intervalu [0,1]
    r = (p .- vmin) ./ (vmax .- vmin)
    # Přepočet na indexy v mřížce
    i_f = r[1] * (mesh.nx - 1) + 1
    j_f = r[2] * (mesh.ny - 1) + 1
    k_f = r[3] * (mesh.nz - 1) + 1
    i0 = clamp(floor(Int, i_f), 1, mesh.nx - 1)
    j0 = clamp(floor(Int, j_f), 1, mesh.ny - 1)
    k0 = clamp(floor(Int, k_f), 1, mesh.nz - 1)
    i1 = i0 + 1
    j1 = j0 + 1
    k1 = k0 + 1
    # Lokální váhy
    xd = i_f - i0
    yd = j_f - j0
    zd = k_f - k0
    # Získání SDF hodnot na osmi rozích buňky
    c000 = mesh.SDF[i0, j0, k0]
    c100 = mesh.SDF[i1, j0, k0]
    c010 = mesh.SDF[i0, j1, k0]
    c110 = mesh.SDF[i1, j1, k0]
    c001 = mesh.SDF[i0, j0, k1]
    c101 = mesh.SDF[i1, j0, k1]
    c011 = mesh.SDF[i0, j1, k1]
    c111 = mesh.SDF[i1, j1, k1]
    # Trilineární interpolace
    c00 = c000 * (1 - xd) + c100 * xd
    c01 = c001 * (1 - xd) + c101 * xd
    c10 = c010 * (1 - xd) + c110 * xd
    c11 = c011 * (1 - xd) + c111 * xd
    c0 = c00 * (1 - yd) + c10 * yd
    c1 = c01 * (1 - yd) + c11 * yd
    return c0 * (1 - zd) + c1 * zd
  end

  function quantize(p::SVector{3,Float64}, tol::Float64)
    return (round(p[1] / tol) * tol, round(p[2] / tol) * tol, round(p[3] / tol) * tol)
  end

function test_single_tetrahedron(tet_nodes::Vector{SVector{3,Float64}}, 
                                node_sdf::Vector{Float64}, 
                                grid_tol::Float64=1e-6)
    # Vytvoříme minimální mesh strukturu obsahující pouze testovací element
    test_mesh = BlockMesh()
    
    # Nastavíme pouze potřebné části meshe
    test_mesh.X = tet_nodes
    test_mesh.node_sdf = node_sdf
    test_mesh.grid_tol = grid_tol
    test_mesh.IEN = [[1, 2, 3, 4]]  # Jeden element s indexy 1-4
    
    # Exportujeme původní tetrahedron
    export_mesh_vtk(test_mesh, "original_test_tet.vtu")
    
    # Aplikujeme šablonu
    new_tets = apply_stencil(test_mesh, [1, 2, 3, 4])
    
    if !isempty(new_tets)
        # Vytvoříme nový mesh pro rozdělené tetraedry
        subdivided_mesh = deepcopy(test_mesh)
        subdivided_mesh.IEN = new_tets
        export_mesh_vtk(subdivided_mesh, "subdivided_test_tet.vtu")
        
        println("Původní tetrahedron: ", [1, 2, 3, 4])
        println("Rozdělené tetraedry: ", new_tets)
    else
        println("Šablona nevrátila žádné nové tetraedry")
    end
end

# Příklad použití:
# Definujeme testovací tetrahedron
test_nodes = [
    SVector{3,Float64}(0.0, 0.0, 0.0),  # vrchol 1
    SVector{3,Float64}(1.0, 0.0, 0.0),  # vrchol 2
    SVector{3,Float64}(0.0, 1.0, 0.0),  # vrchol 3
    SVector{3,Float64}(0.0, 0.0, 1.0)   # vrchol 4
]

# Definujeme SDF hodnoty pro každý vrchol
# test_sdf = [
#     -0.5,  # SDF pro vrchol 1
#      0.3,  # SDF pro vrchol 2
#      0.2,  # SDF pro vrchol 3
#     -0.   # SDF pro vrchol 4
# ]
test_sdf = [
    1.,  # SDF pro vrchol 1
     1.,  # SDF pro vrchol 2
     -1.,  # SDF pro vrchol 3
    -1.   # SDF pro vrchol 4
]

# Otestujeme tetrahedron
test_single_tetrahedron(test_nodes, test_sdf)
#nnn
#TODO: Otestovat příklady co mají fungovat jestli opravdu fungují
# (= otestovat tento test :)))exit()
# vypadá to že nefunguje (zamrzá?)
# konfigurace neodpovídají SDF