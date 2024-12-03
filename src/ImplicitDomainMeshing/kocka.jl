
# Imports
using LinearAlgebra
using Statistics
using StaticArrays
using WriteVTK
using JLD2
using ProgressMeter

# Rozšíření původní BlockMesh struktury
mutable struct BlockMesh
  nx::Int
  ny::Int
  nz::Int
  grid::Array{Vector{Float64},3}   # 3D pole vektorů souřadnic uzlů
  SDF::Array{Float64,3}            # SDF hodnoty
  X::Vector{Vector{Float64}}       # Souřadnice použitých uzlů
  IEN::Vector{Vector{Int64}}       # Konektivita elementů
  INE::Vector{Vector{Int64}}       # Inverzní konektivita (Node -> Elements)
  node_sdf::Vector{Float64}        # SDF hodnoty pro použité uzly
  node_map::Dict{Int64,Int64}      # původní_index => nový_index
  boundary_vertices::Set{Int64}    # Množina hraničních uzlů

  function BlockMesh()
    # Načtení dat (ponecháno z původního kódu)
    # @load "src/ImplicitDomainMeshing/data/Z_chapadlo_FineGrid.jld2" fine_grid
    # @load "src/ImplicitDomainMeshing/data/Z_chapadlo_FineSDF.jld2" fine_sdf
    @load "src/ImplicitDomainMeshing/data/Z_block_FineGrid.jld2" fine_grid
    @load "src/ImplicitDomainMeshing/data/Z_block_FineSDF.jld2" fine_sdf

    # Konverze Float32 na Float64 pro konzistenci
    grid = Array{Vector{Float64},3}(undef, size(fine_grid))
    for i in eachindex(fine_grid)
      grid[i] = Vector{Float64}(fine_grid[i])
    end
    sdf = Float64.(fine_sdf)
    nx, ny, nz = size(grid)

    # Inicializace nové instance
    mesh = new(nx, ny, nz)
    mesh.grid = grid
    mesh.SDF = sdf
    mesh.node_map = Dict{Int64,Int64}()
    mesh.X = Vector{Vector{Float64}}()
    mesh.IEN = Vector{Vector{Int64}}()
    mesh.INE = Vector{Vector{Int64}}()
    mesh.node_sdf = Vector{Float64}()
    mesh.boundary_vertices = Set{Int64}()
    return mesh
  end
end

# Pomocná funkce pro získání původního globálního indexu uzlu z i,j,k indexů
function get_original_node_id(mesh::BlockMesh, i::Int, j::Int, k::Int)
  return i + (j - 1) * mesh.nx + (k - 1) * mesh.nx * mesh.ny
end

# Pomocná funkce pro získání nebo vytvoření nového indexu uzlu
function get_or_create_node!(mesh::BlockMesh, i::Int, j::Int, k::Int)
  orig_id = get_original_node_id(mesh, i, j, k)
  if !haskey(mesh.node_map, orig_id)
    # Vytvoření nového uzlu
    push!(mesh.X, mesh.grid[i, j, k])
    push!(mesh.node_sdf, mesh.SDF[i, j, k])
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
    mesh.SDF[i, j, k],        # front-bottom-left
    mesh.SDF[i+1, j, k],      # front-bottom-right
    mesh.SDF[i+1, j+1, k],    # front-top-right
    mesh.SDF[i, j+1, k],      # front-top-left
    mesh.SDF[i, j, k+1],      # back-bottom-left
    mesh.SDF[i+1, j, k+1],    # back-bottom-right
    mesh.SDF[i+1, j+1, k+1],  # back-top-right
    mesh.SDF[i, j+1, k+1]     # back-top-left
  )
end


# Konstanty pro acute tetrahedral lattice definici podle C++ implementace
const LATTICE_NODES = [
  [1, 1, 1],
  [3, 1, 1],
  [3, 3, 1],
  [1, 3, 1],
  [2, 2, 1], # bottom face centre
  [2, 2, 2], # cube centre
  [1, 2, 2], # left  face centre
  [2, 1, 2], # back  face centre
  [3, 2, 2], # right face centre
  [2, 3, 2], # front face centre
  [1, 1, 3],
  [3, 1, 3],
  [3, 3, 3],
  [1, 3, 3],
  [2, 2, 3], # top face centre
]

const LATTICE_TETS = [
  # spotní +
  [5, 4, 1, 6],
  [5, 1, 2, 6],
  [5, 2, 3, 6],
  [5, 3, 4, 6],
  # levý +
  [1, 4, 7, 6],
  [4, 14, 7, 6],
  [14, 11, 7, 6],
  [11, 1, 7, 6],
  # zadní +
  [2, 1, 8, 6],
  [12, 2, 8, 6],
  [11, 12, 8, 6],
  [1, 11, 8, 6],
  # pravý +
  [3, 2, 9, 6],
  [13, 3, 9, 6],
  [12, 13, 9, 6],
  [2, 12, 9, 6],
  # přední +
  [4, 3, 10, 6],
  [14, 4, 10, 6],
  [13, 14, 10, 6],
  [3, 13, 10, 6],
  # horní +
  [12, 11, 15, 6],
  [13, 12, 15, 6],
  [14, 13, 15, 6],
  [11, 14, 15, 6],
]

# Vylepšená funkce pro process_cell! používající acute tetrahedral lattice
function process_cell!(mesh::BlockMesh, i::Int, j::Int, k::Int)
  # Vytvoření lokální mřížky pro buňku
  function create_lattice_node(base_i::Int, base_j::Int, base_k::Int, lattice_pos)
    dx = mesh.grid[2, 1, 1][1] - mesh.grid[1, 1, 1][1]
    pos = mesh.grid[base_i, base_j, base_k] +
          [lattice_pos[1] * dx / 2, lattice_pos[2] * dx / 2, lattice_pos[3] * dx / 2]
    return pos
  end

  # Mapování pro lokální vrcholy (indexujeme od 1)
  local_vertices = Dict{Int,Vector{Float64}}()
  local_sdf = Dict{Int,Float64}()

  # Vytvoření vrcholů podle lattice definice (indexujeme od 1)
  for (idx, node) in enumerate(LATTICE_NODES)
    pos = create_lattice_node(i, j, k, node)
    # Upravíme indexování aby začínalo od 1
    local_vertices[idx] = pos
    local_sdf[idx] = evaluate_sdf(mesh, pos)
  end

  # Kontrola SDF hodnot pro celou buňku
  has_positive = any(>=(0), values(local_sdf))
  has_negative = any(<(0), values(local_sdf))

  if !has_positive
    return
  end

  # Vytvoření nebo získání globálních indexů vrcholů
  global_indices = Dict{Int,Int}()
  for (idx, pos) in local_vertices
    # Přidání vrcholu do globální sítě
    push!(mesh.X, pos)
    push!(mesh.node_sdf, local_sdf[idx])
    global_indices[idx] = length(mesh.X)
  end

  # Vytvoření tetraedrů podle lattice definice
  for tet in LATTICE_TETS
    tet_sdf = [local_sdf[v] for v in tet]
    if any(>=(0), tet_sdf)
      # Použijeme přímo indexy z tet, protože už jsou 1-based
      push!(mesh.IEN, [global_indices[v] for v in tet])
    end
  end
end

# Funkce pro vyčištění nepoužitých uzlů

function cleanup_unused_nodes!(mesh::BlockMesh)
  # Vytvoření množiny použitých uzlů
  used_nodes = Set{Int64}()
  @showprogress desc = "Sběr použitých uzlů" for element in mesh.IEN
    union!(used_nodes, element)
  end

  # Předalokace polí a jednorázové seřazení
  sorted_nodes = sort(collect(used_nodes))
  n_nodes = length(sorted_nodes)
  new_node_map = Dict{Int64,Int64}(sorted_nodes[i] => i for i in 1:n_nodes)
  new_X = Vector{Vector{Float64}}(undef, n_nodes)
  new_node_sdf = Vector{Float64}(undef, n_nodes)

  # Přemapování uzlů a odstranění nepoužitých
  @showprogress desc = "Přemapování uzlů" for i in 1:n_nodes
    new_X[i] = mesh.X[sorted_nodes[i]]
    new_node_sdf[i] = mesh.node_sdf[sorted_nodes[i]]
  end

  # Aktualizace konektivity elementů
  @showprogress desc = "Aktualizace konektivity" for i in 1:length(mesh.IEN)
    mesh.IEN[i] = [new_node_map[old_id] for old_id in mesh.IEN[i]]
  end

  # Aktualizace dat sítě
  mesh.X = new_X
  mesh.node_sdf = new_node_sdf
  mesh.node_map = new_node_map
  @info "Počet uzlů po vyčištění: $(length(mesh.X))"
end

# Funkce pro vytvoření INE
function create_INE!(mesh::BlockMesh)
  # Inicializace INE pro každý uzel prázdným vektorem
  mesh.INE = [Vector{Int64}() for _ in 1:length(mesh.X)]

  # Procházení všech elementů a přiřazení k uzlům
  for (elem_id, element) in enumerate(mesh.IEN)
    for node_id in element
      push!(mesh.INE[node_id], elem_id)
    end
  end
end

##___

# Pomocné funkce pro interpolaci a výpočet gradientu SDF
function trilinear_interpolate(values::Array{Float64,3}, x::Float64, y::Float64, z::Float64)
  # Trilineární interpolace mezi osmi vrcholy krychle
  c00 = values[1, 1, 1] * (1 - x) * (1 - y) * (1 - z) + values[2, 1, 1] * x * (1 - y) * (1 - z)
  c01 = values[1, 1, 2] * (1 - x) * (1 - y) * z + values[2, 1, 2] * x * (1 - y) * z
  c10 = values[1, 2, 1] * (1 - x) * y * (1 - z) + values[2, 2, 1] * x * y * (1 - z)
  c11 = values[1, 2, 2] * (1 - x) * y * z + values[2, 2, 2] * x * y * z

  return c00 + c01 + c10 + c11
end

function compute_sdf_gradient(mesh::BlockMesh, point::Vector{Float64})
  # Nalezení nejbližší buňky v mřížce
  dx = mesh.grid[2, 1, 1][1] - mesh.grid[1, 1, 1][1]
  dy = mesh.grid[1, 2, 1][2] - mesh.grid[1, 1, 1][2]
  dz = mesh.grid[1, 1, 2][3] - mesh.grid[1, 1, 1][3]

  i = floor(Int, (point[1] - mesh.grid[1, 1, 1][1]) / dx) + 1
  j = floor(Int, (point[2] - mesh.grid[1, 1, 1][2]) / dy) + 1
  k = floor(Int, (point[3] - mesh.grid[1, 1, 1][3]) / dz) + 1

  # Omezení na hranice mřížky
  i = clamp(i, 1, mesh.nx - 1)
  j = clamp(j, 1, mesh.ny - 1)
  k = clamp(k, 1, mesh.nz - 1)

  # Centrální diference pro gradient
  grad_x = (mesh.SDF[i+1, j, k] - mesh.SDF[i-1, j, k]) / (2dx)
  grad_y = (mesh.SDF[i, j+1, k] - mesh.SDF[i, j-1, k]) / (2dy)
  grad_z = (mesh.SDF[i, j, k+1] - mesh.SDF[i, j, k-1]) / (2dz)

  return [grad_x, grad_y, grad_z]
end

# Funkce pro výpočet SDF hodnoty v libovolném bodě
function evaluate_sdf(mesh::BlockMesh, point::Vector{Float64})
  dx = mesh.grid[2, 1, 1][1] - mesh.grid[1, 1, 1][1]
  dy = mesh.grid[1, 2, 1][2] - mesh.grid[1, 1, 1][2]
  dz = mesh.grid[1, 1, 2][3] - mesh.grid[1, 1, 1][3]

  i = floor(Int, (point[1] - mesh.grid[1, 1, 1][1]) / dx) + 1
  j = floor(Int, (point[2] - mesh.grid[1, 1, 1][2]) / dy) + 1
  k = floor(Int, (point[3] - mesh.grid[1, 1, 1][3]) / dz) + 1

  # Omezení na hranice mřížky
  i = clamp(i, 1, mesh.nx - 1)
  j = clamp(j, 1, mesh.ny - 1)
  k = clamp(k, 1, mesh.nz - 1)

  # Relativní pozice v buňce
  x_frac = (point[1] - mesh.grid[i, j, k][1]) / dx
  y_frac = (point[2] - mesh.grid[i, j, k][2]) / dy
  z_frac = (point[3] - mesh.grid[i, j, k][3]) / dz

  return trilinear_interpolate(mesh.SDF[i:i+1, j:j+1, k:k+1], x_frac, y_frac, z_frac)
end

#
#___


# Aktualizovaná funkce pro generování sítě
function generate_mesh!(mesh::BlockMesh)
  # Vyčištění existujících dat
  empty!(mesh.X)
  empty!(mesh.IEN)
  empty!(mesh.node_sdf)
  empty!(mesh.node_map)
  empty!(mesh.boundary_vertices)

  # Generování základní sítě (původní implementace)
  @info "Generování základní sítě..."
  for i in 1:mesh.nx-1
    for j in 1:mesh.ny-1
      for k in 1:mesh.nz-1
        process_cell!(mesh, i, j, k)
      end
    end
  end

  # Vyčištění nepoužitých uzlů
  cleanup_unused_nodes!(mesh)

  # Vytvoření inverzní konektivity
  create_INE!(mesh)

  @info "Základní síť: $(length(mesh.X)) uzlů a $(length(mesh.IEN)) tetraedrů"

end

# Dokončení funkce pro export do VTK
function export_to_vtk(mesh::BlockMesh, filename::String)
  # Převod vektorů na matici bodů
  points = zeros(Float64, 3, length(mesh.X))
  for (i, x) in enumerate(mesh.X)
    points[:, i] = x
  end

  # Převod IEN na pole VTK buněk
  cells = [MeshCell(VTKCellTypes.VTK_TETRA, tet) for tet in mesh.IEN]

  # Vytvoření VTK souboru
  vtk = vtk_grid(filename, points, cells)

  # Přidání skalárních dat (SDF hodnoty)
  vtk["sdf"] = mesh.node_sdf

  # Označení hraničních uzlů
  boundary_markers = zeros(Int, length(mesh.X))
  for v in mesh.boundary_vertices
    boundary_markers[v] = 1
  end
  vtk["boundary"] = boundary_markers

  # Uložení souboru
  vtk_save(vtk)
end

# Příklad použití
function main()
  # Vytvoření instance sítě
  mesh = BlockMesh()

  # Generování a optimalizace sítě
  generate_mesh!(mesh)

  # Export do VTK pro vizualizaci
  export_to_vtk(mesh, "optimized_mesh")
end

# Spuštění
main()
