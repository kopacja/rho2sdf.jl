using WriteVTK
using StaticArrays
using Logging
using JLD2
using LinearAlgebra

# Globální tolerance pro slučování uzlů (můžete upravit podle velikosti elementů)
const TOL = 1e-6

# ----------------------------
# Struktura blokové sítě (rozšířená o node_hash)
# ----------------------------
mutable struct BlockMesh
  nx::Int
  ny::Int
  nz::Int
  grid::Array{Vector{Float64},3}       # 3D pole souřadnic uzlů (základní síť)
  SDF::Array{Float64,3}                 # SDF hodnoty
  X::Vector{Vector{Float64}}            # Seznam fyzikálních souřadnic použitých uzlů
  IEN::Vector{Vector{Int64}}            # Konektivita tetraedrů (elementů)
  INE::Vector{Vector{Int64}}            # Inverzní konektivita: pro každý uzel seznam přilehlých elementů
  node_sdf::Vector{Float64}             # SDF hodnoty pro uzly
  node_map::Dict{Int64,Int64}           # mapování původního indexu gridového uzlu -> nový index
  cell_center_map::Dict{Tuple{Int,Int,Int},Int64}  # mapování buňkových středu (Steiner body)
  node_hash::Dict{NTuple{3,Float64},Int64}  # globální slovník pro slučování uzlů

  function BlockMesh()
    @load "src/ImplicitDomainMeshing/data/Z_block_FineGrid.jld2" fine_grid
    @load "src/ImplicitDomainMeshing/data/Z_block_FineSDF.jld2" fine_sdf

    grid = Array{Vector{Float64},3}(undef, size(fine_grid))
    for i in eachindex(fine_grid)
      grid[i] = Vector{Float64}(fine_grid[i])
    end
    sdf = Float64.(fine_sdf)
    nx, ny, nz = size(grid)
    
    mesh = new(nx, ny, nz)
    mesh.grid = grid
    mesh.SDF = sdf
    mesh.node_map = Dict{Int64,Int64}()
    mesh.cell_center_map = Dict{Tuple{Int,Int,Int},Int64}()
    mesh.X = Vector{Vector{Float64}}()
    mesh.IEN = Vector{Vector{Int64}}()
    mesh.INE = Vector{Vector{Int64}}()
    mesh.node_sdf = Vector{Float64}()
    mesh.node_hash = Dict{NTuple{3,Float64},Int64}()
    return mesh
  end
end

# Pomocná funkce pro získání SDF hodnot osmi vrcholů buňky (nezměněno)
function get_cell_sdf_values(mesh::BlockMesh, i::Int, j::Int, k::Int)
  1 <= i < mesh.nx || throw(BoundsError(mesh.SDF, i))
  1 <= j < mesh.ny || throw(BoundsError(mesh.SDF, j))
  1 <= k < mesh.nz || throw(BoundsError(mesh.SDF, k))
  return SVector{8}(
    mesh.SDF[i, j, k],       # front-bottom-left
    mesh.SDF[i+1, j, k],     # front-bottom-right
    mesh.SDF[i+1, j+1, k],   # front-top-right
    mesh.SDF[i, j+1, k],     # front-top-left
    mesh.SDF[i, j, k+1],     # back-bottom-left
    mesh.SDF[i+1, j, k+1],   # back-bottom-right
    mesh.SDF[i+1, j+1, k+1], # back-top-right
    mesh.SDF[i, j+1, k+1]    # back-top-left
  )
end

# ----------------------------
# NOVÉ: Definice Schlafli orthoscheme
# ----------------------------
# Tetraedrická konektivita – používáme pouze osmi vrcholů buňky
const schlafli_tet_connectivity = [
  [1, 2, 3, 7],  # Cesta 1: x, y, z
  [1, 6, 2, 7],  # Cesta 2: x, z, y
  [1, 3, 4, 7],  # Cesta 3: y, x, z
  [1, 4, 8, 7],  # Cesta 4: y, z, x
  [1, 5, 6, 7],  # Cesta 5: z, x, y
  [1, 8, 5, 7]   # Cesta 6: z, y, x
]

# ----------------------------
# Původní A15 schéma
# ----------------------------
include("A15_scheme.jl")  # Obsahuje definice tile_ref a tetra_connectivity pro A15

# ----------------------------
# Kvantizace souřadnic – stejná funkce jako dříve
# ----------------------------
function quantize(p::Vector{Float64}, tol::Float64=TOL)
  return (round(p[1] / tol)*tol, round(p[2] / tol)*tol, round(p[3] / tol)*tol)
end

# ----------------------------
# Funkce diskretizace buňky pomocí A15 schématu (nezměněno)
# ----------------------------
function process_cell_A15!(mesh::BlockMesh, i::Int, j::Int, k::Int)
  sdf_values = get_cell_sdf_values(mesh, i, j, k)
  if !any(x -> x >= 0, sdf_values)
    return
  end

  v000 = mesh.grid[i,   j,   k]
  v111 = mesh.grid[i+1, j+1, k+1]
  vmins = v000
  vmaxs = v111
  
  cell_sdf = get_cell_sdf_values(mesh, i, j, k)
  f000, f100, f110, f010 = cell_sdf[1], cell_sdf[2], cell_sdf[3], cell_sdf[4]
  f001, f101, f111, f011 = cell_sdf[5], cell_sdf[6], cell_sdf[7], cell_sdf[8]
  
  local_mapping = Dict{Int, Int}()
  
  for li in 1:length(tile_ref)
    ξηζ = tile_ref[li] / 4.0   # Normalizované souřadnice
    ξ, η, ζ = ξηζ[1], ξηζ[2], ξηζ[3]
    p = [ vmins[d] + (vmaxs[d] - vmins[d]) * ξηζ[d] for d in 1:3 ]
    p_key = quantize(p)
    if haskey(mesh.node_hash, p_key)
      local_mapping[li] = mesh.node_hash[p_key]
    else
      push!(mesh.X, p)
      N000 = (1-ξ) * (1-η) * (1-ζ)
      N100 = ξ     * (1-η) * (1-ζ)
      N010 = (1-ξ) * η     * (1-ζ)
      N110 = ξ     * η     * (1-ζ)
      N001 = (1-ξ) * (1-η) * ζ
      N101 = ξ     * (1-η) * ζ
      N011 = (1-ξ) * η     * ζ
      N111 = ξ     * η     * ζ
      sdf_interp = N000 * f000 + N100 * f100 + N010 * f010 + N110 * f110 +
                   N001 * f001 + N101 * f101 + N011 * f011 + N111 * f111
      push!(mesh.node_sdf, sdf_interp)
      local_index = length(mesh.X)
      local_mapping[li] = local_index
      mesh.node_hash[p_key] = local_index
    end
  end
  
  for tet in tetra_connectivity
    global_tet = [ local_mapping[li] for li in tet ]
    tet_sdf = [ mesh.node_sdf[idx] for idx in global_tet ]
    if any(x -> x >= 0, tet_sdf)
      push!(mesh.IEN, global_tet)
    end
  end
end

# ----------------------------
# NOVÉ: Funkce diskretizace buňky pomocí Schlafli orthoscheme
# ----------------------------
function process_cell_Schlafli!(mesh::BlockMesh, i::Int, j::Int, k::Int)
  # Získáme SDF hodnoty osmi vrcholů buňky
  sdf_values = get_cell_sdf_values(mesh, i, j, k)
  if !any(x -> x >= 0, sdf_values)
    return
  end
  
  # Vytvoříme lokální mapování přímo pro 8 uzlů buňky
  local_mapping = Dict{Int,Int}()
  cell_nodes = [
    mesh.grid[i,   j,   k],     # Node 1: front-bottom-left
    mesh.grid[i+1, j,   k],       # Node 2: front-bottom-right
    mesh.grid[i+1, j+1, k],       # Node 3: front-top-right
    mesh.grid[i,   j+1, k],       # Node 4: front-top-left
    mesh.grid[i,   j,   k+1],     # Node 5: back-bottom-left
    mesh.grid[i+1, j,   k+1],     # Node 6: back-bottom-right
    mesh.grid[i+1, j+1, k+1],     # Node 7: back-top-right
    mesh.grid[i,   j+1, k+1]      # Node 8: back-top-left
  ]
  
  # Pořadí SDF hodnot odpovídá pořadí vrcholů (viz get_cell_sdf_values)
  for li in 1:8
    p = cell_nodes[li]
    p_key = quantize(p)
    if haskey(mesh.node_hash, p_key)
      local_mapping[li] = mesh.node_hash[p_key]
    else
      push!(mesh.X, p)
      push!(mesh.node_sdf, sdf_values[li])
      local_index = length(mesh.X)
      local_mapping[li] = local_index
      mesh.node_hash[p_key] = local_index
    end
  end
  
  # Sestavíme tetraedry podle schématu Schlafli
  for tet in schlafli_tet_connectivity
    global_tet = [ local_mapping[li] for li in tet ]
    tet_sdf = [ mesh.node_sdf[idx] for idx in global_tet ]
    if any(x -> x >= 0, tet_sdf)
      push!(mesh.IEN, global_tet)
    end
  end
end

# ----------------------------
# Funkce pro sloučení duplicitních uzlů po tvorbě sítě (nezměněno)
# ----------------------------
function merge_duplicate_nodes!(mesh::BlockMesh, tol::Float64=TOL)
  new_nodes = Vector{Vector{Float64}}()
  new_node_sdf = Vector{Float64}()
  node_map = Dict{Int, Int}()
  global_hash = Dict{NTuple{3,Float64}, Int}()
  for i in 1:length(mesh.X)
    p = mesh.X[i]
    p_key = quantize(p, tol)
    if haskey(global_hash, p_key)
      node_map[i] = global_hash[p_key]
    else
      push!(new_nodes, p)
      push!(new_node_sdf, mesh.node_sdf[i])
      new_index = length(new_nodes)
      node_map[i] = new_index
      global_hash[p_key] = new_index
    end
  end
  for i in 1:length(mesh.IEN)
    mesh.IEN[i] = [ node_map[old] for old in mesh.IEN[i] ]
  end
  mesh.X = new_nodes
  mesh.node_sdf = new_node_sdf
  mesh.node_map = node_map
  @info "Po slučování duplicit: $(length(mesh.X)) uzlů"
end

# ----------------------------
# Vyčištění nepoužitých uzlů a přeindexování konektivity (nezměněno)
# ----------------------------
function cleanup_unused_nodes!(mesh::BlockMesh)
  @info "Počet uzlů před vyčištěním: $(length(mesh.X))"
  used_nodes = Set{Int64}()
  for element in mesh.IEN
    union!(used_nodes, element)
  end
  new_node_map = Dict{Int64,Int64}()
  new_X = Vector{Vector{Float64}}()
  new_node_sdf = Vector{Float64}()
  sorted_used = sort(collect(used_nodes))
  for (new_id, old_id) in enumerate(sorted_used)
    new_node_map[old_id] = new_id
    push!(new_X, mesh.X[old_id])
    push!(new_node_sdf, mesh.node_sdf[old_id])
  end
  for i in 1:length(mesh.IEN)
    mesh.IEN[i] = [ new_node_map[old_id] for old_id in mesh.IEN[i] ]
  end
  mesh.X = new_X
  mesh.node_sdf = new_node_sdf
  mesh.node_map = new_node_map
  @info "Počet uzlů po vyčištění: $(length(mesh.X))"
end

# ----------------------------
# Vytvoření inverzní konektivity (nezměněno)
# ----------------------------
function create_INE!(mesh::BlockMesh)
  mesh.INE = [Vector{Int64}() for _ in 1:length(mesh.X)]
  for (elem_id, element) in enumerate(mesh.IEN)
    for node_id in element
      push!(mesh.INE[node_id], elem_id)
    end
  end
  return mesh
end

# ----------------------------
# Upravená funkce generování sítě s možností volby schématu
# ----------------------------
function generate_mesh!(mesh::BlockMesh, scheme::String)
  empty!(mesh.X)
  empty!(mesh.IEN)
  empty!(mesh.node_sdf)
  empty!(mesh.node_map)
  empty!(mesh.cell_center_map)
  empty!(mesh.node_hash)
  
  for i in 1:mesh.nx-1
    for j in 1:mesh.ny-1
      for k in 1:mesh.nz-1
        if scheme == "A15"
          process_cell_A15!(mesh, i, j, k)
        elseif scheme == "Schlafli"
          process_cell_Schlafli!(mesh, i, j, k)
        else
          error("Neznámé schéma: $scheme")
        end
      end
    end
  end
  
  cleanup_unused_nodes!(mesh)
  merge_duplicate_nodes!(mesh, TOL)
  
  create_INE!(mesh)
  @info "Vytvořeno $(length(mesh.X)) uzlů a $(length(mesh.IEN)) tetraedrů"
end

# ----------------------------
# Export do VTK (nezměněno)
# ----------------------------
function export_mesh_vtk(mesh::BlockMesh, filename::String)
  npoints = length(mesh.X)
  points = zeros(Float64, 3, npoints)
  for i in 1:npoints
    points[:, i] = mesh.X[i]
  end
  cells = [MeshCell(VTKCellTypes.VTK_TETRA, tet) for tet in mesh.IEN]
  vtkfile = vtk_grid(filename, points, cells)
  vtk_point_data(vtkfile, mesh.node_sdf, "sdf")
  vtk_save(vtkfile)
end

# ----------------------------
# Hlavní běh – vytvoříme síť a exportujeme ji do souboru
# ----------------------------
mesh = BlockMesh()
# Volba schématu: "A15" nebo "Schlafli"
@time generate_mesh!(mesh, "Schlafli")
export_mesh_vtk(mesh, "block-mesh.vtu")


mesh.INE
mesh.IEN
mesh.X
# mesh.X 282703-element
# mesh.IEN 475215-element