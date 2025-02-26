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

# ----------------------------
# Helper function: Get SDF values at the 8 corners of a cell (unchanged)
# ----------------------------
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
# New: Definition of Schlafli orthoscheme connectivity (unchanged)
# ----------------------------
const schlafli_tet_connectivity = [
  [1, 2, 3, 7],  # Path 1: x, y, z
  [1, 6, 2, 7],  # Path 2: x, z, y
  [1, 3, 4, 7],  # Path 3: y, x, z
  [1, 4, 8, 7],  # Path 4: y, z, x
  [1, 5, 6, 7],  # Path 5: z, x, y
  [1, 8, 5, 7]   # Path 6: z, y, x
]

# ----------------------------
# Include original A15 scheme (contains definitions tile_ref and tetra_connectivity for A15)
# ----------------------------
include("A15_scheme.jl")  # This file provides tile_ref and tetra_connectivity for A15 scheme

# ----------------------------
# Quantization function – unchanged
# ----------------------------
function quantize(p::SVector{3,Float64}, tol::Float64)
  return (round(p[1] / tol) * tol, round(p[2] / tol) * tol, round(p[3] / tol) * tol)
end

# ----------------------------
# Helper function to compute shape functions for a hex8 element
# Accepts local coordinates in [0,1] (normalized tile_ref) and returns shape functions
# computed with standard transformation to [-1,1]
# ----------------------------
function shape_functions(ξηζ::SVector{3,Float64})::SVector{8,Float64}
  # Transform local coordinates from [0,1] to [-1,1]
  ξ = 2 * ξηζ[1] - 1.0
  η = 2 * ξηζ[2] - 1.0
  ζ = 2 * ξηζ[3] - 1.0
  coef = 1 / 8.0
  return @SVector [
    coef * (1 - ξ) * (1 - η) * (1 - ζ),
    coef * (1 + ξ) * (1 - η) * (1 - ζ),
    coef * (1 + ξ) * (1 + η) * (1 - ζ),
    coef * (1 - ξ) * (1 + η) * (1 - ζ),
    coef * (1 - ξ) * (1 - η) * (1 + ζ),
    coef * (1 + ξ) * (1 - η) * (1 + ζ),
    coef * (1 + ξ) * (1 + η) * (1 + ζ),
    coef * (1 - ξ) * (1 + η) * (1 + ζ)
  ]
end

# ----------------------------
# Pomocná funkce: Vyhodnocení SDF pomocí trilineární interpolace
# ----------------------------
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

# ----------------------------
# Function for discretizing a cell using A15 scheme (unchanged logic, but with refactored shape function and slight optimizations)
# ----------------------------
function process_cell_A15!(mesh::BlockMesh, i::Int, j::Int, k::Int)
  tol = mesh.grid_tol
  sdf_values = get_cell_sdf_values(mesh, i, j, k)
  if !any(x -> x >= 0, sdf_values)
    return
  end

  # Retrieve min and max corners of the cell
  v000 = mesh.grid[i, j, k]
  v111 = mesh.grid[i+1, j+1, k+1]
  vmins = v000
  vmaxs = v111

  # Precompute differences for coordinate interpolation
  Δ = vmaxs .- vmins
  local_mapping = Dict{Int,Int}()

  @inbounds for li in 1:length(tile_ref)
    # tile_ref is assumed to be defined in A15_scheme.jl and normalized (in [0,1]*4 originally)
    local_coord = SVector{3,Float64}(tile_ref[li] ./ 4.0)  # normalized coordinates in [0,1]
    # Compute physical point using linear interpolation
    p = vmins .+ Δ .* local_coord
    p = SVector{3,Float64}(p)  # ensure static vector type
    p_key = quantize(p, tol)
    if haskey(mesh.node_hash, p_key)
      local_mapping[li] = mesh.node_hash[p_key]
    else
      push!(mesh.X, p)

      sdf_aprox = eval_sdf(mesh, p)
      push!(mesh.node_sdf, sdf_aprox)
      local_index = length(mesh.X)
      local_mapping[li] = local_index
      mesh.node_hash[p_key] = local_index
    end
  end

  # Process tetrahedral connectivity from A15 scheme
  @inbounds for tet in tetra_connectivity
    global_tet = [local_mapping[li] for li in tet]
    tet_sdf = [mesh.node_sdf[idx] for idx in global_tet]
    if any(x -> x >= 0, tet_sdf)
      push!(mesh.IEN, global_tet)
    end
  end
end

# ----------------------------
# Function for discretizing a cell using Schlafli orthoscheme (unchanged logic, only minor type annotation changes)
# ----------------------------
function process_cell_Schlafli!(mesh::BlockMesh, i::Int, j::Int, k::Int)
  # Get SDF values at the 8 corners of the cell
  sdf_values = get_cell_sdf_values(mesh, i, j, k)
  if !any(x -> x >= 0, sdf_values)
    return
  end

  local_mapping = Dict{Int,Int}()
  # Define cell nodes as SVectors from grid
  cell_nodes = [
    mesh.grid[i, j, k],     # Node 1: front-bottom-left
    mesh.grid[i+1, j, k],       # Node 2: front-bottom-right
    mesh.grid[i+1, j+1, k],       # Node 3: front-top-right
    mesh.grid[i, j+1, k],       # Node 4: front-top-left
    mesh.grid[i, j, k+1],     # Node 5: back-bottom-left
    mesh.grid[i+1, j, k+1],     # Node 6: back-bottom-right
    mesh.grid[i+1, j+1, k+1],     # Node 7: back-top-right
    mesh.grid[i, j+1, k+1]      # Node 8: back-top-left
  ]

  @inbounds for li in 1:8
    p = cell_nodes[li]
    p_key = quantize(p, tol)
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

  # Construct tetrahedra according to Schlafli scheme
  @inbounds for tet in schlafli_tet_connectivity
    global_tet = [local_mapping[li] for li in tet]
    tet_sdf = [mesh.node_sdf[idx] for idx in global_tet]
    if any(x -> x >= 0, tet_sdf)
      push!(mesh.IEN, global_tet)
    end
  end
end

# ----------------------------
# Merge duplicate nodes after mesh generation (unchanged logic, only type annotations updated)
# ----------------------------
function merge_duplicate_nodes!(mesh::BlockMesh)
  tol = mesh.grid_tol
  new_nodes = Vector{SVector{3,Float64}}()
  new_node_sdf = Vector{Float64}()
  node_map = Dict{Int,Int}()
  global_hash = Dict{NTuple{3,Float64},Int}()
  @inbounds for i in 1:length(mesh.X)
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
  @inbounds for i in 1:length(mesh.IEN)
    mesh.IEN[i] = [node_map[old] for old in mesh.IEN[i]]
  end
  mesh.X = new_nodes
  mesh.node_sdf = new_node_sdf
  mesh.node_map = node_map
  @info "After merging duplicates: $(length(mesh.X)) nodes"
end

# ----------------------------
# Cleanup unused nodes and reindex connectivity (unchanged logic)
# ----------------------------
function cleanup_unused_nodes!(mesh::BlockMesh)
  @info "Number of nodes before cleanup: $(length(mesh.X))"
  used_nodes = Set{Int64}()
  @inbounds for element in mesh.IEN
    union!(used_nodes, element)
  end
  new_node_map = Dict{Int64,Int64}()
  new_coords = Vector{SVector{3,Float64}}()
  new_node_sdf = Vector{Float64}()
  sorted_used = sort(collect(used_nodes))
  for (new_id, old_id) in enumerate(sorted_used)
    new_node_map[old_id] = new_id
    push!(new_coords, mesh.X[old_id])
    push!(new_node_sdf, mesh.node_sdf[old_id])
  end
  @inbounds for i in 1:length(mesh.IEN)
    mesh.IEN[i] = [new_node_map[old_id] for old_id in mesh.IEN[i]]
  end
  mesh.X = new_coords
  mesh.node_sdf = new_node_sdf
  mesh.node_map = new_node_map
  @info "Number of nodes after cleanup: $(length(mesh.X))"
end

# ----------------------------
# Create inverse connectivity (unchanged logic)
# ----------------------------
function create_INE!(mesh::BlockMesh)
  mesh.INE = [Vector{Int64}() for _ in 1:length(mesh.X)]
  @inbounds for (elem_id, element) in enumerate(mesh.IEN)
    for node_id in element
      push!(mesh.INE[node_id], elem_id)
    end
  end
  return mesh
end

# ----------------------------
# Modified mesh generation function with scheme selection
# ----------------------------
function generate_mesh!(mesh::BlockMesh, scheme::String)
  empty!(mesh.X)
  empty!(mesh.IEN)
  empty!(mesh.node_sdf)
  empty!(mesh.node_map)
  empty!(mesh.cell_center_map)
  empty!(mesh.node_hash)

  @inbounds for i in 1:mesh.nx-1
    for j in 1:mesh.ny-1
      for k in 1:mesh.nz-1
        if scheme == "A15"
          process_cell_A15!(mesh, i, j, k)
        elseif scheme == "Schlafli"
          process_cell_Schlafli!(mesh, i, j, k)
        else
          error("Unknown scheme: $scheme")
        end
      end
    end
  end

  cleanup_unused_nodes!(mesh)
  merge_duplicate_nodes!(mesh)

  create_INE!(mesh)
  @info "Mesh created: $(length(mesh.X)) nodes and $(length(mesh.IEN)) tetrahedra"
end

# ----------------------------
# Export mesh to VTK (unchanged logic, minor type adjustments)
# ----------------------------
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

# ----------------------------
# Pomocná funkce: Aproximace gradientu SDF v bodě p pomocí centrálních diferencí
# ----------------------------
function approximate_gradient(mesh::BlockMesh, p::SVector{3,Float64}; h::Float64=1e-3)
  dx = SVector{3,Float64}(h, 0.0, 0.0)
  dy = SVector{3,Float64}(0.0, h, 0.0)
  dz = SVector{3,Float64}(0.0, 0.0, h)
  df_dx = (eval_sdf(mesh, p + dx) - eval_sdf(mesh, p - dx)) / (2 * h)
  df_dy = (eval_sdf(mesh, p + dy) - eval_sdf(mesh, p - dy)) / (2 * h)
  df_dz = (eval_sdf(mesh, p + dz) - eval_sdf(mesh, p - dz)) / (2 * h)
  return SVector{3,Float64}(df_dx, df_dy, df_dz)
end

# Funkce pro výpočet délky nejdelší hrany mezi uzly ve všech tetraedrech
# TODO: Stačí vzít počáteční diskretizaci a vzít pouze jeden element (pro zrychlení)
function longest_edge(mesh::BlockMesh)
  # Pre-allocate maximum length with type stability
  max_length = zero(eltype(mesh.X[1]))

  # Vectorize operations by using broadcast
  for tet in mesh.IEN
    # Use array views for better memory efficiency
    nodes = @view mesh.X[tet]
    # Calculate all edge lengths at once using comprehension
    lengths = [norm(nodes[i] - nodes[j]) for i in 1:3 for j in i+1:4]
    # Use built-in maximum function
    max_length = max(max_length, maximum(lengths))
  end
  return max_length
end


# Pomocná funkce – odhad gradientu SDF v okolí uzlu
# First, let's modify compute_gradient to work with both node indices and positions
function compute_gradient(mesh::BlockMesh, p::SVector{3,Float64}; δ::Float64=1e-3)
  # Pre-allocate unit vectors as static vectors for better performance
  unit_vectors = (
    SVector{3,Float64}(1.0, 0.0, 0.0),
    SVector{3,Float64}(0.0, 1.0, 0.0),
    SVector{3,Float64}(0.0, 0.0, 1.0)
  )

  # Use tuple comprehension for better compile-time optimization
  grad = ntuple(3) do d
    unit_vec = unit_vectors[d]
    # Compute central difference
    (eval_sdf(mesh, p + δ * unit_vec) - eval_sdf(mesh, p - δ * unit_vec)) / (2δ)
  end

  return SVector{3,Float64}(grad)
end

# Add method for node index input
function compute_gradient(mesh::BlockMesh, node_index::Int; δ::Float64=1e-3)
  return compute_gradient(mesh, mesh.X[node_index]; δ=δ)
end

# Update warp_node_to_isocontour! to handle positions directly
function warp_node_to_isocontour!(mesh::BlockMesh, node_index::Int, max_iter)
  tol = mesh.grid_tol
  current_position = mesh.X[node_index]

  for iter in 1:max_iter
    f = eval_sdf(mesh, current_position)

    # Early return if we're close enough to the isocontour
    abs2(f) < tol * tol && break

    grad = compute_gradient(mesh, current_position)
    norm_grad_squared = sum(abs2, grad)

    # Early return if gradient is too small
    norm_grad_squared < 1e-16 && break

    # Newton step
    dp = (f / norm_grad_squared) * grad
    current_position -= dp
  end
  current_sdf = eval_sdf(mesh, current_position)
  if abs(current_sdf) < tol
    mesh.node_sdf[node_index] = 0.
  else
    println("current_sdf: ", current_sdf)
    mesh.node_sdf[node_index] = current_sdf
  end

  mesh.X[node_index] = current_position
end

# Hlavní funkce pro warping uzlů – ordered warping
#
# Nejprve se upraví uzly s negativní SDF hodnotou (uvnitř izopovrchu) a poté uzly s kladnou hodnotou.
# Uzly se posunují směrem k nulové hladině SDF (izopovrchu) a práh posunu se počítá jako
# threshold_sdf = 0.2 * (délka nejdelší hrany tetraedru).
function warp!(mesh::BlockMesh, max_iter::Int=160)
  # Vypočítat nejdelší hranu a následně threshold pro posun
  max_edge = longest_edge(mesh)
  threshold_sdf = 0.5 * max_edge

  @info "Warping: max edge = $max_edge, threshold_sdf = $threshold_sdf"

  # První průchod: uzly s kladnou SDF hodnotou (uvnitř)
  for i in 1:length(mesh.X)
    sdf = mesh.node_sdf[i]
    if sdf > 0 && abs(sdf) < threshold_sdf
      warp_node_to_isocontour!(mesh, i, max_iter)
    end
  end

  # Druhý průchod: uzly s negativní SDF hodnotou (vně)
  for i in 1:length(mesh.X)
    sdf = mesh.node_sdf[i]
    if sdf < 0 && abs(sdf) < threshold_sdf
      warp_node_to_isocontour!(mesh, i, max_iter)
    end
  end
end

# ---------------------------------------------------
# Funkce: Aktualizace topologie meshe (mesh.X, mesh.IEN, mesh.INE)
# ---------------------------------------------------
function update_connectivity!(mesh::BlockMesh)
  @info "Aktualizuji topologii meshe: cleanup, slučování duplicit a tvorba inverzní konektivity..."
  cleanup_unused_nodes!(mesh)        # Přepočítá mesh.X, mesh.node_sdf a reindexuje mesh.IEN a mesh.node_map
  merge_duplicate_nodes!(mesh)    # Sloučí duplicitní uzly a upraví konektivitu v mesh.IEN
  create_INE!(mesh)                    # Vytvoří inverzní konektivitu (mesh.INE)
  @info "Aktualizace topologie dokončena: $(length(mesh.X)) uzlů, $(length(mesh.IEN)) tetraedrů."
end

# === Kompletní šablony podle literatury ===
"""
  apply_complete_stencil(mesh, tet)

Funkce `apply_complete_stencil` vezme daný tetraedr (se čtyřmi uzly a jejich SDF hodnotami)
a na základě signatury (počtu kladných vs. záporných) aplikuje příslušnou triangulační šablonu.
Šablony jsou definovány následovně:

1. **np == 1 (1 kladný, 3 záporné):**  
   Vypočítají se průsečíky na hranách spojujících kladný vrchol s každým záporným.
   Výsledkem je jeden nový tetraedr:  
     `[P, I1, I2, I3]`

2. **np == 3 (3 kladné, 1 záporný):**  
   Vypočítají se průsečíky na hranách spojujících záporný vrchol s každým kladným.
   Výsledkem je rozdělení do tří tetraedrů:  
     ```
     [P1, P2, I(P1,P_neg), I(P2,P_neg)]
     [P2, P3, I(P2,P_neg), I(P3,P_neg)]
     [P3, P1, I(P3,P_neg), I(P1,P_neg)]
     ```
   (v naší implementaci mírně odlišné pořadí, viz níže)

3. **np == 2 (ambiguózní případ):**  
   Nejprve se spočítají průsečíky na všech hranách spojujících kladné a záporné vrcholy.
   Těmito průsečíky vzniká čtyřúhelník. Podle literatury existují dvě varianty:
   
   - **Variant A (standardní):**  
     Po výběru optimální diagonály (na základě minimální cost a “barev” hran)
     se čtyřúhelník reuspořádá a aplikuje se následující rozdělení do tří tetraedrů:
       ```
       T1: [P1, I1, I2, I3]
       T2: [P2, I1, I3, I4]
       T3: [P1, P2, I1, I3]
       ```
       
   - **Variant B (pro uzly na hranici):**  
     Pokud některý z vrcholů leží na hranici, použije se jednodušší šablona rozdělující tetraedr na dvě části:
       ```
       T1: [P1, I1, I2, I3]
       T2: [P2, I1, I3, I4]
       ```
       
   Funkce rozhodne mezi variantami na základě testu `is_on_boundary`.

4. **np == 4 (všechny vrcholy kladné):**  
   Tetraedr se ponechá beze změny.
  
V případě, že všechny hodnoty jsou záporné, tetraedr se vynechá.
"""

# Funkce, která projde všechny tetraedry v konektivitě a aplikuje kompletní stencils
function slice_ambiguous_tetrahedra!(mesh::BlockMesh)
  new_IEN = Vector{Vector{Int64}}()
  for tet in mesh.IEN
    # new_tets = apply_complete_stencil(mesh, tet)
    new_tets = apply_stencil(mesh, tet)
    for nt in new_tets
      t_sdf = [mesh.node_sdf[i] for i in nt]
      if any(x -> x >= 0, t_sdf)
        push!(new_IEN, nt)
      end
    end
  end
  mesh.IEN = new_IEN
  @info "Po kompletním řezu tetraedrů: $(length(mesh.IEN)) tetraedrů"
end



# ----------------------------
# Main execution – create mesh and export it to file
# ----------------------------
mesh = BlockMesh()
# Choose scheme: "A15" or "Schlafli"
@time generate_mesh!(mesh, "A15")

@time warp!(mesh)
# @time warp!(mesh, 0.8)

# Nejprve aktualizujeme topologii meshe
update_connectivity!(mesh)

export_mesh_vtk(mesh, "block-mesh_warped.vtu")

slice_ambiguous_tetrahedra!(mesh)

update_connectivity!(mesh)

export_mesh_vtk(mesh, "block-mesh.vtu")

#________________________
function run_all()
  mesh = BlockMesh()
# Choose scheme: "A15" or "Schlafli"
@time generate_mesh!(mesh, "A15")

@time warp!(mesh)
# @time warp!(mesh, 0.8)

# Nejprve aktualizujeme topologii meshe
update_connectivity!(mesh)

export_mesh_vtk(mesh, "block-mesh_warped.vtu")

slice_ambiguous_tetrahedra!(mesh)

update_connectivity!(mesh)

export_mesh_vtk(mesh, "block-mesh.vtu")
end
run_all()
