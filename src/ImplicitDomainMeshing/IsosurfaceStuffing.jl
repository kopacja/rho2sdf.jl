# Imports
using LinearAlgebra
using Statistics
using StaticArrays
using WriteVTK
using JLD2

# Parametry pro optimalizaci sítě
struct MeshOptimizationParams
  warp_threshold::Float64     # Práh pro warp vertices (default 0.3 z C++)
  max_iterations::Int         # Maximální počet iterací optimalizace
  min_improvement_ratio::Float64 # Minimální poměr zlepšení pro pokračování

  function MeshOptimizationParams(;
    warp_threshold=0.3,
    max_iterations=20,
    min_improvement_ratio=1.0
  )
    new(warp_threshold, max_iterations, min_improvement_ratio)
  end
end

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
  [2, 1, 5], [4, 1, 5], [6, 1, 5],  # +1 ke všem souřadnicím
  [1, 1, 3], [3, 2, 3], [5, 1, 3],
  [2, 1, 1], [4, 1, 1], [6, 1, 1],
  [1, 3, 6], [3, 3, 5], [5, 3, 6],
  [1, 3, 4], [3, 4, 3], [5, 3, 4],
  [1, 3, 2], [3, 3, 1], [5, 3, 2],
  [2, 5, 5], [4, 5, 5], [6, 5, 5],
  [1, 5, 3], [3, 6, 3], [5, 5, 3],
  [2, 5, 1], [4, 5, 1], [6, 5, 1]
]

const LATTICE_TETS = [
  [13, 5, 1, 4], [1, 13, 10, 11], [5, 13, 1, 11], [2, 5, 1, 11], [14, 13, 5, 11],
  [15, 2, 3, 12], [2, 15, 5, 11], [15, 2, 5, 6], [15, 14, 5, 11], [2, 15, 3, 6],
  [2, 12, 15, 11], [7, 8, 17, 5], [7, 5, 16, 4], [16, 5, 13, 4], [13, 14, 5, 16],
  [17, 16, 7, 5], [14, 16, 17, 5], [5, 18, 8, 17], [18, 5, 8, 6], [18, 14, 17, 5],
  [5, 15, 6, 18], [9, 6, 18, 8], [14, 18, 15, 5], [23, 20, 19, 14], [19, 13, 11, 10],
  [14, 20, 19, 11], [14, 13, 19, 22], [13, 14, 19, 11], [19, 23, 14, 22], [20, 14, 15, 11],
  [20, 12, 11, 15], [12, 20, 21, 15], [15, 20, 21, 24], [20, 23, 24, 14], [14, 20, 15, 24],
  [14, 13, 22, 16], [25, 23, 22, 14], [16, 25, 14, 17], [25, 16, 14, 22], [25, 26, 23, 14],
  [26, 25, 17, 14], [26, 23, 14, 24], [14, 18, 17, 26], [26, 14, 18, 24], [14, 15, 18, 24],
  [24, 27, 18, 26]
]

# Vylepšená funkce pro process_cell! používající acute tetrahedral lattice
function process_cell!(mesh::BlockMesh, i::Int, j::Int, k::Int)
  # Vytvoření lokální mřížky pro buňku
  function create_lattice_node(base_i::Int, base_j::Int, base_k::Int, lattice_pos)
    dx = mesh.grid[2, 1, 1][1] - mesh.grid[1, 1, 1][1]
    pos = mesh.grid[base_i, base_j, base_k] +
          [lattice_pos[1] * dx / 6, lattice_pos[2] * dx / 6, lattice_pos[3] * dx / 6]
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
    # Vytvoření množiny použitých uzlů a současně vytvoření seřazeného pole
    used_nodes_set = Set{Int64}()
    max_node_id = 0
    
    # První průchod - zjištění použitých uzlů a maximálního ID
    for element in mesh.IEN
        union!(used_nodes_set, element)
        max_node_id = max(max_node_id, maximum(element))
    end
    
    # Vytvoření mapování starých na nové indexy v jednom průchodu
    new_node_map = zeros(Int, max_node_id)  # Předalokace pole pro mapování
    new_id = 1
    
    # Alokace nových polí s přesnou velikostí
    n_used = length(used_nodes_set)
    new_X = Vector{Vector{Float64}}(undef, n_used)
    new_node_sdf = Vector{Float64}(undef, n_used)
    
    # Vytvoření mapování a současně kopírování dat
    for old_id in 1:max_node_id
        if old_id in used_nodes_set
            new_node_map[old_id] = new_id
            new_X[new_id] = mesh.X[old_id]
            new_node_sdf[new_id] = mesh.node_sdf[old_id]
            new_id += 1
        end
    end
    
    # Aktualizace konektivity elementů pomocí předpočítaného mapování
    @inbounds for i in 1:length(mesh.IEN)
        mesh.IEN[i] = new_node_map[mesh.IEN[i]]
    end
    
    # Aktualizace dat sítě
    mesh.X = new_X
    mesh.node_sdf = new_node_sdf
    
    # Převod pole mapování na Dict pouze pro výstup (pokud je potřeba)
    mesh.node_map = Dict(i => new_node_map[i] for i in 1:max_node_id if new_node_map[i] > 0)
    
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

# Opravená funkce pro výpočet dihedrálních úhlů tetraedru
function compute_dihedral_angles(tet_vertices::Vector{Vector{Float64}})
  # Výpočet normál stěn s lepší numerickou stabilitou
  function face_normal(v1::Vector{Float64}, v2::Vector{Float64}, v3::Vector{Float64})
    e1 = v2 .- v1
    e2 = v3 .- v1
    n = cross(e1, e2)
    norm_n = norm(n)
    return norm_n > 1e-15 ? n ./ norm_n : zeros(Float64, 3)
  end

  # Definice stěn (trojúhelníků) tetraedru
  faces = [
    (1, 2, 3), (1, 2, 4),
    (1, 3, 4), (2, 3, 4)
  ]

  # Výpočet normál stěn
  normals = Vector{Vector{Float64}}(undef, 4)
  for (i, (a, b, c)) in enumerate(faces)
    normals[i] = face_normal(
      tet_vertices[a],
      tet_vertices[b],
      tet_vertices[c]
    )
  end

  # Výpočet dihedrálních úhlů s ošetřením numerických nepřesností
  angles = Vector{Float64}(undef, 6)
  idx = 1
  for i in 1:3
    for j in (i+1):4
      # Ošetření vstupní hodnoty pro acos
      cos_angle = -dot(normals[i], normals[j])
      cos_angle = clamp(cos_angle, -1.0, 1.0)
      angles[idx] = acos(cos_angle)
      idx += 1
    end
  end

  return angles
end

# Upravená funkce pro výpočet kvality tetraedru
function compute_tet_quality(tet_vertices::Vector{Vector{Float64}})
  # Nejdřív ověříme, zda tetrahedron není degenerovaný
  function tet_volume(v1, v2, v3, v4)
    a = v2 .- v1
    b = v3 .- v1
    c = v4 .- v1
    return abs(dot(cross(a, b), c)) / 6.0
  end

  volume = tet_volume(tet_vertices...)
  if volume < 1e-15
    return 0.0
  end

  # Výpočet dihedrálních úhlů
  angles = compute_dihedral_angles(tet_vertices)

  # Kvalita je definována jako minimum ze sinů dihedrálních úhlů
  # Přidáme malou hodnotu pro numerickou stabilitu
  return maximum([1e-15, minimum(sin.(angles))])
end


# Vylepšená funkce pro projekci na izokonturu
function project_to_isosurface!(mesh::BlockMesh, vertex_idx::Int)
  max_iterations = 10
  tolerance = 1e-6
  point = copy(mesh.X[vertex_idx])

  for iter in 1:max_iterations
    sdf_val = evaluate_sdf(mesh, point)
    if abs(sdf_val) < tolerance
      break
    end

    grad = compute_sdf_gradient(mesh, point)
    grad_norm = norm(grad)

    if grad_norm < 1e-10
      break
    end

    # Newton step s omezením délky kroku
    step = -(sdf_val / grad_norm^2) * grad
    max_step = 0.1 * (mesh.grid[2, 1, 1][1] - mesh.grid[1, 1, 1][1])  # 10% velikosti buňky
    if norm(step) > max_step
      step = (max_step / norm(step)) * step
    end

    point += step
  end

  mesh.X[vertex_idx] = point
  mesh.node_sdf[vertex_idx] = evaluate_sdf(mesh, point)
end


# Vylepšená funkce pro warp_vertices s lepší projekcí na izokonturu
function warp_vertices!(mesh::BlockMesh, params::MeshOptimizationParams)
  warped_vertices = Set{Int}()

  # Nejdřív identifikujeme všechny hrany které kříží izokonturu
  for (elem_idx, tet) in enumerate(mesh.IEN)
    for i in 1:4
      for j in (i+1):4
        v1, v2 = tet[i], tet[j]
        sdf1, sdf2 = mesh.node_sdf[v1], mesh.node_sdf[v2]

        if (sdf1 < 0) ≠ (sdf2 < 0)  # Hrana protíná izokonturu
          alpha = sdf1 / (sdf1 - sdf2)

          # Kontrola warping thresholdu pro oba vrcholy
          if alpha < params.warp_threshold
            push!(warped_vertices, v1)
          elseif alpha > 1 - params.warp_threshold
            push!(warped_vertices, v2)
          end
        end
      end
    end
  end

  # Nyní provedeme projekci identifikovaných vrcholů
  for v in warped_vertices
    old_pos = copy(mesh.X[v])
    project_to_isosurface!(mesh, v)

    # Ověření kvality okolních tetraedrů po projekci
    incident_tets = [
      [mesh.X[i] for i in mesh.IEN[elem_idx]]
      for elem_idx in mesh.INE[v]
    ]

    if any(tet -> compute_tet_quality(tet) < 1e-6, incident_tets)
      # Vrátíme původní pozici pokud projekce způsobila degeneraci
      mesh.X[v] = old_pos
    end
  end
end

# Hlavní funkce pro optimalizaci sítě
function optimize_mesh!(mesh::BlockMesh, params::MeshOptimizationParams)
  # Identifikace hraničních vrcholů
  empty!(mesh.boundary_vertices)
  for (node_idx, node_sdf) in enumerate(mesh.node_sdf)
    for elem_idx in mesh.INE[node_idx]
      for other_node in mesh.IEN[elem_idx]
        if sign(node_sdf) != sign(mesh.node_sdf[other_node])
          push!(mesh.boundary_vertices, node_idx)
          break
        end
      end
      if node_idx in mesh.boundary_vertices
        break
      end
    end
  end

  # Warp vertices
  warp_vertices!(mesh, params)

  # Optimalizace kvality sítě
  for iteration in 1:params.max_iterations
    max_improvement = 0.0

    for vertex_idx in mesh.boundary_vertices
      # Výpočet původní kvality
      incident_tets = [
        [mesh.X[i] for i in mesh.IEN[elem_idx]]
        for elem_idx in mesh.INE[vertex_idx]
      ]
      old_quality = minimum(compute_tet_quality.(incident_tets))

      # Pokus o optimalizaci pozice vrcholu
      original_pos = copy(mesh.X[vertex_idx])
      project_to_isosurface!(mesh, vertex_idx)

      # Výpočet nové kvality
      new_incident_tets = [
        [mesh.X[i] for i in mesh.IEN[elem_idx]]
        for elem_idx in mesh.INE[vertex_idx]
      ]
      new_quality = minimum(compute_tet_quality.(new_incident_tets))

      # Pokud se kvalita nezlepšila, vrátíme původní pozici
      if new_quality <= old_quality
        mesh.X[vertex_idx] = original_pos
      else
        max_improvement = max(max_improvement, new_quality / old_quality)
      end
    end

    # Kontrola konvergence
    if max_improvement < params.min_improvement_ratio
      break
    end
  end
end

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

  @info "Cleanup unused nodes..."
  # Vyčištění nepoužitých uzlů
  cleanup_unused_nodes!(mesh)

  @info "Create INE..."
  # Vytvoření inverzní konektivity
  create_INE!(mesh)

  @info "Základní síť: $(length(mesh.X)) uzlů a $(length(mesh.IEN)) tetraedrů"

  # Optimalizace sítě
  @info "Optimalizace sítě..."
  params = MeshOptimizationParams()
  optimize_mesh!(mesh, params)

  @info "Optimalizace dokončena"
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

# Funkce pro výpočet statistik kvality sítě
function compute_mesh_quality_stats(mesh::BlockMesh)
  quality_metrics = Float64[]
  dihedral_angles = Float64[]

  for tet in mesh.IEN
    vertices = [mesh.X[i] for i in tet]

    # Výpočet kvality tetraedru
    quality = compute_tet_quality(vertices)
    push!(quality_metrics, quality)

    # Sběr dihedrálních úhlů
    angles = compute_dihedral_angles(vertices)
    append!(dihedral_angles, angles)
  end

  # Převod úhlů na stupně
  dihedral_angles_deg = rad2deg.(dihedral_angles)

  return Dict(
    "min_quality" => minimum(quality_metrics),
    "max_quality" => maximum(quality_metrics),
    "mean_quality" => mean(quality_metrics),
    "min_dihedral" => minimum(dihedral_angles_deg),
    "max_dihedral" => maximum(dihedral_angles_deg),
    "mean_dihedral" => mean(dihedral_angles_deg)
  )
end

# Funkce pro validaci sítě
function validate_mesh(mesh::BlockMesh)
  issues = String[]

  # Kontrola orientace tetraedrů
  function compute_tet_volume(tet_vertices)
    v1 = tet_vertices[2] - tet_vertices[1]
    v2 = tet_vertices[3] - tet_vertices[1]
    v3 = tet_vertices[4] - tet_vertices[1]
    return dot(cross(v1, v2), v3) / 6.0
  end

  for (i, tet) in enumerate(mesh.IEN)
    vertices = [mesh.X[v] for v in tet]

    # Kontrola objemu
    volume = compute_tet_volume(vertices)
    if volume <= 0
      push!(issues, "Tetrahedron $i has negative or zero volume: $volume")
    end

    # Kontrola dihedrálních úhlů
    angles = rad2deg.(compute_dihedral_angles(vertices))
    min_angle = minimum(angles)
    max_angle = maximum(angles)

    if min_angle < 10.0
      push!(issues, "Tetrahedron $i has small dihedral angle: $min_angle degrees")
    end
    if max_angle > 170.0
      push!(issues, "Tetrahedron $i has large dihedral angle: $max_angle degrees")
    end
  end

  return issues
end

# Pomocná funkce pro výpis statistik sítě
function print_mesh_statistics(mesh::BlockMesh)
  stats = compute_mesh_quality_stats(mesh)

  println("\nMesh Statistics:")
  println("----------------")
  println("Number of vertices: ", length(mesh.X))
  println("Number of tetrahedra: ", length(mesh.IEN))
  println("Number of boundary vertices: ", length(mesh.boundary_vertices))
  println("\nQuality Metrics:")
  println("Minimum quality: ", round(stats["min_quality"], digits=4))
  println("Maximum quality: ", round(stats["max_quality"], digits=4))
  println("Mean quality: ", round(stats["mean_quality"], digits=4))
  println("\nDihedral Angles:")
  println("Minimum angle: ", round(stats["min_dihedral"], digits=2), "°")
  println("Maximum angle: ", round(stats["max_dihedral"], digits=2), "°")
  println("Mean angle: ", round(stats["mean_dihedral"], digits=2), "°")

  # Validace sítě
  issues = validate_mesh(mesh)
  if !isempty(issues)
    println("\nMesh Issues:")
    for issue in issues
      println("- ", issue)
    end
  end
end

# Příklad použití
function main()
  # Vytvoření instance sítě
  mesh = BlockMesh()

  # Generování a optimalizace sítě
  generate_mesh!(mesh)

  # Výpis statistik
  print_mesh_statistics(mesh)

  # Export do VTK pro vizualizaci
  export_to_vtk(mesh, "optimized_mesh")
end

# Spuštění
main()
