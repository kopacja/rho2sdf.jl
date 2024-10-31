using LinearAlgebra
using WriteVTK
using StaticArrays
using JLD2
using JuMP
import Ipopt

mutable struct BlockMesh
  nx::Int
  ny::Int
  nz::Int
  grid::Array{Vector{Float64},3}   # 3D pole vektorů souřadnic uzlů
  SDF::Array{Float64,3}             # SDF hodnoty
  X::Vector{Vector{Float64}}        # Souřadnice použitých uzlů
  IEN::Vector{Vector{Int64}}        # Konektivita elementů
  INE::Vector{Vector{Int64}}       # Inverzní konektivita (Node Element)
  node_sdf::Vector{Float64}         # SDF hodnoty pro použité uzly
  node_map::Dict{Int64,Int64}      # původní_index => nový_index

  function BlockMesh()
    # Načtení dat
    @load "data/Z_block_FineGrid.jld2" fine_grid
    # @load "src/ImplicitDomainMeshing/data/Z_block_FineGrid.jld2" fine_grid
    @load "data/Z_block_FineSDF.jld2" fine_sdf
    # @load "src/ImplicitDomainMeshing/data/Z_block_FineSDF.jld2" fine_sdf

    # Konverze Float32 na Float64 pro konzistenci
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
    mesh.X = Vector{Vector{Float64}}()
    mesh.IEN = Vector{Vector{Int64}}()
    mesh.INE = Vector{Vector{Int64}}()
    mesh.node_sdf = Vector{Float64}()
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

# Zpracování jedné buňky a vytvoření tetraedrů
function process_cell!(mesh::BlockMesh, i::Int, j::Int, k::Int)
  # Získání SDF hodnot
  sdf_values = get_cell_sdf_values(mesh, i, j, k)

  # Kontrola, zda buňka protíná nulovou hladinu (má alespoň jeden kladný a jeden záporný uzel)
  has_positive = any(>=(0), sdf_values)
  has_negative = any(<(0), sdf_values)

  # Pokud buňka nemá žádný kladný uzel, přeskočíme ji
  if !has_positive
    return
  end

  # Získání/vytvoření indexů vrcholů buňky
  node_ids = [
    get_or_create_node!(mesh, i, j, k),      # front-bottom-left
    get_or_create_node!(mesh, i + 1, j, k),    # front-bottom-right
    get_or_create_node!(mesh, i + 1, j + 1, k),  # front-top-right
    get_or_create_node!(mesh, i, j + 1, k),    # front-top-left
    get_or_create_node!(mesh, i, j, k + 1),    # back-bottom-left
    get_or_create_node!(mesh, i + 1, j, k + 1),  # back-bottom-right
    get_or_create_node!(mesh, i + 1, j + 1, k + 1),# back-top-right
    get_or_create_node!(mesh, i, j + 1, k + 1)   # back-top-left
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
    # Získání SDF hodnot pro vrcholy tetraedru
    tet_sdf_values = [sdf_values[i] for i in tet]

    # Přidání tetraedru pokud:
    # 1. Všechny uzly jsou kladné (vnitřní element)
    # 2. Element je přechodový a má alespoň jeden kladný uzel
    if any(>=(0), tet_sdf_values)
      push!(mesh.IEN, [node_ids[i] for i in tet])
    end
  end
end

# Funkce pro vyčištění nepoužitých uzlů
function cleanup_unused_nodes!(mesh::BlockMesh)
  # Vytvoření množiny použitých uzlů
  used_nodes = Set{Int64}()
  for element in mesh.IEN
    union!(used_nodes, element)
  end

  # Vytvoření nového mapování pro použité uzly
  new_node_map = Dict{Int64,Int64}()
  new_X = Vector{Vector{Float64}}()
  new_node_sdf = Vector{Float64}()

  # Přemapování uzlů a odstranění nepoužitých
  for new_id in 1:length(used_nodes)
    old_id = sort(collect(used_nodes))[new_id]
    new_node_map[old_id] = new_id
    push!(new_X, mesh.X[old_id])
    push!(new_node_sdf, mesh.node_sdf[old_id])
  end

  # Aktualizace konektivity elementů
  for i in 1:length(mesh.IEN)
    mesh.IEN[i] = [new_node_map[old_id] for old_id in mesh.IEN[i]]
  end

  # Aktualizace dat sítě
  mesh.X = new_X
  mesh.node_sdf = new_node_sdf
  mesh.node_map = new_node_map

  @info "Počet uzlů po vyčištění: $(length(mesh.X))"
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

  # Vyčištění nepoužitých uzlů
  cleanup_unused_nodes!(mesh)

  create_INE!(mesh)  # Vytvoří INE v mesh struktuře

  @info "Vytvořeno $(length(mesh.X)) uzlů a $(length(mesh.IEN)) tetraedrů"
end


function find_regular_elements(mesh::BlockMesh, node_coords::Vector{Float64})
  # Najít všechny elementy (krychle) obsahující daný uzel
  connected_elements = Vector{NamedTuple{(:indices, :coords, :sdf_values),
    Tuple{Vector{Tuple{Int,Int,Int}},
      Vector{Vector{Float64}},
      Vector{Float64}}}}()

  # Najít i,j,k indexy uzlu v pravidelné síti
  # Procházíme síť a hledáme uzel se stejnými souřadnicemi
  found_indices = Vector{Tuple{Int,Int,Int}}()
  for i in 1:mesh.nx, j in 1:mesh.ny, k in 1:mesh.nz
    if mesh.grid[i, j, k] ≈ node_coords
      push!(found_indices, (i, j, k))
    end
  end

  # Pro každý nalezený výskyt uzlu najít připojené elementy
  for (i, j, k) in found_indices
    # Uzel může být součástí až 8 elementů (kromě krajních uzlů)
    possible_elements = [
      (i - 1, j - 1, k - 1), # element kde je uzel v pozici (right-top-back)
      (i - 1, j - 1, k),   # element kde je uzel v pozici (right-top-front)
      (i - 1, j, k - 1),   # element kde je uzel v pozici (right-bottom-back)
      (i - 1, j, k),     # element kde je uzel v pozici (right-bottom-front)
      (i, j - 1, k - 1),   # element kde je uzel v pozici (left-top-back)
      (i, j - 1, k),     # element kde je uzel v pozici (left-top-front)
      (i, j, k - 1),     # element kde je uzel v pozici (left-bottom-back)
      (i, j, k)        # element kde je uzel v pozici (left-bottom-front)
    ]

    # Zkontrolovat, které z možných elementů skutečně existují
    for (ei, ej, ek) in possible_elements
      if 1 <= ei < mesh.nx && 1 <= ej < mesh.ny && 1 <= ek < mesh.nz
        # Získat souřadnice všech 8 uzlů elementu
        element_coords = [
          mesh.grid[ei, ej, ek],    # front-bottom-left
          mesh.grid[ei+1, ej, ek],    # front-bottom-right
          mesh.grid[ei+1, ej+1, ek],    # front-top-right
          mesh.grid[ei, ej+1, ek],    # front-top-left
          mesh.grid[ei, ej, ek+1],  # back-bottom-left
          mesh.grid[ei+1, ej, ek+1],  # back-bottom-right
          mesh.grid[ei+1, ej+1, ek+1],  # back-top-right
          mesh.grid[ei, ej+1, ek+1]   # back-top-left
        ]

        # Získat SDF hodnoty všech 8 uzlů elementu
        element_sdf = [
          mesh.SDF[ei, ej, ek],    # front-bottom-left
          mesh.SDF[ei+1, ej, ek],    # front-bottom-right
          mesh.SDF[ei+1, ej+1, ek],    # front-top-right
          mesh.SDF[ei, ej+1, ek],    # front-top-left
          mesh.SDF[ei, ej, ek+1],  # back-bottom-left
          mesh.SDF[ei+1, ej, ek+1],  # back-bottom-right
          mesh.SDF[ei+1, ej+1, ek+1],  # back-top-right
          mesh.SDF[ei, ej+1, ek+1]   # back-top-left
        ]

        # Přidat element do výsledku
        push!(connected_elements, (
          indices=[(ei, ej, ek), (ei + 1, ej, ek), (ei + 1, ej + 1, ek), (ei, ej + 1, ek),
            (ei, ej, ek + 1), (ei + 1, ej, ek + 1), (ei + 1, ej + 1, ek + 1), (ei, ej + 1, ek + 1)],
          coords=element_coords,
          sdf_values=element_sdf
        ))
      end
    end
  end

  return connected_elements
end


function compute_coords(
  x::Vector,
  Xₑ::Matrix,
  ρₑ::Vector,
  # starting_points::Vector{Tuple{Float64, Float64, Float64}}
)
  starting_points = [
    (0.0, 0.0, 0.0),
    # (-0.5, -0.5, -0.5),
    # (0.5, -0.5, -0.5),
    # (0.5, 0.5, -0.5),
    # (-0.5, 0.5, -0.5),
    # (-0.5, -0.5, 0.5),
    # (0.5, -0.5, 0.5),
    # (0.5, 0.5, 0.5),
    # (-0.5, 0.5, 0.5),
  ]

  best_solution = nothing
  best_objective = Inf

  for (ξ₁_start, ξ₂_start, ξ₃_start) in starting_points
    model = Model(Ipopt.Optimizer)
    set_silent(model)

    set_optimizer_attribute(model, "tol", 1e-6)
    set_optimizer_attribute(model, "max_iter", 50)
    set_optimizer_attribute(model, "acceptable_tol", 1e-6)

    @variable(model, ξ₁, lower_bound = -1.0, upper_bound = 1.0, start = ξ₁_start)
    @variable(model, ξ₂, lower_bound = -1.0, upper_bound = 1.0, start = ξ₂_start)
    @variable(model, ξ₃, lower_bound = -1.0, upper_bound = 1.0, start = ξ₃_start)

    N8 = [
      -1 / 8 * (ξ₁ - 1) * (ξ₂ - 1) * (ξ₃ - 1),
      1 / 8 * (ξ₁ + 1) * (ξ₂ - 1) * (ξ₃ - 1),
      -1 / 8 * (ξ₁ + 1) * (ξ₂ + 1) * (ξ₃ - 1),
      1 / 8 * (ξ₁ - 1) * (ξ₂ + 1) * (ξ₃ - 1),
      1 / 8 * (ξ₁ - 1) * (ξ₂ - 1) * (ξ₃ + 1),
      -1 / 8 * (ξ₁ + 1) * (ξ₂ - 1) * (ξ₃ + 1),
      1 / 8 * (ξ₁ + 1) * (ξ₂ + 1) * (ξ₃ + 1),
      -1 / 8 * (ξ₁ - 1) * (ξ₂ + 1) * (ξ₃ + 1)
    ]
    @NLobjective(model, Min, sum((x[i] - sum(Xₑ[i, k] * N8[k] for k in 1:length(N8)))^2 for i in 1:length(x)))
    @NLconstraint(model, sum(ρₑ[k] * N8[k] for k in 1:length(N8)) == 0.0)
    JuMP.optimize!(model)

    current_objective = objective_value(model)
    current_solution = value.([ξ₁, ξ₂, ξ₃])

    if current_objective < best_objective
      best_solution = current_solution
      best_objective = current_objective
    end
  end

  return best_solution
end


function euclidean_distance(p1, p2)
  return sqrt(sum((p1 .- p2) .^ 2))
end

function copy_nonzero_vectors!(source, destination)
  @assert length(source) == length(destination) "Vektory musí mít stejnou délku"

  for i in 1:length(source)
    if !all(iszero, source[i])  # Pokud zdrojový vektor není nulový
      destination[i] = copy(source[i])  # Překopíruj ho do cílového vektoru
    end
    # Pokud je zdrojový vektor nulový, ponechej původní hodnotu v destination
  end
  return destination
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

  return mesh
end

# Pomocná funkce pro získání elementů připojených k uzlu
function get_connected_elements(mesh::BlockMesh, node_id::Int)
  return mesh.INE[node_id]
end


mesh = BlockMesh()
@info "Generování sítě..."
generate_mesh!(mesh)

function project_nodes_to_isocontour!(mesh::BlockMesh)
  non = length(mesh.X)
  X_new = fill([0.0, 0.0, 0.0], non)

  # Iterate through all nodes
  for i in 1:non
    if mesh.node_sdf[i] < 0.0
      x = mesh.X[i]
      # Find all connected elements
      elements = find_regular_elements(mesh, x)

      if isempty(elements)
        @warn "No elements found for node $i at position $(x)"
        continue
      end

      # Find crossing elements (those intersecting the zero level)
      crossing_indices = findall(e -> any(>=(0), e.sdf_values) && any(<(0), e.sdf_values), elements)

      if isempty(crossing_indices)
        @warn "No crossing elements found for node $i at position $(x)"
        continue
      end

      nce = length(crossing_indices)
      x_possibilities = fill([0.0, 0.0, 0.0], nce)

      # Process each crossing element
      for j in 1:nce
        element = elements[crossing_indices[j]]
        sdf_e = element.sdf_values
        Xₑ = reduce(hcat, element.coords)

        # Compute projected coordinates
        try
          ξ₁, ξ₂, ξ₃ = compute_coords(x, Xₑ, sdf_e)

          # Compute shape functions
          N = [
            -1 / 8 * (ξ₁ - 1) * (ξ₂ - 1) * (ξ₃ - 1),
            1 / 8 * (ξ₁ + 1) * (ξ₂ - 1) * (ξ₃ - 1),
            -1 / 8 * (ξ₁ + 1) * (ξ₂ + 1) * (ξ₃ - 1),
            1 / 8 * (ξ₁ - 1) * (ξ₂ + 1) * (ξ₃ - 1),
            1 / 8 * (ξ₁ - 1) * (ξ₂ - 1) * (ξ₃ + 1),
            -1 / 8 * (ξ₁ + 1) * (ξ₂ - 1) * (ξ₃ + 1),
            1 / 8 * (ξ₁ + 1) * (ξ₂ + 1) * (ξ₃ + 1),
            -1 / 8 * (ξ₁ - 1) * (ξ₂ + 1) * (ξ₃ + 1)
          ]

          x_new = Xₑ * N
          x_possibilities[j] = x_new
        catch e
          @warn "Failed to compute coordinates for node $i in element $j: $e"
          continue
        end
      end

      # Filter out zero vectors from failed computations
      valid_points = filter(p -> !all(iszero, p), x_possibilities)

      if !isempty(valid_points)
        distances = [euclidean_distance(point, x) for point in valid_points]
        closest_idx = argmin(distances)
        X_new[i] = valid_points[closest_idx]
      else
        @warn "No valid projections found for node $i at position $(x)"
      end
    end
  end

  # Update mesh coordinates
  copy_nonzero_vectors!(X_new, mesh.X)

  return mesh
end

project_nodes_to_isocontour!(mesh)

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


export_vtk(mesh, "test_geom")
