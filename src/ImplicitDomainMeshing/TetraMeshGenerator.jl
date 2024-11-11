using LinearAlgebra
using WriteVTK
using Statistics
using StaticArrays
using JLD2
using NLopt
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
    # @load "data/Z_block_FineGrid.jld2" fine_grid
    @load "data/Z_chapadlo_FineGrid.jld2" fine_grid
    # @load "src/ImplicitDomainMeshing/data/Z_block_FineGrid.jld2" fine_grid
    # @load "data/Z_block_FineSDF.jld2" fine_sdf
    @load "data/Z_chapadlo_FineSDF.jld2" fine_sdf
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

  # Tetrahedra definitions according to Schläfli orthoscheme
  tet_connectivity = [
    [1, 2, 3, 7],  # Path 1: x, y, z
    [1, 2, 6, 7],  # Path 2: x, z, y
    [1, 4, 3, 7],  # Path 3: y, x, z
    [1, 4, 8, 7],  # Path 4: y, z, x
    [1, 5, 6, 7],  # Path 5: z, x, y
    [1, 5, 8, 7]   # Path 6: z, y, x
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

function find_elements_by_negative_sdf(mesh::BlockMesh, connected_elements::Vector{Int64})
  # Inicializace polí pro ukládání ID elementů
  tet_three_negative = Int64[]  # Pro elementy se 3 zápornými uzly
  tet_two_negative = Int64[]    # Pro elementy se 2 zápornými uzly

  # Procházení připojených elementů
  for elem_id in connected_elements
    # Získání SDF hodnot pro všechny uzly elementu
    tet_sdf_values = mesh.node_sdf[mesh.IEN[elem_id]]

    # Spočítání počtu záporných hodnot
    negative_count = count(x -> x < 0.0, tet_sdf_values)

    # Klasifikace elementu podle počtu záporných uzlů
    if negative_count == 3
      push!(tet_three_negative, elem_id)
    elseif negative_count == 2
      push!(tet_two_negative, elem_id)
    end
  end

  return tet_three_negative, tet_two_negative
end


function find_local_coordinates(
  Xₑ,
  xₙ
)
  starting_points::Vector{Tuple{Float64,Float64,Float64}} = [
    (0.0, 0.0, 0.0),
    (-0.5, -0.5, -0.5),
    (0.5, -0.5, -0.5),
    (0.5, 0.5, -0.5),
    (-0.5, 0.5, -0.5),
    (-0.5, -0.5, 0.5),
    (0.5, -0.5, 0.5),
    (0.5, 0.5, 0.5),
    (-0.5, 0.5, 0.5)
  ]

  function objective(ξ::Vector, grad::Vector)
    ξ₁, ξ₂, ξ₃ = ξ

    # Compute N(ξ)
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

    d¹N_dξ¹ = zeros(Float64, 8, 3)
    d¹N_dξ¹[1, :] = [
      -0.125 * (ξ₂ - 1) * (ξ₃ - 1)
      -0.125 * (ξ₁ - 1) * (ξ₃ - 1)
      -0.125 * (ξ₁ - 1) * (ξ₂ - 1)
    ]

    d¹N_dξ¹[2, :] = [
      0.125 * (ξ₂ - 1) * (ξ₃ - 1)
      0.125 * (ξ₁ + 1) * (ξ₃ - 1)
      0.125 * (ξ₁ + 1) * (ξ₂ - 1)
    ]

    d¹N_dξ¹[3, :] = [
      -0.125 * (ξ₂ + 1) * (ξ₃ - 1)
      -0.125 * (ξ₁ + 1) * (ξ₃ - 1)
      -0.125 * (ξ₁ + 1) * (ξ₂ + 1)
    ]

    d¹N_dξ¹[4, :] = [
      0.125 * (ξ₂ + 1) * (ξ₃ - 1)
      0.125 * (ξ₁ - 1) * (ξ₃ - 1)
      0.125 * (ξ₁ - 1) * (ξ₂ + 1)
    ]

    d¹N_dξ¹[5, :] = [
      0.125 * (ξ₂ - 1) * (ξ₃ + 1)
      0.125 * (ξ₁ - 1) * (ξ₃ + 1)
      0.125 * (ξ₁ - 1) * (ξ₂ - 1)
    ]

    d¹N_dξ¹[6, :] = [
      -0.125 * (ξ₂ - 1) * (ξ₃ + 1)
      -0.125 * (ξ₁ + 1) * (ξ₃ + 1)
      -0.125 * (ξ₁ + 1) * (ξ₂ - 1)
    ]

    d¹N_dξ¹[7, :] = [
      0.125 * (ξ₂ + 1) * (ξ₃ + 1)
      0.125 * (ξ₁ + 1) * (ξ₃ + 1)
      0.125 * (ξ₁ + 1) * (ξ₂ + 1)
    ]

    d¹N_dξ¹[8, :] = [
      -0.125 * (ξ₂ + 1) * (ξ₃ + 1)
      -0.125 * (ξ₁ - 1) * (ξ₃ + 1)
      -0.125 * (ξ₁ - 1) * (ξ₂ + 1)
    ]

    # Compute x(ξ)
    x = Xₑ * N  # x is a 3-element vector

    # Compute residual R = x - xₙ
    R = x - xₙ  # R is a 3-element vector

    # Compute objective function value f
    f = dot(R, R)  # Equivalent to sum(R[i]^2 for i in 1:3)

    if length(grad) > 0
      # Compute the derivatives of N with respect to ξ
      dN_dξ = zeros(Float64, 8, 3)  # 8 shape functions x 3 variables

      # Assign your computed derivatives to dN_dξ
      dN_dξ .= d¹N_dξ¹  # Ensure dN_dξ is correctly populated

      # Compute Jacobian J = Xₑ * dN_dξ
      J = Xₑ * dN_dξ  # J is 3x3

      # Compute gradient grad = 2 * Jᵗ * R
      grad .= 2 * (J' * R)
    end

    return f
  end

  best_solution = nothing
  best_objective = Inf

  for start_point in starting_points
    #   opt = Opt(:LN_COBYLA, 3)
    opt = Opt(:LD_LBFGS, 3)
    opt.lower_bounds = [-5.0, -5.0, -5.0]
    opt.upper_bounds = [5.0, 5.0, 5.0]
    opt.xtol_rel = 1e-6
    opt.maxeval = 500
    opt.min_objective = objective
    opt.maxtime = 1.0

    (minf, minx, ret) = NLopt.optimize(opt, collect(start_point))

    if minf < best_objective
      best_objective = minf
      best_solution = minx
    end
  end

  if best_solution === nothing
    return (false, [10.0, 10.0, 10.0])
  else
    return (true, best_solution)
  end
end


mesh = BlockMesh()
@info "Generování sítě..."
generate_mesh!(mesh)
# mesh.INE[5]

# connected_elements = get_connected_elements(mesh, 1) # id elements
# for j in 1:length(connected_elements)
#   mesh.node_sdf[mesh.IEN[j]] # hodnoty sdf pro uzly j-tého tetra elementu, který je součástí itého uzlu

function project_nodes_to_isocontour!(mesh::BlockMesh)
  non = length(mesh.X)
  X_new = fill([0.0, 0.0, 0.0], non)

  # Počítadla pro sledování zpracovaných případů
  processed = Dict{String,Int}()
  processed["three_negative"] = 0
  processed["two_negative"] = 0
  processed["one_negative"] = 0

  # Iterate through all nodes
  for i in 1:non
    if mesh.node_sdf[i] < 0.0
      connected_elements = get_connected_elements(mesh, i)
      x = mesh.X[i]
      elements = find_regular_elements(mesh, x)

      if isempty(elements)
        continue
      end

      # Najít elementy protínající nulovou hladinu
      crossing_indices = findall(e -> any(>=(0), e.sdf_values) && any(<(0), e.sdf_values), elements)

      if isempty(crossing_indices)
        continue
      end

      tet_three_negative, tet_two_negative = find_elements_by_negative_sdf(mesh, connected_elements)

      # Přidáno logování pro debugging
      if !isempty(tet_three_negative)
        #INFO: 3 uzly:
        processed["three_negative"] += 1
        tet_id = first(tet_three_negative)

        # Získání uzlů tetrahedra
        tet_nodes = mesh.IEN[tet_id]
        tet_sdf_values = mesh.node_sdf[tet_nodes]

        # Nalezení uzlu s kladnou SDF hodnotou
        positive_node_idx = findfirst(x -> x >= 0.0, tet_sdf_values)
        if positive_node_idx === nothing
          @warn "Nenalezen uzel s kladnou SDF hodnotou v tetrahedru $tet_id"
          continue
        end

        # Nalezení pravidelného elementu
        regular_element = nothing
        for element in elements
          contains_all_nodes = true
          for node in mesh.X[tet_nodes]
            if !any(e_node -> all(isapprox.(node, e_node)), element.coords)
              contains_all_nodes = false
              break
            end
          end
          if contains_all_nodes
            regular_element = element
            break
          end
        end

        if regular_element === nothing
          @warn "Nenalezen pravidelný element pro tetrahedr $tet_id"
          continue
        end

        # Bod s kladnou SDF hodnotou
        positive_point = mesh.X[tet_nodes[positive_node_idx]]

        # Iterativní hledání bodu na nulové hladině
        max_iterations = 30
        tolerance = 1e-4
        current_point = x
        Xₑ = reduce(hcat, regular_element.coords)
        sdf_e = regular_element.sdf_values

        for iter in 1:max_iterations
          # Výpočet směrového vektoru
          direction = positive_point - current_point
          # direction = - positive_point + current_point

          # Výpočet nového bodu v polovině vzdálenosti
          new_point = current_point + 0.5 * direction

          # Nalezení lokálních souřadnic pro nový bod
          try
            (_, local_coords) = find_local_coordinates(Xₑ, new_point)
            ξ₁, ξ₂, ξ₃ = local_coords

            # Výpočet tvarových funkcí
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

            # Výpočet SDF hodnoty v novém bodě
            sdf_new = dot(sdf_e, N)

            # Kontrola konvergence
            if abs(sdf_new) < tolerance
              X_new[i] = new_point
              break
            end

            # Aktualizace bodů pro další iteraci
            if sdf_new < 0
              current_point = new_point
            else
              positive_point = new_point
            end

            # Pokud jsme v poslední iteraci a stále jsme nenašli bod
            if iter == max_iterations
              X_new[i] = new_point  # Použijeme poslední vypočtený bod
              @warn "Nedosažena požadovaná přesnost pro uzel $i po $max_iterations iteracích"
            end

          catch e
            @warn "Chyba při výpočtu lokálních souřadnic pro uzel $i: $e"
            X_new[i] = current_point  # Použijeme poslední platný bod
            break
          end
        end
      elseif !isempty(tet_two_negative)
        #INFO: 2 uzly:
        processed["two_negative"] += 1
        println("Zpracovávám konfiguraci se dvěma zápornými uzly pro uzel $i")
        tet_id = first(tet_two_negative)

        # Získání uzlů tetrahedru
        tet_nodes = mesh.IEN[tet_id]
        tet_sdf_values = mesh.node_sdf[tet_nodes]

        # Najít uzly s kladnou a zápornou SDF
        positive_indices = findall(x -> x >= 0.0, tet_sdf_values)
        negative_indices = findall(x -> x < 0.0, tet_sdf_values)

        # Použít průměr interpolovaných bodů
        x_new = zeros(3)
        count = 0

        for neg_idx in negative_indices
          for pos_idx in positive_indices
            neg_point = mesh.X[tet_nodes[neg_idx]]
            pos_point = mesh.X[tet_nodes[pos_idx]]
            neg_sdf = tet_sdf_values[neg_idx]
            pos_sdf = tet_sdf_values[pos_idx]

            # Lineární interpolace
            t = neg_sdf / (neg_sdf - pos_sdf)
            interpolated_point = neg_point + t * (pos_point - neg_point)

            x_new += interpolated_point
            count += 1
          end
        end

        X_new[i] = x_new / count
      else
        processed["one_negative"] += 1
        #INFO: 1 uzel:
        println("Zpracovávám konfiguraci s jedním záporným uzlem pro uzel $i")
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
  end

  # Výpis statistik
  println("\nStatistiky zpracování:")
  println("Zpracováno případů s třemi zápornými uzly: $(processed["three_negative"])")
  println("Zpracováno případů s dvěma zápornými uzly: $(processed["two_negative"])")
  println("Zpracováno případů s jedním záporným uzlem: $(processed["one_negative"])")


  # Update mesh coordinates
  copy_nonzero_vectors!(X_new, mesh.X)

  return mesh
end

project_nodes_to_isocontour!(mesh)

# Výpočet objemu tetrahedru
function tetrahedron_volume(vertices::Vector{Vector{Float64}})
    # Matrix pro výpočet objemu
    a = vertices[2] - vertices[1]
    b = vertices[3] - vertices[1]
    c = vertices[4] - vertices[1]
    
    # Objem = 1/6 * |det(a b c)|
    return abs(dot(a, cross(b, c))) / 6.0
end

# Výpočet plochy trojúhelníku
function triangle_area(v1::Vector{Float64}, v2::Vector{Float64}, v3::Vector{Float64})
    # Výpočet plochy pomocí vektorového součinu
    return norm(cross(v2 - v1, v3 - v1)) / 2.0
end

# Výpočet délky hrany
function edge_length(v1::Vector{Float64}, v2::Vector{Float64})
    return norm(v2 - v1)
end

# Výpočet kvality elementu pomocí poměru objemu k součtu ploch stěn
function volume_surface_quality(vertices::Vector{Vector{Float64}})
    # Výpočet objemu
    volume = tetrahedron_volume(vertices)
    
    # Výpočet ploch všech stěn
    areas = [
        triangle_area(vertices[1], vertices[2], vertices[3]),  # Face 123
        triangle_area(vertices[1], vertices[2], vertices[4]),  # Face 124
        triangle_area(vertices[1], vertices[3], vertices[4]),  # Face 134
        triangle_area(vertices[2], vertices[3], vertices[4])   # Face 234
    ]
    
    total_surface = sum(areas)
    
    # Normalizovaný poměr (0 = degenerovaný, 1 = pravidelný tetrahedr)
    actual_ratio = (volume / total_surface) * 36.0  # 36 je normalizační faktor
    
    return max(0.0, min(1.0, actual_ratio))
end

# Výpočet kvality elementu pomocí poměru nejkratší a nejdelší hrany
function aspect_ratio_quality(vertices::Vector{Vector{Float64}})
    # Výpočet délky všech hran
    edges = [
        edge_length(vertices[1], vertices[2]),  # Edge 12
        edge_length(vertices[1], vertices[3]),  # Edge 13
        edge_length(vertices[1], vertices[4]),  # Edge 14
        edge_length(vertices[2], vertices[3]),  # Edge 23
        edge_length(vertices[2], vertices[4]),  # Edge 24
        edge_length(vertices[3], vertices[4])   # Edge 34
    ]
    
    min_edge = minimum(edges)
    max_edge = maximum(edges)
    
    # Poměr nejkratší a nejdelší hrany (0 = degenerovaný, 1 = všechny hrany stejně dlouhé)
    return min_edge / max_edge
end

# Výpočet celkové kvality elementu jako kombinace obou metrik
function calculate_element_quality(vertices::Vector{Vector{Float64}})
    vol_surf_q = volume_surface_quality(vertices)
    aspect_q = aspect_ratio_quality(vertices)
    
    # Vážený průměr obou metrik (můžete upravit váhy podle potřeby)
    return 0.5 * (vol_surf_q + aspect_q)
end

# Výpočet kvality všech elementů v síti
function calculate_mesh_quality(mesh::BlockMesh)
    qualities = zeros(Float64, length(mesh.IEN))
    
    for (i, element) in enumerate(mesh.IEN)
        # Získání souřadnic vrcholů elementu
        vertices = [mesh.X[node_id] for node_id in element]
        qualities[i] = calculate_element_quality(vertices)
    end
    
    return qualities
end

# Upravená funkce pro export do VTK
function export_vtk_with_quality(mesh::BlockMesh, filename::String)
    if !endswith(filename, ".vtu")
        filename = filename * ".vtu"
    end

    # Výpočet kvality elementů
    element_qualities = calculate_mesh_quality(mesh)
    
    # Převod uzlů do formátu pro VTK (3×N matice)
    points = reduce(hcat, mesh.X)

    # Vytvoření VTK buněk
    cells = [MeshCell(VTKCellTypes.VTK_TETRA, tet) for tet in mesh.IEN]

    # Export
    vtk = vtk_grid(filename, points, cells)
    vtk["SDF"] = mesh.node_sdf
    vtk["element_quality"] = element_qualities  # Přidání kvality elementů
    
    # Statistiky kvality
    min_quality = minimum(element_qualities)
    max_quality = maximum(element_qualities)
    avg_quality = mean(element_qualities)
    
    @info "Statistiky kvality sítě:"
    @info "  Minimální kvalita: $(round(min_quality, digits=3))"
    @info "  Průměrná kvalita: $(round(avg_quality, digits=3))"
    @info "  Maximální kvalita: $(round(max_quality, digits=3))"
    
    vtk_save(vtk)
    @info "VTK soubor uložen jako: $(abspath(filename))"
end

export_vtk_with_quality(mesh, "test_geom_quality")

