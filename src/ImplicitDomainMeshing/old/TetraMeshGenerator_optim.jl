using LinearAlgebra
using WriteVTK
using Statistics
using StaticArrays
using JLD2
using NLopt
using JuMP
using BenchmarkTools
using Base.Threads
import Ipopt

using Rho2sdf
using Rho2sdf.SignedDistances

# Optimized structure for cell data using StaticArrays
struct CellData
  coords::SVector{8,SVector{3,Float64}}
  sdf_values::SVector{8,Float64}
end

# Optimized BlockMesh structure with constructor
mutable struct BlockMesh
  nx::Int
  ny::Int
  nz::Int
  grid::Array{Vector{Float64},3}
  SDF::Array{Float64,3}
  X::Vector{Vector{Float64}}
  IEN::Vector{Vector{Int64}}
  INE::Vector{Vector{Int64}}
  node_sdf::Vector{Float64}
  node_map::Dict{Int64,Int64}

  function BlockMesh()
    # Load data from JLD2 files
    @load "src/ImplicitDomainMeshing/data/Z_block_FineGrid.jld2" fine_grid
    @load "src/ImplicitDomainMeshing/data/Z_block_FineSDF.jld2" fine_sdf

    # Convert Float32 to Float64 for consistency
    grid = Array{Vector{Float64},3}(undef, size(fine_grid))
    for i in eachindex(fine_grid)
      grid[i] = Vector{Float64}(fine_grid[i])
    end

    sdf = Float64.(fine_sdf)
    nx, ny, nz = size(grid)

    # Initialize other fields
    X = Vector{Vector{Float64}}()
    IEN = Vector{Vector{Int64}}()
    INE = Vector{Vector{Int64}}()
    node_sdf = Vector{Float64}()
    node_map = Dict{Int64,Int64}()

    new(nx, ny, nz, grid, sdf, X, IEN, INE, node_sdf, node_map)
  end
end

# Rest of the optimized code remains the same
function get_cell_data(mesh::BlockMesh, i::Int, j::Int, k::Int)
  coords = SVector{8,SVector{3,Float64}}(
    SVector{3,Float64}(mesh.grid[i, j, k]...),
    SVector{3,Float64}(mesh.grid[i+1, j, k]...),
    SVector{3,Float64}(mesh.grid[i+1, j+1, k]...),
    SVector{3,Float64}(mesh.grid[i, j+1, k]...),
    SVector{3,Float64}(mesh.grid[i, j, k+1]...),
    SVector{3,Float64}(mesh.grid[i+1, j, k+1]...),
    SVector{3,Float64}(mesh.grid[i+1, j+1, k+1]...),
    SVector{3,Float64}(mesh.grid[i, j+1, k+1]...)
  )

  sdf_values = SVector{8,Float64}(
    mesh.SDF[i, j, k],
    mesh.SDF[i+1, j, k],
    mesh.SDF[i+1, j+1, k],
    mesh.SDF[i, j+1, k],
    mesh.SDF[i, j, k+1],
    mesh.SDF[i+1, j, k+1],
    mesh.SDF[i+1, j+1, k+1],
    mesh.SDF[i, j+1, k+1]
  )

  return CellData(coords, sdf_values)
end

# Thread synchronization
const NodeLock = ReentrantLock()


function get_or_create_node_atomic!(mesh::BlockMesh, orig_id::Int, coords::SVector{3,Float64}, sdf::Float64)
  if !haskey(mesh.node_map, orig_id)
    lock(NodeLock) do
      if !haskey(mesh.node_map, orig_id)
        push!(mesh.X, Vector(coords))
        push!(mesh.node_sdf, sdf)
        mesh.node_map[orig_id] = length(mesh.X)
      end
    end
  end
  return mesh.node_map[orig_id]
end

# Pre-computed tetrahedra patterns
const TETRAHEDRA_PATTERNS = [
  SVector{4,Int}(1, 2, 3, 7),
  SVector{4,Int}(1, 2, 6, 7),
  SVector{4,Int}(1, 4, 3, 7),
  SVector{4,Int}(1, 4, 8, 7),
  SVector{4,Int}(1, 5, 6, 7),
  SVector{4,Int}(1, 5, 8, 7)
]

# Thread-safe container for collecting tetrahedra
struct TetrahedraCollector
  elements::Vector{Vector{Int64}}
  lock::ReentrantLock
end

TetrahedraCollector() = TetrahedraCollector(Vector{Vector{Int64}}(), ReentrantLock())

function add_tetrahedron!(collector::TetrahedraCollector, element::Vector{Int64})
  lock(collector.lock) do
    push!(collector.elements, element)
  end
end

function cleanup_unused_nodes!(mesh::BlockMesh)
  used_nodes = Set{Int64}()
  for element in mesh.IEN
    union!(used_nodes, element)
  end

  # Vytvoření nového mapování pro použité uzly
  new_node_map = Dict{Int64,Int64}()
  new_X = Vector{Vector{Float64}}()
  new_node_sdf = Vector{Float64}()

  # Přemapování uzlů a odstranění nepoužitých
  for (new_id, old_id) in enumerate(sort(collect(used_nodes)))
    new_node_map[old_id] = new_id
    push!(new_X, mesh.X[old_id])
    push!(new_node_sdf, mesh.node_sdf[old_id])
  end

  # Aktualizace konektivity elementů
  for i in 1:length(mesh.IEN)
    mesh.IEN[i] = [new_node_map[old_id] for old_id in mesh.IEN[i]]
  end

  mesh.X = new_X
  mesh.node_sdf = new_node_sdf
  mesh.node_map = new_node_map
end

# Efektivní deduplikace uzlů po vygenerování sítě
function deduplicate_nodes!(mesh::BlockMesh)
    # Vytvoření slovníku pro mapování pozic na indexy
    pos_to_idx = Dict{NTuple{3,Float64},Int}()
    new_indices = Int[]
    new_X = Vector{Vector{Float64}}()
    new_node_sdf = Vector{Float64}()
    
    # První průchod - identifikace unikátních uzlů
    for (i, pos) in enumerate(mesh.X)
        pos_tuple = (pos[1], pos[2], pos[3])
        if haskey(pos_to_idx, pos_tuple)
            # Pro duplicitní uzel použij existující index
            push!(new_indices, pos_to_idx[pos_tuple])
        else
            # Pro nový uzel vytvoř nový index
            new_idx = length(new_X) + 1
            pos_to_idx[pos_tuple] = new_idx
            push!(new_indices, new_idx)
            push!(new_X, pos)
            push!(new_node_sdf, mesh.node_sdf[i])
        end
    end
    
    # Aktualizace IEN s novými indexy
    for i in 1:length(mesh.IEN)
        mesh.IEN[i] = [new_indices[idx] for idx in mesh.IEN[i]]
    end
    
    # Aktualizace meshe
    mesh.X = new_X
    mesh.node_sdf = new_node_sdf
    
    # Přepočítání node_map
    old_node_map = mesh.node_map
    mesh.node_map = Dict{Int64,Int64}()
    for (orig_id, old_idx) in old_node_map
        mesh.node_map[orig_id] = new_indices[old_idx]
    end
    
    return length(mesh.X)
end

# Optimized mesh generation function
function generate_mesh!(mesh::BlockMesh)
  empty!(mesh.X)
  empty!(mesh.IEN)
  empty!(mesh.node_sdf)
  empty!(mesh.node_map)

  # Create thread-local collectors
  n_threads = Threads.nthreads()
  collectors = [TetrahedraCollector() for _ in 1:n_threads]

  # Parallel processing of cells
  Threads.@threads for linear_idx in 1:((mesh.nx-1)*(mesh.ny-1)*(mesh.nz-1))
    # Convert linear index to 3D indices
    i = (linear_idx - 1) % (mesh.nx - 1) + 1
    j = ((linear_idx - 1) ÷ (mesh.nx - 1)) % (mesh.ny - 1) + 1
    k = ((linear_idx - 1) ÷ ((mesh.nx - 1) * (mesh.ny - 1))) + 1

    # Get thread-local collector
    collector = collectors[Threads.threadid()]

    # Process cell
    cell_data = get_cell_data(mesh, i, j, k)

    # Skip if no positive values
    any(>=(0), cell_data.sdf_values) || continue

    # Create nodes and get their IDs
    orig_id_base = i + (j - 1) * mesh.nx + (k - 1) * mesh.nx * mesh.ny
    node_ids = Vector{Int64}(undef, 8)

    for n in 1:8
      orig_id = orig_id_base +
                (n - 1) ÷ 4 * (mesh.nx * mesh.ny) +
                ((n - 1) % 4) ÷ 2 * mesh.nx +
                (n - 1) % 2
      node_ids[n] = get_or_create_node_atomic!(
        mesh,
        orig_id,
        cell_data.coords[n],
        cell_data.sdf_values[n]
      )
    end

    # Generate tetrahedra
    for pattern in TETRAHEDRA_PATTERNS
      tet_sdf = SVector{4,Float64}(
        cell_data.sdf_values[pattern[1]],
        cell_data.sdf_values[pattern[2]],
        cell_data.sdf_values[pattern[3]],
        cell_data.sdf_values[pattern[4]]
      )

      if any(>=(0), tet_sdf)
        add_tetrahedron!(collector,
          [node_ids[pattern[1]],
            node_ids[pattern[2]],
            node_ids[pattern[3]],
            node_ids[pattern[4]]]
        )
      end
    end
  end

  # Combine results from all threads
  total_elements = sum(length(c.elements) for c in collectors)
  sizehint!(mesh.IEN, total_elements)
  for collector in collectors
    append!(mesh.IEN, collector.elements)
  end

  # Vyčištění nepoužitých uzlů
  cleanup_unused_nodes!(mesh)
  deduplicate_nodes!(mesh)

  # Create INE mapping
  create_INE!(mesh)

  @info "Generated $(length(mesh.X)) nodes and $(length(mesh.IEN)) tetrahedra using $(Threads.nthreads()) threads"
end

# Optimized INE creation
function create_INE!(mesh::BlockMesh)
  mesh.INE = [Vector{Int64}() for _ in 1:length(mesh.X)]

  # Pre-allocate approximate size for each node's element list
  avg_elements_per_node = 8  # Typical average for tetrahedral meshes
  for vec in mesh.INE
    sizehint!(vec, avg_elements_per_node)
  end

  # Parallel construction of INE using atomic operations
  element_locks = [ReentrantLock() for _ in 1:length(mesh.X)]

  Threads.@threads for elem_id in 1:length(mesh.IEN)
    for node_id in mesh.IEN[elem_id]
      lock(element_locks[node_id]) do
        push!(mesh.INE[node_id], elem_id)
      end
    end
  end

  return mesh
end

#WARNING:
@time mesh = BlockMesh()
@info "Generování sítě..."
@time generate_mesh!(mesh)

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
      else
        @warn println("Neznámá konfigurace, pro výpočet posunu uzlu.")
        println("Node ID: ", i)
      end
    end
  end

  # Update mesh coordinates
  copy_nonzero_vectors!(X_new, mesh.X)

  return mesh
end

println("Nodes to isocontour")
# @time project_nodes_to_isocontour!(mesh)

#INFO: Výpočet kvality:

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

println("export data:")
@time export_vtk_with_quality(mesh, "test_geom_quality")

