using WriteVTK
using StaticArrays
using Logging
using JLD2
using LinearAlgebra

# Rozšířená definice struktury BlockMesh – nyní obsahuje i mapování pro středy buňky
mutable struct BlockMesh
  nx::Int
  ny::Int
  nz::Int
  grid::Array{Vector{Float64},3}       # 3D pole souřadnic uzlů
  SDF::Array{Float64,3}                 # SDF hodnoty
  X::Vector{Vector{Float64}}            # Seznam souřadnic použitých uzlů
  IEN::Vector{Vector{Int64}}            # Konektivita tetraedrů (elementů)
  INE::Vector{Vector{Int64}}            # Inverzní konektivita: pro každý uzel seznam přilehlých elementů
  node_sdf::Vector{Float64}             # SDF hodnoty pro uzly
  node_map::Dict{Int64,Int64}           # mapování původního indexu gridového uzlu -> nový index
  cell_center_map::Dict{Tuple{Int,Int,Int},Int64}  # mapování buňkových středu (Steiner body)

  function BlockMesh()
    # Načtení dat – cestu k souborům upravte dle potřeby
    @load "src/ImplicitDomainMeshing/data/Z_block_FineGrid.jld2" fine_grid
    @load "src/ImplicitDomainMeshing/data/Z_block_FineSDF.jld2" fine_sdf

    # Konverze dat pro konzistenci (Float32 -> Float64)
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
    return mesh
  end
end

# Pomocná funkce pro získání původního globálního indexu uzlu z indexů i,j,k
function get_original_node_id(mesh::BlockMesh, i::Int, j::Int, k::Int)
  return i + (j - 1)*mesh.nx + (k - 1)*mesh.nx*mesh.ny
end

# Vrací (nebo vytvoří) nový index pro gridový uzel (vrchol buňky)
function get_or_create_node!(mesh::BlockMesh, i::Int, j::Int, k::Int)
  orig_id = get_original_node_id(mesh, i, j, k)
  if !haskey(mesh.node_map, orig_id)
    push!(mesh.X, mesh.grid[i, j, k])
    push!(mesh.node_sdf, mesh.SDF[i, j, k])
    mesh.node_map[orig_id] = length(mesh.X)
  end
  return mesh.node_map[orig_id]
end

# Vrací (nebo vytvoří) uzel ve středu buňky (Steiner bod) pro buňku s indexy (i,j,k)
function get_or_create_cell_center!(mesh::BlockMesh, i::Int, j::Int, k::Int)
  key = (i,j,k)
  if haskey(mesh.cell_center_map, key)
    return mesh.cell_center_map[key]
  end

  # Získání souřadnic všech 8 vrcholů buňky
  pts = [
    mesh.grid[i, j, k],
    mesh.grid[i+1, j, k],
    mesh.grid[i, j+1, k],
    mesh.grid[i+1, j+1, k],
    mesh.grid[i, j, k+1],
    mesh.grid[i+1, j, k+1],
    mesh.grid[i, j+1, k+1],
    mesh.grid[i+1, j+1, k+1]
  ]
  center = zeros(length(pts[1]))
  for pt in pts
    center .+= pt
  end
  center ./= 8.0

  # Výpočet průměrné hodnoty SDF v buňce
  sdf_vals = [
    mesh.SDF[i, j, k],
    mesh.SDF[i+1, j, k],
    mesh.SDF[i, j+1, k],
    mesh.SDF[i+1, j+1, k],
    mesh.SDF[i, j, k+1],
    mesh.SDF[i+1, j, k+1],
    mesh.SDF[i, j+1, k+1],
    mesh.SDF[i+1, j+1, k+1]
  ]
  avg_sdf = sum(sdf_vals) / 8.0

  push!(mesh.X, center)
  push!(mesh.node_sdf, avg_sdf)
  new_index = length(mesh.X)
  mesh.cell_center_map[key] = new_index
  return new_index
end

# Získá SDF hodnoty osmi vrcholů buňky
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

# lattice discretizace
function process_cell!(mesh::BlockMesh, i::Int, j::Int, k::Int)
    # Získáme SDF hodnoty všech 8 vrcholů buňky
    sdf_values = get_cell_sdf_values(mesh, i, j, k)
    # Pokud buňka nemá žádný uzel s f ≥ 0, přeskočíme ji
    if !any(x -> x >= 0, sdf_values)
      return
    end
  
    # Získáme indexy vrcholů buňky
    v1 = get_or_create_node!(mesh, i, j, k)         # front-bottom-left
    v2 = get_or_create_node!(mesh, i+1, j, k)         # front-bottom-right
    v3 = get_or_create_node!(mesh, i+1, j+1, k)       # front-top-right
    v4 = get_or_create_node!(mesh, i, j+1, k)         # front-top-left
    v5 = get_or_create_node!(mesh, i, j, k+1)         # back-bottom-left
    v6 = get_or_create_node!(mesh, i+1, j, k+1)       # back-bottom-right
    v7 = get_or_create_node!(mesh, i+1, j+1, k+1)     # back-top-right
    v8 = get_or_create_node!(mesh, i, j+1, k+1)       # back-top-left
  
    # Získáme středu buňky (Steiner bod)
    vc = get_or_create_cell_center!(mesh, i, j, k)
  
    # Pomocná funkce pro výpočet orientovaného objemu tetraedru
    tet_volume = function(a::Vector{Float64}, b::Vector{Float64}, c::Vector{Float64}, d::Vector{Float64})
        return dot(d .- a, cross(b .- a, c .- a)) / 6.0
    end
  
    # Definujeme obličejové čtyřúhelníky buňky
    faces = [
      [v1, v2, v3, v4],  # front
      [v5, v6, v7, v8],  # back
      [v1, v4, v8, v5],  # left
      [v2, v3, v7, v6],  # right
      [v4, v3, v7, v8],  # top
      [v1, v2, v6, v5]   # bottom
    ]
  
    for face in faces
      # Rozdělení čtyřúhelníku na dvě trojúhelníkové oblasti
      triangles = [
        [face[1], face[2], face[3]],
        [face[1], face[3], face[4]]
      ]
      for tri in triangles
        # Vytvoříme tetraedr spojením středu buňky a vrcholů trojúhelníku
        tet = [vc, tri[1], tri[2], tri[3]]
        # Získáme souřadnice jednotlivých vrcholů
        a = mesh.X[tet[1]]
        b = mesh.X[tet[2]]
        c = mesh.X[tet[3]]
        d = mesh.X[tet[4]]
        # Pokud je orientovaný objem záporný nebo nulový, zaměníme poslední dva vrcholy
        if tet_volume(a, b, c, d) <= 0
          tet = [vc, tri[1], tri[3], tri[2]]
        end
        # Přidáme tetraedr do konektivity, pokud alespoň jeden z jeho uzlů má f ≥ 0
        tet_sdf = [mesh.node_sdf[idx] for idx in tet]
        if any(x -> x >= 0, tet_sdf)
          push!(mesh.IEN, tet)
        end
      end
    end
  end
  

# Vyčištění nepoužitých uzlů a přeindexování konektivity
function cleanup_unused_nodes!(mesh::BlockMesh)
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
    mesh.IEN[i] = [new_node_map[old_id] for old_id in mesh.IEN[i]]
  end
  mesh.X = new_X
  mesh.node_sdf = new_node_sdf
  mesh.node_map = new_node_map
  @info "Počet uzlů po vyčištění: $(length(mesh.X))"
end

# Vytvoří INE – inverzní konektivitu
function create_INE!(mesh::BlockMesh)
  mesh.INE = [Vector{Int64}() for _ in 1:length(mesh.X)]
  for (elem_id, element) in enumerate(mesh.IEN)
    for node_id in element
      push!(mesh.INE[node_id], elem_id)
    end
  end
  return mesh
end

# Generování sítě
function generate_mesh!(mesh::BlockMesh)
  empty!(mesh.X)
  empty!(mesh.IEN)
  empty!(mesh.node_sdf)
  empty!(mesh.node_map)
  empty!(mesh.cell_center_map)
  for i in 1:mesh.nx-1
    for j in 1:mesh.ny-1
      for k in 1:mesh.nz-1
        process_cell!(mesh, i, j, k)
      end
    end
  end
  cleanup_unused_nodes!(mesh)
  create_INE!(mesh)
  @info "Vytvořeno $(length(mesh.X)) uzlů a $(length(mesh.IEN)) tetraedrů"
end

function export_mesh_vtk(mesh::BlockMesh, filename::String)
    # First, prepare the point coordinates
    # We need to reshape our collection of 3D points into a (3, npoints) array
    # where each column represents a point's coordinates
    npoints = length(mesh.X)
    points = zeros(Float64, 3, npoints)
    for i in 1:npoints
        points[:, i] = mesh.X[i]  # Assign each point's coordinates as a column
    end
    
    # Convert our tetrahedra connectivity list to the format expected by WriteVTK
    # VTK_TETRA = 10 is the VTK cell type for tetrahedra
    cells = [MeshCell(VTKCellTypes.VTK_TETRA, tet) for tet in mesh.IEN]
    
    # Create the VTK grid with our points and cells
    vtkfile = vtk_grid(filename, points, cells)
    
    # Add point data - we can include the SDF values as a scalar field
    vtk_point_data(vtkfile, mesh.node_sdf, "sdf")
    
    # Add cell data if needed
    # vtk_cell_data(vtkfile, cell_data, "cell_field_name")
    
    # Write the file
    vtk_save(vtkfile)
end


# Hlavní běh – vytvoříme síť a exportujeme ji do souboru
mesh = BlockMesh()
generate_mesh!(mesh)
export_mesh_vtk(mesh, "final_mesh.vtu")
