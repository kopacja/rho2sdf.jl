using LinearAlgebra
using WriteVTK
using StaticArrays: SVector

# Struktura reprezentující síť
struct TetMesh
  vertices::Vector{Vector{Float64}}         # Souřadnice uzlů tetraedrů
  tets::Vector{Vector{Int}}                # Konektivita tetrahedrů
  hex_vertices::Vector{Vector{Float64}}     # 8 vrcholů původního hexahedru
end

# Definice acute tetrahedral lattice pro jeden element (27)
const LATTICE_NODES = [
  [2, 1, 5],
  [4, 1, 5],
  [6, 1, 5],
  [1, 1, 3],
  [3, 2, 3],
  [5, 1, 3],
  [2, 1, 1],
  [4, 1, 1],
  [6, 1, 1],
  [1, 3, 6],
  [3, 3, 5],
  [5, 3, 6],
  [1, 3, 4],
  [3, 4, 3],
  [5, 3, 4],
  [1, 3, 2],
  [3, 3, 1],
  [5, 3, 2],
  [2, 5, 5],
  [4, 5, 5],
  [6, 5, 5],
  [1, 5, 3],
  [3, 6, 3],
  [5, 5, 3],
  [2, 5, 1],
  [4, 5, 1],
  [6, 5, 1]
  # [3, 3, 2],
  # [4, 4, 2],
  # [5, 5, 2],
  # [2, 2, 4],
  # [3, 1, 4],
  # [5, 1, 4]
]

const LATTICE_TETS = [
  [13, 5, 1, 4],
  [1, 13, 10, 11],
  [5, 13, 1, 11],
  [2, 5, 1, 11],
  [14, 13, 5, 11],

  # [15, 2, 3, 12],
  # [2, 15, 5, 11],
  # [15, 2, 5, 6],
  # [15, 14, 5, 11],
  # [2, 15, 3, 6],
  # [2, 12, 15, 11],
  #
  # [7, 8, 17, 5],
  # [7, 5, 16, 4],
  # [16, 5, 13, 4],
  # [13, 14, 5, 16],
  # [17, 16, 7, 5],
  # [14, 16, 17, 5],
  #
  # [5, 18, 8, 17],
  # [18, 5, 8, 6],
  # [18, 14, 17, 5],
  # [5, 15, 6, 18],
  # [9, 6, 18, 8],
  # [14, 18, 15, 5],
  #
  # [23, 20, 19, 14],
  # [19, 13, 11, 10],
  # [14, 20, 19, 11],
  # [14, 13, 19, 22],
  # [13, 14, 19, 11],
  # [19, 23, 14, 22], # ok
  #
  # [20, 14, 15, 11],
  # [20, 12, 11, 15],
  # [12, 20, 21, 15],
  # [15, 20, 21, 24],
  # [20, 23, 24, 14],
  # [14, 20, 15, 24],
  #
  # [14, 13, 22, 16],
  # [25, 23, 22, 14],
  # [16, 25, 14, 17],
  # [25, 16, 14, 22],
  # [25, 26, 23, 14],
  # [26, 25, 17, 14],
  #
  # [26, 23, 14, 24],
  # [14, 18, 17, 26],
  # [26, 14, 18, 24],
  # [14, 15, 18, 24],
  # [24, 27, 18, 26]
]

function create_hex_vertices(origin::Vector{Float64}, size::Float64)
  # Vytvoření 8 vrcholů hexahedru
  vertices = Vector{Vector{Float64}}()

  # Pořadí vrcholů podle VTK konvence:
  # 0: (0,0,0), 1: (1,0,0), 2: (1,1,0), 3: (0,1,0)
  # 4: (0,0,1), 5: (1,0,1), 6: (1,1,1), 7: (0,1,1)
  push!(vertices, origin)                                    # 1
  push!(vertices, origin + [size, 0, 0])                    # 2
  push!(vertices, origin + [size, size, 0])                 # 3
  push!(vertices, origin + [0, size, 0])                    # 4
  push!(vertices, origin + [0, 0, size])                    # 5
  push!(vertices, origin + [size, 0, size])                 # 6
  push!(vertices, origin + [size, size, size])              # 7
  push!(vertices, origin + [0, size, size])                 # 8

  return vertices
end

function create_single_hex_element(origin::Vector{Float64}, size::Float64)
  vertices = Vector{Vector{Float64}}()
  local_vertices = Dict{Int,Vector{Float64}}()

  # # Vytvoření všech uzlů podle lattice
  # for (idx, node) in enumerate(LATTICE_NODES)
  #   x = origin[1] + (node[1] - 1) * size / 5
  #   y = origin[2] + (node[2] - 1) * size / 5
  #   z = origin[3] + (node[3] - 1) * size / 5
  #   local_vertices[idx] = [x, y, z]
  # end

  # Vytvoření všech uzlů podle lattice
  for (idx, node) in enumerate(LATTICE_NODES)
    x = origin[1] + (node[1] - 1) * size / 5
    y = origin[2] + (node[2] - 1) * size / 5
    z = origin[3] + (node[3] - 1) * size / 5
    local_vertices[idx] = [x, y, z]
  end


  # Převedení lokálních vrcholů do globálního pole
  # vertices = [v for (_, v) in sort(local_vertices)]
  vertices = [local_vertices[i] for i in eachindex(LATTICE_NODES)]

  # Vytvoření tetraedrů - všechny jsou validní pro jednoduchý element
  tets = copy(LATTICE_TETS)  # Použijeme všechny tetraedry

  hex_vertices = create_hex_vertices(origin, size)

  return TetMesh(vertices, tets, hex_vertices)
end

# Pomocná funkce pro konverzi bodů do VTK formátu
function points_to_vtk_format(points::Vector{Vector{Float64}})
  return [SVector{3,Float64}(v...) for v in points]
end


function export_to_vtk(mesh::TetMesh, filename::String)
    # Převod bodů do formátu pro VTK
    points_tet = points_to_vtk_format(mesh.vertices)
    points_hex = points_to_vtk_format(mesh.hex_vertices)

    # Tetrahedra
    cells_tet = [MeshCell(VTKCellTypes.VTK_TETRA, tet) for tet in mesh.tets]
    vtk = vtk_grid(filename * "_tets", points_tet, cells_tet)
    vtk_save(vtk)

    # Hexahedron
    cells_hex = [MeshCell(VTKCellTypes.VTK_HEXAHEDRON, collect(1:8))]
    vtk_hex = vtk_grid(filename * "_hex", points_hex, cells_hex)
    vtk_save(vtk_hex)

    # Vrcholy hexahedru
    vtk_vertices = vtk_grid(filename * "_hex_vertices", points_hex)
    vtk_save(vtk_vertices)
    
    # Vrcholy tetraedrální sítě
    vtk_tet_vertices = vtk_grid(filename * "_tet_vertices", points_tet)
    vtk_save(vtk_tet_vertices)
end

# Příklad použití
function main()
  # Vytvoření jednoho elementu se středem v počátku a velikostí 1.0
  origin = [0.0, 0.0, 0.0]
  size = 1.0
  mesh = create_single_hex_element(origin, size)

  # Export do VTK pro vizualizaci
  export_to_vtk(mesh, "hex_discretization")

  println("Počet uzlů tetraedrů: ", length(mesh.vertices))
  println("Počet tetraedrů: ", length(mesh.tets))
  println("Počet vrcholů hexahedru: ", length(mesh.hex_vertices))
end

main()
