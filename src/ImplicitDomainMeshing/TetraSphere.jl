using LinearAlgebra
using WriteVTK
using StaticArrays

# Reprezentace čtyřstěnu s použitím StaticArrays pro lepší výkon
struct TetrahedronElement
  vertices::SMatrix{4,3,Float64}  # 4×3 statická matice souřadnic vrcholů
  sdf_values::SVector{4,Float64}  # 4 SDF hodnoty pro vrcholy

  # Konstruktor s kontrolou dimenzí a konverzí typů
  function TetrahedronElement(vertices::AbstractMatrix{Float64}, sdf_values::AbstractVector{Float64})
    size(vertices, 1) == 4 || throw(ArgumentError("Vertices must have 4 rows"))
    size(vertices, 2) == 3 || throw(ArgumentError("Vertices must have 3 columns"))
    length(sdf_values) == 4 || throw(ArgumentError("SDF values must have length 4"))

    # Konverze na statické typy
    return new(SMatrix{4,3}(vertices), SVector{4}(sdf_values))
  end
end

# Hlavní třída pro generování sítě koule
mutable struct SphereMesh
  nx::Int
  ny::Int
  nz::Int
  x::Vector{Float64}
  y::Vector{Float64}
  z::Vector{Float64}
  SDF::Array{Float64,3}
  tetrahedra::Vector{TetrahedronElement}

  # Konstruktor s validací vstupů
  function SphereMesh(nx::Int=30, ny::Int=30, nz::Int=30)
    nx > 1 || throw(ArgumentError("nx must be > 1"))
    ny > 1 || throw(ArgumentError("ny must be > 1"))
    nz > 1 || throw(ArgumentError("nz must be > 1"))

    mesh = new(nx, ny, nz)
    mesh.tetrahedra = TetrahedronElement[]
    setup_grid!(mesh)
    return mesh
  end
end

# Inicializace výpočetní mřížky
function setup_grid!(mesh::SphereMesh)
  # Použití stejných rozsahů jako v Pythonu pro konzistenci
  mesh.x = collect(range(-1, 1, length=mesh.nx))
  mesh.y = collect(range(-1, 1, length=mesh.ny))
  mesh.z = collect(range(-1, 1, length=mesh.nz))

  # Pre-alokace pole pro SDF
  mesh.SDF = Array{Float64}(undef, mesh.nx, mesh.ny, mesh.nz)
  radius = 0.5

  # Vektorizovaný výpočet SDF
  for i in 1:mesh.nx, j in 1:mesh.ny, k in 1:mesh.nz
    mesh.SDF[i, j, k] = radius - sqrt(mesh.x[i]^2 + mesh.y[j]^2 + mesh.z[k]^2)
  end
end

# Bezpečné získání SDF hodnot pro vrcholy buňky
function get_cell_sdf_values(mesh::SphereMesh, i::Int, j::Int, k::Int)
  # Kontrola hranic
  1 <= i < mesh.nx || throw(BoundsError(mesh.SDF, i))
  1 <= j < mesh.ny || throw(BoundsError(mesh.SDF, j))
  1 <= k < mesh.nz || throw(BoundsError(mesh.SDF, k))

  return SVector{8}(
    mesh.SDF[i, j, k],        # 0: front-bottom-left
    mesh.SDF[i+1, j, k],      # 1: front-bottom-right
    mesh.SDF[i+1, j+1, k],    # 2: front-top-right
    mesh.SDF[i, j+1, k],      # 3: front-top-left
    mesh.SDF[i, j, k+1],      # 4: back-bottom-left
    mesh.SDF[i+1, j, k+1],    # 5: back-bottom-right
    mesh.SDF[i+1, j+1, k+1],  # 6: back-top-right
    mesh.SDF[i, j+1, k+1]     # 7: back-top-left
  )
end

# Klasifikace buňky podle SDF hodnot
function classify_cell(sdf_values::AbstractVector{Float64})
  if all(>=(0), sdf_values)  # Optimalizovaná verze all(sdf_values .>= 0)
    return :full
  elseif all(<(0), sdf_values)  # Optimalizovaná verze all(sdf_values .< 0)
    return :empty
  else
    return :transition
  end
end

# Lineární interpolace pro nalezení nulového bodu
function interpolate_zero_crossing(
  p1::AbstractVector{Float64},
  p2::AbstractVector{Float64},
  sdf1::Float64,
  sdf2::Float64
)
  # Ochrana proti dělení nulou
  denom = sdf1 - sdf2
  if abs(denom) < eps(Float64)
    return p1
  end

  t = clamp(sdf1 / denom, 0.0, 1.0)  # Zajištění stability interpolace
  return @. p1 + t * (p2 - p1)  # Vektorizovaná operace
end

# Úprava vrcholů čtyřstěnu na povrchu koule
function adjust_tetrahedron_vertices(tet::TetrahedronElement)
  vertices_mat = Matrix(tet.vertices)  # Konverze na běžnou matici pro modifikace

  for i in 1:4
    if tet.sdf_values[i] < 0
      vertex = view(vertices_mat, i, :)
      distances = [norm(view(vertices_mat, j, :) - vertex) for j in 1:4]
      positive_indices = findall(>=(0), tet.sdf_values)

      if !isempty(positive_indices)
        _, closest_idx = findmin(view(distances, positive_indices))
        closest_positive_idx = positive_indices[closest_idx]

        p_positive = view(vertices_mat, closest_positive_idx, :)
        sdf_positive = tet.sdf_values[closest_positive_idx]

        vertices_mat[i, :] .= interpolate_zero_crossing(
          Vector(vertex), Vector(p_positive),
          tet.sdf_values[i], sdf_positive
        )
      end
    end
  end

  return TetrahedronElement(vertices_mat, tet.sdf_values)
end

# Vytvoření pravidelné tetrahedrální sítě pro buňku
function tetrahedralize_regular_cell(mesh::SphereMesh, i::Int, j::Int, k::Int)
  # Pre-alokace matice vrcholů
  vertices = Matrix{Float64}(undef, 8, 3)

  # Naplnění matice vrcholů
  vertices[1, :] .= [mesh.x[i], mesh.y[j], mesh.z[k]]    # front-bottom-left
  vertices[2, :] .= [mesh.x[i+1], mesh.y[j], mesh.z[k]]    # front-bottom-right
  vertices[3, :] .= [mesh.x[i+1], mesh.y[j+1], mesh.z[k]]    # front-top-right
  vertices[4, :] .= [mesh.x[i], mesh.y[j+1], mesh.z[k]]    # front-top-left
  vertices[5, :] .= [mesh.x[i], mesh.y[j], mesh.z[k+1]]  # back-bottom-left
  vertices[6, :] .= [mesh.x[i+1], mesh.y[j], mesh.z[k+1]]  # back-bottom-right
  vertices[7, :] .= [mesh.x[i+1], mesh.y[j+1], mesh.z[k+1]]  # back-top-right
  vertices[8, :] .= [mesh.x[i], mesh.y[j+1], mesh.z[k+1]]  # back-top-left

  # Definice čtyřstěnů (1-based indexing)
  tets = SMatrix{6,4,Int}([
    1 2 4 6;  # První čtyřstěn
    2 3 4 6;  # Druhý čtyřstěn
    1 4 6 8;  # Třetí čtyřstěn
    1 6 5 8;  # Čtvrtý čtyřstěn
    3 4 7 6;  # Pátý čtyřstěn
    4 6 7 8   # Šestý čtyřstěn
  ])

  return vertices, tets
end

# Zpracování jedné buňky
function process_cell(mesh::SphereMesh, i::Int, j::Int, k::Int)
  vertices, tets = tetrahedralize_regular_cell(mesh, i, j, k)
  sdf_values = get_cell_sdf_values(mesh, i, j, k)
  cell_type = classify_cell(sdf_values)

  result_tets = TetrahedronElement[]

  if cell_type != :empty
    for tet_idx in 1:size(tets, 1)
      tet_indices = view(tets, tet_idx, :)
      tet_vertices = vertices[tet_indices, :]
      # Konverze na Vector pro konstruktor
      tet_sdf = Vector(sdf_values[tet_indices])

      tet = TetrahedronElement(tet_vertices, tet_sdf)

      if any(<(0), tet_sdf) && any(>=(0), tet_sdf)
        tet = adjust_tetrahedron_vertices(tet)
      end

      if !all(<(0), tet_sdf)
        push!(result_tets, tet)
      end
    end
  end

  return result_tets
end


# Generování celé sítě
function generate_mesh!(mesh::SphereMesh)
  empty!(mesh.tetrahedra)
  sizehint!(mesh.tetrahedra, mesh.nx * mesh.ny * mesh.nz * 6)  # Pre-alokace

  # Paralelní zpracování buněk by bylo možné přidat zde
  for i in 1:mesh.nx-1
    for j in 1:mesh.ny-1
      for k in 1:mesh.nz-1
        append!(mesh.tetrahedra, process_cell(mesh, i, j, k))
      end
    end
  end
end

# Export do VTK formátu
function export_vtk(mesh::SphereMesh, filename::String)
    # Zajistíme, že filename končí .vtu
    if !endswith(filename, ".vtu")
        filename = filename * ".vtu"
    end
    
    # Pre-alokace polí pro body a jejich hodnoty
    points_dict = Dict{NTuple{3,Float64}, Int}()
    points = Matrix{Float64}(undef, 3, 0)  # Inicializace prázdné matice s správnými dimenzemi
    sdf_values = Float64[]
    
    # Seznam čtyřstěnů pro VTK
    cells = Vector{MeshCell{VTKCellTypes.VTKCellType, Vector{Int}}}()
    sizehint!(cells, length(mesh.tetrahedra))
    
    # Zpracování všech tetraedrů
    for tet in mesh.tetrahedra
        tet_points = zeros(Int, 4)
        
        for i in 1:4
            vertex = Tuple(tet.vertices[i,:])  # Konverze na tuple pro použití jako klíč
            
            if !haskey(points_dict, vertex)
                # Přidání nového bodu do matice bodů
                points = hcat(points, [vertex[1], vertex[2], vertex[3]])
                push!(sdf_values, tet.sdf_values[i])
                points_dict[vertex] = size(points, 2)  # Nový index je počet sloupců
            end
            
            tet_points[i] = points_dict[vertex]
        end
        
        push!(cells, MeshCell(VTKCellTypes.VTK_TETRA, tet_points))
    end
    
    # Kontrola, zda máme nějaké body k exportu
    if size(points, 2) == 0
        @warn "No points to export!"
        return
    end
    
    # Vytvoření VTK souboru
    vtk = vtk_grid(filename, points, cells)
    vtk["SDF"] = sdf_values
    vtk["radius"] = fill(0.5, length(sdf_values))
    vtk_save(vtk)
    
    @info "VTK soubor uložen jako: $(abspath(filename))"
    @info "Exported $(size(points, 2)) points and $(length(cells)) tetrahedra"
end

# Hlavní funkce pro test
function main()
  mesh = SphereMesh(30, 30, 30)
  @info "Generování sítě..."
  generate_mesh!(mesh)
  @info "Export do VTK..."
  export_vtk(mesh, "sphere_mesh")
  @info "Hotovo! Počet tetraedrů: $(length(mesh.tetrahedra))"
end

# Spuštění testu
main()
