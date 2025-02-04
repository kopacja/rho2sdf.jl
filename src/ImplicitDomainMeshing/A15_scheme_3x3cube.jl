using StaticArrays
using WriteVTK

# Definice datové struktury pro blokovou síť
struct BlockMesh
    X::Vector{SVector{3,Float64}}      # Seznam souřadnic uzlů
    IEN::Vector{Vector{Int}}            # Konektivita tetraedrických elementů (každý čtyřprvkový vektor obsahuje indexy uzlů)
    node_sdf::Vector{Float64}           # Skalární pole – zde pouze nulové hodnoty (lze rozšířit)
end

# Export geometrie do VTK (inspirace z předloženého kódu)
function export_mesh_vtk(mesh::BlockMesh, filename::String)
    npoints = length(mesh.X)
    points = zeros(Float64, 3, npoints)
    for i in 1:npoints
        points[:, i] = mesh.X[i]
    end

    # VTK typ tetraedru je 10; pro každý tetraedr vytvoříme MeshCell
    cells = [MeshCell(VTKCellTypes.VTK_TETRA, tet) for tet in mesh.IEN]
    vtkfile = vtk_grid(filename, points, cells)
    vtk_point_data(vtkfile, mesh.node_sdf, "sdf")
    vtk_save(vtkfile)
end

# -------------------------------
# Definice A15 schématu
# -------------------------------

# tile_ref obsahuje lokální souřadnice uzlů (pozn.: hodnoty mimo interiér elementu zůstávají zachovány)
const tile_ref = [
    SVector(1.0, 0.0, 0.0),   # 0
    SVector(2.0, 2.0, 0.0),   # 1
    SVector(1.0, 4.0, 0.0),   # 2
    SVector(3.0, 4.0, 0.0),   # 3
    SVector(1.0, 0.0, 4.0),   # 4
    SVector(3.0, 0.0, 4.0),   # 5
    SVector(2.0, 1.0, 2.0),   # 6
    SVector(0.0, 2.0, 1.0),   # 7
    SVector(0.0, 2.0, 3.0),   # 8
    SVector(2.0, 2.0, 4.0),   # 9
    SVector(1.0, 4.0, 4.0),   # 10
    SVector(3.0, 4.0, 4.0),   # 11
    SVector(0.0, 4.0, 2.0),   # 12
    SVector(2.0, 3.0, 2.0),   # 13
    SVector(2.0, 5.0, 2.0),   # 14
    SVector(4.0, 2.0, 1.0),   # 15
    SVector(4.0, 2.0, 3.0),   # 16
    SVector(5.0, 4.0, 4.0),   # 17
    SVector(4.0, 4.0, 2.0),   # 18
    SVector(0.0, 2.0, 5.0),   # 19
    SVector(4.0, 2.0, 5.0),   # 20
    SVector(0.0, 0.0, 2.0),   # 21
    SVector(5.0, 0.0, 4.0),   # 22
    SVector(4.0, 0.0, 2.0),   # 23
    SVector(3.0, 0.0, 0.0),   # 24
    SVector(5.0, 0.0, 0.0),   # 25
    SVector(5.0, 4.0, 0.0)    # 26
]

# Definice konektivity tetraedrů (v literatuře jsou indexy 0-based, převedeme je na 1-based při zpracování)
const tetra_connectivity = [
    [3,4,15,14], # 0
    [3,15,13,14],
    [6,21,17,10],
    [6,17,23,24],
    [12,14,17,10],
    [1,25,2,7],
    [21,12,18,17],
    [4,14,16,19],
    [4,15,14,19],
    [14,16,2,4],
    [1,7,8,22],
    [7,16,25,2],
    [9,7,5,22],
    [8,7,9,22],
    [14,12,17,19],
    [12,21,10,17],
    [14,9,13,8],
    [8,3,14,2],
    [14,17,16,19],
    [17,14,7,10],
    [16,7,25,24],
    [17,6,7,24],
    [11,15,14,13],
    [4,3,2,14],
    [9,14,7,8],
    [9,14,11,10],
    [2,8,1,7],
    [14,8,2,7],
    [12,15,19,14],
    [19,27,16,4], # zde se vyskytuje index 27 (0-based), který po posunu dává 28 – předpokládáme, že jde o drobnou nepřesnost v literatuře
    [17,12,18,19],
    [14,9,11,13],
    [14,12,11,10],
    [3,8,14,13],
    [6,17,7,10],
    [5,9,20,10],
    [14,9,7,10],
    [11,9,10,20],
    [7,9,5,10],
    [17,6,23,21],
    [6,7,5,10],
    [15,12,11,14],
    [16,14,2,7],
    [7,17,24,16],
    [26,24,16,25],
    [14,16,17,7]
]

"""
    remove_duplicate_nodes(nodes, connectivity; tol=1e-8)

Funkce projde všechny uzly (nodes) a najde ty, které se shodují (s ohledem na toleranci `tol`).
Vrací nový seznam unikátních uzlů a aktualizovanou konektivitu (connectivity), kde jsou duplicitní
uzly sloučeny.
"""
function remove_duplicate_nodes(nodes::Vector{SVector{3,Float64}}, connectivity::Vector{Vector{Int}}; tol=1e-8)
    unique_dict = Dict{Tuple{Float64,Float64,Float64}, Int}()
    new_nodes = SVector{3,Float64}[]
    old_to_new = Dict{Int,Int}()
    
    for (i, node) in enumerate(nodes)
        # Zaokrouhlení souřadnic s ohledem na toleranci
        rnode = (round(node[1]/tol)*tol, round(node[2]/tol)*tol, round(node[3]/tol)*tol)
        if haskey(unique_dict, rnode)
            old_to_new[i] = unique_dict[rnode]
        else
            new_index = length(new_nodes) + 1
            push!(new_nodes, node)
            unique_dict[rnode] = new_index
            old_to_new[i] = new_index
        end
    end
    
    # Aktualizace konektivity: každý starý index nahradíme novým
    new_connectivity = Vector{Vector{Int}}()
    for tet in connectivity
         new_tet = [old_to_new[idx] for idx in tet]
         push!(new_connectivity, new_tet)
    end
    return new_nodes, new_connectivity
end

# -------------------------------
# Sestavení globální sítě pro více elementů
# -------------------------------
"""
    build_mesh(n_x, n_y, n_z)

Sestaví globální síť z osmiuzlových elementů podle A15 schématu.
Elementy jsou umístěny v mřížce (n_x × n_y × n_z), přičemž každý element
je translací základního tile_ref. Po vytvoření lokálních konektivit se provede
konverze z 0-based indexace na 1-based (pomocí `li + 1`) a následně odstranění duplicitních uzlů.
"""
function build_mesh(n_x::Int, n_y::Int, n_z::Int)
    all_nodes = SVector{3,Float64}[]
    global_IEN = Vector{Vector{Int}}()
    
    # Předpokládaná velikost jednoho elementu je [0,5] ve všech směrech.
    element_size = 4.0
    
    # Projdeme všechny elementy v 3D mřížce
    for k in 0:(n_z-1)
       for j in 0:(n_y-1)
          for i in 0:(n_x-1)
             # Translací posuneme základní element
             T = SVector(i * element_size, j * element_size, k * element_size)
             # Vytvoření lokálního seznamu uzlů pro aktuální element
             local_mapping = Dict{Int, Int}()
             for li in 1:length(tile_ref)
                node = tile_ref[li] + T
                new_index = length(all_nodes) + 1
                push!(all_nodes, node)
                local_mapping[li] = new_index
             end
             # Vytvoření konektivity tetraedrických elementů pro aktuální element.
             # KONVERZE: Původní konektivity jsou v 0-based indexaci, proto pro každý lokální index použijeme li + 1.
             for tet in tetra_connectivity
                global_tet = [local_mapping[li] for li in tet]
                push!(global_IEN, global_tet)
             end
          end
       end
    end

    # Odstranění duplicitních uzlů (když se elementy dotýkají, budou mít společné uzly)
    new_nodes, new_connectivity = remove_duplicate_nodes(all_nodes, global_IEN)
    node_sdf = zeros(Float64, length(new_nodes))
    return BlockMesh(new_nodes, new_connectivity, node_sdf)
end

# -------------------------------
# Hlavní funkce – test diskretizace
# -------------------------------
# function main()
    # Test: 3×3×3 osmiuzlových elementů, které se dotýkají a tvoří jednu velkou krychli
    mesh = build_mesh(3, 3, 3)
    export_mesh_vtk(mesh, "mesh_3x3cube")
    println("Mesh bylo exportováno do souboru mesh.vtk")
# end

# main()

mesh.x

3^3*27
# Jakým způsobem dochází k mazání duplicitních uzlů
