

using StaticArrays
using LinearAlgebra
using WriteVTK
using Logging

# Pomocná funkce – trilineární interpolace bodu (r,s,t) do fyzikálních souřadnic buňky
# Vstupní buňkou jsou její 8 uzlů ve standardním pořadí:
# A = (0,0,0), B = (1,0,0), C = (1,1,0), D = (0,1,0),
# E = (0,0,1), F = (1,0,1), G = (1,1,1), H = (0,1,1)
function trilinear_map(r::Float64, s::Float64, t::Float64,
                         A::Vector{Float64}, B::Vector{Float64},
                         C::Vector{Float64}, D::Vector{Float64},
                         E::Vector{Float64}, F::Vector{Float64},
                         G::Vector{Float64}, H::Vector{Float64})
    one_r = 1 - r
    one_s = 1 - s
    one_t = 1 - t
    # Váhové koeficienty pro jednotlivé rohy
    wA = one_r * one_s * one_t
    wB = r * one_s * one_t
    wC = r * s * one_t
    wD = one_r * s * one_t
    wE = one_r * one_s * t
    wF = r * one_s * t
    wG = r * s * t
    wH = one_r * s * t
    x = wA .* A .+ wB .* B .+ wC .* C .+ wD .* D .+ wE .* E .+ wF .* F .+ wG .* G .+ wH .* H
    return x
end

# Testovací funkce: diskretizace jednoho osmiuzlového elementu (krychle) pomocí modifikovaného A15 vzoru
function test_single_cell_A15_modified()
    println("Test diskretizace modifikovaným A15 vzorem podle [1].")

    # Definice krychle: mřížka 2x2x2 (jediný element s jednotkovou geometrií)
    A = [0.0, 0.0, 0.0]
    B = [1.0, 0.0, 0.0]
    C = [1.0, 1.0, 0.0]
    D = [0.0, 1.0, 0.0]
    E = [0.0, 0.0, 1.0]
    F = [1.0, 0.0, 1.0]
    G = [1.0, 1.0, 1.0]
    H = [0.0, 1.0, 1.0]

    # Definice referenční geometrie modifikované A15 dlaždice (podle Appendix A z [1])
    tile_ref = [
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
        SVector(5.0, 4.0, 0.0)    # 26 ok - kontrolováno
    ]

    # Definice konektivity (46 tetraedrů) – indexy jsou 0-based v literatuře, proto je převedeme na 1-based.
    tetra_connectivity = [
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
        [19,27,16,4], # 29
        [17,12,18,19], # 30
        [14,9,11,13], # 31
        [14,12,11,10], # 32
        [3,8,14,13], # 33
        [6,17,7,10], # 34
        [5,9,20,10], # 35
        [14,9,7,10],
        [11,9,10,20],
        [7,9,5,10],
        [17,6,23,21],
        [6,7,5,10],
        [15,12,11,14],
        [16,14,2,7],
        [7,17,24,16],
        [26,24,16,25],
        [14,16,17,7] # 45
    ]

    # Převedeme referenční dlaždicové souřadnice do parametrických (r,s,t)
    # V referenční dlaždici jsou x ∈ [0,5], y ∈ [0,4], z ∈ [0,4] 
    # ⇒ r = x/5, s = y/4, t = z/4.
    tile_vertices = Vector{Vector{Float64}}(undef, length(tile_ref))
    for i in 1:length(tile_ref)
        x_tile, y_tile, z_tile = tile_ref[i]  # Opraveno: místo "tile_ref[i]..." používáme "tile_ref[i]"
        r = x_tile / 5.0
        s = y_tile / 4.0
        t = z_tile / 4.0
        tile_vertices[i] = trilinear_map(r, s, t, A, B, C, D, E, F, G, H)
    end

    # Pomocná funkce pro orientovaný objem tetraedru
    tet_volume(a, b, c, d) = dot(d .- a, cross(b .- a, c .- a)) / 6.0

    tetra_volumes = Vector{Float64}()
    for (tet_idx, tet) in enumerate(tetra_connectivity)
        a = tile_vertices[tet[1]]
        b = tile_vertices[tet[2]]
        c = tile_vertices[tet[3]]
        d = tile_vertices[tet[4]]
        vol = tet_volume(a, b, c, d)
        println("Tetraedr $tet_idx s uzly $tet má objem: $vol")
        push!(tetra_volumes, vol)
    end

    # Export do VTK
    npoints = length(tile_vertices)
    points = zeros(Float64, 3, npoints)
    for i in 1:npoints
        points[:, i] = tile_vertices[i]
    end
    cells = [MeshCell(VTKCellTypes.VTK_TETRA, tet) for tet in tetra_connectivity]
    vtkfile = vtk_grid("A15_modified_cell.vtu", points, cells)
    vtk_point_data(vtkfile, ones(npoints), "dummy_scalar")
    vtk_save(vtkfile)
    println("Export do VTK dokončen: A15_modified_cell.vtu")
    return tetra_volumes
end


# Spuštění testu
@time tetra_volumes = test_single_cell_A15_modified()
