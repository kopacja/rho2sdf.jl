using LinearAlgebra
using Statistics
using GeometryBasics
using MarchingCubes
using Combinatorics

ρₜ = 0.5 # Threshold density (isosurface level)

Xₑ = [
                [-1.0, -1.0, -1.0],
                [1.0, -1.0, -1.0],
                [1.0, 1.0, -1.0],
                [-1.0, 1.0, -1.0],
                [-1.0, -1.0, 1.0],
                [1.0, -1.0, 1.0],
                [1.0, 1.0, 1.0],
                [-1.0, 1.0, 1.0],
            ]
            IEN = [[1, 2, 3, 4, 5, 6, 7, 8]]
            ρₑ = [0.0, 0.0, 0.0, 0.0, 1, 1, 1, 1]

edges = ((1,2), (2,3), (3,4), (4,1), 
        (5,6), (6,7), (7,8), (8,5), 
        (1,5), (2,6), (3,7), (4,8))

        
function ProjectionIntoIsocontourVertices(
            mesh::Mesh,
            ρₑ::Vector{Float64},
            ρₜ::Float64,
            Xₑ::Matrix)

    edges = mesh.edges
    nes = mesh.nes
    NodeOnEdge = zeros(Float64, (length(edges), 3))

    i = 0
    for edge in edges
        i = i + 1
        vrt1 = edge[1]
        vrt2 = edge[2]
        ρ_min, min_idx = findmin([ρₑ[vrt1], ρₑ[vrt2]])
        ρ_max, max_idx = findmax([ρₑ[vrt1], ρₑ[vrt2]])

        a_min = [vrt1, vrt2][min_idx]
        a_max = [vrt1, vrt2][max_idx]
        if (ρ_min <= ρₜ && ρ_max >= ρₜ)

            ratio = (ρₜ - ρ_min) / (ρ_max - ρ_min)
            xₚ = Xₑ[:, a_min] + ratio .* (Xₑ[:, a_max] - Xₑ[:, a_min])
            NodeOnEdge[i, :] = xₚ
        end
    end
    ISE = ((1, 2, 3, 4),
        (1, 9, 5, 10),
        (2, 11, 6, 10),
        (3, 12, 7, 11),
        (4, 12, 8, 9),
        (5, 6, 7, 8))

    vector_of_vector_pairs = [Vector{Float64}[] for _ in 1:nes]


    for i in 1:nes
        for j in 1:4
            vector_of_vector_pairs[i][j] = NodeOnEdge[ISE[i][j],:]
        end
    end
    return vector_of_vector_pairs
end


(dρ_dΞ, d²ρ_dΞ², d³ρ_dΞ³) = ρ_derivatives(ρₑ, [0.0, 0.0, 0.0])
norm_dρ_dΞ = norm(dρ_dΞ)
n = dρ_dΞ / norm_dρ_dΞ

for i in 1:nes
    inside, xp = ProjOnIsoEdge(vector_of_vector_pairs[i], x)
    if inside 
        dist_tmp = sign(dot(x - xₚ, n)) * norm(x - xₚ)
        if (abs(dist_tmp) < abs(dist[v]))
            dist[v] = dist_tmp
            xp[:, v] = xₚ
        end
    end
end



function ProjOnIsoEdge(pairs::Vector, x::Vector{Float64})
    a = pairs[1]
    b = pairs[2]
    p = x
    
    # Calculate the directional vector of the line segment AB
    ab = b - a
    # Calculate the vector AP
    ap = p - a
    
    # Project P onto AB to find the projection vector AP_proj
    dp = dot(ap, ab)
    ab2 = dot(ab, ab)
    P_proj = dp / ab2

    proj = a + P_proj * ab
    inside = 0 <= P_proj <= 1
    # The projection is inside the segment if the scalar is between 0 and 1 (inclusive)
    return inside, proj
end

vector_of_vector_pairs = [Vector{Float64}[] for _ in 1:nes]
NodeOnEdge = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0; -1.0 -1.0 0.0; 1.0 -1.0 0.0; 1.0 1.0 0.0; -1.0 1.0 0.0]

    # println("NodeOnEdge", NodeOnEdge)
    # println("vector_of_vector_pairs", vector_of_vector_pairs)
    for i in 1:nes
        for j in 1:4
            # if mean(abs.(NodeOnEdge[ISE[i][j])) != 0.
            if NodeOnEdge[ISE[i][j]] != [0., 0., 0.]
                vector_of_vector_pairs[i][j] = NodeOnEdge[ISE[i][j],:]
            end
        end
    end

    vector_of_vector_pairs[1] = [0., 0., 0.]
    vector_of_vector_pairs[1]
    NodeOnEdge

    ρₜ = 0.5 # Threshold density (isosurface level)

    Xₑ = [
                    [-1.0, -1.0, -1.0],
                    [1.0, -1.0, -1.0],
                    [1.0, 1.0, -1.0],
                    [-1.0, 1.0, -1.0],
                    [-1.0, -1.0, 1.0],
                    [1.0, -1.0, 1.0],
                    [1.0, 1.0, 1.0],
                    [-1.0, 1.0, 1.0],
                ]

        ISE = ((1, 2, 3, 4),
                (1, 9, 5, 10),
                (2, 11, 6, 10),
                (3, 12, 7, 11),
                (4, 12, 8, 9),
                (5, 6, 7, 8))

                IEN = [[1, 2, 3, 4, 5, 6, 7, 8]]
                ρₑ = [0.0, 0.0, 0.0, 0.0, 1, 1, 1, 1]
    
    edges = ((1,2), (2,3), (3,4), (4,1), 
            (5,6), (6,7), (7,8), (8,5), 
            (1,5), (2,6), (3,7), (4,8))

    # ρₑ = [0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0]
    # ρₑ = [1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
    ρₑ = [1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0]

    nsd = 3

ρ = [ρₑ[1], ρₑ[2], ρₑ[4], ρₑ[3], ρₑ[5], ρₑ[6], ρₑ[8], ρₑ[7]] .-ρₜ   
Y = reshape(ρ, (2, 2, 2))

mc1 = MC(Y, Int; x = [-1., 1], y = [-1., 1], z = [-1., 1]) # cube inicialization (MC not working with general mesh, but with grid)
    
# Perform MC:
march(mc1)
triangles = mc1.triangles
vertices = mc1.vertices

"""
- zjistit kterou hranu mi to protíná
- dopočítat ratio abych zjistil přesnou pozici protnutí 

- Co vlastně potřebuji jsou ty trojúhelníky které se nedotýkají
  
"""
# Every triangle is defined by three nodes (ID node) 
# Im looking if the triangles touching eachother (if they have same/s)
function have_common_component(vec1, vec2)
    return !isempty(intersect(vec1, vec2))
end

number_of_tri = length(triangles) # number of all triangle in the element
common_tri = Int[]                # inicialization for commons triangles
common_tri2 = Int[]               # -||- in one element there can be max 2 faces with multiple triangles
tri_prev = triangles[1]           # Assuming there's at least one triangle
is_common = false                 # for seperating triangle sets

# ID of triangles forming sets
for i in 2:number_of_tri
    tri = triangles[i]
    if !is_common && have_common_component(tri_prev, tri) # první plocha
        common_tri = vcat(common_tri, i-1, i)
    elseif is_common && have_common_component(tri_prev, tri) # druhá plocha (jsou maximálně dvě v elementu)
        common_tri2 = vcat(common_tri2, i-1, i)
    else
        is_common = true
    end
    tri_prev = tri
end

common_tri = unique!(common_tri)    # removing duplicates
common_tri2 = unique!(common_tri2)  # removing duplicates
# ok funguje dobře (testováno)

# Combinations of possible connestion of vertices (ID of vertices)
# for commons triangles
common_tri_uniq_vec1 = unique(reduce(vcat, triangles[common_tri]))
common_comb1 = collect(combinations(common_tri_uniq_vec1,2))

common_tri_uniq_vec2 = unique(reduce(vcat, triangles[common_tri2]))
common_comb2 = collect(combinations(common_tri_uniq_vec2,2))

## Find triangles that are alone
v = collect(1:number_of_tri)
uncommon_tri = setdiff(v, vcat(common_tri, common_tri2))
# Combinations of possible connestion of vertices 
uncommon_comb = Int[]
for i in eachindex(uncommon_tri)
    pairs = collect(combinations(collect(triangles[uncommon_tri[i]]),2))
    uncommon_comb = vcat(uncommon_comb, pairs)
end

# All possible connestions:
vert_connections = vcat(uncommon_comb, common_comb1, common_comb2) ## RESULT (ID of vertices)
# Does nodes correlate with my element notation?

# Compute real position of vertices (MC works with cubes, I have general 8 nodes elements)
number_of_vert = length(vertices)
vertices_coords = [zeros(Float64, nsd) for _ in 1:number_of_vert]
ID_edges = zeros(Int, number_of_vert)

# To be able to compute position, I need to know on which edge am I:
for i in 1:number_of_vert
    vert = vertices[i]
    # Cases (which edge?)
    if vert[2] == -1. && vert[3] == -1. # edge 1
        edge_ID = 1
        node₁, node₂ = edges[edge_ID]
        ratio = (vert[1] + 1) / 2
        vert_coords =  Xₑ[node₁] .+ (Xₑ[node₂] .- Xₑ[node₁]).*ratio

    elseif vert[1] == 1. && vert[3] == -1. # edge 2
        edge_ID = 2
        node₁, node₂ = edges[edge_ID] 
        ratio = (vert[2] + 1) / 2
        vert_coords =  Xₑ[node₁] .+ (Xₑ[node₂] .- Xₑ[node₁]).*ratio

    elseif vert[2] == 1. && vert[3] == -1. # edge 3
        edge_ID = 3
        node₁, node₂ = edges[edge_ID] 
        ratio = (vert[1] + 1) / 2
        vert_coords =  Xₑ[node₁] .+ (Xₑ[node₂] .- Xₑ[node₁]).*ratio

    elseif vert[1] == -1. && vert[3] == -1. # edge 4
        edge_ID = 4
        node₁, node₂ = edges[edge_ID] 
        ratio = (vert[2] + 1) / 2
        vert_coords =  Xₑ[node₁] .+ (Xₑ[node₂] .- Xₑ[node₁]).*ratio

    elseif vert[2] == -1. && vert[3] == 1. # edge 5
        edge_ID = 5
        node₁, node₂ = edges[edge_ID] 
        ratio = (vert[1] + 1) / 2
        vert_coords =  Xₑ[node₁] .+ (Xₑ[node₂] .- Xₑ[node₁]).*ratio

    elseif vert[1] == 1. && vert[3] == 1. # edge 6
        edge_ID = 6
        node₁, node₂ = edges[edge_ID] 
        ratio = (vert[2] + 1) / 2
        vert_coords =  Xₑ[node₁] .+ (Xₑ[node₂] .- Xₑ[node₁]).*ratio

    elseif vert[2] == 1. && vert[3] == 1. # edge 7
        edge_ID = 7
        node₁, node₂ = edges[edge_ID] 
        ratio = (vert[1] + 1) / 2
        vert_coords =  Xₑ[node₁] .+ (Xₑ[node₂] .- Xₑ[node₁]).*ratio

    elseif vert[1] == -1. && vert[3] == 1. # edge 8
        edge_ID = 8
        node₁, node₂ = edges[edge_ID] 
        ratio = (vert[2] + 1) / 2
        vert_coords =  Xₑ[node₁] .+ (Xₑ[node₂] .- Xₑ[node₁]).*ratio

    elseif vert[1] == -1. && vert[2] == -1. # edge 9
        edge_ID = 9
        node₁, node₂ = edges[edge_ID] 
        ratio = (vert[3] + 1) / 2
        vert_coords =  Xₑ[node₁] .+ (Xₑ[node₂] .- Xₑ[node₁]).*ratio

    elseif vert[1] == 1. && vert[2] == -1. # edge 10
        edge_ID = 10
        node₁, node₂ = edges[edge_ID] 
        ratio = (vert[3] + 1) / 2
        vert_coords =  Xₑ[node₁] .+ (Xₑ[node₂] .- Xₑ[node₁]).*ratio

    elseif vert[1] == 1. && vert[2] == 1. # edge 11
        edge_ID = 11
        node₁, node₂ = edges[edge_ID] 
        ratio = (vert[3] + 1) / 2
        vert_coords =  Xₑ[node₁] .+ (Xₑ[node₂] .- Xₑ[node₁]).*ratio

    elseif vert[1] == -1. && vert[2] == 1. # edge 12
        edge_ID = 12
        node₁, node₂ = edges[edge_ID] 
        ratio = (vert[3] + 1) / 2
        vert_coords =  Xₑ[node₁] .+ (Xₑ[node₂] .- Xₑ[node₁]).*ratio
    else
        println("chyba při nalezení hrany")

    end

    # Checking if the result makes sense
    if Xₑ[node₁] <= vert_coords <= Xₑ[node₂] || Xₑ[node₂] <= vert_coords <= Xₑ[node₁] 
        # println("ok") 
    else
        println("chyba")
    end
    # Putting it all together:
    vertices_coords[i] = vert_coords
    ID_edges[i] = edge_ID

end
# RESULT:
# vertices_coords
# ID_edges
#
# Dále budu pro každý element procházet všechny možné konektivity hran (budu se dívat do tabulky segment-edge)
# pokud konektivita dává smysl vypočítám hranu na kterou se bude promítat vzdálenostní funkce
