using LinearAlgebra
using Statistics
using GeometryBasics
using MarchingCubes
using Combinatorics
using BenchmarkTools

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

    ISE = ((1, 2, 3, 4),
        (1, 9, 5, 10),
        (2, 11, 6, 10),
        (3, 12, 7, 11),
        (4, 12, 8, 9),
        (5, 6, 7, 8))

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
have_common_component(vec1, vec2) = any(x -> x in vec2, vec1)

number_of_tri = length(triangles) # number of all triangle in the element
common_tri = Int[]                # inicialization for commons triangles
common_tri2 = Int[]               # -||- in one element there can be max 2 faces with multiple triangles
tri_prev = triangles[1]           # Assuming there's at least one triangle

# ID of triangles forming sets
for i in 2:number_of_tri
    tri = triangles[i]
    if have_common_component(tri_prev, tri)
        if isempty(common_tri) || last(common_tri) == i-1
            push!(common_tri, i-1, i)
        else
            push!(common_tri2, i-1, i)
        end
    end
    tri_prev = tri
end
tri_prev
common_tri
common_tri2
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
# vert_connections
# vertices_coords
# ID_edges
#
# Dále budu pro každý element procházet všechny možné konektivity hran (budu se dívat do tabulky segment-edge)
# pokud konektivita dává smysl vypočítám hranu na kterou se bude promítat vzdálenostní funkce


function find_matching_ISE_indices(vert_connections, ISE)
    edges2faceID = zeros(Int, length(vert_connections)) # Ensure it's Int to store indices
    
    for (i, vc) in enumerate(vert_connections)
        vc_set = Set(vc) # Convert to Set for efficient intersection calculation
        for (j, ise) in enumerate(ISE)
            # Convert ise to Set only once per iteration for efficiency
            ise_set = Set(ise)
            if length(intersect(ise_set, vc_set)) == 2
                edges2faceID[i] = j
                break # Break early since we found a match
            end
        end
    end
    
    return edges2faceID
end

edges2faceID = find_matching_ISE_indices(vert_connections, ISE)
ID_edges
vertices_coords
vert_connections

x = [2.1, 2.2, 3.3]

for i in eachindex(edges2faceID)
    if edges2faceID[i] != 0



##________________________

# function compute_ratio_and_edge(vert, edges_conditions)
#     for (edge_ID, conditions) in enumerate(edges_conditions)
#         x_cond, y_cond, z_cond, ratio_index = conditions
#         if vert[x_cond[1]] == x_cond[2] && vert[y_cond[1]] == y_cond[2] && vert[z_cond[1]] == z_cond[2]
#             ratio = (vert[ratio_index] + 1) / 2
#             return edge_ID, ratio
#         end
#     end
#     println("Error finding edge")
#     return nothing, nothing
# end

# function compute_vertices_positions(vertices, edges, Xₑ)
#     nsd = size(Xₑ[1], 1) # Assuming Xₑ is a collection of coordinate arrays for nodes
#     number_of_vert = length(vertices)
#     vertices_coords = [zeros(Float64, nsd) for _ in 1:number_of_vert]
#     ID_edges = zeros(Int, number_of_vert)

#     edges_conditions = [
#         ((2, -1.0), (3, -1.0), (1, nothing), 1), # Edge 1
#         ((1,  1.0), (3, -1.0), (2, nothing), 2), # Edge 2
#         ((2,  1.0), (3, -1.0), (1, nothing), 2), # Edge 3
#         ((1, -1.0), (3, -1.0), (2, nothing), 2), # Edge 4
#         ((2, -1.0), (3,  1.0), (1, nothing), 2), # Edge 5
#         ((1,  1.0), (3,  1.0), (2, nothing), 2), # Edge 6
#         ((2,  1.0), (3,  1.0), (1, nothing), 2), # Edge 7
#         ((1, -1.0), (3,  1.0), (2, nothing), 2), # Edge 8
#         ((1, -1.0), (2, -1.0), (3, nothing), 2), # Edge 9
#         ((1,  1.0), (2, -1.0), (3, nothing), 2), # Edge 10
#         ((1,  1.0), (2,  1.0), (3, nothing), 2), # Edge 11
#         ((1, -1.0), (2,  1.0), (3, nothing), 2), # Edge 12
#         # Add other edges conditions following the pattern
#     ]

#     for i in 1:number_of_vert
#         vert = vertices[i]
#         edge_ID, ratio = compute_ratio_and_edge(vert, edges_conditions)
#         if edge_ID !== nothing
#             node₁, node₂ = edges[edge_ID]
#             vert_coords = Xₑ[node₁] .+ (Xₑ[node₂] .- Xₑ[node₁]).*ratio
#             vertices_coords[i] = vert_coords
#             ID_edges[i] = edge_ID
#         end
#     end

#     return vertices_coords, ID_edges
# end

# # Example usage (assuming `vertices`, `edges`, and `Xₑ` are defined)
# vertices_coords, ID_edges = compute_vertices_positions(vertices, edges, Xₑ)

#_________________________

        


