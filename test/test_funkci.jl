using LinearAlgebra
using Statistics
using GeometryBasics
using MarchingCubes
using Combinatorics
using BenchmarkTools


    ρₜ = 0.5 # Threshold density (isosurface level)

    Xₑ = [
                    [-4.1, -1.0, -1.0],
                    [1.0, -1.0, -1.0],
                    [1.0, 1.0, -1.0],
                    [-1.0, 1.0, -1.0],
                    [-1.0, -1.0, 1.0],
                    [1.0, -1.0, 1.0],
                    [1.0, 1.0, 1.0],
                    [-1.0, 1.0, 1.0],
                ] * 0.5

        ISE = ((1, 2, 3, 4),
                (1, 9, 5, 10),
                (2, 11, 6, 10),
                (3, 12, 7, 11),
                (4, 12, 8, 9),
                (5, 6, 7, 8))

                IEN = [[1, 2, 3, 4, 5, 6, 7, 8]]
    
    edges = ((1,2), (2,3), (3,4), (4,1), 
            (5,6), (6,7), (7,8), (8,5), 
            (1,5), (2,6), (3,7), (4,8))

    # ρₑ = [0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0]
    # ρₑ = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0] # 1, 7
    ρₑ = [1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0] # dve roviny 4 uzlové
    # ρₑ = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1] # 8

    nsd = 3

    ISE = ((1, 2, 3, 4),
        (1, 9, 5, 10),
        (2, 11, 6, 10),
        (3, 12, 7, 11),
        (4, 12, 8, 9),
        (5, 6, 7, 8))

function MC_OnCube(ρₑ::Vector{Float64},
                   ρₜ::Float64)

    # Sort to match my elements definitions
    ρ = [ρₑ[1], ρₑ[2], ρₑ[4], ρₑ[3], ρₑ[5], ρₑ[6], ρₑ[8], ρₑ[7]] .-ρₜ   
    Y = reshape(ρ, (2, 2, 2)) # Array for MC algorithm
    
    mc = MC(Y, Int; x = [-1., 1], y = [-1., 1], z = [-1., 1]) # cube inicialization (MC not working with general mesh, but with grid)
        
    # Perform MC:
    march(mc)
    return (mc.triangles, mc.vertices)
end

(triangles, vertices) =  MC_OnCube(ρₑ, ρₜ)
triangles
vertices

"""
- zjistit kterou hranu mi to protíná
- dopočítat ratio abych zjistil přesnou pozici protnutí 

- Co vlastně potřebuji jsou ty trojúhelníky které se nedotýkají
  
"""
# Every triangle is defined by three nodes (ID node) 
# Im looking if the triangles touching eachother (if they have same/s)
have_common_component(vec1, vec2) = any(x -> x in vec2, vec1)

# Adjusted function to handle empty input
function get_combinations(triangle_indices)
    if isempty(triangle_indices)
        return Int[]
    else
        unique_vertices = unique(reduce(vcat, triangles[triangle_indices]))
        return collect(combinations(unique_vertices, 2))
    end
end


function PossibleEdgeConnections(
    triangles::Vector,
)
    # ISE = mesh.ISE

    # Assuming `triangles` is an array of arrays, where each sub-array contains the vertex IDs of a triangle
    number_of_tri = length(triangles)
    common_tri = Int[]
    common_tri2 = Int[]
    tri_prev = triangles[1]
    
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
    # Since we're always adding pairs, the need for `unique!` is eliminated
    # However, if triangles can repeat in `triangles`, consider reintroducing `unique!`
    
    common_comb1 = get_combinations(common_tri)
    common_comb2 = get_combinations(common_tri2)
    
    # The rest of the code remains the same since it already handles empty cases for `uncommon_comb`
    all_tri_indices = 1:number_of_tri
    uncommon_tri = setdiff(all_tri_indices, union(common_tri, common_tri2))
    
    # Calculate uncommon combinations, handling the case where uncommon_comb might be empty
    uncommon_comb = [collect(combinations(triangles[i], 2)) for i in uncommon_tri]
    uncommon_comb_flat = isempty(uncommon_comb) ? Int[] : reduce(vcat, uncommon_comb) 
    
    # Combine all connections, already handling empty cases correctly
    vert_connections = vcat(uncommon_comb_flat, common_comb1, common_comb2)
    return vert_connections
end

vert_connections = PossibleEdgeConnections(triangles, vertices, ISE)



function EdgeIdDetection(
    vertices::Vector,
    #mesh::Mesh,
    Xₑ::Vector,
    edges::NTuple,
    nsd::Int,
    )

    # edges = mesh.edges
    # nsd = mesh.nsd

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
    return (vertices_coords, ID_edges)
end

(vertices_coords, ID_edges) = EdgeIdDetection(vertices, Xₑ, edges, nsd)

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


function IsocontourEdgesForElement(ρₑ::Vector, ρₜ::Float64, edges, Xₑ, ISE, nsd)
    # ISE = mesh.ISE
    # nsd = mesh.nsd

    # Cutting triangles single cube element based on nodal densities (MS algorithm):
    (triangles, vertices) =  MC_OnCube(ρₑ, ρₜ) # POTŘEBUJI

    # Possible connections of vertices to forming isocontour edges:
    vert_connections = PossibleEdgeConnections(triangles) # POTŘEBUJI

    # Transform vertices coords from cube to real element:
    (vertices_coords, ID_edges) = EdgeIdDetection(vertices, Xₑ, edges, nsd) # POTŘEBUJI vertices_coords


    # Check if the vertices connections forming meaningful pairs and if so assigned ID of face
    edges2faceID = find_matching_ISE_indices(vert_connections, ISE)

    real_vert_connections = [vert_connections[i] for i in 1:length(edges2faceID) if edges2faceID[i] != 0]

    # return (edges2faceID, ID_edges, vertices_coords, vert_connections)
    return (vertices_coords, real_vert_connections)

    #WARNING: Pozor na inverzní řešení. V tom případě prohodit hustotu.
end

(vertices_coords, real_vert_connections) = IsocontourEdgesForElement(ρₑ, ρₜ, edges, Xₑ, ISE, nsd)



    a, b = real_vert_connections[1]
    coords_a = vertices_coords[a]
    coords_b = vertices_coords[b]
bb = [coords_a, coords_b]

bb[1]
