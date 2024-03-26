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

x = [0., 4., -0.1]

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

    ρₑ = [0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0]
    # ρₑ = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0] # 1, 7
    # ρₑ = [1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0] # dve roviny 4 uzlové
    # ρₑ = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1] # 8
    # ρₑ = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0] # 1-7
    # ρₑ = [0.0, 0.0, 0.0, 0.0, 1, 1, 1, 1]
    # ρₑ = [0.0, 0.0, 0.0, 0.0, 0.8, 0.8, 1, 1]

    nsd = 3

        # ρ = [ρₑ[1], ρₑ[2], ρₑ[4], ρₑ[3], ρₑ[5], ρₑ[6], ρₑ[8], ρₑ[7]] .-ρₜ   
        # Y = reshape(ρ, (2, 2, 2)) # Array for MC algorithm
        
        # mc = MC(Y, Int; x = [-1., 1], y = [-1., 1], z = [-1., 1]) # cube inicialization (MC not working with general mesh, but with grid)
        # fieldnames(typeof(mc))
        # # Perform MC:
        # march(mc)

        # mc.cube
        # triangles = mc.triangles
        # vertices = mc.vertices
        # mc.normals
        # mc.nrm
        # mc.normal_sign
        # mc.x

        # triangle = triangles[1]


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


function Normal2Tri(
    vertices::Vector{Vector{Float64}}, # More specific type for vertices
    triangle,
)
    # Extracting vertices based on triangle indices
    A = vertices[triangle[1]]
    B = vertices[triangle[2]]
    C = vertices[triangle[3]]

    # Calculating the normal
    N = cross(B - A, C - A)
    N_unit = normalize(N)

    centre = (A + B + C) ./ 3

    return (N_unit, centre)
end

function findClosestPoint(Xe::Vector{Vector{Float64}}, center::Vector{Float64}, normal::Vector{Float64}, X::Vector{Float64})
    direction = dot(normal, X - center) >= 0 ? 1 : -1  # Určení směru hledání
    min_distance = Inf
    closest_index = 0
    distances = zeros(length(Xe))

    for (index, point) in enumerate(Xe)
        distance = dot(normal, point - center)
        distances[index] = distance # Vypočet vzdálenosti ve správném směru
        if direction > 0
            # Hledání minimální kladné vzdálenosti
            if distance >= 0 && distance < min_distance
                min_distance = distance
                closest_index = index
            end
        else
            # Hledání minimální záporné vzdálenosti (v absolutní hodnotě)
            if distance <= 0 && abs(distance) < abs(min_distance)
                min_distance = distance  # Zde ukládáme zápornou vzdálenost
                closest_index = index
            end
        end
    end

    return (closest_index, distances)
end

function SignDetermination(
    closest_index::Int,
    ρₑ::Vector{Float64},
    ρₜ::Float64
)
    S = ρₑ[closest_index] >= ρₜ ? 1 : -1
    return S
end

function Sign4Triangles(
    vertices_coords::Vector{Vector{Float64}}, # More specific type
    triangles,
    Xe::Vector{Vector{Float64}},
    x::Vector{Float64},
)

    # normals = Vector{Vector{Float64}}(undef, length(triangles))
    # centres = Vector{Vector{Float64}}(undef, length(triangles))
    sign4Triangles = Vector{Int}(undef, length(triangles))

    for i in 1:length(triangles)
        triangle = triangles[i]
        (normal, center) = Normal2Tri(vertices_coords, triangle)

        (closest_index, distances) = findClosestPoint(Xe, center, normal, x)
        sign4Triangles[i] = SignDetermination(closest_index, ρₑ, ρₜ)
    end
    #TODO: VÝSTUPEM JE PŘIŘAZENÍ ZNAMÉNKA KE KAŽDÉMU TROJÚHELNÍKU

    return sign4Triangles
end

sign4Triangles = Sign4Triangless(vertices_coords, triangles, Xₑ, x)

function findUniqueEdgesAndSigns(triangles, sign4Triangles)
    edgesDict = Dict{Tuple{Int, Int}, Tuple{Int, Int}}()

    for (i, triangle) in enumerate(triangles)
        # Convert SVector to a regular array for sorting purposes
        sortedTriangle = sort(collect(triangle))

        # Assign the current triangle's sign
        currentSign = sign4Triangles[i]

        # Create edges as sorted pairs
        edges = [(sortedTriangle[1], sortedTriangle[2]), (sortedTriangle[2], sortedTriangle[3]), (sortedTriangle[3], sortedTriangle[1])]

        for edge in edges
            if haskey(edgesDict, edge)
                # If edge already exists, just update the count
                edgesDict[edge] = (edgesDict[edge][1] + 1, edgesDict[edge][2])
            else
                # Else, add the new edge with its count set to 1 and its sign
                edgesDict[edge] = (1, currentSign)
            end
        end
    end

    uniqueEdges = Vector{Tuple{Int, Int}}()
    correspondingSigns = Vector{Int}()

    for (edge, data) in edgesDict
        if data[1] == 1
            push!(uniqueEdges, edge)
            push!(correspondingSigns, data[2])
        end
    end

    return uniqueEdges, correspondingSigns
end

# Example usage

(uniqueEdges, correspondingSigns) = findUniqueEdgesAndSigns(triangles, sign4Triangles)
uniqueEdges
correspondingSigns


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

edges2faceID = find_matching_ISE_indices(uniqueEdges, ISE)

function check_for_zero(vector::Vector{Int})
    if 0 in vector
        println("Error: The vector contains a zero value.")
    end
end

# Example usage:
# check_for_zero(edges2faceID)


function IsocontourEdgesForElement(ρₑ::Vector, ρₜ::Float64, edges, Xₑ, ISE, nsd, x)
    # ISE = mesh.ISE
    # nsd = mesh.nsd

    # Cutting triangles single cube element based on nodal densities (MS algorithm):
    (triangles, vertices) =  MC_OnCube(ρₑ, ρₜ) # POTŘEBUJI

    # Transform vertices coords from cube to real element:
    (vertices_coords, ID_edges) = EdgeIdDetection(vertices, Xₑ, edges, nsd) # POTŘEBUJI vertices_coords

    sign4Triangles = Sign4Triangles(vertices_coords, triangles, Xₑ, x)

    # (uniqueEdges, correspondingNormals, correspondingCentres) = findUniqueEdgesAndNormals(triangles, normals, centres)
    (uniqueEdges, correspondingSigns) = findUniqueEdgesAndSigns(triangles, sign4Triangles)

    # Check if the vertices connections forming meaningful pairs and if so assigned ID of face
    edges2faceID = find_matching_ISE_indices(uniqueEdges, ISE)

    real_vert_connections = [uniqueEdges[i] for i in 1:length(edges2faceID) if edges2faceID[i] != 0]
    real_Signs = [correspondingSigns[i] for i in 1:length(edges2faceID) if edges2faceID[i] != 0]

    check_for_zero(edges2faceID)

    # return (edges2faceID, ID_edges, vertices_coords, vert_connections)
    return (vertices_coords, real_vert_connections, real_Signs)

    #WARNING: Pozor na inverzní řešení. V tom případě prohodit hustotu.
end

(vertices_coords, real_vert_connections, real_Signs) = IsocontourEdgesForElement(ρₑ, ρₜ, edges, Xₑ, ISE, nsd, x)
vertices_coords
real_vert_connections
correspondingNormals
correspondingCentres

edges2faceID[2] = 0

real_vert_connections = [uniqueEdges[i] for i in 1:length(edges2faceID) if edges2faceID[i] != 0]
real_Signs = [correspondingSigns[i] for i in 1:length(edges2faceID) if edges2faceID[i] != 0]


function findClosestPoint(Xe::Vector{Vector{Float64}}, center::Vector{Float64}, normal::Vector{Float64}, X::Vector{Float64})
    direction = dot(normal, X - center) >= 0 ? 1 : -1  # Určení směru hledání
    min_distance = Inf
    closest_index = 0
    distances = zeros(length(Xe))

    for (index, point) in enumerate(Xe)
        distance = dot(normal, point - center)
        distances[index] = distance # Vypočet vzdálenosti ve správném směru
        if direction > 0
            # Hledání minimální kladné vzdálenosti
            if distance >= 0 && distance < min_distance
                min_distance = distance
                closest_index = index
            end
        else
            # Hledání minimální záporné vzdálenosti (v absolutní hodnotě)
            if distance <= 0 && abs(distance) < abs(min_distance)
                min_distance = distance  # Zde ukládáme zápornou vzdálenost
                closest_index = index
            end
        end
    end

    return (closest_index, distances)
end

(closest_index, distances) = findClosestPoint(Xe, center, normal)


function SignDetermination(closest_index, ρₑ, ρₜ)
    if ρₑ[closest_index] >= ρₜ 
        S = 1
    else
        S = -1
    end
    return S
end

SignDetermination(closest_index, ρₑ, ρₜ)

0 <= 1
point_vector = Vector{Float64}([-0.5, 0.5, 0.5])
[-0.5, 0.5, 0.5]

# Data
Xe = [ [-0.5, -0.5, -0.5],
       [0.5, -0.5, -0.5],
       [0.5, 0.5, -0.5],
       [-0.5, 0.5, -0.5],
       [-0.5, -0.5, 0.5],
       [0.5, -0.5, 0.5],
       [0.5, 0.5, 0.5],
       [-0.5, 0.5, 0.5] ]

center = [0., 0., 0.]
normal = [0.5773502691896258, 0.5773502691896258, 0.5773502691896258]
normal = [0., 0., -1]

X = [0., 0., -0.1]

# Najděte index nejbližšího bodu
(closest_index, distances) = findClosestPoint(Xe, center, normal, X)

println("Index nejbližšího bodu je: ", closest_index)

Xe = [-1.0 1.0 1.0 -1.0 -1.0 1.0 1.0 -1.0; -1.0 -1.0 1.0 1.0 -1.0 -1.0 1.0 1.0; -1.0 -1.0 -1.0 -1.0 1.0 1.0 1.0 1.0]

function findClosestPoint(
    Xe::Matrix{Float64},
    center::Vector{Float64},
    normal::Vector{Float64},
    x::Vector{Float64}
)
    println(size(Xe))
    println(Xe)
    direction = dot(normal, x - center) >= 0 ? 1 : -1  # Určení směru hledání
    min_distance = Inf
    closest_index = 0
    distances = zeros(size(Xe, 2))   # Vytvoření pole vzdáleností pro každý bod
    #WARNING: Xe je špatně použité

    for index in 1:size(Xe, 2)
        point = Xe[:, index]  # Získání bodu jako vektoru
        distance = dot(normal, point - center)
        distances[index] = distance  # Vypočet a uložení vzdálenosti
        if direction > 0
            # Hledání minimální kladné vzdálenosti
            if distance >= 0 && distance < min_distance
                min_distance = distance
                closest_index = index
            end
        else
            # Hledání minimální záporné vzdálenosti (v absolutní hodnotě)
            if distance <= 0 && abs(distance) < abs(min_distance)
                min_distance = distance  # Zde ukládáme zápornou vzdálenost
                closest_index = index
            end
        end
    end

    return (closest_index, distances)
end
(closest_index, distances) = findClosestPoint(Xe, [0., 1., 1.], [0., 1., 0.1], x)
distances = zeros(size(Xe, 2))  # Vytvoření pole vzdáleností pro každý bod
Xe[1, :] 