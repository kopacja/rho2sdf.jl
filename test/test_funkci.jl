using LinearAlgebra
using Statistics

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