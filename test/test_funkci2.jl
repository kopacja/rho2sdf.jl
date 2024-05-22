using MarchingCubes
using LinearAlgebra
using StaticArrays

ρₜ = 0.5 # Threshold density (isosurface level)

Xₑᵥ = [
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

edges = ((1, 2), (2, 3), (3, 4), (4, 1),
  (5, 6), (6, 7), (7, 8), (8, 5),
  (1, 5), (2, 6), (3, 7), (4, 8))

ρₑ = [0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0]

nsd = 3

function convert_to_matrix(Xₑ::Vector)
  matrix = hcat(Xₑ...)
  return matrix
end
Xₑ = convert_to_matrix(Xₑᵥ)

function MC_OnCube(ρₑ::Vector{Float64},
  ρₜ::Float64)

  # Sort to match my elements definitions
  ρ = [ρₑ[1], ρₑ[2], ρₑ[4], ρₑ[3], ρₑ[5], ρₑ[6], ρₑ[8], ρₑ[7]] .- ρₜ
  Y = reshape(ρ, (2, 2, 2)) # Array for MC algorithm

  mc = MC(Y, Int; x=[-1.0, 1], y=[-1.0, 1], z=[-1.0, 1]) # cube inicialization (MC not working with general mesh, but with grid)

  # Perform MC:
  march(mc)
  return (mc.triangles, mc.vertices)
end

(triangles, vertices) = MC_OnCube(ρₑ, ρₜ)


function EdgeIdDetection1(
  vertices,
  Xₑ::Matrix,
  edges::NTuple,
  nsd::Int,
)

  number_of_vert = length(vertices)
  vertices_coords = [zeros(Float64, nsd) for _ in 1:number_of_vert]
  ID_edges = zeros(Int, number_of_vert)

  # To be able to compute position, I need to know on which edge am I:
  for i in 1:number_of_vert
    vert = vertices[i]
    # Cases (which edge?)
    if vert[2] == -1.0 && vert[3] == -1.0 # edge 1
      edge_ID = 1
      node₁, node₂ = edges[edge_ID]
      ratio = (vert[1] + 1) / 2
      vert_coords = Xₑ[:, node₁] .+ (Xₑ[:, node₂] .- Xₑ[:, node₁]) .* ratio

    elseif vert[1] == 1.0 && vert[3] == -1.0 # edge 2
      edge_ID = 2
      node₁, node₂ = edges[edge_ID]
      ratio = (vert[2] + 1) / 2
      vert_coords = Xₑ[:, node₁] .+ (Xₑ[:, node₂] .- Xₑ[:, node₁]) .* ratio

    elseif vert[2] == 1.0 && vert[3] == -1.0 # edge 3
      edge_ID = 3
      node₁, node₂ = edges[edge_ID]
      ratio = (vert[1] + 1) / 2
      vert_coords = Xₑ[:, node₁] .+ (Xₑ[:, node₂] .- Xₑ[:, node₁]) .* ratio

    elseif vert[1] == -1.0 && vert[3] == -1.0 # edge 4
      edge_ID = 4
      node₁, node₂ = edges[edge_ID]
      ratio = (vert[2] + 1) / 2
      vert_coords = Xₑ[:, node₁] .+ (Xₑ[:, node₂] .- Xₑ[:, node₁]) .* ratio

    elseif vert[2] == -1.0 && vert[3] == 1.0 # edge 5
      edge_ID = 5
      node₁, node₂ = edges[edge_ID]
      ratio = (vert[1] + 1) / 2
      vert_coords = Xₑ[:, node₁] .+ (Xₑ[:, node₂] .- Xₑ[:, node₁]) .* ratio

    elseif vert[1] == 1.0 && vert[3] == 1.0 # edge 6
      edge_ID = 6
      node₁, node₂ = edges[edge_ID]
      ratio = (vert[2] + 1) / 2
      vert_coords = Xₑ[:, node₁] .+ (Xₑ[:, node₂] .- Xₑ[:, node₁]) .* ratio

    elseif vert[2] == 1.0 && vert[3] == 1.0 # edge 7
      edge_ID = 7
      node₁, node₂ = edges[edge_ID]
      ratio = (vert[1] + 1) / 2
      vert_coords = Xₑ[:, node₁] .+ (Xₑ[:, node₂] .- Xₑ[:, node₁]) .* ratio

    elseif vert[1] == -1.0 && vert[3] == 1.0 # edge 8
      edge_ID = 8
      node₁, node₂ = edges[edge_ID]
      ratio = (vert[2] + 1) / 2
      vert_coords = Xₑ[:, node₁] .+ (Xₑ[:, node₂] .- Xₑ[:, node₁]) .* ratio

    elseif vert[1] == -1.0 && vert[2] == -1.0 # edge 9
      edge_ID = 9
      node₁, node₂ = edges[edge_ID]
      ratio = (vert[3] + 1) / 2
      vert_coords = Xₑ[:, node₁] .+ (Xₑ[:, node₂] .- Xₑ[:, node₁]) .* ratio

    elseif vert[1] == 1.0 && vert[2] == -1.0 # edge 10
      edge_ID = 10
      node₁, node₂ = edges[edge_ID]
      ratio = (vert[3] + 1) / 2
      vert_coords = Xₑ[:, node₁] .+ (Xₑ[:, node₂] .- Xₑ[:, node₁]) .* ratio

    elseif vert[1] == 1.0 && vert[2] == 1.0 # edge 11
      edge_ID = 11
      node₁, node₂ = edges[edge_ID]
      ratio = (vert[3] + 1) / 2
      vert_coords = Xₑ[:, node₁] .+ (Xₑ[:, node₂] .- Xₑ[:, node₁]) .* ratio

    elseif vert[1] == -1.0 && vert[2] == 1.0 # edge 12
      edge_ID = 12
      node₁, node₂ = edges[edge_ID]
      ratio = (vert[3] + 1) / 2
      vert_coords = Xₑ[:, node₁] .+ (Xₑ[:, node₂] .- Xₑ[:, node₁]) .* ratio
    else
      println("chyba při nalezení hrany")
      continue
    end

    # Checking if the result makes sense
    if Xₑ[:, node₁] <= vert_coords <= Xₑ[:, node₂] || Xₑ[:, node₂] <= vert_coords <= Xₑ[:, node₁]
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

(vertices_coords, ID_edges) = EdgeIdDetection1(vertices, Xₑ, edges, nsd)
vertices

triangles
vertices_coords
ID_edges


# Function to assign coordinates to each node
function assign_coordinates(triangles::Vector{SVector{3, Int}}, coordinates::Vector{Vector{Float64}})
    # Initialize the result as a vector of vectors of vectors
    result = Vector{Vector{Vector{Float64}}}(undef, length(triangles))

    # Iterate over each triangle
    for (i, triangle) in enumerate(triangles)
        # Initialize the coordinate vector for this triangle
        triangle_coords = Vector{Vector{Float64}}(undef, 3)
        
        # Map each node ID to its corresponding coordinates
        for (j, node_id) in enumerate(triangle)
            triangle_coords[j] = coordinates[node_id]
        end

        # Assign the coordinates to the result vector
        result[i] = triangle_coords
    end

    return result
end

# Call the function
result = assign_coordinates(triangles, vertices_coords)
