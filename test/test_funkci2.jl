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

vertices
VertCoord = []
for i in 1:length(vertices)
    VertCoord = Xₑ * ShapeFunctionTemp(vertices[i])

end

VertCoord = Vector{Vector{Float64}}(undef, length(vertices))  # Preallocate VertCoord

for i in 1:length(vertices)
    VertCoord[i] = Xₑ * ShapeFunctionTemp(vertices[i])  # Directly assign the result to each index
end

size(veces)
function ShapeFunctionTemp(Ξ)
    ξ₁ = Ξ[1]
    ξ₂ = Ξ[2]
    ξ₃ = Ξ[3]

    # N = Array{Float64}(undef, 8)
    N = zeros(Float64, 8)
    N[1] = -1 / 8 * (ξ₁ - 1) * (ξ₂ - 1) * (ξ₃ - 1)
    N[2] =  1 / 8 * (ξ₁ + 1) * (ξ₂ - 1) * (ξ₃ - 1)
    N[3] = -1 / 8 * (ξ₁ + 1) * (ξ₂ + 1) * (ξ₃ - 1)
    N[4] =  1 / 8 * (ξ₁ - 1) * (ξ₂ + 1) * (ξ₃ - 1)

    N[5] =  1 / 8 * (ξ₁ - 1) * (ξ₂ - 1) * (ξ₃ + 1)
    N[6] = -1 / 8 * (ξ₁ + 1) * (ξ₂ - 1) * (ξ₃ + 1)
    N[7] =  1 / 8 * (ξ₁ + 1) * (ξ₂ + 1) * (ξ₃ + 1)
    N[8] = -1 / 8 * (ξ₁ - 1) * (ξ₂ + 1) * (ξ₃ + 1)
    return N
end




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

VertCoord == vertices_coords

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
length(result)
triangles


ρ = [ρₑ[1], ρₑ[2], ρₑ[4], ρₑ[3], ρₑ[5], ρₑ[6], ρₑ[8], ρₑ[7]] .- ρₜ

triangles == SVector([5, 6, 7],
[2, 3, 1],
[3, 2, 4])


Xt = [35.95015 37.24444 36.597295; 56.52723 53.90267 55.21495; -86.08696 -86.08696 -88.04348]
a = Xt[1,:]

X = [a, a, a]



using LinearAlgebra

# Definice tvarových funkcí a jejich derivací
function shape_functions(ξ)
    ξ₁, ξ₂, ξ₃ = ξ
    N = zeros(Float64, 8)
    N[1] = -1 / 8 * (ξ₁ - 1) * (ξ₂ - 1) * (ξ₃ - 1)
    N[2] = 1 / 8 * (ξ₁ + 1) * (ξ₂ - 1) * (ξ₃ - 1)
    N[3] = -1 / 8 * (ξ₁ + 1) * (ξ₂ + 1) * (ξ₃ - 1)
    N[4] = 1 / 8 * (ξ₁ - 1) * (ξ₂ + 1) * (ξ₃ - 1)
    N[5] = 1 / 8 * (ξ₁ - 1) * (ξ₂ - 1) * (ξ₃ + 1)
    N[6] = -1 / 8 * (ξ₁ + 1) * (ξ₂ - 1) * (ξ₃ + 1)
    N[7] = 1 / 8 * (ξ₁ + 1) * (ξ₂ + 1) * (ξ₃ + 1)
    N[8] = -1 / 8 * (ξ₁ - 1) * (ξ₂ + 1) * (ξ₃ + 1)
    return N
end

function shape_function_derivatives(ξ)
    ξ₁, ξ₂, ξ₃ = ξ
    dN_dξ = zeros(Float64, 8, 3)
    # Derivatives with respect to ξ₁
    dN_dξ[1, 1] = -1 / 8 * (ξ₂ - 1) * (ξ₃ - 1)
    dN_dξ[2, 1] = 1 / 8 * (ξ₂ - 1) * (ξ₃ - 1)
    dN_dξ[3, 1] = -1 / 8 * (ξ₂ + 1) * (ξ₃ - 1)
    dN_dξ[4, 1] = 1 / 8 * (ξ₂ + 1) * (ξ₃ - 1)
    dN_dξ[5, 1] = 1 / 8 * (ξ₂ - 1) * (ξ₃ + 1)
    dN_dξ[6, 1] = -1 / 8 * (ξ₂ - 1) * (ξ₃ + 1)
    dN_dξ[7, 1] = 1 / 8 * (ξ₂ + 1) * (ξ₃ + 1)
    dN_dξ[8, 1] = -1 / 8 * (ξ₂ + 1) * (ξ₃ + 1)
    # Derivatives with respect to ξ₂
    dN_dξ[1, 2] = -1 / 8 * (ξ₁ - 1) * (ξ₃ - 1)
    dN_dξ[2, 2] = 1 / 8 * (ξ₁ + 1) * (ξ₃ - 1)
    dN_dξ[3, 2] = -1 / 8 * (ξ₁ + 1) * (ξ₃ + 1)
    dN_dξ[4, 2] = 1 / 8 * (ξ₁ - 1) * (ξ₃ + 1)
    dN_dξ[5, 2] = 1 / 8 * (ξ₁ - 1) * (ξ₃ + 1)
    dN_dξ[6, 2] = -1 / 8 * (ξ₁ + 1) * (ξ₃ + 1)
    dN_dξ[7, 2] = 1 / 8 * (ξ₁ + 1) * (ξ₃ + 1)
    dN_dξ[8, 2] = -1 / 8 * (ξ₁ - 1) * (ξ₃ + 1)
    # Derivatives with respect to ξ₃
    dN_dξ[1, 3] = -1 / 8 * (ξ₁ - 1) * (ξ₂ - 1)
    dN_dξ[2, 3] = 1 / 8 * (ξ₁ + 1) * (ξ₂ - 1)
    dN_dξ[3, 3] = -1 / 8 * (ξ₁ + 1) * (ξ₂ + 1)
    dN_dξ[4, 3] = 1 / 8 * (ξ₁ - 1) * (ξ₂ + 1)
    dN_dξ[5, 3] = 1 / 8 * (ξ₁ - 1) * (ξ₂ - 1)
    dN_dξ[6, 3] = -1 / 8 * (ξ₁ + 1) * (ξ₂ - 1)
    dN_dξ[7, 3] = 1 / 8 * (ξ₁ + 1) * (ξ₂ + 1)
    dN_dξ[8, 3] = -1 / 8 * (ξ₁ - 1) * (ξ₂ + 1)
    return dN_dξ
end

# Newton-Raphson metoda
function find_local_coordinates(global_coords, node_coords)
    ξ = [0.0, 0.0, 0.0]  # Initial guess
    tolerance = 1e-8
    max_iterations = 100
    X = node_coords[:, 1]
    Y = node_coords[:, 2]
    Z = node_coords[:, 3]
    
    for iter in 1:max_iterations
        N = shape_functions(ξ)
        dN_dξ = shape_function_derivatives(ξ)
        
        # Compute global coordinates from shape functions
        x = dot(N, X)
        y = dot(N, Y)
        z = dot(N, Z)
        
        # Residuals
        R = [x - global_coords[1], y - global_coords[2], z - global_coords[3]]
        
        # Check convergence
        if norm(R) < tolerance
            return ξ
        end
        
        # Compute Jacobian
        J = [sum(dN_dξ[:, 1] .* X), sum(dN_dξ[:, 2] .* X), sum(dN_dξ[:, 3] .* X);
             sum(dN_dξ[:, 1] .* Y), sum(dN_dξ[:, 2] .* Y), sum(dN_dξ[:, 3] .* Y);
             sum(dN_dξ[:, 1] .* Z), sum(dN_dξ[:, 2] .* Z), sum(dN_dξ[:, 3] .* Z)]
        
        # Update ξ
        Δξ = - J \ R
        ξ += Δξ
    end
    
    error("Newton-Raphson method did not converge")
end

# Příklad použití
global_coords = [1.0, 1.0, 1.0]  # Zadejte globální souřadnice bodu
node_coords = [
    0.0 0.0 0.0;
    2.0 0.0 0.0;
    2.0 2.0 0.0;
    0.0 2.0 0.0;
    0.0 0.0 2.0;
    2.0 0.0 2.0;
    2.0 2.0 2.0;
    0.0 2.0 2.0
]

local_coords = find_local_coordinates(global_coords, node_coords)
println("Lokální souřadnice: $local_coords")
# ___________________________________________

using LinearAlgebra

# Definice tvarových funkcí a jejich derivací
function shape_functions(ξ)
    ξ₁, ξ₂, ξ₃ = ξ
    N = zeros(Float64, 8)
    N[1] = -1 / 8 * (ξ₁ - 1) * (ξ₂ - 1) * (ξ₃ - 1)
    N[2] = 1 / 8 * (ξ₁ + 1) * (ξ₂ - 1) * (ξ₃ - 1)
    N[3] = -1 / 8 * (ξ₁ + 1) * (ξ₂ + 1) * (ξ₃ - 1)
    N[4] = 1 / 8 * (ξ₁ - 1) * (ξ₂ + 1) * (ξ₃ - 1)
    N[5] = 1 / 8 * (ξ₁ - 1) * (ξ₂ - 1) * (ξ₃ + 1)
    N[6] = -1 / 8 * (ξ₁ + 1) * (ξ₂ - 1) * (ξ₃ + 1)
    N[7] = 1 / 8 * (ξ₁ + 1) * (ξ₂ + 1) * (ξ₃ + 1)
    N[8] = -1 / 8 * (ξ₁ - 1) * (ξ₂ + 1) * (ξ₃ + 1)
    return N
end

function shape_function_derivatives(ξ)
    ξ₁, ξ₂, ξ₃ = ξ
    dN_dξ = zeros(Float64, 8, 3)
    # Derivatives with respect to ξ₁
    dN_dξ[1, 1] = -1 / 8 * (ξ₂ - 1) * (ξ₃ - 1)
    dN_dξ[2, 1] = 1 / 8 * (ξ₂ - 1) * (ξ₃ - 1)
    dN_dξ[3, 1] = -1 / 8 * (ξ₂ + 1) * (ξ₃ - 1)
    dN_dξ[4, 1] = 1 / 8 * (ξ₂ + 1) * (ξ₃ - 1)
    dN_dξ[5, 1] = 1 / 8 * (ξ₂ - 1) * (ξ₃ + 1)
    dN_dξ[6, 1] = -1 / 8 * (ξ₂ - 1) * (ξ₃ + 1)
    dN_dξ[7, 1] = 1 / 8 * (ξ₂ + 1) * (ξ₃ + 1)
    dN_dξ[8, 1] = -1 / 8 * (ξ₂ + 1) * (ξ₃ + 1)
    # Derivatives with respect to ξ₂
    dN_dξ[1, 2] = -1 / 8 * (ξ₁ - 1) * (ξ₃ - 1)
    dN_dξ[2, 2] = 1 / 8 * (ξ₁ + 1) * (ξ₃ - 1)
    dN_dξ[3, 2] = -1 / 8 * (ξ₁ + 1) * (ξ₃ - 1)
    dN_dξ[4, 2] = 1 / 8 * (ξ₁ - 1) * (ξ₃ - 1)
    dN_dξ[5, 2] = 1 / 8 * (ξ₁ - 1) * (ξ₃ + 1)
    dN_dξ[6, 2] = -1 / 8 * (ξ₁ + 1) * (ξ₃ + 1)
    dN_dξ[7, 2] = 1 / 8 * (ξ₁ + 1) * (ξ₃ + 1)
    dN_dξ[8, 2] = -1 / 8 * (ξ₁ - 1) * (ξ₃ + 1)
    # Derivatives with respect to ξ₃
    dN_dξ[1, 3] = -1 / 8 * (ξ₁ - 1) * (ξ₂ - 1)
    dN_dξ[2, 3] = 1 / 8 * (ξ₁ + 1) * (ξ₂ - 1)
    dN_dξ[3, 3] = -1 / 8 * (ξ₁ + 1) * (ξ₂ + 1)
    dN_dξ[4, 3] = 1 / 8 * (ξ₁ - 1) * (ξ₂ + 1)
    dN_dξ[5, 3] = 1 / 8 * (ξ₁ - 1) * (ξ₂ - 1)
    dN_dξ[6, 3] = -1 / 8 * (ξ₁ + 1) * (ξ₂ - 1)
    dN_dξ[7, 3] = 1 / 8 * (ξ₁ + 1) * (ξ₂ + 1)
    dN_dξ[8, 3] = -1 / 8 * (ξ₁ - 1) * (ξ₂ + 1)
    return dN_dξ
end

# Newton-Raphson metoda
function find_local_coordinates(global_coords, node_coords)
    ξ = [0.0, 0.0, 0.0]  # Initial guess
    tolerance = 1e-8
    max_iterations = 100
    X = node_coords[:, 1]
    Y = node_coords[:, 2]
    Z = node_coords[:, 3]
    
    for iter in 1:max_iterations
        N = shape_functions(ξ)
        dN_dξ = shape_function_derivatives(ξ)
        
        # Compute global coordinates from shape functions
        x = dot(N, X)
        y = dot(N, Y)
        z = dot(N, Z)
        
        # Residuals
        R = [x - global_coords[1], y - global_coords[2], z - global_coords[3]]
        
        # Check convergence
        if norm(R) < tolerance
            return ξ
        end
        
        # Compute Jacobian
        J = [sum(dN_dξ[:, 1] .* X) sum(dN_dξ[:, 2] .* X) sum(dN_dξ[:, 3] .* X);
             sum(dN_dξ[:, 1] .* Y) sum(dN_dξ[:, 2] .* Y) sum(dN_dξ[:, 3] .* Y);
             sum(dN_dξ[:, 1] .* Z) sum(dN_dξ[:, 2] .* Z) sum(dN_dξ[:, 3] .* Z)]
        
        # Update ξ
        Δξ = - J \ R
        ξ += Δξ
    end
    
    error("Newton-Raphson method did not converge")
end

# Příklad použití
global_coords = [1.0, 1.5, 1.0]  # Zadejte globální souřadnice bodu
node_coords = [
    0.0 0.0 0.0;
    2.0 0.0 0.0;
    2.0 2.0 0.0;
    0.0 2.0 0.0;
    0.0 0.0 2.0;
    2.0 0.0 2.0;
    2.0 2.0 2.0;
    0.0 2.0 2.0
]

local_coords = find_local_coordinates(global_coords, node_coords)
println("Lokální souřadnice: $local_coords")


using JuMP
using Ipopt

# Definice tvarových funkcí
function shape_functions(ξ₁, ξ₂, ξ₃)
  N = Vector{JuMP.NonlinearExpr}(undef, 8)
  N[1] = -1 / 8 * (ξ₁ - 1) * (ξ₂ - 1) * (ξ₃ - 1)
  N[2] = 1 / 8 * (ξ₁ + 1) * (ξ₂ - 1) * (ξ₃ - 1)
  N[3] = -1 / 8 * (ξ₁ + 1) * (ξ₂ + 1) * (ξ₃ - 1)
  N[4] = 1 / 8 * (ξ₁ - 1) * (ξ₂ + 1) * (ξ₃ - 1)
  N[5] = 1 / 8 * (ξ₁ - 1) * (ξ₂ - 1) * (ξ₃ + 1)
  N[6] = -1 / 8 * (ξ₁ + 1) * (ξ₂ - 1) * (ξ₃ + 1)
  N[7] = 1 / 8 * (ξ₁ + 1) * (ξ₂ + 1) * (ξ₃ + 1)
  N[8] = -1 / 8 * (ξ₁ - 1) * (ξ₂ + 1) * (ξ₃ + 1)
  return N
end

# Funkce pro nalezení lokálních souřadnic
function find_local_coordinates_Iopt(global_coords, node_coords)
    model = Model(Ipopt.Optimizer)
    
    # Definice proměnných ξ₁, ξ₂, ξ₃
    @variable(model, ξ₁, start = 0.0)
    @variable(model, ξ₂, start = 0.0)
    @variable(model, ξ₃, start = 0.0)
    
    # Definice globálních souřadnic uzlů
    X = node_coords[:, 1]
    Y = node_coords[:, 2]
    Z = node_coords[:, 3]
    
    # Tvarové funkce
    N = shape_functions(ξ₁, ξ₂, ξ₃)
    
    # Globální souřadnice jako funkce lokálních souřadnic
    x = sum(N[i] * X[i] for i in 1:8)
    y = sum(N[i] * Y[i] for i in 1:8)
    z = sum(N[i] * Z[i] for i in 1:8)
    
    # Přidání omezení (residua by měla být nulová)
    @NLconstraint(model, x == global_coords[1])
    @NLconstraint(model, y == global_coords[2])
    @NLconstraint(model, z == global_coords[3])
    
    # Optimalizační problém (minimalizujeme nulovou funkci)
    @objective(model, Min, 0)
    
    # Řešení problému
    optimize!(model)
    
    # Výsledné lokální souřadnice
    return value(ξ₁), value(ξ₂), value(ξ₃)
end

# Příklad použití
global_coords = [1.5, 1.1, 0.5]  # Zadejte globální souřadnice bodu
node_coords = [
    0.0 0.0 0.0;
    2.0 0.0 0.0;
    2.0 2.0 0.0;
    0.0 2.0 0.0;
    0.0 0.0 2.0;
    2.0 0.0 2.0;
    2.0 2.0 2.0;
    0.0 2.0 2.0
]

local_coords = @time find_local_coordinates_Iopt(global_coords, node_coords)
println("Lokální souřadnice: $local_coords")
local_coords = @time find_local_coordinates(global_coords, node_coords)