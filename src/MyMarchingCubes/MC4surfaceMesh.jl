# function convert_to_matrix(Xₑ::Vector)
#   matrix = hcat(Xₑ...)
#   return matrix
# end
# Xₑ = convert_to_matrix(Xₑᵥ)

function MC_OnCube(ρₑ::Vector{Float64},
  ρₜ::Float64)
  # println("rho ele: ", ρₑ)

  # Sort to match my elements definitions
  ρ = [ρₑ[1], ρₑ[2], ρₑ[4], ρₑ[3], ρₑ[5], ρₑ[6], ρₑ[8], ρₑ[7]] .- ρₜ
  Y = reshape(ρ, (2, 2, 2)) # Array for MC algorithm

  mc = MC(Y, Int; x=[-1.0, 1.0], y=[-1.0, 1.0], z=[-1.0, 1.0]) # cube inicialization (MC not working with general mesh, but with grid)

  # Perform MC:
  march(mc)
  # println("triangles: ", mc.triangles)
  # println("vertices: ", mc.vertices)
  # exit()
  return (mc.triangles, mc.vertices)
end

# (triangles, vertices) = MC_OnCube(ρₑ, ρₜ)

function ShapeFunctionTemp(Ξ)
  ξ₁ = Ξ[1]
  ξ₂ = Ξ[2]
  ξ₃ = Ξ[3]

  # N = Array{Float64}(undef, 8)
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


# Function to assign coordinates to each node
function assign_coordinates(
  triangles,
  coordinates::Vector{Vector{Float64}}
)

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

# result = assign_coordinates(triangles, vertices_coords)


function MC4trianglesInElement(
  ρₑ::Vector,
  ρₜ::Float64,
  Xₑ::Matrix,
)

  # Cutting triangles single cube element based on nodal densities (MS algorithm):
  (triangles, vertices) = MC_OnCube(ρₑ, ρₜ) # POTŘEBUJI

  # Transform vertices coords from cube to real element:
  VertCoord = Vector{Vector{Float64}}(undef, length(vertices))  # Preallocate VertCoord

  for i in 1:length(vertices)
    VertCoord[i] = Xₑ * ShapeFunctionTemp(vertices[i])  # Directly assign the result to each index
  end

  Coords4Triangles = assign_coordinates(triangles, VertCoord)

  return (Coords4Triangles, triangles)
end

function MC_SurfaceTriangularMesh(
  ρₙ::Vector{Float64},
  ρₜ::Float64,
  mesh::MGMesh,
)

  X = mesh.X
  IEN = mesh.IEN # ID element -> nodes
  nel = mesh.nel # number of elements

  X_new = Vector{Vector{Float64}}()
  # push!(X_new, vec(X[:, 1]))

  IEN_new = Vector{Vector{Int64}}()

  for el = 1:nel
    ρₑ = ρₙ[IEN[:, el]]

    ρₑ_min = minimum(ρₑ)
    ρₑ_max = maximum(ρₑ)
    if (ρₑ_max > ρₜ) && (ρₑ_min < ρₜ) # The boundary (isocontour) goes through the element

      Xₑ = X[:, IEN[:, el]]
      (Coords4Triangles, triangles) = MC4trianglesInElement(ρₑ, ρₜ, Xₑ)

      not = length(Coords4Triangles) # number of trinagles in the element

      for i in 1:not
        Xt = Coords4Triangles[i]
        IEN_el = zeros(Int64, 3)

        for j in 1:3
          a = findfirst(x -> norm(x - Xt[j]) < 1.0e-5, X_new)
          if (a === nothing)
            push!(X_new, Xt[j])
            IEN_el[j] = length(X_new)
          else
            IEN_el[j] = a[1]
          end
        end
        push!(IEN_new, IEN_el)  # This should be outside the inner loop
      end
    end
  end # el

  return (X_new, IEN_new)
end
