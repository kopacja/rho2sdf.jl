# function convert_to_matrix(Xₑ::Vector)
#   matrix = hcat(Xₑ...)
#   return matrix
# end
# Xₑ = convert_to_matrix(Xₑᵥ)

function MC_OnCube(ρₑ::Vector{Float64},
  ρₜ::Float64)

  # Sort to match my elements definitions
  ρ = [ρₑ[1], ρₑ[2], ρₑ[4], ρₑ[3], ρₑ[5], ρₑ[6], ρₑ[8], ρₑ[7]] .- ρₜ
  Y = reshape(ρ, (2, 2, 2)) # Array for MC algorithm

  mc = MC(Y, Int; x=[-1.0, 1.0], y=[-1.0, 1.0], z=[-1.0, 1.0]) # cube inicialization (MC not working with general mesh, but with grid)

  # Perform MC:
  march(mc)
  return (mc.triangles, mc.vertices)
end

# (triangles, vertices) = MC_OnCube(ρₑ, ρₜ)


function EdgeIdDetection(
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

# (vertices_coords, ID_edges) = EdgeIdDetection1(vertices, Xₑ, edges, nsd)

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
  mesh::MGMesh,
  Xₑ::Matrix,
)

  # ISE = mesh.ISE
  nsd = mesh.nsd
  edges = mesh.edges

  # Cutting triangles single cube element based on nodal densities (MS algorithm):
  (triangles, vertices) = MC_OnCube(ρₑ, ρₜ) # POTŘEBUJI

  # Transform vertices coords from cube to real element:
  (vertices_coords, ID_edges) = EdgeIdDetection(vertices, Xₑ, edges, nsd) # POTŘEBUJI vertices_coords

  Coords4Triangles = assign_coordinates(triangles, vertices_coords)

  return (Coords4Triangles, triangles)
end

function MC_SurfaceTriangularMesh(
  ρₙ::Vector{Float64},
  ρₜ::Float64,
  mesh::MGMesh,
)

  X = mesh.X
  IEN = mesh.IEN # ID element -> nodes
  INE = mesh.INE # ID node -> ID elements
  ISN = mesh.ISN # connectivity face - edges
  nsd = mesh.nsd # number of special dimension
  nel = mesh.nel # number of elements
  nes = mesh.nes # number of element segments (faces) 6
  nsn = mesh.nsn # number of segment nodes (kolik má stěna uzlů) 4
  sfce = mesh.sfce

  X_new = Vector{Vector{Float64}}()
  # push!(X_new, vec(X[:, 1]))

  IEN_new = Vector{Vector{Int64}}()

  for el = 1:nel
    ρₑ = ρₙ[IEN[:, el]]

    ρₑ_min = minimum(ρₑ)
    ρₑ_max = maximum(ρₑ)
    if (ρₑ_max > ρₜ) && (ρₑ_min < ρₜ) # The boundary (isocontour) goes through the element

      Xₑ = X[:, IEN[:, el]]
      (Coords4Triangles, triangles) = MC4trianglesInElement(ρₑ, ρₜ, mesh, Xₑ)

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
