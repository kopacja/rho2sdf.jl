function MeshInformations(data::Any)

  rho = vec(data["rho"])
  mesh = data["msh"]
  X = mesh["X"] # Matrix
  X = [X[:, i] for i in axes(X, 2)]
  IEN = convert(Array{Int64}, mesh["IEN"] .+ 1) # Matrix
  IEN = [IEN[:, i] for i in axes(IEN, 2)] # Vector{Vector{Int64}}

  return (X, IEN, rho)
end

abstract type AbstractMesh end

mutable struct Mesh <: AbstractMesh
  X::Matrix{Float64} # vector of nodes positions
  IEN::Matrix{Int64} # ID element -> ID nodes
  INE::Vector{Vector{Int64}} # ID node -> ID elements
  ISN::Vector{Vector{Int64}} # connectivity face - edges
  sfce::Function # shape functions
  nsd::Int64 # number of spacial dimensions
  nnp::Int64 # number of all nodes
  nen::Int64 # number of element nodes
  nel::Int64 # number of all elements
  nes::Int64 # number of element segments (faces)
  nsn::Int64 # number of face nodes
  edges::NTuple #
  ISE::NTuple #
  rho::Vector{Float64}
  V_domain::Float64
  V_frac::Float64
  ρₜ::Float64
  # INN::Vector{Vector{Int64}}

  function Mesh(X::Vector{Vector{Float64}}, IEN::Vector{Vector{Int64}}, ρₙ::Vector{Float64}, sfce::Function)
    (V_domain, V_frac) = calculate_mesh_volume(X, IEN, ρₙ)

    X = reduce(hcat, X)
    IEN = reduce(hcat, IEN)
    INE = nodeToElementConnectivity(X, IEN)
    ISN = [
      [1, 4, 3, 2], # nodes define face 1
      [1, 2, 6, 5],
      [2, 3, 7, 6],
      [3, 4, 8, 7],
      [4, 1, 5, 8],
      [5, 6, 7, 8],
    ]
    nsd = size(X, 1)
    nnp = size(X, 2)
    nen = size(IEN, 1)
    nel = size(IEN, 2)
    nes = length(ISN)
    nsn = length(ISN[1])
    # INN = [
    #      [4, 2, 5], # Neighbor nodes ID for the corresponding node
    #      [1, 3, 6],
    #      [2, 4, 7],
    #      [3, 1, 8],
    #      [8, 6, 1],
    #      [5, 7, 2],
    #      [6, 8, 3],
    #      [7, 5, 4],
    # ]                               
    edges = ((1, 2), (2, 3), (3, 4), (4, 1),
      (5, 6), (6, 7), (7, 8), (8, 5),
      (1, 5), (2, 6), (3, 7), (4, 8))

    ISE = ((1, 2, 3, 4), # ID edges for each segment
      (1, 9, 5, 10),
      (2, 11, 6, 10),
      (3, 12, 7, 11),
      (4, 12, 8, 9),
      (5, 6, 7, 8))

    # initial_rho_t = @isdefined(ρₜ) ? ρₜ : 0.0
    ρₜ = 0.

    return new(X, IEN, INE, ISN, sfce, nsd, nnp, nen, nel, nes, nsn, edges, ISE, ρₙ, V_domain, V_frac, ρₜ) #INN)
  end
end

function nodeToElementConnectivity(
  X::Matrix,
  IEN::Matrix,
)
  INE = [Vector{Int64}() for _ = 1:size(X, 2)]
  for el = 1:size(IEN, 2)
    for i = 1:size(IEN, 1)
      push!(INE[IEN[i, el]], el)
    end
  end
  return INE
end

mutable struct TriangularMesh <: AbstractMesh
  X::Matrix{Float64} # vector of nodes positions
  IEN::Matrix{Int64} # ID element -> ID nodes
  INE::Vector{Vector{Int64}} # ID node -> ID elements
  nsd::Int64 # number of spacial dimensions (3)
  nnp::Int64 # number of all nodes
  nen::Int64 # number of element nodes (3)
  nel::Int64 # number of all elements

  function TriangularMesh(X::Vector{Vector{Float64}}, IEN::Vector{Vector{Int64}})
    X = reduce(hcat, X)
    IEN = reduce(hcat, IEN)
    INE = nodeToElementConnectivity(X, IEN)
    nsd = size(X, 1)
    nnp = size(X, 2)
    nen = size(IEN, 1)
    nel = size(IEN, 2)

    return new(X, IEN, INE, nsd, nnp, nen, nel)
  end
end
