using Rho2sdf.ElementTypes

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

mutable struct Mesh{T<:AbstractElement} <: AbstractMesh
  element_type::Type{T}           # Element type for dispatch
  X::Matrix{Float64}              # vector of nodes positions
  IEN::Matrix{Int64}              # ID element -> ID nodes
  INE::Vector{Vector{Int64}}      # ID node -> ID elements
  ISN::Vector{Vector{Int64}}      # connectivity face - edges
  sfce::Function                  # shape functions
  nsd::Int64                      # number of spacial dimensions
  nnp::Int64                      # number of all nodes
  nen::Int64                      # number of element nodes
  nel::Int64                      # number of all elements
  nes::Int64                      # number of element segments (faces)
  nsn::Int64                      # number of face nodes
  edges::NTuple                   # edge definitions
  ISE::NTuple                     # edge-to-segment mapping
  rho::Vector{Float64}            # element densities
  V_domain::Float64               # total domain volume
  V_frac::Float64                 # volume fraction
  ρₜ::Float64                     # threshold density

  function Mesh(X::Vector{Vector{Float64}}, IEN::Vector{Vector{Int64}}, ρₙ::Vector{Float64}, 
                sfce::Function; element_type::Type{T}=HEX8) where {T<:AbstractElement}
    
    # Calculate mesh volume
    (V_domain, V_frac) = calculate_mesh_volume(X, IEN, ρₙ, element_type)

    # Convert to matrices
    X_mat = reduce(hcat, X)
    IEN_mat = reduce(hcat, IEN)
    
    # Create inverse connectivity
    INE = nodeToElementConnectivity(X_mat, IEN_mat)
    
    # Get element topology information
    nen, nes, nsn, ISN, edges, ISE = get_element_topology(element_type)
    
    # Basic mesh properties
    nsd = size(X_mat, 1)
    nnp = size(X_mat, 2)
    nel = size(IEN_mat, 2)
    
    # Validate element consistency
    if size(IEN_mat, 1) != nen
        error("Element connectivity size ($(size(IEN_mat, 1))) doesn't match element type nodes ($nen)")
    end

    ρₜ = 0.0  # Initialize threshold density

    return new{T}(element_type, X_mat, IEN_mat, INE, ISN, sfce, nsd, nnp, nen, nel, 
                  nes, nsn, edges, ISE, ρₙ, V_domain, V_frac, ρₜ)
  end
end

function nodeToElementConnectivity(X::Matrix, IEN::Matrix)
  INE = [Vector{Int64}() for _ = 1:size(X, 2)]
  for el = 1:size(IEN, 2)
    for i = 1:size(IEN, 1)
      push!(INE[IEN[i, el]], el)
    end
  end
  return INE
end

mutable struct TriangularMesh <: AbstractMesh
  X::Matrix{Float64}              # vector of nodes positions
  IEN::Matrix{Int64}              # ID element -> ID nodes
  INE::Vector{Vector{Int64}}      # ID node -> ID elements
  nsd::Int64                      # number of spacial dimensions (3)
  nnp::Int64                      # number of all nodes
  nen::Int64                      # number of element nodes (3)
  nel::Int64                      # number of all elements

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

# Convenience constructors for specific element types
function HexMesh(X::Vector{Vector{Float64}}, IEN::Vector{Vector{Int64}}, ρₙ::Vector{Float64}, sfce::Function)
    return Mesh(X, IEN, ρₙ, sfce; element_type=HEX8)
end

function TetMesh(X::Vector{Vector{Float64}}, IEN::Vector{Vector{Int64}}, ρₙ::Vector{Float64}, sfce::Function)
    return Mesh(X, IEN, ρₙ, sfce; element_type=TET4)
end

function get_vector_format(mesh::Mesh)
    X_vec = [mesh.X[:, i] for i in 1:size(mesh.X, 2)]
    IEN_vec = [mesh.IEN[:, i] for i in 1:size(mesh.IEN, 2)]
    return X_vec, IEN_vec
end
