function MeshInformations(data::Any)

    rho = vec(data["rho"])
    mesh = data["msh"]
    X = mesh["X"] # Matrix
    X = [X[:, i] for i in axes(X,2)]
    IEN = convert(Array{Int64}, mesh["IEN"] .+ 1) # Matrix
    IEN = [IEN[:, i] for i in axes(IEN, 2)] # Vector{Vector{Int64}}

    return (X, IEN, rho)
end


mutable struct Mesh
    X::Matrix{Float64} # vector of nodes positions
    IEN::Matrix{Int64} # ID element -> ID nodes
    INE::Vector{Vector{Int64}} # ID node -> ID elements
    ISN::Vector{Vector{Int64}} # connectivity face - edges
    nsd::Int64 # number of spacial dimensions
    nnp::Int64 # number of all nodes
    nen::Int64 # number of element nodes
    nel::Int64 # number of all elements
    nes::Int64 # number of element segments (faces)
    nsn::Int64 # number of face nodes

    function Mesh(X::Vector{Vector{Float64}}, IEN::Vector{Vector{Int64}})
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
        return new(X, IEN, INE, ISN, nsd, nnp, nen, nel, nes, nsn)
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


