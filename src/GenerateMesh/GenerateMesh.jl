module GenerateMesh


export MeshInformations, elementToNodalValues, Mesh, PartOfModel

using Statistics
using LinearAlgebra

include("MeshAnalysis.jl")
include("PartOfModel.jl")
include("NodalDensities.jl")


end

