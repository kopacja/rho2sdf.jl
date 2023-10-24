module GenerateMesh


export MeshInformations, elementToNodalValues, Mesh

using Statistics
using LinearAlgebra

include("MeshAnalysis.jl")
# include("PartOfModel.jl")
include("NodalDensities.jl")


end

