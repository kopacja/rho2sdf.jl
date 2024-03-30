module MeshGrid


export MeshInformations, elementToNodalValues, Mesh, PartOfModel, Grid, LinkedList, getMesh_AABB, generateGridPoints, extractSurfaceTriangularMesh, calculateMiniAABB_grid, ModiffElementalDensities, generateConnectivityArray

using Statistics
using LinearAlgebra

include("MeshInformations.jl")
include("PartOfModel.jl")
include("NodalDensities.jl")
include("SurfaceTriangularMesh.jl")
include("Grid.jl")


end

