module SignedDistances

export evalSignedDistancesOnTriangularMesh, evalSignedDistances, computePseudoNormals, barycentricCoordinates, calculate_triangle_edges, update_distance!, SelectProjectedNodes

using Base.Threads
using Einsum
using Statistics
using LinearAlgebra
using DelimitedFiles
using JuMP
import Ipopt

using Rho2sdf.ShapeFunctions
using Rho2sdf.MeshGrid
# using Rho2sdf.MyMarchingCubes
using Rho2sdf

include("Derivatives.jl")
include("PseudoNormals.jl")
include("sdfOnTriangularMesh.jl")
include("sdfOnDensityField_new.jl")
# include("sdfOnDensityField.jl")



end
