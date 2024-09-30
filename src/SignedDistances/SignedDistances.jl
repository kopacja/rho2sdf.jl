module SignedDistances

export evalSignedDistancesOnTriangularMesh, evalSignedDistances, computePseudoNormals, barycentricCoordinates, calculate_triangle_edges, update_distance!, SelectProjectedNodes, find_local_coordinates

using Base.Threads
using Einsum
using Statistics
using LinearAlgebra
using DelimitedFiles
using ProgressMeter
using LazySets
using JLD2
using NLopt

using Rho2sdf.ShapeFunctions
using Rho2sdf.MeshGrid
using Rho2sdf.MyMarchingCubes
using Rho2sdf

include("Derivatives.jl")
include("PseudoNormals.jl")
include("sdfOnTriangularMesh.jl")
# include("sdfOnDensityField_clean.jl")
include("sdfOnDensityField_parallel_NLopt.jl")
# include("sdfOnDensityField_new.jl")
# include("sdfOnDensityField.jl")



end
