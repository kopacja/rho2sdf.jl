module SignedDistances

export evalSignedDistancesOnTriangularMesh, evalDistances, computePseudoNormals, barycentricCoordinates, calculate_triangle_edges, update_distance!, SelectProjectedNodes, find_local_coordinates, Sign_Detection, compute_coords

using Base.Threads
using Einsum
using Statistics
using StaticArrays
using LinearAlgebra
using DelimitedFiles
using ProgressMeter
using NLopt
using BenchmarkTools

using Rho2sdf.ShapeFunctions
using Rho2sdf.TerminalUtils
using Rho2sdf.MeshGrid
using Rho2sdf

# Compute local coords:
include("FindLocalCoordinates.jl")
include("ComputeCoordsOnIso.jl")

include("Derivatives.jl")
include("PseudoNormals.jl")
include("sdfOnTriangularMesh.jl")
include("sdfOnDensityField.jl")
include("SignDetection.jl")


end
