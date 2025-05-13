module SignedDistances

export evalDistances, calculate_triangle_edges, SelectProjectedNodes, Sign_Detection

using Base.Threads
using Einsum
using Statistics
using StaticArrays
using LinearAlgebra
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

include("PseudoNormals.jl")
include("TriangularMeshUtils.jl")
include("sdfOnDensityField.jl")
include("SignDetection.jl")


end
