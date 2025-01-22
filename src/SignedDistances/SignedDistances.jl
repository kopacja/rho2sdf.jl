module SignedDistances

export evalSignedDistancesOnTriangularMesh, evalDistances, computePseudoNormals, barycentricCoordinates, calculate_triangle_edges, update_distance!, SelectProjectedNodes, find_local_coordinates, Sign_Detection, compute_coords

using Base.Threads
using Einsum
using Statistics
using LinearAlgebra
using DelimitedFiles
using ProgressMeter
using LazySets
using JLD2
using NLopt
using JuMP
import Ipopt

using Rho2sdf.ShapeFunctions
using Rho2sdf.MeshGrid
using Rho2sdf

include("Derivatives.jl")
include("PseudoNormals.jl")
include("sdfOnTriangularMesh.jl")
include("sdfOnDensityField.jl")
include("SignDetection.jl")



end
