module Rho2sdf

export MeshInformations

using LinearAlgebra
using Statistics
using DelimitedFiles
using Einsum
using BenchmarkTools
using JLD2

# Tools for monitoring the computation process:
include("TerminalUtils/TerminalUtils.jl")
using .TerminalUtils

# Element type system
include("ElementTypes/ElementTypes.jl")
using .ElementTypes
export AbstractElement, HEX8, TET4

# Data import capabilities
include("DataImport/DataImport.jl")
using .DataImport
export import_vtu_mesh, validate_vtu_mesh

# Updated shape functions with type dispatch
include("ShapeFunctions/ShapeFunctions.jl")
using .ShapeFunctions
export shape_functions, compute_shape_and_derivatives!

include("MeshGrid/MeshGrid.jl")
using .MeshGrid

include("SignedDistances/SignedDistances.jl")
using .SignedDistances

include("PrimitiveGeometries/PrimitiveGeometries.jl")
using .PrimitiveGeometries

include("DataExport/DataExport.jl")
using .DataExport

include("SdfSmoothing/SdfSmoothing.jl")
using .SdfSmoothing

include("Visualizations/Visualizations.jl")
using .Visualizations

include("RhoToSDF.jl")
export Rho2sdfOptions, rho2sdf, rho2sdf_hex8, rho2sdf_tet4

end # module Rho2sdf
