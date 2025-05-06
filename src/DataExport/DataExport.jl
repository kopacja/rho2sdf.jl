module DataExport
   
export SelectProjectedNodes, exportToVTU, exportStructuredPointsToVTK, InputDataToVTU, exportSdfToVTI, export_sdf_results

using Statistics
using LinearAlgebra
using WriteVTK
using JLD2
using Rho2sdf
using Rho2sdf.MeshGrid

# Exporting input data for post-processing:
include("InputDataToVTU.jl")

# Post-processing utilities for filtering and analyzing projected nodes
include("DataPostProcess.jl")

# Exporting mesh data to VTU format (unstructured grid visualization)
include("ExportToVTU.jl")

# Exporting structured point data to legacy VTK format
include("ExportToVTK.jl")

# Exporting SDF grid data to VTI format (image data visualization)
include("ExportToVTI.jl")

# Exporting complete SDF results with metadata
include("ExportSdfResults.jl")

end
