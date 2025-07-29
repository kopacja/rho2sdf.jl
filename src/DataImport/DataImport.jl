# src/DataImport/DataImport.jl

module DataImport

export import_vtu_mesh, validate_vtu_mesh

using ReadVTK
using Rho2sdf.TerminalUtils

include("VTUImport.jl")

end
