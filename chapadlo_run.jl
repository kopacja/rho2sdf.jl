using Rho2sdf
using Rho2sdf.PrimitiveGeometries
using Rho2sdf.ShapeFunctions
using Rho2sdf.MeshGrid
using Rho2sdf.SignedDistances
using Rho2sdf.DataExport
using Rho2sdf.SdfSmoothing
using Rho2sdf.Visualizations
using MAT
using JLD2
using LinearAlgebra

## Inputs:
taskName = "chapadlo"

## Read FEM mesh:
data = matread("test/" * taskName * ".mat")
(X, IEN, rho) = MeshInformations(data)

## Generate FEM mesh structure:
mesh = MeshGrid.Mesh(X, IEN, rho, hex8_shape)

## Grid:
sdf_grid = MeshGrid.interactive_sdf_grid_setup(mesh)
# sdf_grid = MeshGrid.noninteractive_sdf_grid_setup(mesh)
points = MeshGrid.generateGridPoints(sdf_grid) # uzly pravidelné mřížky

## Map elemental densities to the nodes:
ρₙ = MeshGrid.DenseInNodes(mesh, rho) # LSQ
ρₜ = find_threshold_for_volume(mesh, ρₙ)

VTK_CODE = 12 # https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html
Rho2sdf.exportToVTU(taskName * "_nodal_densities.vtu", X, IEN, VTK_CODE, ρₙ)

## SDF from densities:
(dists, xp) = SignedDistances.evalDistances(mesh, sdf_grid, points, ρₙ, ρₜ)
signs = SignedDistances.Sign_Detection(mesh, sdf_grid, points, ρₙ, ρₜ)
sdf_dists = dists .* signs

## Export to VTK:
B = round(sdf_grid.cell_size, digits=4)
exportSdfToVTI(taskName * "_SDF_B-" * string(B) * ".vti", sdf_grid, sdf_dists, "distance")

# RBF smoothing:
is_interp = false
smooth = 1
(fine_sdf, fine_grid) = RBFs_smoothing(mesh, sdf_dists, sdf_grid, is_interp, smooth, taskName) # interpolation == true, aproximation == false, smooth
export_sdf_results(fine_sdf, fine_grid, sdf_grid, taskName, smooth, is_interp)

# @save "Z_$(taskName)_cele_SDF_B-$(B).jld2" sdf_dists
# @save "Z_$(taskName)_cele_Grid_B-$(B).jld2" sdf_grid
# @save "Z_$(taskName)_cele_Points_B-$(B).jld2" points
# @save "Z_$(taskName)_cele_Mesh_B-$(B).jld2" mesh
