using Rho2sdf
using Rho2sdf.PrimitiveGeometries
using Rho2sdf.ShapeFunctions
using Rho2sdf.MeshGrid
using Rho2sdf.SignedDistances
using Rho2sdf.DataExport
using MAT
using JLD2
using LinearAlgebra

## Inputs:
taskName = "chapadlo"
ρₜ = 0.5 # Threshold density (isosurface level)

## Read FEM mesh:
data = matread("test/" * taskName * ".mat")
(X, IEN, rho) = MeshInformations(data)

## Generate FEM mesh structure:
mesh = MeshGrid.Mesh(X, IEN, C3D8_SFaD)

## Grid:
# sdf_grid = MeshGrid.interactive_sdf_grid_setup(mesh)
sdf_grid = MeshGrid.noninteractive_sdf_grid_setup(mesh, 2.3)
points = MeshGrid.generateGridPoints(sdf_grid) # uzly pravidelné mřížky

## Map elemental densities to the nodes:
ρₙ = MeshGrid.DenseInNodes(mesh, rho) # LSQ

VTK_CODE = 12 # https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html
Rho2sdf.exportToVTU(taskName * "_nodal_densities.vtu", X, IEN, VTK_CODE, ρₙ)

## SDF from densities:
(dists, xp) = SignedDistances.evalDistances(mesh, sdf_grid, points, ρₙ, ρₜ)
signs = SignedDistances.Sign_Detection(mesh, sdf_grid, points, ρₙ, ρₜ)
sdf_dists = dists .* signs


## Export to VTK:
B = sdf_grid.cell_size
Rho2sdf.exportStructuredPointsToVTK(taskName * "_sdf_CellSize-" * (round(B, digits=4)) * ".vtk", sdf_grid, sdf_dists, "distance")

# RBF smoothing:
RBFs_smoothing(sdf_dists, sdf_grid, false, 2, taskName) # interpolation == true, aproximation == false, smooth

@save "Z_$(taskName)_cele_SDF_CellSize$(round(B, digits=4)).jld2" sdf_dists
@save "Z_$(taskName)_cele_Grid_CellSize$(round(B, digits=4)).jld2" sdf_grid
@save "Z_$(taskName)_cele_Points_CellSize$(round(B, digits=4)).jld2" points
@save "Z_$(taskName)_cele_Mesh_CellSize$(round(B, digits=4)).jld2" mesh

