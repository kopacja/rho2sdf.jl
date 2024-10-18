using Rho2sdf
using Rho2sdf.PrimitiveGeometries
using Rho2sdf.ShapeFunctions
using Rho2sdf.MeshGrid
using Rho2sdf.SignedDistances
using Rho2sdf.DataExport
using Rho2sdf.SdfSmoothing
using MAT
using JLD2
using LinearAlgebra

## Inputs:
taskName = "chapadlo"
N = 60  # Number of cells along the longest side
ρₜ = 0.5 # Threshold density (isosurface level)

## Read FEM mesh:
data = matread("test/chapadlo.mat")
(X, IEN, rho) = MeshGrid.MeshInformations(data)

## Generate FEM mesh structure:
mesh = MeshGrid.Mesh(X, IEN, C3D8_SFaD)

## Map elemental densities to the nodes:
ρₙ = MeshGrid.DenseInNodes(mesh, rho) # LSQ

VTK_CODE = 12 # https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html
Rho2sdf.exportToVTU(taskName * "part_nodal_densities.vtu", X, IEN, VTK_CODE, ρₙ)

## Grid:
X_min, X_max = MeshGrid.getMesh_AABB(mesh.X)
X_min[2] = 0
X_min[3] = 50
sdf_grid = MeshGrid.Grid(X_min, X_max, N, 3) # cartesian grid
points = MeshGrid.generateGridPoints(sdf_grid) # uzly pravidelné mřížky

## SDF from densities:
(dists, xp) = SignedDistances.evalDistances(mesh, sdf_grid, points, ρₙ, ρₜ)
signs = SignedDistances.Sign_Detection(mesh, sdf_grid, points, ρₙ, ρₜ)
sdf_dists = dists .* signs

## Export to VTK:
Rho2sdf.exportStructuredPointsToVTK(taskName * "_polovina_sdf_new.vtk", sdf_grid, sdf_dists, "distance")

@save "$(taskName)_polovina_SDF.jld2" sdf_dists
@save "$(taskName)_polovina_Grid.jld2" sdf_grid
@save "$(taskName)_polovina_Points.jld2" points
@save "$(taskName)_polovina_Mesh.jld2" mesh
