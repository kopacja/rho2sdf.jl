using Base: print_without_params
using Rho2sdf
using Rho2sdf.TerminalUtils
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

taskName = "block"

N = 10  # Number of cells along the longest side
# ρₜ = 0.5 # Threshold density (isosurface level)

(X, IEN, rho) = PrimitiveGeometries.selectPrimitiveGeometry("block", [2, 1, 1])
# ρₙ = [0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.5, 0.5]

mesh = MeshGrid.Mesh(X, IEN, rho, hex8_shape)
# ρₙ = MeshGrid.DenseInNodes(mesh, rho) # LSQ

# Modif ρₙ:
ρₙ = [0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 0.0, 0.0, 0.5, 0.5]
ρₜ = find_threshold_for_volume(mesh, ρₙ)

VTK_CODE = 12 # https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html
Rho2sdf.exportToVTU(taskName * "_nodal_densities.vtu", X, IEN, VTK_CODE, ρₙ)

## Grid:
X_min, X_max = MeshGrid.getMesh_AABB(mesh.X)
sdf_grid = MeshGrid.Grid(X_min, X_max, N, 3) # cartesian grid
points = MeshGrid.generateGridPoints(sdf_grid) # uzly pravidelné mřížky

## SDF from densities:
(dists, xp) = SignedDistances.evalDistances(mesh, sdf_grid, points, ρₙ, ρₜ)
signs = SignedDistances.Sign_Detection(mesh, sdf_grid, points, ρₙ, ρₜ)
sdf_dists = dists .* signs

## Export to VTK:
Rho2sdf.exportStructuredPointsToVTK(taskName * "_SDF.vtk", sdf_grid, sdf_dists, "distance")

# RBF smoothing:
fine_LSF = RBFs_smoothing(sdf_dists, sdf_grid, false, 2, taskName) # interpolation == true, aproximation == false, smooth

fig = visualize_stable_isosurface(fine_LSF)
display(fig)

