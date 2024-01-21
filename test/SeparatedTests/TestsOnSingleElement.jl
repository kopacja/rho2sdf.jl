module TestsOnSingleElement


using Test
using Rho2sdf
using Rho2sdf.PrimitiveGeometries
using Rho2sdf.ShapeFunctions
using Rho2sdf.MeshGrid
using Rho2sdf.SignedDistances
using MAT
using SymPy
using LinearAlgebra
using JLD

    
taskName = "block"
(X, IEN, rho) = PrimitiveGeometries.selectPrimitiveGeometry("block", [1,1,1])

# generate mesh struct:
mesh = MeshGrid.Mesh(X, IEN)

# rho in nodes:
# ρₙ = MeshGrid.DenseInNodes(mesh, rho) # LSQ
ρₙ = [0., 0., 0., 0., 1., 1., 1., 1.]

## Grid:
X_min, X_max = MeshGrid.getMesh_AABB(mesh.X) # vec, vec
    
N = 2  #Number of cells along the longest side (along some axis)
sdf_grid = MeshGrid.Grid(X_min, X_max, N) # cartesian grid
    

## SDF from densities:
ρₜ = 0.5
sdf_dists = SignedDistances.evalSignedDistances(mesh, sdf_grid, ρₙ , ρₜ)

## Data export to VTK:
# Rho2sdf.DataProcessing.exportStructuredPointsToVTK(taskName*"_sdf.vtk", sdf_grid, sdf_dists, "distance")
Rho2sdf.exportStructuredPointsToVTK("sdf2-test-" *taskName*"_sdf.vtk", sdf_grid, sdf_dists, "distance")



end
