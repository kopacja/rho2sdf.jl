module TestOnPrimitiveGeometry

using Test
using Rho2sdf
using Rho2sdf.PrimitiveGeometries
using Rho2sdf.ShapeFunctions
using Rho2sdf.MeshGrid
using Rho2sdf.SignedDistances
using Rho2sdf.MarchingCubes
using Rho2sdf.DataExport
using MAT
using SymPy
using LinearAlgebra
using JLD

     
# taskName = "cube" 
# taskName = "block" 
taskName = "sphere" 
    
# (X, IEN, rho) = PrimitiveGeometries.selectPrimitiveGeometry("cube", 4)
# (X, IEN, rho) = PrimitiveGeometries.selectPrimitiveGeometry("block", [9,1,1])
(X, IEN, rho) = PrimitiveGeometries.selectPrimitiveGeometry("sphere", 9)

mesh = MeshGrid.Mesh(X, IEN, C3D8_SFaD)
ρₙ = MeshGrid.DenseInNodes(mesh, rho) # LSQ

## Grid:
X_min, X_max = MeshGrid.getMesh_AABB(mesh.X) # vec, vec
    
N = 6  #Number of cells along the longest side (along some axis)
sdf_grid = MeshGrid.Grid(X_min, X_max, N, 3) # cartesian grid
    
## SDF from densities:
ρₜ = 0.5
(sdf_dists, xp) = SignedDistances.evalSignedDistances(mesh, sdf_grid, ρₙ, ρₜ)

## Data export to VTK:
# DataExport.exportStructuredPointsToVTK(taskName*"_sdf.vtk", sdf_grid, sdf_dists, "distance")
DataExport.exportStructuredPointsToVTK("primitive-geometry-" *taskName*"_sdf.vtk", sdf_grid, sdf_dists, "distance")

end
