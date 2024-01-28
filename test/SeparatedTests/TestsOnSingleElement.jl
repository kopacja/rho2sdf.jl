module TestsOnSingleElement


using Test
using Rho2sdf
using Rho2sdf.PrimitiveGeometries
using Rho2sdf.ShapeFunctions
using Rho2sdf.MeshGrid
using Rho2sdf.SignedDistances
using Rho2sdf.DataExport
using MAT
using SymPy
using LinearAlgebra
using JLD
using DelimitedFiles

taskName = "block"
(X, IEN, rho) = PrimitiveGeometries.selectPrimitiveGeometry("block", [1,1,1])

# generate mesh struct:
mesh = MeshGrid.Mesh(X, IEN)

# rho in nodes:
# ρₙ = MeshGrid.DenseInNodes(mesh, rho) # LSQ
# ρₙ = [0.3, 0.8, 0., 0., 1., 1., 1., 1.]
ρₙ = [0., 0., 0., 0., 1., 1., 1., 1.]

## Grid:
X_min, X_max = MeshGrid.getMesh_AABB(mesh.X) # vec, vec
    
N = 2  #Number of cells along the longest side (along some axis)
sdf_grid = MeshGrid.Grid(X_min, X_max, N) # cartesian grid

sdf_grid.AABB_min = [0.0, 0.0, 0.0]
sdf_grid.AABB_max = [1.0, 1.0, 1.0]
sdf_grid.N = [1.0, 1.0, 1.0]
sdf_grid.cell_size = 0.5
sdf_grid.ngp = 1.

## SDF from densities:
ρₜ = 0.5
sdf_dists, xp = SignedDistances.evalSignedDistances(mesh, sdf_grid, ρₙ , ρₜ)

points = MeshGrid.generateGridPoints(sdf_grid) # uzly pravidelné mřížky
Xg, Xp, mean_PD, max_PD = DataExport.SelectProjectedNodes(mesh, sdf_grid, xp, points)

nnp = Int(length(Xg)/2)
IEN = [[i; i + nnp] for i = 1:nnp]

nnp₂ = Int(length(Xg)/2)
IEN₂ = [[i; i + nnp₂] for i = 1:nnp₂]

X_combined = [Xg; Xp] 
# X_combined_couples = [X Xp]

DataExport.exportToVTU("xp.vtu", X_combined, IEN)
DataExport.exportToVTU("Xg.vtu", Xg, IEN₂)
DataExport.exportToVTU("Xp.vtu", Xp, IEN₂)

open("xp.csv", "w") do io
    writedlm(io, ['x' 'y' 'z'], ',')
    writedlm(io, xp', ',')
end





## Data export to VTK:
DataExport.exportStructuredPointsToVTK(taskName*"_sdf.vtk", sdf_grid, sdf_dists, "distance")
# Rho2sdf.exportStructuredPointsToVTK("sdf2-test-" *taskName*"_sdf.vtk", sdf_grid, sdf_dists, "distance")



end
