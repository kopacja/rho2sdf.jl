using Rho2sdf
using Rho2sdf.PrimitiveGeometries
using Rho2sdf.ShapeFunctions
using Rho2sdf.MeshGrid
using Rho2sdf.SignedDistances
using Rho2sdf.DataExport
using Rho2sdf.SdfSmoothing
using Rho2sdf.ImplicitDomainMeshing
using MAT
using JLD2
using LinearAlgebra

## Inputs:
# taskName = "chapadlo"
taskName = "cantilever_beam_vfrac_04"

## Read FEM mesh:
data = matread("test/" * taskName * ".mat")
(X, IEN, rho) = MeshInformations(data)
IEN = [subvector .- 1 for subvector in IEN]  # This will work

## Generate FEM mesh structure:
mesh = MeshGrid.Mesh(X, IEN, rho, hex8_shape)

InputDataToVTU(mesh, (taskName * "_Raw"))

## Grid:
sdf_grid = MeshGrid.interactive_sdf_grid_setup(mesh)
# sdf_grid = MeshGrid.noninteractive_sdf_grid_setup(mesh, 2.3)
points = MeshGrid.generateGridPoints(sdf_grid) # uzly pravidelné mřížky

## Map elemental densities to the nodes:
ρₙ = MeshGrid.DenseInNodes(mesh, rho) # LSQ
# ρₜ = find_threshold_for_volume(mesh, ρₙ)
ρₜ = 0.5

VTK_CODE = 12 # https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html
Rho2sdf.exportToVTU(taskName * "_nodal_densities.vtu", X, IEN, VTK_CODE, ρₙ)

## SDF from densities:
(dists, xp) = SignedDistances.evalDistances(mesh, sdf_grid, points, ρₙ, ρₜ)
signs = SignedDistances.Sign_Detection(mesh, sdf_grid, points, ρₙ, ρₜ)
sdf_dists = dists .* signs


## Export to VTK:
B = round(sdf_grid.cell_size, digits=4)
Rho2sdf.exportStructuredPointsToVTK(taskName * "_SDF_B-" * string(B) * ".vtk", sdf_grid, sdf_dists, "distance")

# RBF smoothing:
(fine_sdf, fine_grid) = RBFs_smoothing(sdf_dists, sdf_grid, false, 1, taskName) # interpolation == true, aproximation == false, smooth

# Definice rovin
plane_definitions = [
    PlaneDefinition([-1.0, 0.0, 0.0], [0.0, 10.0, 0.0], Square(30.)),
    PlaneDefinition([1.0, 0.0, 0.0], [60.0, 2.0, 2.0], Square(5.))
]

tetMesh = GenerateTetMesh(fine_sdf, fine_grid, "A15", taskName, plane_definitions)

slice_mesh_with_plane!(tetMesh, "x", 0.6, export_file="sliced_mesh.vtu")

# @save "Z_$(taskName)_cele_SDF_B-$(B).jld2" sdf_dists
# @save "Z_$(taskName)_cele_Grid_B-$(B).jld2" sdf_grid
# @save "Z_$(taskName)_cele_Points_B-$(B).jld2" points
# @save "Z_$(taskName)_cele_Mesh_B-$(B).jld2" mesh

