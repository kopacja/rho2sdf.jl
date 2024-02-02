module TestOnLegGripper

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
    
    # Data from Matlab:
    taskName = "chapadlo"
    # data = matread(taskName * ".mat")
    data = matread("test/" * taskName * ".mat")
    
    part_name = "test/elementy_trubky.txt"

    (X, IEN, rho) = MeshGrid.MeshInformations(data)

    # input data propertis (mesh, density)
    mesh = MeshGrid.Mesh(X, IEN, C3D8_SFaD)
    (mesh, rho) = MeshGrid.PartOfModel(mesh, rho, part_name)
    # rho = MeshGrid.ModiffElementalDensities(mesh, rho) # change density along the object

    ρₙ = MeshGrid.DenseInNodes(mesh, rho) # LSQ
    ## Grid:
    X_min, X_max = MeshGrid.getMesh_AABB(mesh.X) # vec, vec
    
    N = 50  #Number of divisions along the longest side (along some axis)
    sdf_grid = MeshGrid.Grid(X_min, X_max, N, 5) # cartesian grid
    
    ## SDF from densities:
    ρₜ = 0.5
    (sdf_dists, xp) = SignedDistances.evalSignedDistances(mesh, sdf_grid, ρₙ, ρₜ)

    ## Data export to VTK:
    Rho2sdf.exportStructuredPointsToVTK("sdf2-leg-" * taskName *"_sdf.vtk", sdf_grid, sdf_dists, "distance")


end
