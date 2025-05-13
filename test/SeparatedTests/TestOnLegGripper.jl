module TestOnLegGripper

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
    
    # Data from Matlab:
    taskName = "chapadlo"
    data = matread(taskName * ".mat")
    # data = matread("test/" * taskName * ".mat")
    

    # part_name = "2elementy.txt"
    part_name = "elementy_trubky.txt"
    # part_name = "test/elementy_trubky.txt"

    (X, IEN, rho) = MeshGrid.MeshInformations(data)

    # input data propertis (mesh, density)
    mesh = MeshGrid.Mesh(X, IEN, hex8_shape)
    rho = MeshGrid.ModiffElementalDensities(mesh, rho) # change density along the object

    ρₙ = MeshGrid.DenseInNodes(mesh, rho) # LSQ
    ## Grid:
    X_min, X_max = MeshGrid.getMesh_AABB(mesh.X) # vec, vec
    

    N = 18  #Number of divisions along the longest side (along some axis)
    sdf_grid = MeshGrid.Grid(X_min, X_max, N, 3) # cartesian grid

    
    ## SDF from densities:
    ρₜ = 0.5
    (sdf_dists, xp) = SignedDistances.evalSignedDistances(mesh, sdf_grid, ρₙ, ρₜ)
    # points = generateGridPoints(sdf_grid)
    # Xg, Xp, mean_PD, max_PD = SignedDistances.SelectProjectedNodes(mesh, sdf_grid, xp, points) # PD (=projected distance)
    #
    # nnp = Int(length(Xg)/2)
    # IEN = [[i; i + nnp] for i = 1:nnp]
    #
    # nnp₂ = Int(length(Xg)/2)
    # IEN₂ = [[i; i + nnp₂] for i = 1:nnp₂]
    #
    # # X_combined = [Xg; Xp] 
    # # X_combined_couples = [X Xp]
    #
    # nnp = length(Xg)
    # IEN = [[i; i + nnp] for i = 1:nnp]
    #
    # Xx = vec([Xg Xp])
    #
    # DataExport.exportToVTU("lines.vtu", X, IEN, 3)
    #
    # DataExport.exportToVTU("xp.vtu", X, IEN, 1)
    #
    # # Rho2sdf.exportToVTU("xp.vtu", X_combined, IEN)
    # DataExport.exportToVTU("Xg.vtu", Xg, IEN₂, 1)
    # DataExport.exportToVTU("Xp.vtu", Xp, IEN₂, 1)


    ## Data export to VTK:
    Rho2sdf.exportStructuredPointsToVTK("sdf2-leg-" * taskName *".vtk", sdf_grid, sdf_dists, "distance")


end
