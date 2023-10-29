# using REPLVim; @async REPLVim.serve()
using Test
using Rho2sdf
using Rho2sdf.ShapeFunctions
using Rho2sdf.MeshGrid
using Rho2sdf.Derivatives
using MAT
using SymPy
using LinearAlgebra
using JLD

@testset "Rho2sdf.jl" begin

    taskName = "chapadlo"

    # # Data from Matlab:
    data = matread(taskName * ".mat")
    # data = matread("test/" * taskName * ".mat")
    part_name = "elementy_trubky.txt"
    # part_name = "test/elementy_trubky.txt"
    
    (X, IEN, rho) = MeshGrid.MeshInformations(data)

    # # input data propertis (mesh, density)
    mesh = MeshGrid.Mesh(X, IEN)
    (mesh, rho) = MeshGrid.PartOfModel(mesh, rho, part_name)
    
    ρₙ = MeshGrid.DenseInNodes(mesh, rho) # LSQ
    # ρₙ = MeshGrid.elementToNodalValues(mesh, rho) # average
    # exit()


    ## Face triangular mesh:
    # mesh = Rho2sdf.extractSurfaceTriangularMesh(mesh, ρₙ) # přepsání X a IEN (pro trojúhelníky) smazat rho

    # save("taskName" * "_triangular_mesh.jld", "mesh", mesh)
    # mesh = load("taskName" * "_triangular_mesh.jld", "mesh") # načtení chapadla (stl)

    # X = [mesh.X[:,i] for i in 1:size(mesh.X,2)]
    # IEN = [mesh.IEN[:,i] for i in 1:size(mesh.IEN,2)]
    # Rho2sdf.exportToVTU("triChapadlo.vtu", X, IEN)

    ## Grid:
    X_min, X_max = MeshGrid.getMesh_AABB(mesh.X) # vec, vec
    N = [300, 300, 300]
    sdf_grid = MeshGrid.Grid(X_min, X_max, N) # cartesian grid
   
    ## SFD from triangular mesh:
    # sdf_dists = Rho2sdf.evalSignedDiscancesOnTriangularMesh(mesh, sdf_grid) # Vector{Float64}
    
    ## SDF from densities:
    ρₜ = 0.5
    sdf_dists = Rho2sdf.evalSignedDiscances(mesh, sdf_grid, ρₙ , ρₜ)

    exit()
    ## Data export to VTK:
    # Rho2sdf.DataProcessing.exportStructuredPointsToVTK(taskName*"_sdf.vtk", sdf_grid, sdf_dists, "distance")
    Rho2sdf.exportStructuredPointsToVTK("trubka_" *taskName*"_sdf.vtk", sdf_grid, sdf_dists, "distance")



end
