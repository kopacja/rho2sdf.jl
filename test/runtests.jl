# using REPLVim; @async REPLVim.serve()
using Test
using Rho2sdf
using Rho2sdf.ShapeFunctions
using Rho2sdf.GenerateMesh
using MAT
using SymPy
using LinearAlgebra
using JLD

@testset "Rho2sdf.jl" begin

    taskName = "chapadlo"

    # # Data from Matlab:
    data = matread(taskName * ".mat")
    # # data = matread("test/" * taskName * ".mat")
    # rho = vec(data["rho"])
    # mesh = data["msh"]
    # X = mesh["X"] # Matrix
    # X = [X[:, i] for i in axes(X,2)]
    # IEN = convert(Array{Int64}, mesh["IEN"] .+ 1) # Matrix
    # IEN = [IEN[:, i] for i in axes(IEN, 2)] # Vector{Vector{Int64}}
    (X, IEN, rho) = GenerateMesh.MeshInformations(data)

    # # input data propertis (mesh, density)
    mesh = Rho2sdf.Mesh(X, IEN)
    # (mesh, rho) = GenerateMesh.MeshInformations(data)
    ρₙ = Rho2sdf.elementToNodalValues(mesh, rho) # nodal values calculation (AVERAGE!! -> least squares)
    # ρₙ = GenerateMesh.elementToNodalValues(mesh, rho) # nodal values calculation (AVERAGE!! -> least squares)


    ## Face triangular mesh:
    mesh = Rho2sdf.extractSurfaceTriangularMesh(mesh, ρₙ) # přepsání X a IEN (pro trojúhelníky) smazat rho

    save("taskName" * "_triangular_mesh.jld", "mesh", mesh)
    mesh = load("taskName" * "_triangular_mesh.jld", "mesh") # načtení chapadla (stl)

    X = [mesh.X[:,i] for i in 1:size(mesh.X,2)]
    IEN = [mesh.IEN[:,i] for i in 1:size(mesh.IEN,2)]
    Rho2sdf.exportToVTU("triChapadlo.vtu", X, IEN)

    ## Grid:
    X_min, X_max = Rho2sdf.getMesh_AABB(mesh.X) # vec, vec
    N = [300, 300, 300]
    sdf_grid = Rho2sdf.Grid(X_min, X_max, N) # cartesian grid
   
    ## SFD from triangular mesh:
    sdf_dists = Rho2sdf.evalSignedDiscancesOnTriangularMesh(mesh, sdf_grid) # Vector{Float64}
    
    ## SDF from densities:
    # ρₜ = 0.5
    # sdf_dists = Rho2sdf.evalSignedDiscances(mesh, sdf_grid, ρₙ , ρₜ)

    ## Data export to VTK:
    # Rho2sdf.DataProcessing.exportStructuredPointsToVTK(taskName*"_sdf.vtk", sdf_grid, sdf_dists, "distance")



end
