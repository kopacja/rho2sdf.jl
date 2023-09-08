# using REPLVim; @async REPLVim.serve()
using Test
using Rho2sdf
using Rho2sdf.ShapeFunctions
using MAT
using SymPy
using LinearAlgebra
using JLD

@testset "Rho2sdf.jl" begin

    taskName = "chapadlo"

    # Data from Matlab:
    data = matread(taskName * ".mat")
    # data = matread("test/" * taskName * ".mat")
    rho = vec(data["rho"])
    mesh = data["msh"]
    X = mesh["X"] # Matrix
    X = [X[:, i] for i in axes(X,2)]
    IEN = convert(Array{Int64}, mesh["IEN"] .+ 1) # Matrix
    IEN = [IEN[:, i] for i in axes(IEN, 2)] # Vector{Vector{Int64}}

    # input data propertis (mesh, density)
    mesh = Rho2sdf.Mesh(X, IEN)
    ρₙ = Rho2sdf.elementToNodalValues(mesh, rho) # nodal values calculation (AVERAGE!! -> least squares)

    # face triangular mesh from
    # mesh = Rho2sdf.extractSurfaceTriangularMesh(mesh, ρₙ) # přepsání X a IEN (pro trojúhelníky) smazat rho

    # save("taskName" * "_triangular_mesh.jld", "mesh", mesh)
    # mesh = load("taskName" * "_triangular_mesh.jld", "mesh") # načtení chapadla (stl)

    # X = [mesh.X[:,i] for i in 1:size(mesh.X,2)]
    # IEN = [mesh.IEN[:,i] for i in 1:size(mesh.IEN,2)]
    # Rho2sdf.exportToVTU("triChapadlo.vtu", X, IEN)

    X_min, X_max = Rho2sdf.getMesh_AABB(mesh.X) # vec, vec
    N = [300, 300, 300]
    sdf_grid = Rho2sdf.Grid(X_min, X_max, N) # cartesian grid

    ρₜ = 0.5

    # sdf_dists = Rho2sdf.evalSignedDiscancesOnTriangularMesh(mesh, sdf_grid) # Vector{Float64}
    sdf_dists = Rho2sdf.evalSignedDiscances(mesh, sdf_grid, ρₙ , ρₜ)

    Rho2sdf.exportStructuredPointsToVTK(taskName*"_sdf.vtk", sdf_grid, sdf_dists, "distance")






#     H = ["-1/8*(ξ₁-1)*(ξ₂-1)*(ξ₃-1)",
#          "1/8*(ξ₁+1)*(ξ₂-1)*(ξ₃-1)",
#          "-1/8*(ξ₁+1)*(ξ₂+1)*(ξ₃-1)",
#          "1/8*(ξ₁-1)*(ξ₂+1)*(ξ₃-1)",
#
#          "1/8*(ξ₁-1)*(ξ₂-1)*(ξ₃+1)",
#          "-1/8*(ξ₁+1)*(ξ₂-1)*(ξ₃+1)",
#          "1/8*(ξ₁+1)*(ξ₂+1)*(ξ₃+1)",
#          "-1/8*(ξ₁-1)*(ξ₂+1)*(ξ₃+1)"]
#
# ξ₁,ξ₂,ξ₃ = Sym("ξ₁, ξ₂, ξ₃")
#
# println("d¹N_dξ¹[6] = [")
# println(diff(-1/8*(ξ₁+1)*(ξ₂-1)*(ξ₃+1), ξ₁))
# println(diff(-1/8*(ξ₁+1)*(ξ₂-1)*(ξ₃+1), ξ₂))
# println(diff(-1/8*(ξ₁+1)*(ξ₂-1)*(ξ₃+1), ξ₃),"]")
#
# println("d²N_dξ²[6] = [")
# println(diff(-1/8*(ξ₁+1)*(ξ₂-1)*(ξ₃+1), ξ₁, ξ₁),"\t",
#         diff(-1/8*(ξ₁+1)*(ξ₂-1)*(ξ₃+1), ξ₁, ξ₂),"\t",
#         diff(-1/8*(ξ₁+1)*(ξ₂-1)*(ξ₃+1), ξ₁, ξ₃))
# println(diff(-1/8*(ξ₁+1)*(ξ₂-1)*(ξ₃+1), ξ₂, ξ₁),"\t",
#         diff(-1/8*(ξ₁+1)*(ξ₂-1)*(ξ₃+1), ξ₂, ξ₂),"\t",
#         diff(-1/8*(ξ₁+1)*(ξ₂-1)*(ξ₃+1), ξ₂, ξ₃))
# println(diff(-1/8*(ξ₁+1)*(ξ₂-1)*(ξ₃+1), ξ₃, ξ₁),"\t",
#         diff(-1/8*(ξ₁+1)*(ξ₂-1)*(ξ₃+1), ξ₃, ξ₂),"\t",
#         diff(-1/8*(ξ₁+1)*(ξ₂-1)*(ξ₃+1), ξ₃, ξ₃),"]")
#
# println("d³N_dξ³[6] = [")
# println(diff(-1/8*(ξ₁+1)*(ξ₂-1)*(ξ₃+1), ξ₁, ξ₁, ξ₁),"\t",
#         diff(-1/8*(ξ₁+1)*(ξ₂-1)*(ξ₃+1), ξ₁, ξ₂, ξ₁),"\t",
#         diff(-1/8*(ξ₁+1)*(ξ₂-1)*(ξ₃+1), ξ₁, ξ₃, ξ₁),";")
# println(diff(-1/8*(ξ₁+1)*(ξ₂-1)*(ξ₃+1), ξ₂, ξ₁, ξ₁),"\t",
#         diff(-1/8*(ξ₁+1)*(ξ₂-1)*(ξ₃+1), ξ₂, ξ₂, ξ₁),"\t",
#         diff(-1/8*(ξ₁+1)*(ξ₂-1)*(ξ₃+1), ξ₂, ξ₃, ξ₁),";")
# println(diff(-1/8*(ξ₁+1)*(ξ₂-1)*(ξ₃+1), ξ₃, ξ₁, ξ₁),"\t",
#         diff(-1/8*(ξ₁+1)*(ξ₂-1)*(ξ₃+1), ξ₃, ξ₂, ξ₁),"\t",
#         diff(-1/8*(ξ₁+1)*(ξ₂-1)*(ξ₃+1), ξ₃, ξ₃, ξ₁),";;;")
#
# println(diff(-1/8*(ξ₁+1)*(ξ₂-1)*(ξ₃+1), ξ₁, ξ₁, ξ₂),"\t",
#         diff(-1/8*(ξ₁+1)*(ξ₂-1)*(ξ₃+1), ξ₁, ξ₂, ξ₂),"\t",
#         diff(-1/8*(ξ₁+1)*(ξ₂-1)*(ξ₃+1), ξ₁, ξ₃, ξ₂),";")
# println(diff(-1/8*(ξ₁+1)*(ξ₂-1)*(ξ₃+1), ξ₂, ξ₁, ξ₂),"\t",
#         diff(-1/8*(ξ₁+1)*(ξ₂-1)*(ξ₃+1), ξ₂, ξ₂, ξ₂),"\t",
#         diff(-1/8*(ξ₁+1)*(ξ₂-1)*(ξ₃+1), ξ₂, ξ₃, ξ₂),";")
# println(diff(-1/8*(ξ₁+1)*(ξ₂-1)*(ξ₃+1), ξ₃, ξ₁, ξ₂),"\t",
#         diff(-1/8*(ξ₁+1)*(ξ₂-1)*(ξ₃+1), ξ₃, ξ₂, ξ₂),"\t",
#         diff(-1/8*(ξ₁+1)*(ξ₂-1)*(ξ₃+1), ξ₃, ξ₃, ξ₂),";;;")
#
#
# println(diff(-1/8*(ξ₁+1)*(ξ₂-1)*(ξ₃+1), ξ₁, ξ₁, ξ₃),"\t",
#         diff(-1/8*(ξ₁+1)*(ξ₂-1)*(ξ₃+1), ξ₁, ξ₂, ξ₃),"\t",
#         diff(-1/8*(ξ₁+1)*(ξ₂-1)*(ξ₃+1), ξ₁, ξ₃, ξ₃),";")
# println(diff(-1/8*(ξ₁+1)*(ξ₂-1)*(ξ₃+1), ξ₂, ξ₁, ξ₃),"\t",
#         diff(-1/8*(ξ₁+1)*(ξ₂-1)*(ξ₃+1), ξ₂, ξ₂, ξ₃),"\t",
#         diff(-1/8*(ξ₁+1)*(ξ₂-1)*(ξ₃+1), ξ₂, ξ₃, ξ₃),";")
# println(diff(-1/8*(ξ₁+1)*(ξ₂-1)*(ξ₃+1), ξ₃, ξ₁, ξ₃),"\t",
#         diff(-1/8*(ξ₁+1)*(ξ₂-1)*(ξ₃+1), ξ₃, ξ₂, ξ₃),"\t",
#         diff(-1/8*(ξ₁+1)*(ξ₂-1)*(ξ₃+1), ξ₃, ξ₃, ξ₃),"]")



end
