# using REPLVim; @async REPLVim.serve()
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

@testset "Rho2sdf.jl" begin

    # @time @testset "PrimitiveGeometriesTest" begin include("PrimitiveGeometriesTest/runtests.jl") end
    # @time @testset "MeshGridTest" begin include("MeshGridTest/runtests.jl") end
    #     @time @testset "SignedDistancesTest" begin include("SignedDistancesTest/runtests.jl") end
    # end
    # exit()
    #     
    # # Data from Matlab:
    # taskName = "chapadlo"

    RUN_PLANE = false
    RUN_SPHERE = false
    RUN_CHAPADLO = true

    if (RUN_PLANE)
        @testset "Plane" begin

            taskName = "plane"
            N = 1  # Number of cells along the longest side
            ρₜ = 0.5 # Threshold density (isosurface level)

            X = [
                [-1.0, -1.0, -1.0],
                [1.0, -1.0, -1.0],
                [1.0, 1.0, -1.0],
                [-1.0, 1.0, -1.0],
                [-1.0, -1.0, 1.0],
                [1.0, -1.0, 1.0],
                [1.0, 1.0, 1.0],
                [-1.0, 1.0, 1.0],
            ]
            IEN = [[1, 2, 3, 4, 5, 6, 7, 8]]
            ρₙ = [0.0, 0.0, 0.0, 0.0, 1, 1, 1, 1]

            ## Generate FEM mesh structure:
            mesh = MeshGrid.Mesh(X, IEN)

            VTK_CODE = 12 # https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html
            Rho2sdf.exportToVTU(taskName * "_nodal_densities.vtu", X, IEN, VTK_CODE, ρₙ)

            ## Grid:
            X_min, X_max = MeshGrid.getMesh_AABB(mesh.X)
            sdf_grid = MeshGrid.Grid(X_min, X_max, N, 0) # cartesian grid

            ## SDF from densities:
            sdf_dists = SignedDistances.evalSignedDistances(mesh, sdf_grid, ρₙ, ρₜ)

            ## Export to VTK:
            Rho2sdf.exportStructuredPointsToVTK(taskName * "_sdf.vtk", sdf_grid, sdf_dists, "distance")
        end
    end

    if (RUN_SPHERE)
        @testset "Sphere" begin
            ## Inputs:
            taskName = "sphere"
            N = 3  # Number of cells along the longest side
            ρₜ = 0.5 # Threshold density (isosurface level)

            ## Read FEM mesh:
            data = matread(taskName * ".mat")
            (X, IEN, rho) = MeshGrid.MeshInformations(data)

            ## Generate FEM mesh structure:
            mesh = MeshGrid.Mesh(X, IEN)

            ## Map elemental densities to the nodes:
            ρₙ = MeshGrid.DenseInNodes(mesh, rho) # LSQ
            #ρₙ = MeshGrid.elementToNodalValues(mesh, rho) # average

            VTK_CODE = 12 # https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html
            Rho2sdf.exportToVTU(taskName * "_nodal_densities.vtu", X, IEN, VTK_CODE, ρₙ)

            ## Grid:
            X_min, X_max = MeshGrid.getMesh_AABB(mesh.X)
            sdf_grid = MeshGrid.Grid(X_min, X_max, N, 0) # cartesian grid

            ## SDF from densities:
            sdf_dists = SignedDistances.evalSignedDistances(mesh, sdf_grid, ρₙ, ρₜ)

            ## Export to VTK:
            Rho2sdf.exportStructuredPointsToVTK(taskName * "_sdf.vtk", sdf_grid, sdf_dists, "distance")
        end
    end

    if (RUN_CHAPADLO)
        @testset "Chapadlo" begin
            ## Inputs:
            taskName = "chapadlo"
            N = 50  # Number of cells along the longest side
            ρₜ = 0.5 # Threshold density (isosurface level)

            ## Read FEM mesh:
            data = matread(taskName * ".mat")
            (X, IEN, rho) = MeshGrid.MeshInformations(data)
            #Z,idx_Z = findall(x->X[3,i] > 50 for i in [1:size(X,2)])

            ## Generate FEM mesh structure:
            mesh = MeshGrid.Mesh(X, IEN)

            ## Map elemental densities to the nodes:
            ρₙ = MeshGrid.DenseInNodes(mesh, rho) # LSQ
            #ρₙ = MeshGrid.elementToNodalValues(mesh, rho) # average

            VTK_CODE = 12 # https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html
            Rho2sdf.exportToVTU(taskName * "_nodal_densities.vtu", X, IEN, VTK_CODE, ρₙ)

            ## Grid:
            X_min, X_max = MeshGrid.getMesh_AABB(mesh.X)
            X_min[2] = 0
            X_min[3] = 50
            sdf_grid = MeshGrid.Grid(X_min, X_max, N, 0) # cartesian grid

            ## SDF from densities:
            sdf_dists = SignedDistances.evalSignedDistances(mesh, sdf_grid, ρₙ, ρₜ)

            ## Export to VTK:
            Rho2sdf.exportStructuredPointsToVTK(taskName * "_sdf.vtk", sdf_grid, sdf_dists, "distance")
        end
    end
end