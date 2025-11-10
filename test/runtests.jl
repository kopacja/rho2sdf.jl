# using Base: print_without_params
using Test
using Rho2sdf
using Rho2sdf.ElementTypes
using Rho2sdf.TerminalUtils
using Rho2sdf.PrimitiveGeometries
using Rho2sdf.ShapeFunctions
using Rho2sdf.MeshGrid
using Rho2sdf.SignedDistances
using Rho2sdf.DataExport
using Rho2sdf.SdfSmoothing
using MAT
using JLD2
using LinearAlgebra
using BenchmarkTools

@testset "Rho2sdf.jl" begin
  RUN_1HEX_EL = false
  RUN_ROOF = false
  RUN_SPHERE = false
  RUN_BEAM = true
  # Separated tests:
  RUN_CUBE_MESHES = false
  RUN_COMP_TESTS = false
  RUN_CONV_TESTS = false

  if (RUN_CUBE_MESHES)
    @time @testset "SphereInCube-Meshes" begin
      include("PrimitiveGeometriesTest/SphereInCube-Meshes.jl")
    end
  end

  if (RUN_COMP_TESTS)
    @time @testset "HexBlockSdfTest" begin
      include("HexBlockSdfTest.jl")
    end
    @time @testset "HexSphereSdfTest" begin
      include("HexSphereSdfTest.jl")
    end
  end

  if (RUN_CONV_TESTS)
    @time @testset "SphereConvergenceTest" begin
      include("ConvergenceTests/SphereConvergenceTest.jl")
    end
    @time @testset "CubeConvergenceTest" begin
      include("ConvergenceTests/CubeConvergenceTest.jl")
    end
  end

  if (RUN_1HEX_EL)
    @testset "1hex_el" begin
      taskName = "1hex_el"
      N = 15    # Number of cells along the longest side
      ρₜ = 0.5  # Threshold density (isosurface level)
      rho = [0.0, 0.0] # Elementy density (something for Mesh struct)

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
      ρₙ = [1.0, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 1.0]

      # Generate FEM mesh structure:
      shape_func = coords -> shape_functions(HEX8, coords)
      mesh = Mesh(X, IEN, rho, shape_func; element_type = HEX8)

      VTK_CODE = mesh.element_type == HEX8 ? 12 : 10 # hex/tet element
      exportToVTU(taskName * "_nodal_densities.vtu", X, IEN, VTK_CODE, ρₙ)

      # Grid:
      X_min, X_max = getMesh_AABB(mesh.X)
      sdf_grid = Grid(X_min, X_max, N, 3)
      points = generateGridPoints(sdf_grid)

      # SDF from densities:
      (dists, xp) = evalDistances(mesh, sdf_grid, points, ρₙ, ρₜ)
      signs = Sign_Detection(mesh, sdf_grid, points, ρₙ, ρₜ)
      sdf_dists = dists .* signs

      # Export SDF to vti - Paraview
      exportSdfToVTI(taskName * "_SDF.vti", sdf_grid, sdf_dists, "distance")
    end
  end

  if (RUN_ROOF)
    @testset "Roof" begin
      taskName = "Roof"
      N = 20    # Number of cells along the longest side
      ρₜ = 0.5  # Threshold density (isosurface level)

      (X, IEN, rho) = PrimitiveGeometries.selectPrimitiveGeometry("block", [2, 1, 1])

      # Generate FEM mesh structure:
      shape_func = coords -> shape_functions(HEX8, coords)
      mesh = Mesh(X, IEN, rho, shape_func; element_type = HEX8)

      # Nodal densities for roof-like isocotour:
      ρₙ = [0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 0.0, 0.0, 0.5, 0.5]

      VTK_CODE = mesh.element_type == HEX8 ? 12 : 10
      exportToVTU(taskName * "_nodal_densities.vtu", X, IEN, VTK_CODE, ρₙ)

      # Grid:
      X_min, X_max = getMesh_AABB(mesh.X)
      sdf_grid = Grid(X_min, X_max, N, 3)
      points = generateGridPoints(sdf_grid)

      # SDF from densities:
      (dists, xp) = evalDistances(
        mesh,
        sdf_grid,
        points,
        ρₙ,
        ρₜ,
        taskName = taskName,
        plot_projection_points_and_lines = true
      )
      signs = Sign_Detection(mesh, sdf_grid, points, ρₙ, ρₜ)
      sdf_dists = dists .* signs

      # Export SDF to vti - Paraview
      exportSdfToVTI(taskName * "_SDF.vti", sdf_grid, sdf_dists, "distance")

      # RBFs smoothing:
      is_interp = true
      smooth = 1
      (fine_sdf, fine_grid) =
        RBFs_smoothing(mesh, sdf_dists, sdf_grid, is_interp, smooth, taskName)
      export_sdf_results(fine_sdf, fine_grid, sdf_grid, taskName, smooth, is_interp)
    end
  end

  if (RUN_SPHERE)
    @testset "Sphere" begin
      taskName = "sphere"
      N = 10    # Number of cells along the longest side

      # Load data from mat format (Matlab):
      data = matread(taskName * ".mat")
      (X, IEN, rho) = MeshInformations(data)

      # Generate FEM mesh structure:
      shape_func = coords -> shape_functions(HEX8, coords)
      mesh = Mesh(X, IEN, rho, shape_func; element_type = HEX8)

      # Export input data to Paraview (density field):
      InputDataToVTU(mesh, taskName * "-input_data")

      ρₙ = DenseInNodes(mesh, rho) # Compute nodal densities
      ρₜ = @time find_threshold_for_volume(mesh, ρₙ) # Find threshold density (to maintain folume frac)

      # Export nodal densities to Paraview:
      VTK_CODE = mesh.element_type == HEX8 ? 12 : 10
      exportToVTU(taskName * "_nodal_densities.vtu", X, IEN, VTK_CODE, ρₙ)

      # Grid:
      X_min, X_max = getMesh_AABB(mesh.X)
      sdf_grid = Grid(X_min, X_max, N, 3)
      points = generateGridPoints(sdf_grid)

      # SDF from densities:
      (dists, xp) = evalDistances(mesh, sdf_grid, points, ρₙ, ρₜ)
      signs = @time Sign_Detection(mesh, sdf_grid, points, ρₙ, ρₜ)
      sdf_dists = dists .* signs

      # Export SDF to vti - Paraview
      exportSdfToVTI(taskName * "_SDF.vti", sdf_grid, sdf_dists, "distance")

      # RBFs smoothing:
      is_interp = true
      smooth = 1
      (fine_sdf, fine_grid) =
        RBFs_smoothing(mesh, sdf_dists, sdf_grid, is_interp, smooth, taskName)
      export_sdf_results(fine_sdf, fine_grid, sdf_grid, taskName, smooth, is_interp)
    end
  end

  if (RUN_BEAM)
    @testset "Cantilever_beam_hex_vfrac04" begin
      taskName = "cantilever_beam_vfrac_04"

      # Read FEM mesh:
      data = matread("$(taskName).mat")
      (X, IEN, rho) = MeshInformations(data)
      IEN = [subvector .- 1 for subvector in IEN]  # Data correction

      # Custom options
      options = Rho2sdfOptions(
        export_input_data = true,         # Export input (raw) data to Paraview
        threshold_density = 0.518555,     # Value for isocotour (0, 1)
        sdf_grid_setup = :automatic,      # Automatic/manual grid setup
        export_nodal_densities = true,    # Export nodal field to Paraview
        export_raw_sdf = true,            # Export not smoothed SDF to Paraview
        rbf_interp = true,                # Interp/approx SDF values using RBFs
        rbf_grid = :same                  # Same/fine grid for RBFs interp/approx
      )

      result = rho2sdf("beam", X, IEN, rho, options = options)
    end
  end
end
