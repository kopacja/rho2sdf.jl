using Base: print_without_params
# using REPLVim; @async REPLVim.serve()
using Test
using Rho2sdf
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

  # @time @testset "PrimitiveGeometriesTest" begin include("PrimitiveGeometriesTest/runtests.jl") end
  # @time @testset "MeshGridTest" begin include("MeshGridTest/runtests.jl") end
  # @time @testset "SignedDistancesTest" begin include("SignedDistancesTest/runtests.jl") end
  # @time @testset "SeparatedTests" begin include("SeparatedTests/runtests.jl") end
  #
  ### Tests on geometries: ###
  # @testset "TestOnLegGripper" begin
  # include("SeparatedTests/TestOnLegGripper.jl")
  # end
  # @testset "TestOnLegGripperSTL" begin include("SeparatedTests/TestOnLegGripperSTL.jl") end
  # @testset "TestOnGripperMC" begin include("SeparatedTests/TestOnGripperMC.jl") end
  # @testset "TestOnPrimitiveGeometry" begin include("SeparatedTests/TestOnPrimitiveGeometry.jl") end
  # end
  # exit()
  #     
  # # Data from Matlab:
  # taskName = "chapadlo"

  RUN_PLANE = false
  RUN_BLOCK = true
  RUN_SPHERE = true
  RUN_CHAPADLO = true
  RUN_CHAPADLO_cele = false

  if (RUN_PLANE)
    @testset "Plane" begin

      taskName = "plane"
      N = 15  # Number of cells along the longest side
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
      # ρₙ = [0.0, 0.0, 0.0, 0.0, 1, 1, 1, 1]
      # ρₙ = [0.0, 0.0, 0.0, 0.0, 1, 1, 1, 0.5]
      # ρₙ = [1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0] # dve roviny 4 uzlové
      # ρₙ = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1]
      # ρₙ = [0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 1]
      ρₙ = [1.0, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 1]
      # ρₙ = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0] # 1, 7

      ## Generate FEM mesh structure:
      mesh = MeshGrid.Mesh(X, IEN, rho, hex8_shape)

      VTK_CODE = 12 # https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html
      Rho2sdf.exportToVTU(taskName * "_nodal_densities.vtu", X, IEN, VTK_CODE, ρₙ)

      ## Grid:
      X_min, X_max = MeshGrid.getMesh_AABB(mesh.X)
      sdf_grid = MeshGrid.Grid(X_min, X_max, N, 3) # cartesian grid
      points = MeshGrid.generateGridPoints(sdf_grid) # uzly pravidelné mřížky

      ## SDF from densities:
      (dists, xp) = SignedDistances.evalDistances(mesh, sdf_grid, points, ρₙ, ρₜ)
      signs = SignedDistances.Sign_Detection(mesh, sdf_grid, points, ρₙ, ρₜ)
      sdf_dists = dists .* signs

      ## Export to VTK:
      Rho2sdf.exportStructuredPointsToVTK(taskName * "_SDF.vtk", sdf_grid, sdf_dists, "distance")

      # @save "Z_$(taskName)_xp.jld2" xp
      # @save "Z_$(taskName)_Mesh.jld2" mesh
      # @save "Z_$(taskName)_Grid.jld2" sdf_grid
      # @save "Z_$(taskName)_Points.jld2" points
      # @save "Z_$(taskName)_SDF.jld2" sdf_dists

    end
  end

  if (RUN_BLOCK)
    @testset "Block" begin
      ## Inputs:
      taskName = "block"

      N = 10  # Number of cells along the longest side
      ρₜ = 0.5 # Threshold density (isosurface level)

      (X, IEN, rho) = PrimitiveGeometries.selectPrimitiveGeometry("block", [2, 1, 1])
      # ρₙ = [0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.5, 0.5]

      mesh = MeshGrid.Mesh(X, IEN, rho, hex8_shape)
      # ρₙ = MeshGrid.DenseInNodes(mesh, rho) # LSQ

      # Modif ρₙ:
      ρₙ = [0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 0.0, 0.0, 0.5, 0.5]
      
      # ρₜ= find_threshold_for_volume(mesh, ρₙ)

      VTK_CODE = 12 # https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html
      Rho2sdf.exportToVTU(taskName * "_nodal_densities.vtu", X, IEN, VTK_CODE, ρₙ)

      ## Grid:
      X_min, X_max = MeshGrid.getMesh_AABB(mesh.X)
      sdf_grid = MeshGrid.Grid(X_min, X_max, N, 3) # cartesian grid
      points = MeshGrid.generateGridPoints(sdf_grid) # uzly pravidelné mřížky

      ## SDF from densities:
      (dists, xp) = SignedDistances.evalDistances(mesh, sdf_grid, points, ρₙ, ρₜ)
      signs = SignedDistances.Sign_Detection(mesh, sdf_grid, points, ρₙ, ρₜ)
      sdf_dists = dists .* signs

      ## Export to VTK:
      exportSdfToVTI(taskName * "_SDF.vti", sdf_grid, sdf_dists, "distance")

      # RBF smoothing:
      is_interp = true
      smooth = 1
      (fine_sdf, fine_grid) = RBFs_smoothing(mesh, sdf_dists, sdf_grid, is_interp, smooth, taskName) # interpolation == true, aproximation == false, smooth
      export_sdf_results(fine_sdf, fine_grid, sdf_grid, taskName, smooth, is_interp)
     
      # fig = visualize_stable_isosurface(fine_LSF)
      # display(fig)

      # @save "Z_$(taskName)_xp.jld2" xp
      # @save "Z_$(taskName)_Mesh.jld2" mesh
      # @save "Z_$(taskName)_Grid.jld2" sdf_grid
      # @save "Z_$(taskName)_Points.jld2" points
      # @save "Z_$(taskName)_SDF.jld2" sdf_dists
      # @save "Z_$(taskName)_rho.jld2" ρₙ

    end
  end


  if (RUN_SPHERE)
    @testset "Sphere" begin
      ## Inputs:
      taskName = "sphere"

      N = 10  # Number of cells along the longest side
      # ρₜ = 0.5 # Threshold density (isosurface level)
      ## Read FEM mesh:
      data = matread(taskName * ".mat")
      # data = matread("test/" * taskName * ".mat")
      (X, IEN, rho) = MeshGrid.MeshInformations(data)

      ## Generate FEM mesh structure:
      mesh = MeshGrid.Mesh(X, IEN, rho, hex8_shape)

      ## Map elemental densities to the nodes:
      ρₙ = MeshGrid.DenseInNodes(mesh, rho) # LSQ
      #ρₙ = MeshGrid.elementToNodalValues(mesh, rho) # average

      ρₜ= @time find_threshold_for_volume(mesh, ρₙ)

      VTK_CODE = 12 # https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html
      Rho2sdf.exportToVTU(taskName * "_nodal_densities.vtu", X, IEN, VTK_CODE, ρₙ)

      ## Grid:
      X_min, X_max = MeshGrid.getMesh_AABB(mesh.X)
      sdf_grid = MeshGrid.Grid(X_min, X_max, N, 3) # cartesian grid
      points = MeshGrid.generateGridPoints(sdf_grid) # uzly pravidelné mřížky

      ## SDF from densities:
      (dists, xp) = SignedDistances.evalDistances(mesh, sdf_grid, points, ρₙ, ρₜ)
      signs = @time SignedDistances.Sign_Detection(mesh, sdf_grid, points, ρₙ, ρₜ)
      sdf_dists = dists .* signs
      
      ## Export to VTK:
      exportSdfToVTI(taskName * "_SDF.vti", sdf_grid, sdf_dists, "distance")

      # RBF smoothing:
      is_interp = true
      smooth = 1
      (fine_sdf, fine_grid) = RBFs_smoothing(mesh, sdf_dists, sdf_grid, is_interp, smooth, taskName) # interpolation == true, aproximation == false, smooth
      export_sdf_results(fine_sdf, fine_grid, sdf_grid, taskName, smooth, is_interp)
     
      # @save "Z_$(taskName)_xp.jld2" xp
      # @save "Z_$(taskName)_Mesh.jld2" mesh
      # @save "Z_$(taskName)_Grid.jld2" sdf_grid
      # @save "Z_$(taskName)_Points.jld2" points
      # @save "Z_$(taskName)_SDF.jld2" sdf_dists
      # @save "Z_$(taskName)_rho.jld2" ρₙ

    end
  end

  if (RUN_BEAM)
    @testset "Chapadlo" begin
      ## Inputs:
      taskName = "cantilever_beam_vfrac_04"
      
      ## Read FEM mesh:
      data = matread(taskName * ".mat")
      (X, IEN, rho) = MeshInformations(data)
      IEN = [subvector .- 1 for subvector in IEN]  # This will work

      # Custom options
      options = Rho2sdfOptions(
          threshold_density=0.5,
          sdf_grid_setup=:noninteractive,
          export_nodal_densities=false,
          export_raw_sdf=false,
          rbf_interp=true,
          rbf_grid=:normal
      )

      result = rho2sdf("beam", X, IEN, rho, options=options)

    end
  end

  if (RUN_CHAPADLO)
    @testset "Chapadlo" begin
      ## Inputs:
      taskName = "chapadlo"
      N = 60  # Number of cells along the longest side
      # ρₜ = 0.5 # Threshold density (isosurface level)

      ## Read FEM mesh:
      data = matread(taskName * ".mat")
      (X, IEN, rho) = MeshGrid.MeshInformations(data)
      #Z,idx_Z = findall(x->X[3,i] > 50 for i in [1:size(X,2)])

      ## Generate FEM mesh structure:
      mesh = MeshGrid.Mesh(X, IEN, rho, hex8_shape)

      ## Map elemental densities to the nodes:
      ρₙ = MeshGrid.DenseInNodes(mesh, rho) # LSQ
      #ρₙ = MeshGrid.elementToNodalValues(mesh, rho) # average

      ρₜ= @time find_threshold_for_volume(mesh, ρₙ)

      VTK_CODE = 12 # https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html
      Rho2sdf.exportToVTU(taskName * "_nodal_densities.vtu", X, IEN, VTK_CODE, ρₙ)

      ## Grid:
      X_min, X_max = MeshGrid.getMesh_AABB(mesh.X)
      X_min[2] = 0
      X_min[3] = 50
      sdf_grid = MeshGrid.Grid(X_min, X_max, N, 3) # cartesian grid
      points = MeshGrid.generateGridPoints(sdf_grid) # uzly pravidelné mřížky

      ## SDF from densities:
      (dists, xp) = SignedDistances.evalDistances(mesh, sdf_grid, points, ρₙ, ρₜ)
      signs = SignedDistances.Sign_Detection(mesh, sdf_grid, points, ρₙ, ρₜ)
      sdf_dists = dists .* signs
      
      ## Export to VTK:
      exportSdfToVTI(taskName * "_SDF.vti", sdf_grid, sdf_dists, "distance")

      # RBF smoothing:
      is_interp = true
      smooth = 1
      (fine_sdf, fine_grid) = RBFs_smoothing(mesh, sdf_dists, sdf_grid, is_interp, smooth, taskName) # interpolation == true, aproximation == false, smooth
      export_sdf_results(fine_sdf, fine_grid, sdf_grid, taskName, smooth, is_interp)
    
      # @save "Z_$(taskName)_xp.jld2" xp
      # @save "Z_$(taskName)_Mesh.jld2" mesh
      # @save "Z_$(taskName)_Grid.jld2" sdf_grid
      # @save "Z_$(taskName)_Points.jld2" points
      # @save "Z_$(taskName)_SDF.jld2" sdf_dists
      # @save "Z_$(taskName)_rho.jld2" ρₙ

    end
  end

  if (RUN_CHAPADLO_cele)
    @testset "Chapadlo" begin
      ## Inputs:
      taskName = "chapadlo"
      # ρₜ = 0.5 # Threshold density (isosurface level)

      ## Read FEM mesh:
      data = matread(taskName * ".mat")
      (X, IEN, rho) = MeshGrid.MeshInformations(data)

      ## Generate FEM mesh structure:
      mesh = MeshGrid.Mesh(X, IEN, rho, hex8_shape)

      ## Grid:
      # sdf_grid = MeshGrid.interactive_sdf_grid_setup(mesh)
      sdf_grid = MeshGrid.noninteractive_sdf_grid_setup(mesh, 2.0)
      points = MeshGrid.generateGridPoints(sdf_grid) # uzly pravidelné mřížky

      ## Map elemental densities to the nodes:
      ρₙ = MeshGrid.DenseInNodes(mesh, rho) # LSQ
      
      ρₜ= find_threshold_for_volume(mesh, ρₙ)

      VTK_CODE = 12 # https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html
      Rho2sdf.exportToVTU(taskName * "_nodal_densities.vtu", X, IEN, VTK_CODE, ρₙ)

      ## SDF from densities:
      (dists, xp) = SignedDistances.evalDistances(mesh, sdf_grid, points, ρₙ, ρₜ)
      signs = SignedDistances.Sign_Detection(mesh, sdf_grid, points, ρₙ, ρₜ)
      sdf_dists = dists .* signs

      ## Export to VTK:
      B = round(sdf_grid.cell_size, digits=4)
      exportSdfToVTI(taskName * "_SDF_CellSize-" * string(B) * ".vti", sdf_grid, sdf_dists, "distance")

      # RBF smoothing:
      is_interp = true
      smooth = 1
      (fine_sdf, fine_grid) = RBFs_smoothing(mesh, sdf_dists, sdf_grid, is_interp, smooth, taskName) # interpolation == true, aproximation == false, smooth
      export_sdf_results(fine_sdf, fine_grid, sdf_grid, taskName, smooth, is_interp)

      # @save "Z_$(taskName)_cele_xp_B-$(B).jld2" xp
      # @save "Z_$(taskName)_cele_SDF_B-$(B).jld2" sdf_dists
      # @save "Z_$(taskName)_cele_Grid_B-$(B).jld2" sdf_grid
      # @save "Z_$(taskName)_cele_Points_B-$(B).jld2" points
      # @save "Z_$(taskName)_cele_Mesh_B-$(B).jld2" mesh
      #
    end
  end
end
