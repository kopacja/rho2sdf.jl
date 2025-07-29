RUN_CUBE_HEX      = true
RUN_CUBE_HEX_REF  = true
RUN_CUBE_TET      = true
RUN_CUBE_TET_REF  = true


if (RUN_CUBE_HEX)
  @testset "Sphere_In_Cube_Hex-Manual_Steps" begin
      ## Inputs:
      taskName = "sphere_in_cube_hex"
      include("SimpleCube.jl")
      # N = 3  # Number of cells along the longest side (coarse)
      # N = 10  # Number of cells along the longest side (optimal)
      N = 80  # Number of cells along the longest side (fine)
     
      ## Generate FEM mesh structure with new type system:
      mesh, ρₙ = create_test_cube_with_linear_density()
      
      # Verify element type
      print_info("Using element type: $(mesh.element_type) with $(get_num_nodes(mesh.element_type)) nodes per element")
      
      ## Find threshold density:
      # ρₜ = @time find_threshold_for_volume(mesh, ρₙ)
      ρₜ = 0.5
      
      ## Export nodal densities:
      X_vec, IEN_vec = get_vector_format(mesh)
      VTK_CODE = 12 # VTK_HEXAHEDRON for HEX8 elements
      exportToVTU(taskName * "_nodal_densities.vtu", X_vec, IEN_vec, VTK_CODE, ρₙ)
      
      ## Grid setup:
      X_min, X_max = getMesh_AABB(mesh.X)
      sdf_grid = Grid(X_min, X_max, N, 3) # cartesian grid
      points = generateGridPoints(sdf_grid) # grid nodes
      
      ## SDF from densities (updated functions):
      (dists, xp) = evalDistances(mesh, sdf_grid, points, ρₙ, ρₜ)
      signs = @time Sign_Detection(mesh, sdf_grid, points, ρₙ, ρₜ)
      sdf_dists = dists .* signs
      
      ## Export raw SDF:
      exportSdfToVTI("$(taskName)_SDF_N-$(N).vti", sdf_grid, sdf_dists, "distance")
   end
end

if (RUN_CUBE_HEX_REF)
  @testset "Sphere_In_Cube_Ref_Hex-Manual_Steps" begin
      ## Inputs:
      taskName = "sphere_in_cube_ref_hex"
      include("CubeWithRefinedBottome.jl")
      # N = 3  # Number of cells along the longest side (coarse)
      # N = 10  # Number of cells along the longest side (optimal)
      N = 80  # Number of cells along the longest side (fine)
     
      ## Generate FEM mesh structure with new type system:
      mesh, ρₙ = create_test_cube_with_refined_bottom()
      
      # Verify element type
      print_info("Using element type: $(mesh.element_type) with $(get_num_nodes(mesh.element_type)) nodes per element")
      
      ## Find threshold density:
      # ρₜ = @time find_threshold_for_volume(mesh, ρₙ)
      ρₜ = 0.5
      
      ## Export nodal densities:
      X_vec, IEN_vec = get_vector_format(mesh)
      VTK_CODE = 12 # VTK_HEXAHEDRON for HEX8 elements
      exportToVTU(taskName * "_nodal_densities.vtu", X_vec, IEN_vec, VTK_CODE, ρₙ)
      
      ## Grid setup:
      X_min, X_max = getMesh_AABB(mesh.X)
      sdf_grid = Grid(X_min, X_max, N, 3) # cartesian grid
      points = generateGridPoints(sdf_grid) # grid nodes
      
      ## SDF from densities (updated functions):
      (dists, xp) = evalDistances(mesh, sdf_grid, points, ρₙ, ρₜ)
      signs = @time Sign_Detection(mesh, sdf_grid, points, ρₙ, ρₜ)
      sdf_dists = dists .* signs
      
      ## Export raw SDF:
      exportSdfToVTI("$(taskName)_SDF_N-$(N).vti", sdf_grid, sdf_dists, "distance")
   end
end

if (RUN_CUBE_TET)
  @testset "Sphere_In_Cube_Tet-Manual_Steps" begin
      ## Inputs:
      taskName = "sphere_in_cube_tet"
      include("SimpleCubeWithSchlafli.jl")
      # N = 3  # Number of cells along the longest side (coarse)
      # N = 10  # Number of cells along the longest side (optimal)
      N = 80  # Number of cells along the longest side (fine)
     
      ## Generate FEM mesh structure with new type system:
      print_info("\n1. Generate tet mesh...")
      mesh, ρₙ = create_test_cube_with_schlafli_tetrahedra()
      
      # Verify element type
      print_info("Using element type: $(mesh.element_type) with $(get_num_nodes(mesh.element_type)) nodes per element")
      
      ## Find threshold density:
      # ρₜ = @time find_threshold_for_volume(mesh, ρₙ)
      ρₜ = 0.5
      
      ## Export nodal densities:
      print_info("\n2. Export nodal values...")
      X_vec, IEN_vec = get_vector_format(mesh)
      VTK_CODE = 10 # VTK_TETRA for TET4 elements
      exportToVTU(taskName * "_nodal_densities.vtu", X_vec, IEN_vec, VTK_CODE, ρₙ)
      
      ## Grid setup:
      print_info("\n3. Grid setup...")
      X_min, X_max = getMesh_AABB(mesh.X)
      sdf_grid = Grid(X_min, X_max, N, 3) # cartesian grid
      points = generateGridPoints(sdf_grid) # grid nodes
      
      ## SDF from densities (updated functions):
      print_info("\n4. Compute SDF...")
      (dists, xp) = evalDistances(mesh, sdf_grid, points, ρₙ, ρₜ)
      signs = @time Sign_Detection(mesh, sdf_grid, points, ρₙ, ρₜ)
      sdf_dists = dists .* signs
      
      ## Export raw SDF:
      print_info("\n5. Export SDF to VTI...")
      exportSdfToVTI("$(taskName)_SDF_N-$(N).vti", sdf_grid, sdf_dists, "distance")
   end
end


if (RUN_CUBE_TET_REF)
  @testset "Sphere_In_Cube_Ref_Tet-Manual_Steps" begin
      ## Inputs:
      taskName = "sphere_in_cube_ref_tet"
      include("CubeWithRefinedBottomeSchlafli.jl")
      # N = 3  # Number of cells along the longest side (coarse)
      # N = 10  # Number of cells along the longest side (optimal)
      N = 80  # Number of cells along the longest side (fine)
     
      ## Generate FEM mesh structure with new type system:
      print_info("\n1. Generate tet mesh...")
      mesh, ρₙ = create_test_cube_with_refined_bottom_schlafli()
      
      # Verify element type
      print_info("Using element type: $(mesh.element_type) with $(get_num_nodes(mesh.element_type)) nodes per element")
      
      ## Find threshold density:
      # ρₜ = @time find_threshold_for_volume(mesh, ρₙ)
      ρₜ = 0.5
      
      ## Export nodal densities:
      print_info("\n2. Export nodal values...")
      X_vec, IEN_vec = get_vector_format(mesh)
      VTK_CODE = 10 # VTK_TETRA for TET4 elements
      exportToVTU(taskName * "_nodal_densities.vtu", X_vec, IEN_vec, VTK_CODE, ρₙ)
      
      ## Grid setup:
      print_info("\n3. Grid setup...")
      X_min, X_max = getMesh_AABB(mesh.X)
      sdf_grid = Grid(X_min, X_max, N, 3) # cartesian grid
      points = generateGridPoints(sdf_grid) # grid nodes
      
      ## SDF from densities (updated functions):
      print_info("\n4. Compute SDF...")
      (dists, xp) = evalDistances(mesh, sdf_grid, points, ρₙ, ρₜ)
      signs = @time Sign_Detection(mesh, sdf_grid, points, ρₙ, ρₜ)
      sdf_dists = dists .* signs
      
      ## Export raw SDF:
      print_info("\n5. Export SDF to VTI...")
      exportSdfToVTI("$(taskName)_SDF_N-$(N).vti", sdf_grid, sdf_dists, "distance")
   end
end

