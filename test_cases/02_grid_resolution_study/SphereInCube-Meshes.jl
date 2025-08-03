using LinearAlgebra
using Statistics

if (RUN_CUBE_HEX)
   ## Inputs:
   taskName = "sphere_in_cube_hex"
   include("SimpleCube.jl")
   
   ## Generate FEM mesh structure with new type system:
   mesh, ρₙ = create_test_cube_with_linear_density()
   
   # Verify element type
   print_info("Using element type: $(mesh.element_type) with $(get_num_nodes(mesh.element_type)) nodes per element")
   
   ## Threshold density:
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

   print_success("01_SDF_VALIDATION - (simple cube N = $(N)) done")
end

if (RUN_CUBE_HEX_REF)
   ## Inputs:
   taskName = "sphere_in_cube_ref_hex"
   include("CubeWithRefinedBottome.jl")
   
   ## Generate FEM mesh structure with new type system:
   mesh, ρₙ = create_test_cube_with_refined_bottom()
   
   # Verify element type
   print_info("Using element type: $(mesh.element_type) with $(get_num_nodes(mesh.element_type)) nodes per element")
   
   ## Threshold density:
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
   
   print_success("01_SDF_VALIDATION - (cube with refined bottome N = $(N)) done")
end

if (RUN_CUBE_TET)
   ## Inputs:
   taskName = "sphere_in_cube_tet"
   include("SimpleCubeWithSchlafli.jl")
   
   ## Generate FEM mesh structure with new type system:
   print_info("\n1. Generate tet mesh...")
   mesh, ρₙ = create_test_cube_with_schlafli_tetrahedra()
   
   # Verify element type
   print_info("Using element type: $(mesh.element_type) with $(get_num_nodes(mesh.element_type)) nodes per element")
   
   ## Threshold density:
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
   
   print_success("01_SDF_VALIDATION - (simple cube with schlafli N = $(N)) done")
end


if (RUN_CUBE_TET_REF)
   ## Inputs:
   taskName = "sphere_in_cube_ref_tet"
   include("CubeWithRefinedBottomeSchlafli.jl")

   ## Generate FEM mesh structure with new type system:
   print_info("\n1. Generate tet mesh...")
   mesh, ρₙ = create_test_cube_with_refined_bottom_schlafli()
   
   # Verify element type
   print_info("Using element type: $(mesh.element_type) with $(get_num_nodes(mesh.element_type)) nodes per element")
   
   ## Threshold density:
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
   
   print_success("01_SDF_VALIDATION - (cube with refined bottome schlafli cube N = $(N)) done")
end

