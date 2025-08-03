taskName = "roof_geometry"
N = 20    # Number of cells along the longest side
ρₜ = 0.5  # Threshold density (isosurface level)

(X, IEN, rho) = PrimitiveGeometries.selectPrimitiveGeometry("block", [2, 1, 1])

# Generate FEM mesh structure:
shape_func = coords -> shape_functions(HEX8, coords)
mesh = Mesh(X, IEN, rho, shape_func; element_type=HEX8)

# Nodal densities for roof-like isocotour:
ρₙ = [0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 0.0, 0.0, 0.5, 0.5]

VTK_CODE = mesh.element_type == HEX8 ? 12 : 10
exportToVTU(taskName * "_nodal_densities.vtu", X, IEN, VTK_CODE, ρₙ)

# Grid:
X_min, X_max = getMesh_AABB(mesh.X)
sdf_grid = Grid(X_min, X_max, N, 2)
points = generateGridPoints(sdf_grid)

# SDF from densities:
(dists, xp) = evalDistances(mesh, sdf_grid, points, ρₙ, ρₜ, taskName=taskName, plot_projection_points_and_lines=true)
signs = Sign_Detection(mesh, sdf_grid, points, ρₙ, ρₜ)
sdf_dists = dists .* signs

# Export SDF to vti - Paraview
exportSdfToVTI(taskName * "_SDF.vti", sdf_grid, sdf_dists, "distance")

print_success("01_SDF_VALIDATION - (roof_geometry) done")
