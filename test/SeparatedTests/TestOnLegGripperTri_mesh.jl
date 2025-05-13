module TestOnLegGripperSTL

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
    # data = matread(taskName * ".mat")
    data = matread("test/" * taskName * ".mat")
    
    # part_name = "elementy_trubky.txt"
    part_name = "test/elementy_trubky.txt"

    (X, IEN, rho) = MeshGrid.MeshInformations(data)

    # input data propertis (mesh, density)
    mesh = MeshGrid.Mesh(X, IEN, hex8_shape)
    rho = MeshGrid.ModiffElementalDensities(mesh, rho) # change density along the object

    ρₙ = MeshGrid.DenseInNodes(mesh, rho) # LSQ

    ## Face triangular mesh:
    tri_mesh = Rho2sdf.extractSurfaceTriangularMesh(mesh) 

    EN = Rho2sdf.NodePosition3D(tri_mesh)
    EN.x
    EN.y
    EN.z

    x₁ = [EN.x[:, 1663][1], EN.y[:, 1663][1], EN.z[:, 1663][1]]
    # x₂ = [EN.x[:, 1663][2], EN.y[:, 1663][2], EN.z[:, 1663][2]]
    # x₃  = [EN.x[:, 1663][3], EN.y[:, 1663][3], EN.z[:, 1663][3]]
    x₂  = [EN.x[:, 1663][3], EN.y[:, 1663][3], EN.z[:, 1663][3]]
    x₃ = [EN.x[:, 1663][2], EN.y[:, 1663][2], EN.z[:, 1663][2]]


function find_triangle_position(EN::Rho2sdf.MeshGrid.NodalCoordinatesInElement, vertices)
    # Preprocess input vertices to sort each vertex set for unordered comparison
    sorted_vertices = sort!(vec(vertices))

    for i in 1:size(EN.x, 2)  # Iterate over each triangle (column)
        # Extract and sort the vertices of the current triangle for comparison
        current_vertices = sort!(vec([EN.x[:, i] EN.y[:, i] EN.z[:, i]]))

        if all(current_vertices .== sorted_vertices)
            return i
        end
    end
    return -1  # Triangle not found
end

position = find_triangle_position(EN, [x₃ x₁ x₂])

#______________
# Output the result

    
    ## Grid:
    X_min, X_max = MeshGrid.getMesh_AABB(mesh.X) # vec, vec
    
    N = 18  #Number of divisions along the longest side (along some axis)
    sdf_grid = MeshGrid.Grid(X_min, X_max, N, 3) # cartesian grid
    
    ## SFD from triangular mesh:
    sdf_dists = SignedDistances.evalSignedDistancesOnTriangularMesh(mesh, sdf_grid) # Vector{Float64}
    print("done")

    ## Data export to VTK:
    Rho2sdf.exportStructuredPointsToVTK("STL-sdf-leg-" * taskName *".vtk", sdf_grid, sdf_dists, "distance")


end
