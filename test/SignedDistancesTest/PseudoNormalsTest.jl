module PseudoNormalsTest


using Test
using LinearAlgebra
using Statistics
using Rho2sdf
using Rho2sdf.ShapeFunctions
using Rho2sdf.PrimitiveGeometries
using Rho2sdf.MeshGrid
using Rho2sdf.SignedDistances
using MAT

# Data:
# taskName = "chapadlo"  # Define the task name, used for data file naming
# data = matread(taskName * ".mat")  # Read the Matlab data file with the specified task name
# (X, IEN, rho) = MeshGrid.MeshInformations(data)  # Extract mesh information from the data
(X, IEN, rho) = PrimitiveGeometries.selectPrimitiveGeometry("sphere", 6)
mesh = MeshGrid.Mesh(X, IEN, rho, hex8_shape)
mesh = Rho2sdf.extractSurfaceTriangularMesh(mesh) 


function computePseudoNormals1(mesh::TriangularMesh)
    verticesPseudoNormals = [zeros(mesh.nsd) for _ = 1:mesh.nnp]
    edgesPseudoNormals = [[zeros(mesh.nsd) for _ = 1:3] for _ = 1:mesh.nel]

    for el = 1:mesh.nel
        As = mesh.IEN[:, el]
        Xt = mesh.X[:, As]
        n = cross(Xt[:, 2] - Xt[:, 1], Xt[:, 3] - Xt[:, 1])
        n = n / norm(n)

        for i = 1:3
            j = mod(i, 3) + 1
            el_adj = setdiff(intersect(mesh.INE[As[i]], mesh.INE[As[j]]), el)
            As_adj = mesh.IEN[:, el_adj]
            Xt_adj = mesh.X[:, As_adj]
            n_adj =
                cross(Xt_adj[:, 2] - Xt_adj[:, 1], Xt_adj[:, 3] - Xt_adj[:, 1])
            n_adj = π * n_adj / norm(n_adj)
            edgesPseudoNormals[el][i] = n + n_adj
        end
    end

    for A = 1:mesh.nnp
        pseudo_n = zeros(3)
        for k = 1:length(mesh.INE[A])
            el = mesh.INE[A][k]
            As = mesh.IEN[:, el]
            Xt = mesh.X[:, As]
            i = findfirst(x -> x == A, As)
            B = setdiff([1, 2, 3], i)
            e₁ = Xt[:, B[1]] - Xt[:, i]
            e₂ = Xt[:, B[2]] - Xt[:, i]
            α = acos(dot(e₁ / norm(e₁), e₂ / norm(e₂)))
            n = cross(Xt[:, 2] - Xt[:, 1], Xt[:, 3] - Xt[:, 1])
            n = n / norm(n)
            pseudo_n += α * n
        end
        verticesPseudoNormals[A] = pseudo_n
    end

    return verticesPseudoNormals, edgesPseudoNormals
end

verticesPseudoNormals, edgesPseudoNormals = SignedDistances.computePseudoNormals(mesh)
verticesPseudoNormals1, edgesPseudoNormals1 = computePseudoNormals1(mesh)

@test verticesPseudoNormals == verticesPseudoNormals1
@test norm(edgesPseudoNormals .- edgesPseudoNormals1) < 10e-6

end
