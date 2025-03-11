# Normalized normal vector of the triangle.
function computeTriangleNormal(Xt)
    n = cross(Xt[:, 2] - Xt[:, 1], Xt[:, 3] - Xt[:, 1])
    return n / norm(n)
end

# Calculates the pseudo normals for edges in the mesh.
# This function iterates over each element in the mesh and computes the pseudo normals
# for its edges based on adjacent elements.
function computeEdgePseudoNormals(mesh)
    edgesPseudoNormals = [[zeros(mesh.nsd) for _ = 1:3] for _ = 1:mesh.nel]

    # for el = 1:mesh.nel # cycle over the number of all elements
    @threads for el = 1:mesh.nel # cycle over the number of all elements
        # Get the vertices of the current element
        As = mesh.IEN[:, el] # ID of nodes for each triangle
        Xt = mesh.X[:, As]   # cooradinates of nodes
        n = computeTriangleNormal(Xt)

        for i = 1:3
            j = mod(i, 3) + 1
            # Find adjacent element
            el_adj = setdiff(intersect(mesh.INE[As[i]], mesh.INE[As[j]]), el)
            As_adj = mesh.IEN[:, el_adj]
            Xt_adj = mesh.X[:, As_adj]
            n_adj = computeTriangleNormal(Xt_adj) * π
            edgesPseudoNormals[el][i] = n + n_adj
        end
    end

    # Array of pseudo normal vectors for each edge of each element in the mesh.
    return edgesPseudoNormals
end

# Computes the pseudo normals for vertices in the mesh. Iterates over each vertex and calculates the pseudo normal based on the adjacent elements.
function computeVertexPseudoNormals(mesh)
    verticesPseudoNormals = [zeros(mesh.nsd) for _ = 1:mesh.nnp]

    # for A = 1:mesh.nnp # cycle over the number of all nodes
    @threads for A = 1:mesh.nnp # cycle over the number of all nodes
        pseudo_n = zeros(mesh.nsd)
        for el in mesh.INE[A]
            As = mesh.IEN[:, el]
            Xt = mesh.X[:, As]
            i = findfirst(==(A), As)
            B = setdiff([1, 2, 3], i)
            e₁, e₂ = Xt[:, B[1]] - Xt[:, i], Xt[:, B[2]] - Xt[:, i]
            α = acos(dot(e₁ / norm(e₁), e₂ / norm(e₂)))
            n = computeTriangleNormal(Xt)
            pseudo_n += α * n
        end
        verticesPseudoNormals[A] = pseudo_n
    end

    #  Array of pseudo normal vectors for each vertex in the mesh:
    return verticesPseudoNormals
end

## Main Function: computePseudoNormals ##
# Computation of pseudo normals for both edges and vertices in a mesh.
function computePseudoNormals(mesh::TriangularMesh)
    edgesPseudoNormals = computeEdgePseudoNormals(mesh)
    verticesPseudoNormals = computeVertexPseudoNormals(mesh)
    return verticesPseudoNormals, edgesPseudoNormals
end

