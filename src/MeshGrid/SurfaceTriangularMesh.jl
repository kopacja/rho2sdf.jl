

function extractSurfaceTriangularMesh(mesh::Mesh)
    X = mesh.X
    IEN = mesh.IEN # ID element -> nodes
    INE = mesh.INE # ID node -> ID elements
    ISN = mesh.ISN # connectivity face - edges
    nel = mesh.nel # number of elements
    nes = mesh.nes # number of element segments (faces) 6
    nsn = mesh.nsn # number of segment nodes (kolik má stěna uzlů) 4

    X_new = Vector{Vector{Float64}}()
    push!(X_new, vec(X[:, 1]))

    IEN_new = Vector{Vector{Int64}}()

    for el = 1:nel
        # ρₑ = ρₙ[IEN[:, el]]
        commonEls = []
        for sg = 1:nes # 1:6 je face součástí pouze jednoho elementu?
            commonEls = INE[IEN[ISN[sg][1], el]]
            for a = 2:nsn # 2:4
                idx = findall(in(INE[IEN[ISN[sg][a], el]]), commonEls) # mají uzly (jedné stěny) společný pouze jeden element
                commonEls = commonEls[idx]
            end

            if (length(commonEls) == 1) # is a part of the outer boundary of the body

                Xs = X[:, IEN[ISN[sg], el]]
                Xc = mean(Xs, dims = 2) # center of the outside face

                for a = 1:nsn

                    IEN_el = zeros(Int64, 3)

                    x₁ = vec(Xs[:, a])
                    x₂ = vec(Xs[:, (a%nsn)+1])
                    x₃ = vec(Xc)

                    Xt = [x₁, x₂, x₃]
                    # Xt = reduce(hcat, Xt)

                    for i = 1:3
                        a = findfirst(x -> norm(x - Xt[i]) < 1.0e-5, X_new)
                        if (a === nothing)
                            push!(X_new, Xt[i])
                            IEN_el[i] = length(X_new)
                        else
                            IEN_el[i] = a[1]
                        end
                    end
                    push!(IEN_new, IEN_el)
                end # a = 1:nsn
            end # if (length(commonEls) == 1)
        end # sg
    end # el

    return TriangularMesh(X_new, IEN_new)
end

function find_triangle_position(EN::NodalCoordinatesInElement, vertices::Matrix)
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
