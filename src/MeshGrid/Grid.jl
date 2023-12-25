
struct Grid
    AABB_min::Vector{Float64}
    AABB_max::Vector{Float64}
    N::Vector{Int64}
    cell_size::Vector{Float64}
    ngp::Int64 # Number of Grid Points

    function Grid(
        AABB_min::Vector{Float64},
        AABB_max::Vector{Float64},
        N::Vector{Int64},
        margineCells::Int64 = 3,
    )

        cell_size = (AABB_max .- AABB_min) ./ N

        AABB_min = AABB_min .- margineCells * cell_size
        AABB_max = AABB_max .+ margineCells * cell_size

        AABB_size = AABB_max .- AABB_min

        N = N .+ 2 * margineCells
        ngp = prod(N .+ 1)

        return new(AABB_min, AABB_max, N, cell_size, ngp)

    end
end

mutable struct LinkedList # rozdělení pravidelné sítě na regiony
    grid::Grid
    head::Vector{Int64}
    next::Vector{Int64}
    N::Vector{Float64}

    function LinkedList(grid::Grid, X::Matrix{Float64})

        np = size(X, 2) # Number of Points
        N = grid.N
        AABB_min = grid.AABB_min
        AABB_max = grid.AABB_max

        head = -1 * ones(Int64, prod(N .+ 1))
        next = -1 * ones(Int64, np)

        I = floor.(N .* (X .- AABB_min) ./ (AABB_max .- AABB_min))
        I =
            Int.(
                I[3, :] .* (N[1] + 1) * (N[2] + 1) .+ I[2, :] .* (N[1] + 1) .+
                I[1, :] .+ 1,
            )


        for i = 1:np
            next[i] = head[I[i]]
            head[I[i]] = i
        end

        return new(grid, head, next, N)
    end
end


function getMesh_AABB(X::Matrix{Float64})
    X_min = vec(minimum(X, dims = 2))
    X_max = vec(maximum(X, dims = 2))
    return X_min, X_max
end

function generateGridPoints(grid::Grid)::Matrix{Float64}

    X = zeros(3, grid.ngp)
    a = 1
    for k = 0:grid.N[3]
        for j = 0:grid.N[2]
            for i = 0:grid.N[1]
                X[:, a] = grid.AABB_min .+ grid.cell_size .* [i, j, k]
                a += 1
            end
        end
    end
    return X
end


