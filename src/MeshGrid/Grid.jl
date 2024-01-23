
# N = [20, 20, 50]
# struct Grid
#     AABB_min::Vector{Float64}  # Minimum coordinates of the Axis-Aligned Bounding Box (AABB)
#     AABB_max::Vector{Float64}  # Maximum coordinates of the AABB
#     N::Vector{Int64}           # Number of divisions along each axis
#     cell_size::Vector{Float64} # Size of each cell in the grid
#     ngp::Int64                 # Total number of grid points

#     # Marginal cells are added around the AABB to provide a buffer zone.
#     function Grid(
#         AABB_min::Vector{Float64},
#         AABB_max::Vector{Float64},
#         N::Vector{Int64},
#         margineCells::Int64 = 3,
#     )

#         cell_size = (AABB_max .- AABB_min) ./ N

#         # Adjusting the AABB with marginal cells
#         AABB_min = AABB_min .- margineCells * cell_size
#         AABB_max = AABB_max .+ margineCells * cell_size

#         AABB_size = AABB_max .- AABB_min

#         # Recalculating grid dimensions and total grid points
#         N = N .+ 2 * margineCells # (2 = both sides)
#         ngp = prod(N .+ 1) # number of grid points

#         return new(AABB_min, AABB_max, N, cell_size, ngp)

#     end
# end

mutable struct Grid
    AABB_min::Vector{Float64}  # Minimum coordinates of the Axis-Aligned Bounding Box (AABB)
    AABB_max::Vector{Float64}  # Maximum coordinates of the AABB
    N::Vector{Int64}           # Number of cells along the longest side (along some axis)
    cell_size::Float64         # Size of each cell in the grid
    ngp::Int64                 # Total number of grid points

    # Marginal cells are added around the AABB to provide a buffer zone.
    function Grid(
        AABB_min::Vector{Float64},
        AABB_max::Vector{Float64},
        N_max::Int64,
        margineCells::Int64 = 3,
    )

        cell_size = maximum(AABB_max .- AABB_min) ./ N_max

        # Adjusting the AABB with marginal cells
        AABB_min = AABB_min .- margineCells * cell_size
        AABB_max = AABB_max .+ margineCells * cell_size

        # N = ceil.((AABB_max .- AABB_min)./cell_size)
        N = ceil.(Int64, (AABB_max .- AABB_min) / cell_size)
        
        AABB_max = AABB_min + N .* cell_size

        # Recalculating grid dimensions and total grid points
        # N = N .+ 2 * margineCells # (2 = both sides)
        ngp = prod(N .+ 1) # number of grid points

        return new(AABB_min, AABB_max, N, cell_size, ngp)

    end
end


# Struct for managing a linked list in the context of a grid. (Useful for spatial hashing or similar applications.)
mutable struct LinkedList # rozdělení pravidelné sítě na regiony
    grid::Grid           # The grid associated with the linked list
    head::Vector{Int64}  # Array representing the head of each list
    next::Vector{Int64}  # Array representing the next element in each list
    N::Vector{Float64}   # Number of cells along each axis of the grid

    # Constructor for LinkedList.
    # Maps points in a 3D space to the corresponding grid cells.
    function LinkedList(grid::Grid, X::Matrix{Float64})

        np = grid.ngp # Number of grid points
        N = grid.N
        AABB_min = grid.AABB_min
        AABB_max = grid.AABB_max

        head = -1 * ones(Int64, prod(N .+ 1))
        next = -1 * ones(Int64, np)

        # Calculate grid indices for each point
        I = floor.(N .* (X .- AABB_min) ./ (AABB_max .- AABB_min))
        Ia = Int.(I[3, :] .* (N[1] + 1) * (N[2] + 1) .+ I[2, :] .* (N[1] + 1) .+ I[1, :] .+ 1)

        # Construct the linked list
        for i = 1:np
            next[i] = head[Ia[i]]
            head[Ia[i]] = i
        end

        return new(grid, head, next, N)
    end
end


# Function to determine the Axis-Aligned Bounding Box (AABB) of a mesh.
function getMesh_AABB(X::Matrix{Float64})
    X_min = vec(minimum(X, dims = 2))
    X_max = vec(maximum(X, dims = 2))
    return X_min, X_max
end

# Generates grid points for a given grid structure.
# Returns a matrix where each column represents the coordinates of a grid point.
function generateGridPoints(grid::Grid)::Matrix{Float64}
    X = zeros(3, grid.ngp)
    a = 1
    if sum(grid.N) == 3
        X = reshape(grid.AABB_min .+ grid.cell_size .* [1, 1, 1], 3 ,1)
    else
        for k = 0:grid.N[3]
            for j = 0:grid.N[2]
                for i = 0:grid.N[1]
                    X[:, a] = grid.AABB_min .+ grid.cell_size .* [i, j, k]
                    a += 1
                end
            end
        end
    end
    return X
end

# Paralel version:
# function generateGridPoints(grid::Grid)::Matrix{Float64}
#     X = zeros(3, grid.ngp)
#     a = Atomic{Int}(1)
#     @threads for k = 0:grid.N[3]
#         for j = 0:grid.N[2]
#             for i = 0:grid.N[1]
#                 idx = atomic_add!(a, 1)
#                 X[:, idx] = grid.AABB_min .+ grid.cell_size .* [i, j, k]
#             end
#         end
#     end
#     return X
# end


function calculateMiniAABB_grid(
    Xt::Matrix{Float64}, # 
    δ::Float64,
    N::Vector{Int64},
    AABB_min::Vector{Float64},
    AABB_max::Vector{Float64},
    nsd::Int)
    
    Xt_min = minimum(Xt, dims = 2) .- δ # bottom left corner coord
    Xt_max = maximum(Xt, dims = 2) .+ δ # top righ corner coord

    # Mini AABB for triangle:
    I_min = floor.(N .* (Xt_min .- AABB_min) ./ (AABB_max .- AABB_min)) # Triangle location (index) within the grid
    I_max = floor.(N .* (Xt_max .- AABB_min) ./ (AABB_max .- AABB_min)) # Triangle location (index) within the grid


    for j = 1:nsd # am I inside AABB?
        if (I_min[j] < 0)
            I_min[j] = 0
        end
        if (I_max[j] >= N[j])
            I_max[j] = N[j] - 1
        end
    end

    Is = Iterators.product( # step range of mini AABB
        I_min[1]:I_max[1],
        I_min[2]:I_max[2],
        I_min[3]:I_max[3],
    )

    return Is
end

