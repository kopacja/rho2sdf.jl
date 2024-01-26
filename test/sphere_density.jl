using Rho2sdf
using LinearAlgebra
using MAT

a = 1.
N = 10  # Number of cells along the longest side (along some axis)

X_min = [-a, -a, -a]
X_max = [a, a, a]

sdf_grid = Rho2sdf.Grid(X_min, X_max, N, 0) # cartesian grid

radius = 0.9 * a / sdf_grid.cell_size

function generate_sphere_density(N::Int64, radius::Float64, cell_size::Float64)
     
    weights = [0.5688888888888889, 0.4786286704993665, 0.4786286704993665, 0.2369268850561891, 0.2369268850561891]

    points = [0.0000000000000000, -0.5384693101056831, 0.5384693101056831, -0.9061798459386640, 0.9061798459386640]

    points_3D  = Iterators.product(points, points, points)
    weights_3D = Iterators.product(weights, weights, weights)

    Is = Iterators.product( # step range of mini AABB
        0:N-1,
        0:N-1,
        0:N-1,
    )

    sdf_density = zeros((N)^3)

    for I in Is
        i = Int(I[3] * N^2 + I[2] * N + I[1] + 1)
        println("i: ", i)
        X = collect(I) # I is tuple and collect change it to vector
        Xₑ = [[X[1],   X[2],   X[3]],
              [X[1]+1, X[2],   X[3]],
              [X[1]+1, X[2]+1, X[3]],
              [X[1],   X[2]+1, X[3]],
              [X[1],   X[2],   X[3]+1],
              [X[1]+1, X[2],   X[3]+1],
              [X[1]+1, X[2]+1, X[3]+1],
              [X[1],   X[2]+1, X[3]+1]]
        Xₑ = reduce(hcat, Xₑ)

        for (weight, point) ∈ zip(weights_3D, points_3D)
            Ξ = collect(point)
            H, d¹N_dξ¹, d²N_dξ², d³N_dξ³ = Rho2sdf.sfce(Ξ)
            x_g = Xₑ * H
            if(norm(x_g - [N, N, N]./2) < radius)
                sdf_density[i]  += weight[1]*weight[2]*weight[3] / 8
            end
        end
    end

    return sdf_density
end

rho = generate_sphere_density(N, radius, sdf_grid.cell_size)

Rho2sdf.exportStructuredPointsToVTK("sphere.vtk", sdf_grid, rho, "density")

function build_X(X_min::Vector{Float64}, X_max::Vector{Float64}, N::Int64, cell_size::Float64)
    X = [fill(0.0, 3) for _ in 1:(N+1)^3 ]
    i = 1
    for z = X_min[3]:cell_size:X_max[3]
        for y = X_min[2]:cell_size:X_max[2]
            for x = X_min[1]:cell_size:X_max[1]
                X[i] = [x,y,z]
                i = i+1
            end
        end
    end
    X = reduce(hcat, X)
    return X
end

X = build_X(X_min, X_max, N, sdf_grid.cell_size)

function build_IEN(N::Int64) # N is a number of vertex in row not number of cells
    IEN = [fill(0, 8) for _ in 1:(N-1)^3 ]
    m = 1
    for k = 1:N-1
        for j = 1:N-1
            for i = 1:N-1
                IEN[m] = [
                    Int((k-1) * N^2 + (j-1) * N + i-1 ),
                    Int((k-1) * N^2 + (j-1) * N + i ),
                    Int((k-1) * N^2 + j * N + i ),
                    Int((k-1) * N^2 + j * N + i-1 ),
                    Int(k * N^2 + (j-1) * N + i-1),
                    Int(k * N^2 + (j-1) * N + i ),
                    Int(k * N^2 + j * N + i ),
                    Int(k * N^2 + j * N + i-1 ),
                ]
                m = m+1
            end
        end
    end
    IEN = reduce(hcat, IEN)
    return IEN
end

IEN = build_IEN(N+1) # because N is now number of nodes no elements

msh = Dict(
    "X" => X,
    "IEN" => IEN
)

file = matopen("test/sphere.mat", "w")

write(file, "msh", msh)
write(file, "rho", rho)

close(file)