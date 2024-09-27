
using Base.Threads
using Statistics
using LinearAlgebra
using ProgressMeter
using JLD2
using NLopt

using Rho2sdf.ShapeFunctions
using Rho2sdf.MeshGrid
using Rho2sdf.MyMarchingCubes
using Rho2sdf
using Rho2sdf.SignedDistances

include("sign_assignment_chat.jl")
include("sign_assignment_claude.jl")

# @load "Z_Chapadlo_xp.jld2" xp
@load "Z_Chapadlo_dist.jld2" dist
@load "Z_Chapadlo_mesh.jld2" mesh
@load "Z_Chapadlo_grid.jld2" grid
@load "Z_Chapadlo_points.jld2" points
@load "Z_Chapadlo_rho.jld2" ρₙ

# Function to create AABB from a set of points
# function compute_aabb(points::Matrix{Float64})
function compute_aabb(points)
  min_bounds = minimum(points, dims=2)
  max_bounds = maximum(points, dims=2)
  return min_bounds, max_bounds
end

# Function to check if a point is inside the AABB
# function is_point_inside_aabb(x::Vector{Float64}, min_bounds, max_bounds)
function is_point_inside_aabb(x, min_bounds, max_bounds)
  return all(min_bounds .<= x) && all(x .<= max_bounds)
end

function find_local_coordinates(
    sfce::Function,
    Xₑ,
    xₙ,
    # Xₑ::Matrix,
    # xₙ::Vector
)
    starting_points = [
        (0.0, 0.0, 0.0),
        (-0.5, -0.5, -0.5),
        (0.5, -0.5, -0.5),
        (0.5, 0.5, -0.5),
        (-0.5, 0.5, -0.5),
        (-0.5, -0.5, 0.5),
        (0.5, -0.5, 0.5),
        (0.5, 0.5, 0.5),
        (-0.5, 0.5, 0.5)
    ]

    function objective(ξ::Vector, grad::Vector)
        ξ₁, ξ₂, ξ₃ = ξ
        N = [
            -1/8 * (ξ₁ - 1) * (ξ₂ - 1) * (ξ₃ - 1),
             1/8 * (ξ₁ + 1) * (ξ₂ - 1) * (ξ₃ - 1),
            -1/8 * (ξ₁ + 1) * (ξ₂ + 1) * (ξ₃ - 1),
             1/8 * (ξ₁ - 1) * (ξ₂ + 1) * (ξ₃ - 1),
             1/8 * (ξ₁ - 1) * (ξ₂ - 1) * (ξ₃ + 1),
            -1/8 * (ξ₁ + 1) * (ξ₂ - 1) * (ξ₃ + 1),
             1/8 * (ξ₁ + 1) * (ξ₂ + 1) * (ξ₃ + 1),
            -1/8 * (ξ₁ - 1) * (ξ₂ + 1) * (ξ₃ + 1)
        ]
        x = Xₑ * N
        R = x - xₙ
        return sum(R[i]^2 for i in 1:3)
    end

    best_solution = nothing
    best_objective = Inf

    for start_point in starting_points
        opt = Opt(:LN_COBYLA, 3)
        opt.lower_bounds = [-10.2, -10.2, -10.2]
        opt.upper_bounds = [10.2, 10.2, 10.2]
        opt.xtol_rel = 1e-6
        opt.maxeval = 50
        opt.min_objective = objective

        (minf, minx, ret) = optimize(opt, collect(start_point))

        if minf < best_objective
            best_objective = minf
            best_solution = minx
        end
    end

    if best_solution === nothing
        return (false, [10.0, 10.0, 10.0])
    else
        return (true, best_solution)
    end
end

X = mesh.X   # vector of nodes positions
IEN = mesh.IEN # ID element -> ID nodes
INE = mesh.INE # ID node -> ID elements
ISN = mesh.ISN # connectivity face - edges
sfce = mesh.sfce # shape function handler
nsd = mesh.nsd # number of spacial dimensions
nel = mesh.nel # number of all elements
nes = mesh.nes # number of element segments (faces)
nsn = mesh.nsn # number of face nodes
# println("number of all elements: ", nel)

ngp = grid.ngp # number of nodes in grid

ρₜ = 0.5
# signs = -1 * ones(ngp)

# signs = compute_sign_chat(points, X, INE, IEN, ρₙ, ρₜ, ngp, sfce, nel)
signs = compute_sign_claude(points, X, INE, IEN, ρₙ, ρₜ, ngp, sfce, nel)


dist = abs.(dist) .* signs

taskName = "kontrola_znamenka"
Rho2sdf.exportStructuredPointsToVTK(taskName * "_sdf.vtk", grid, dist, "distance")

