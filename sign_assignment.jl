
using Base.Threads
using Statistics
using LinearAlgebra
using ProgressMeter
using JLD2
using NLopt
# using Atomics

using Rho2sdf.ShapeFunctions
using Rho2sdf.MeshGrid
using Rho2sdf
using Rho2sdf.SignedDistances

include("sign_assignment_original.jl")

@load "Z_Chapadlo_final_dist.jld2" dist
@load "Z_Chapadlo_final_mesh.jld2" mesh
@load "Z_Chapadlo_final_grid.jld2" grid
@load "Z_Chapadlo_final_points.jld2" points
@load "Z_Chapadlo_final_rho.jld2" ρₙ

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
  xₙ
)
  starting_points::Vector{Tuple{Float64,Float64,Float64}} = [
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

    # Compute N(ξ)
    N = [
      -1 / 8 * (ξ₁ - 1) * (ξ₂ - 1) * (ξ₃ - 1),
      1 / 8 * (ξ₁ + 1) * (ξ₂ - 1) * (ξ₃ - 1),
      -1 / 8 * (ξ₁ + 1) * (ξ₂ + 1) * (ξ₃ - 1),
      1 / 8 * (ξ₁ - 1) * (ξ₂ + 1) * (ξ₃ - 1),
      1 / 8 * (ξ₁ - 1) * (ξ₂ - 1) * (ξ₃ + 1),
      -1 / 8 * (ξ₁ + 1) * (ξ₂ - 1) * (ξ₃ + 1),
      1 / 8 * (ξ₁ + 1) * (ξ₂ + 1) * (ξ₃ + 1),
      -1 / 8 * (ξ₁ - 1) * (ξ₂ + 1) * (ξ₃ + 1)
    ]

    d¹N_dξ¹ = zeros(Float64, 8, 3)
    d¹N_dξ¹[1, :] = [
      -0.125 * (ξ₂ - 1) * (ξ₃ - 1)
      -0.125 * (ξ₁ - 1) * (ξ₃ - 1)
      -0.125 * (ξ₁ - 1) * (ξ₂ - 1)
    ]

    d¹N_dξ¹[2, :] = [
      0.125 * (ξ₂ - 1) * (ξ₃ - 1)
      0.125 * (ξ₁ + 1) * (ξ₃ - 1)
      0.125 * (ξ₁ + 1) * (ξ₂ - 1)
    ]

    d¹N_dξ¹[3, :] = [
      -0.125 * (ξ₂ + 1) * (ξ₃ - 1)
      -0.125 * (ξ₁ + 1) * (ξ₃ - 1)
      -0.125 * (ξ₁ + 1) * (ξ₂ + 1)
    ]

    d¹N_dξ¹[4, :] = [
      0.125 * (ξ₂ + 1) * (ξ₃ - 1)
      0.125 * (ξ₁ - 1) * (ξ₃ - 1)
      0.125 * (ξ₁ - 1) * (ξ₂ + 1)
    ]

    d¹N_dξ¹[5, :] = [
      0.125 * (ξ₂ - 1) * (ξ₃ + 1)
      0.125 * (ξ₁ - 1) * (ξ₃ + 1)
      0.125 * (ξ₁ - 1) * (ξ₂ - 1)
    ]

    d¹N_dξ¹[6, :] = [
      -0.125 * (ξ₂ - 1) * (ξ₃ + 1)
      -0.125 * (ξ₁ + 1) * (ξ₃ + 1)
      -0.125 * (ξ₁ + 1) * (ξ₂ - 1)
    ]

    d¹N_dξ¹[7, :] = [
      0.125 * (ξ₂ + 1) * (ξ₃ + 1)
      0.125 * (ξ₁ + 1) * (ξ₃ + 1)
      0.125 * (ξ₁ + 1) * (ξ₂ + 1)
    ]

    d¹N_dξ¹[8, :] = [
      -0.125 * (ξ₂ + 1) * (ξ₃ + 1)
      -0.125 * (ξ₁ - 1) * (ξ₃ + 1)
      -0.125 * (ξ₁ - 1) * (ξ₂ + 1)
    ]

    # Compute x(ξ)
    x = Xₑ * N  # x is a 3-element vector

    # Compute residual R = x - xₙ
    R = x - xₙ  # R is a 3-element vector

    # Compute objective function value f
    f = dot(R, R)  # Equivalent to sum(R[i]^2 for i in 1:3)

    if length(grad) > 0
      # Compute the derivatives of N with respect to ξ
      dN_dξ = zeros(Float64, 8, 3)  # 8 shape functions x 3 variables

      # Assign your computed derivatives to dN_dξ
      dN_dξ .= d¹N_dξ¹  # Ensure dN_dξ is correctly populated

      # Compute Jacobian J = Xₑ * dN_dξ
      J = Xₑ * dN_dξ  # J is 3x3

      # Compute gradient grad = 2 * Jᵗ * R
      grad .= 2 * (J' * R)
    end

    return f
  end

  best_solution = nothing
  best_objective = Inf

  for start_point in starting_points
    #   opt = Opt(:LN_COBYLA, 3)
    opt = Opt(:LD_LBFGS, 3)
    opt.lower_bounds = [-5.0, -5.0, -5.0]
    opt.upper_bounds = [5.0, 5.0, 5.0]
    opt.xtol_rel = 1e-6
    opt.maxeval = 500
    opt.min_objective = objective
    opt.maxtime = 1.0

    (minf, minx, ret) = NLopt.optimize(opt, collect(start_point))

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

signs = @time compute_sign_original(points, X, INE, IEN, ρₙ, ρₜ, ngp, sfce, nel)

dist = abs.(dist) .* signs

taskName = "kontrola_znamenka-chapadlo"
Rho2sdf.exportStructuredPointsToVTK(taskName * "_sdf.vtk", grid, dist, "distance")

