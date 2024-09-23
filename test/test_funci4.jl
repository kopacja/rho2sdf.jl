using Rho2sdf
using LinearAlgebra
using Rho2sdf.ShapeFunctions
using Rho2sdf.MeshGrid
using BenchmarkTools


# Newton-Raphson method
Xₑ = [
    [-1.0, -1.0, -1.0],
    [1.0, -1.0, -1.0],
    [1.0, 1.0, -1.0],
    [-1.0, 1.0, -1.0],
    [-1.0, -1.0, 1.0],
    [1.0, -1.0, 1.0],
    [1.0, 1.0, 1.0],
    [-1.0, 1.0, 1.0],
  ]
  IEN = [[1, 2, 3, 4, 5, 6, 7, 8]]

  ## Generate FEM mesh structure:
  mesh = MeshGrid.Mesh(X, IEN, C3D8_SFaD)

sfce = C3D8_SFaD

function find_local_coordinates(
    sfce::Function,
    Xₑ::Matrix, # Xe
    xₙ::Vector, # x
  )
  
    ξ = [0.0, 0.0, 0.0]  # Initial guess
    # ξ = [0.4, 0.4, 0.4]  # Initial guess
    tolerance = 1e-8
    max_iterations = 40
  
    for iter in 1:max_iterations
      N, dN_dξ, _, _ = sfce(ξ) # Shape functions and their derivatives
  
      # Compute global coordinates from shape functions
      x = Xₑ * N
      R = x - xₙ
  
      # Check convergence
      if norm(R) < tolerance
        return ξ
      end
      println("norm(R): ", norm(R))
  
      J = Xₑ * dN_dξ
  
      # Update ξ
      # Δξ = -J \ R
      Δξ = J \ -R
      ξ += 0.5 * Δξ
    end
    println("NR method for local coords did not converge")
    println("Xₑ: ", Xₑ)
    println("xₙ: ", xₙ)
  
    return ξ
    # return ([10.0, 10.0, 10.0])
  end
  
  @time local_coords = find_local_coordinates(sfce, Xₑ, xₙ)
  # 0.008642 seconds (20.92 k allocations: 2.618 MiB)

max_local = maximum(abs.(local_coords))
max_local < 1.0001


Xₑ = [29.56773 33.82964 36.74095 32.72983 29.77495 34.01323 36.71943 32.76317; 37.76198 38.29159 40.52183 42.24267 37.15282 37.73094 40.01894 41.611; 57.5 57.5 57.5 57.5 63.0 63.0 63.0 63.0]
xₙ = [36.416666666666664, 38.0, 59.5]

# Xₑ = [29.56773 33.82964 36.74095 32.72983 29.77495 34.01323 36.71943 32.76317; 37.76198 38.29159 40.52183 42.24267 37.15282 37.73094 40.01894 41.611; 57.5 57.5 57.5 57.5 63.0 63.0 63.0 63.0]
# xₙ = [36.416666666666664, 38.0, 61.08333333333333]

# Xₑ = [29.81633 31.84215 29.30687 25.43636 29.9796 31.92058 29.36322 25.52554; 46.92605 50.16105 51.09008 50.05788 46.26536 49.29867 50.13349 49.05251; 57.5 57.5 57.5 57.5 63.0 63.0 63.0 63.0]
# xₙ = [30.08333333333333, 50.666666666666664, 62.666666666666664]

# Xₑ = [29.9796 35.81804 34.55206 31.92058 30.14286 35.91945 34.66492 31.99901; 46.26536 45.72567 48.5848 49.29867 45.60467 45.14916 47.83677 48.43628; 63.0 63.0 63.0 63.0 68.5 68.5 68.5 68.5]
# xₙ = [34.83333333333333, 49.08333333333333, 65.83333333333333]

# Xₑ = [29.9796 35.81804 34.55206 31.92058 30.14286 35.91945 34.66492 31.99901; 46.26536 45.72567 48.5848 49.29867 45.60467 45.14916 47.83677 48.43628; 63.0 63.0 63.0 63.0 68.5 68.5 68.5 68.5]
# xₙ = [34.83333333333333, 49.08333333333333, 67.41666666666666]

# Xₑ = [30.14286 31.99901 29.41957 25.61472 30.30613 32.07744 29.47592 25.7039; 45.60467 48.43628 49.17691 48.04714 44.94399 47.57389 48.22033 47.04177; 68.5 68.5 68.5 68.5 74.0 74.0 74.0 74.0]
# xₙ = [30.08333333333333, 49.08333333333333, 72.16666666666666]

# Xₑ = [30.14286 31.99901 29.41957 25.61472 30.30613 32.07744 29.47592 25.7039; 45.60467 48.43628 49.17691 48.04714 44.94399 47.57389 48.22033 47.04177; 68.5 68.5 68.5 68.5 74.0 74.0 74.0 74.0]
# xₙ = [30.08333333333333, 49.08333333333333, 73.75]

using JuMP, Ipopt

function find_local_coordinates(
  # sfce::Function,
  Xₑ::Matrix, # Xe
  xₙ::Vector  # x
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

  best_solution = nothing
  best_objective = Inf

  for (ξ₁_start, ξ₂_start, ξ₃_start) in starting_points
    model = Model(Ipopt.Optimizer)
    set_silent(model)

    set_optimizer_attribute(model, "tol", 1e-6)
    set_optimizer_attribute(model, "max_iter", 50)
    set_optimizer_attribute(model, "acceptable_tol", 1e-6)

    @variable(model, ξ₁, lower_bound = -10.2, upper_bound = 10.2, start = ξ₁_start)
    @variable(model, ξ₂, lower_bound = -10.2, upper_bound = 10.2, start = ξ₂_start)
    @variable(model, ξ₃, lower_bound = -10.2, upper_bound = 10.2, start = ξ₃_start)


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

    x = Xₑ * N
    R = x - xₙ

    @objective(model, Min, sum(R[i]^2 for i in 1:3))
    optimize!(model)

    solution = [value(ξ₁), value(ξ₂), value(ξ₃)]
    obj_value = objective_value(model)

    if obj_value < best_objective
      best_objective = obj_value
      best_solution = solution
    end
  end

  if best_solution === nothing
    println("Optimization did not converge to a solution")
    return (false, [10.0, 10.0, 10.0])
  else
    return (true, best_solution)
  end
end

@time (computed, best_solution) = find_local_coordinates(Xₑ::Matrix, xₙ::Vector)

using Optim

function find_local_coordinates_Optim(
    Xₑ::Matrix,
    xₙ::Vector
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

    function objective(ξ)
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
        return sum(R.^2)
    end

    best_solution = nothing
    best_objective = Inf

    for start_point in starting_points
        result = optimize(objective, collect(start_point), LBFGS(),
                          Optim.Options(iterations=50, g_tol=1e-6))
        
        if Optim.minimum(result) < best_objective
            best_objective = Optim.minimum(result)
            best_solution = Optim.minimizer(result)
        end
    end

    if best_solution === nothing
        println("Optimalizace nekonvergovala k řešení")
        return (false, [10.0, 10.0, 10.0])
    else
        return (true, best_solution)
    end
end

@time (computed, best_solution) = find_local_coordinates_Optim(Xₑ, xₙ)


using NLsolve
using LinearAlgebra

function find_local_coordinates_NLsolve1(
    Xₑ::Matrix,
    xₙ::Vector
)
    starting_points = [
        [0.0, 0.0, 0.0],
        [-0.5, -0.5, -0.5],
        [0.5, -0.5, -0.5],
        [0.5, 0.5, -0.5],
        [-0.5, 0.5, -0.5],
        [-0.5, -0.5, 0.5],
        [0.5, -0.5, 0.5],
        [0.5, 0.5, 0.5],
        [-0.5, 0.5, 0.5]
    ]

    function f!(F, ξ)
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
        F .= x - xₙ
    end

    best_solution = nothing
    best_residual_norm = Inf

    for (i, start_point) in enumerate(starting_points)
        result = nlsolve(f!, start_point, method = :trust_region, iterations = 1000, ftol = 1e-8)
        
        println("Start point $i: $(start_point)")
        println("Converged: $(converged(result))")
        # println("Residual norm: $(norm(result.residual))")
        println("Solution: $(result.zero)")
        println()

        if converged(result) && norm(result.residual) < best_residual_norm
            best_residual_norm = norm(result.residual)
            best_solution = result.zero
        end
    end

    if best_solution === nothing
        println("NLsolve nekonvergovalo k řešení")
        return (false, [10.0, 10.0, 10.0])
    else
        println("Nejlepší řešení:")
        println("Residual norm: $best_residual_norm")
        println("Solution: $best_solution")
        return (true, best_solution)
    end
end

# Předpokládám, že Xₑ a xₙ jsou definovány před voláním funkce
@time (computed, best_solution) = find_local_coordinates_NLsolve1(Xₑ, xₙ)