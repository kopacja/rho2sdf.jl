using NLopt
using LinearAlgebra
using StaticArrays

function compute_coords_on_iso(
    x::Vector,
    ρₜ::Float64,
    Xₑ::Matrix,
    ρₑ::Vector
)
    # Define the dimension of our optimization problem
    n_vars = 3  # ξ₁, ξ₂, ξ₃

    # Create the optimizer with SLSQP algorithm
    opt = Opt(:LD_SLSQP, n_vars)
    
    # Set optimization parameters
    ftol_rel!(opt, 1e-8)    # Relative tolerance for objective function
    ftol_abs!(opt, 1e-8)    # Absolute tolerance for objective function
    xtol_rel!(opt, 1e-8)    # Relative tolerance for optimization variables
    maxeval!(opt, 1000)     # Maximum number of evaluations
    
    # Set bound constraints for ξᵢ ∈ [-1, 1]
    lower_bounds!(opt, [-1.0, -1.0, -1.0])
    upper_bounds!(opt, [1.0, 1.0, 1.0])

    # Helper function to compute shape functions N8
    function compute_N8(ξ)
        ξ₁, ξ₂, ξ₃ = ξ
        return [
            -1/8 * (ξ₁ - 1) * (ξ₂ - 1) * (ξ₃ - 1),
             1/8 * (ξ₁ + 1) * (ξ₂ - 1) * (ξ₃ - 1),
            -1/8 * (ξ₁ + 1) * (ξ₂ + 1) * (ξ₃ - 1),
             1/8 * (ξ₁ - 1) * (ξ₂ + 1) * (ξ₃ - 1),
             1/8 * (ξ₁ - 1) * (ξ₂ - 1) * (ξ₃ + 1),
            -1/8 * (ξ₁ + 1) * (ξ₂ - 1) * (ξ₃ + 1),
             1/8 * (ξ₁ + 1) * (ξ₂ + 1) * (ξ₃ + 1),
            -1/8 * (ξ₁ - 1) * (ξ₂ + 1) * (ξ₃ + 1)
        ]
    end

    # Objective function
    function objective(ξ::Vector, grad::Vector)
        N8 = compute_N8(ξ)
        
        # If gradient is required
        if length(grad) > 0
            # Compute gradients analytically or use finite differences
            # This is a simplified gradient computation - you might want to implement
            # the analytical gradient for better performance
            ϵ = 1e-8
            for i in 1:length(ξ)
                ξ_plus = copy(ξ)
                ξ_plus[i] += ϵ
                N8_plus = compute_N8(ξ_plus)
                
                obj_plus = sum((x[j] - sum(Xₑ[j,k] * N8_plus[k] for k in 1:8))^2 
                             for j in 1:length(x))
                obj_minus = sum((x[j] - sum(Xₑ[j,k] * N8[k] for k in 1:8))^2 
                              for j in 1:length(x))
                
                grad[i] = (obj_plus - obj_minus) / ϵ
            end
        end
        
        # Return objective value
        return sum((x[i] - sum(Xₑ[i,k] * N8[k] for k in 1:8))^2 
                  for i in 1:length(x))
    end

    # Constraint function (density interpolation)
    function constraint(ξ::Vector, grad::Vector, ρₜ::Float64, ρₑ::Vector)
        N8 = compute_N8(ξ)
        
        # If gradient is required
        if length(grad) > 0
            # Compute constraint gradients
            ϵ = 1e-8
            for i in 1:length(ξ)
                ξ_plus = copy(ξ)
                ξ_plus[i] += ϵ
                N8_plus = compute_N8(ξ_plus)
                
                con_plus = sum(ρₑ[k] * N8_plus[k] for k in 1:8) - ρₜ
                con_minus = sum(ρₑ[k] * N8[k] for k in 1:8) - ρₜ
                
                grad[i] = (con_plus - con_minus) / ϵ
            end
        end
        
        # Return constraint value
        return sum(ρₑ[k] * N8[k] for k in 1:8) - ρₜ
    end

    # Set objective function
    function obj_wrapper(ξ::Vector, grad::Vector)
        return objective(ξ, grad)
    end
    min_objective!(opt, obj_wrapper)

    # Add equality constraint
    function constraint_wrapper(ξ::Vector, grad::Vector)
        return constraint(ξ, grad, ρₜ, ρₑ)
    end
    equality_constraint!(opt, constraint_wrapper, 1e-8)

    # Initial guess (you might want to try multiple starting points)
    initial_ξ = [0.0, 0.0, 0.0]

    # Perform optimization
    (minf, minx, ret) = optimize(opt, initial_ξ)

    # Check convergence
    if ret == :FAILURE
        @warn "Optimization failed to converge"
    end

    return minx
end


"""
Old version based on Ipopt&JuMP that works fine but it is much slower:
"""

# function compute_coords_on_iso(
#   x::Vector,
#   ρₜ::Float64,
#   Xₑ::Matrix,
#   ρₑ::Vector,
#   # starting_points::Vector{Tuple{Float64, Float64, Float64}}
# )
#   starting_points = [
#     (0.0, 0.0, 0.0),
#     # (-0.5, -0.5, -0.5),
#     # (0.5, -0.5, -0.5),
#     # (0.5, 0.5, -0.5),
#     # (-0.5, 0.5, -0.5),
#     # (-0.5, -0.5, 0.5),
#     # (0.5, -0.5, 0.5),
#     # (0.5, 0.5, 0.5),
#     # (-0.5, 0.5, 0.5),
#   ]
#
#   best_solution = nothing
#   best_objective = Inf
#
#   for (ξ₁_start, ξ₂_start, ξ₃_start) in starting_points
#     model = Model(Ipopt.Optimizer)
#     set_silent(model)
#
#     set_optimizer_attribute(model, "tol", 1e-6)
#     set_optimizer_attribute(model, "max_iter", 50)
#     set_optimizer_attribute(model, "acceptable_tol", 1e-6)
#
#     @variable(model, ξ₁, lower_bound = -1.0, upper_bound = 1.0, start = ξ₁_start)
#     @variable(model, ξ₂, lower_bound = -1.0, upper_bound = 1.0, start = ξ₂_start)
#     @variable(model, ξ₃, lower_bound = -1.0, upper_bound = 1.0, start = ξ₃_start)
#
#     N8 = [
#       -1 / 8 * (ξ₁ - 1) * (ξ₂ - 1) * (ξ₃ - 1),
#       1 / 8 * (ξ₁ + 1) * (ξ₂ - 1) * (ξ₃ - 1),
#       -1 / 8 * (ξ₁ + 1) * (ξ₂ + 1) * (ξ₃ - 1),
#       1 / 8 * (ξ₁ - 1) * (ξ₂ + 1) * (ξ₃ - 1),
#       1 / 8 * (ξ₁ - 1) * (ξ₂ - 1) * (ξ₃ + 1),
#       -1 / 8 * (ξ₁ + 1) * (ξ₂ - 1) * (ξ₃ + 1),
#       1 / 8 * (ξ₁ + 1) * (ξ₂ + 1) * (ξ₃ + 1),
#       -1 / 8 * (ξ₁ - 1) * (ξ₂ + 1) * (ξ₃ + 1)
#     ]
#     @NLobjective(model, Min, sum((x[i] - sum(Xₑ[i, k] * N8[k] for k in 1:length(N8)))^2 for i in 1:length(x)))
#     @NLconstraint(model, sum(ρₑ[k] * N8[k] for k in 1:length(N8)) == ρₜ)
#     JuMP.optimize!(model)
#
#     current_objective = objective_value(model)
#     current_solution = value.([ξ₁, ξ₂, ξ₃])
#
#     if current_objective < best_objective
#       best_solution = current_solution
#       best_objective = current_objective
#     end
#   end
#
#   return best_solution
# end

