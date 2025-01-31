
function compute_coords_on_iso(
  x::Vector,
  ρₜ::Float64,
  Xₑ::Matrix,
  ρₑ::Vector
)

  opt = Opt(:LD_SLSQP, 3)
  ftol_rel!(opt, 1e-5)
  ftol_abs!(opt, 1e-5)
  xtol_rel!(opt, 1e-5)
  maxeval!(opt, 1000)

  lower_bounds!(opt, SVector{3}(-1.0, -1.0, -1.0))
  upper_bounds!(opt, SVector{3}(1.0, 1.0, 1.0))

  # Přealokace všech pomocných vektorů
  N = MVector{8,Float64}(undef)
  dN = MMatrix{8,3,Float64}(undef)
  temp_sums = MVector{3,Float64}(undef)

  # Closure pro sdílení přealokovaných vektorů
  let N = N, dN = dN, temp_sums = temp_sums
    function objective(ξ::Vector, grad::Vector)
      compute_shape_functions!(N, dN, ξ[1], ξ[2], ξ[3])

      # Předpočítání projekcí pro všechny dimenze najednou
      @inbounds for i in 1:3
        temp_sums[i] = sum(Xₑ[i, k] * N[k] for k in 1:8)
      end

      # Výpočet residuálů s využitím předpočítaných sum
      rx = x[1] - temp_sums[1]
      ry = x[2] - temp_sums[2]
      rz = x[3] - temp_sums[3]

      if !isempty(grad)
        @inbounds for j in 1:3
          # Předpočítání derivací pro všechny dimenze
          @inbounds for i in 1:3
            temp_sums[i] = sum(Xₑ[i, k] * dN[k, j] for k in 1:8)
          end
          grad[j] = -2.0 * (rx * temp_sums[1] + ry * temp_sums[2] + rz * temp_sums[3])
        end
      end

      return rx * rx + ry * ry + rz * rz
    end

    function constraint(ξ::Vector, grad::Vector)
      compute_shape_functions!(N, dN, ξ[1], ξ[2], ξ[3])

      if !isempty(grad)
        @inbounds for j in 1:3
          grad[j] = sum(ρₑ[k] * dN[k, j] for k in 1:8)
        end
      end

      return sum(ρₑ[k] * N[k] for k in 1:8) - ρₜ
    end

    min_objective!(opt, objective)
    equality_constraint!(opt, constraint, 1e-8)
  end

  # Optimalizace s přednastavenou přesností
  initial_ξ = zeros(SVector{3,Float64})
  minf, minx, ret = optimize(opt, initial_ξ)

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

