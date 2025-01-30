# using StaticArrays
# using NLopt
# using LinearAlgebra

"""
Hlavní funkce pro výpočet lokálních souřadnic kombinující Newton-Raphsonovu metodu
s NLopt optimalizací jako záložním řešením.

Parametry:
- x: Cílový bod v prostoru
- ρₜ: Cílová hustota
- Xₑ: Matice souřadnic uzlů elementu
- ρₑ: Vektor hustot v uzlech elementu

Vrací: Vektor lokálních souřadnic [ξ₁, ξ₂, ξ₃] nebo nothing při selhání
"""
function compute_coords_on_iso(
  x::Vector{Float64},
  ρₜ::Float64,
  Xₑ::Matrix{Float64},
  ρₑ::Vector{Float64}
)
  # Nejprve zkusíme Newton-Raphsonovu metodu
  # result = newton_raphson_solve(x, ρₜ, Xₑ, ρₑ)
  result = nlopt_solve(x, ρₜ, Xₑ, ρₑ)

  # Pokud Newton-Raphson selže nebo dá neplatné řešení, použijeme NLopt
  # if result === nothing || !is_solution_valid(result)
  #   return nlopt_solve(x, ρₜ, Xₑ, ρₑ)
  # end

  return result
end

"""
Implementace Newton-Raphsonovy metody pro nalezení lokálních souřadnic.
"""
function newton_raphson_solve(
  x::Vector{Float64},
  ρₜ::Float64,
  Xₑ::Matrix{Float64},
  ρₑ::Vector{Float64}
)
  N = MVector{8,Float64}(undef)
  d¹N_dξ¹ = MMatrix{8,3,Float64}(undef)

  ξ = [0.0, 0.0, 0.0]
  max_iter = 20  # Více iterací
  tol = 1e-10    # Přísnější tolerance

  for iter in 1:max_iter
    compute_shape_functions!(N, d¹N_dξ¹, ξ[1], ξ[2], ξ[3])

    # Výpočet reziduí s větší váhou pro hustotní podmínku
    R = zeros(4)
    for i in 1:3
      R[i] = x[i] - sum(Xₑ[i, k] * N[k] for k in 1:8)
    end
    R[4] = 1e3 * (sum(ρₑ[k] * N[k] for k in 1:8) - ρₜ)

    # Odpovídající úprava Jacobiánu
    J = zeros(4, 3)
    for i in 1:3
      for j in 1:3
        J[i, j] = -sum(Xₑ[i, k] * d¹N_dξ¹[k, j] for k in 1:8)
      end
    end
    for j in 1:3
      J[4, j] = 1e3 * sum(ρₑ[k] * d¹N_dξ¹[k, j] for k in 1:8)
    end

    Δξ = zeros(3)
    try
      Δξ = J \ (-R)
    catch e
      return nothing
    end

    ξ_new = ξ + Δξ

    if norm(Δξ, Inf) < tol  # Použití inf-normy pro přísnější kontrolu
      if all(-1.0 .<= ξ_new .<= 1.0)
        # Extra kontrola přesnosti řešení
        compute_shape_functions!(N, d¹N_dξ¹, ξ_new[1], ξ_new[2], ξ_new[3])
        density_error = abs(sum(ρₑ[k] * N[k] for k in 1:8) - ρₜ)
        if density_error < 1e-10
          return ξ_new
        end
      end
      return nothing
    end

    ξ = ξ_new
  end

  return nothing
end

"""
Výpočet reziduí a Jacobiánu pro Newton-Raphsonovu metodu.
"""
function compute_residual_and_jacobian(
  x::Vector{Float64},
  ρₜ::Float64,
  Xₑ::Matrix{Float64},
  ρₑ::Vector{Float64},
  ξ::Vector{Float64},
  N::MVector{8,Float64},
  d¹N_dξ¹::MMatrix{8,3,Float64}
)
  # Výpočet shape funkcí a jejich derivací
  compute_shape_functions!(N, d¹N_dξ¹, ξ[1], ξ[2], ξ[3])

  # Výpočet reziduí
  R = zeros(4)  # 3 prostorové + 1 hustotní podmínka

  # Geometrická rezidua
  for i in 1:3
    projected = sum(Xₑ[i, k] * N[k] for k in 1:8)
    R[i] = x[i] - projected
  end

  # Hustotní reziduum
  R[4] = sum(ρₑ[k] * N[k] for k in 1:8) - ρₜ

  # Sestavení Jacobiánu
  J = zeros(4, 3)

  # Geometrická část Jacobiánu
  for i in 1:3
    for j in 1:3
      J[i, j] = -sum(Xₑ[i, k] * d¹N_dξ¹[k, j] for k in 1:8)
    end
  end

  # Hustotní část Jacobiánu
  for j in 1:3
    J[4, j] = sum(ρₑ[k] * d¹N_dξ¹[k, j] for k in 1:8)
  end

  return R, J
end

"""
Záložní řešení pomocí NLopt optimalizace.
"""
function nlopt_solve(
  x::Vector{Float64},
  ρₜ::Float64,
  Xₑ::Matrix{Float64},
  ρₑ::Vector{Float64}
)
  N = MVector{8,Float64}(undef)
  d¹N_dξ¹ = MMatrix{8,3,Float64}(undef)

  # Create local optimizer with SLSQP
  local_opt = Opt(:LD_SLSQP, 3)
  local_opt.xtol_rel = 1e-12
  local_opt.ftol_rel = 1e-12
  local_opt.maxeval = 100

  # Create main optimizer with AUGLAG
  opt = Opt(:AUGLAG, 3)
  opt.local_optimizer = local_opt
  opt.xtol_rel = 1e-12
  opt.ftol_rel = 1e-12
  opt.maxeval = 200

  opt.lower_bounds = [-1.0, -1.0, -1.0]
  opt.upper_bounds = [1.0, 1.0, 1.0]

  # Objective function focused on minimizing distance
  function objective(ξ::Vector{Float64}, grad::Vector{Float64})
    compute_shape_functions!(N, d¹N_dξ¹, ξ[1], ξ[2], ξ[3])

    current_pos = zeros(3)
    for i in 1:3
      current_pos[i] = sum(Xₑ[i, k] * N[k] for k in 1:8)
    end

    return norm(current_pos - x)
  end

  # Density constraint as equality
  function constraint(ξ::Vector{Float64}, grad::Vector{Float64})
    compute_shape_functions!(N, d¹N_dξ¹, ξ[1], ξ[2], ξ[3])
    return sum(ρₑ[k] * N[k] for k in 1:8) - ρₜ
  end

  opt.min_objective = objective
  equality_constraint!(opt, constraint, 1e-12)

  # Generate starting points considering input point position
  function generate_starting_points(x::Vector{Float64})
    # Basic vertices and center point
    points = [
      [0.0, 0.0, 0.0],
      [-0.8, -0.8, -0.8], [0.8, -0.8, -0.8],
      [-0.8, 0.8, -0.8], [0.8, 0.8, -0.8],
      [-0.8, -0.8, 0.8], [0.8, -0.8, 0.8],
      [-0.8, 0.8, 0.8], [0.8, 0.8, 0.8]
    ]

    # Add points based on input position
    if norm(x) > 0
      norm_x = normalize(x)
      push!(points, [0.5 * norm_x[1], 0.5 * norm_x[2], 0.5 * norm_x[3]])
      push!(points, [0.8 * norm_x[1], 0.8 * norm_x[2], 0.8 * norm_x[3]])
    end

    return points
  end

  # Optimization with multiple starting points
  best_result = [0.0, 0.0, 0.0]  # Default result if everything fails
  best_value = Inf
  starting_points = generate_starting_points(x)

  for start_point in starting_points
    try
      (minf, minx, ret) = optimize(opt, start_point)

      # Validate solution
      if minf < best_value && all(-1.0 .<= minx .<= 1.0)
        compute_shape_functions!(N, d¹N_dξ¹, minx[1], minx[2], minx[3])
        density_error = abs(sum(ρₑ[k] * N[k] for k in 1:8) - ρₜ)

        if density_error < 1e-10
          best_value = minf
          best_result = minx
        end
      end
    catch e
      # Log error for debugging if needed
      # println("Optimization failed for starting point: ", start_point)
      continue
    end
  end

  # If no valid solution found, try one last time with relaxed constraints
  if best_value == Inf
    try
      opt.xtol_rel = 1e-8  # Relax tolerances
      opt.ftol_rel = 1e-8
      (minf, minx, ret) = optimize(opt, [0.0, 0.0, 0.0])
      if all(-1.0 .<= minx .<= 1.0)
        best_result = minx
      end
    catch e
      # Keep default result
    end
  end

  return best_result  # Always returns a Vector{Float64}, never nothing
end

# function nlopt_solve(
#   x::Vector{Float64},
#   ρₜ::Float64,
#   Xₑ::Matrix{Float64},
#   ρₑ::Vector{Float64}
# )
#   N = MVector{8,Float64}(undef)
#   d¹N_dξ¹ = MMatrix{8,3,Float64}(undef)
#
#   opt = Opt(:LD_SLSQP, 3)
#   opt.lower_bounds = [-1.0, -1.0, -1.0]
#   opt.upper_bounds = [1.0, 1.0, 1.0]
#
#   # Přísnější tolerance
#   opt.xtol_rel = 1e-10
#   opt.ftol_rel = 1e-10
#   opt.maxeval = 100
#
#   # Vyvážená objektivní funkce
#   function objective(ξ::Vector{Float64}, grad::Vector{Float64})
#     compute_shape_functions!(N, d¹N_dξ¹, ξ[1], ξ[2], ξ[3])
#
#     # Geometrická část
#     geom_error = 0.0
#     for i in 1:3
#       diff = x[i] - sum(Xₑ[i, k] * N[k] for k in 1:8)
#       geom_error += abs(diff)  # Používáme L1 normu místo L2
#     end
#
#     # Hustotní část
#     density_error = abs(sum(ρₑ[k] * N[k] for k in 1:8) - ρₜ)
#
#     # Vyvážení obou částí
#     return geom_error + 1e3 * density_error  # Větší váha pro hustotní podmínku
#   end
#
#   # Hustotní podmínka jako přesná rovnost
#   function constraint(ξ::Vector{Float64}, grad::Vector{Float64})
#     compute_shape_functions!(N, d¹N_dξ¹, ξ[1], ξ[2], ξ[3])
#     return sum(ρₑ[k] * N[k] for k in 1:8) - ρₜ
#   end
#
#   opt.min_objective = objective
#   equality_constraint!(opt, constraint, 1e-10)  # Přísnější tolerance pro hustotní podmínku
#
#   # Vícenásobné počáteční odhady pro lepší pokrytí prostoru řešení
#   best_result = nothing
#   best_value = Inf
#
#   starting_points = [
#     [0.0, 0.0, 0.0],
#     [-0.5, -0.5, -0.5],
#     [0.5, -0.5, -0.5],
#     [0.5, 0.5, -0.5],
#     [-0.5, 0.5, -0.5],
#     [-0.5, -0.5, 0.5],
#     [0.5, -0.5, 0.5],
#     [0.5, 0.5, 0.5],
#     [-0.5, 0.5, 0.5]
#   ]
#
#   for start_point in starting_points
#     try
#       (minf, minx, ret) = optimize(opt, start_point)
#       if minf < best_value && all(-1.0 .<= minx .<= 1.0)
#         best_value = minf
#         best_result = minx
#       end
#     catch e
#       continue
#     end
#   end
#
#   return best_result
# end
#

"""
Pomocná funkce pro kontrolu platnosti řešení.
"""
function is_solution_valid(ξ::Vector{Float64})
  return all(-1.0 .<= ξ .<= 1.0)
end
