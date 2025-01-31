@inline function compute_shape_functions!(
  N::MVector{8,Float64},
  d¹N_dξ¹::MMatrix{8,3,Float64},
  ξ₁::Float64,
  ξ₂::Float64,
  ξ₃::Float64
)
  # Předpočítání základních výrazů - compiler může optimalizovat do registrů
  ξ₁m1 = ξ₁ - 1
  ξ₁p1 = ξ₁ + 1
  ξ₂m1 = ξ₂ - 1
  ξ₂p1 = ξ₂ + 1
  ξ₃m1 = ξ₃ - 1
  ξ₃p1 = ξ₃ + 1

  # Předpočítání společných součinů
  t1 = ξ₁m1 * ξ₂m1  # Pro uzly 1 a 5
  t2 = ξ₁p1 * ξ₂m1  # Pro uzly 2 a 6
  t3 = ξ₁p1 * ξ₂p1  # Pro uzly 3 a 7
  t4 = ξ₁m1 * ξ₂p1  # Pro uzly 4 a 8

  # Konstantní koeficient
  coef = 0.125

  # Výpočet hodnot shape funkcí - přímo do výstupního vektoru
  N[1] = -coef * t1 * ξ₃m1
  N[2] = coef * t2 * ξ₃m1
  N[3] = -coef * t3 * ξ₃m1
  N[4] = coef * t4 * ξ₃m1
  N[5] = coef * t1 * ξ₃p1
  N[6] = -coef * t2 * ξ₃p1
  N[7] = coef * t3 * ξ₃p1
  N[8] = -coef * t4 * ξ₃p1

  # Předpočítání společných výrazů pro derivace
  d1_coef = coef * ξ₃m1
  d1_coefp = coef * ξ₃p1

  # Derivace podle ξ₁
  d¹N_dξ¹[1, 1] = -d1_coef * ξ₂m1
  d¹N_dξ¹[2, 1] = d1_coef * ξ₂m1
  d¹N_dξ¹[3, 1] = -d1_coef * ξ₂p1
  d¹N_dξ¹[4, 1] = d1_coef * ξ₂p1
  d¹N_dξ¹[5, 1] = d1_coefp * ξ₂m1
  d¹N_dξ¹[6, 1] = -d1_coefp * ξ₂m1
  d¹N_dξ¹[7, 1] = d1_coefp * ξ₂p1
  d¹N_dξ¹[8, 1] = -d1_coefp * ξ₂p1

  # Derivace podle ξ₂
  d¹N_dξ¹[1, 2] = -d1_coef * ξ₁m1
  d¹N_dξ¹[2, 2] = d1_coef * ξ₁p1
  d¹N_dξ¹[3, 2] = -d1_coef * ξ₁p1
  d¹N_dξ¹[4, 2] = d1_coef * ξ₁m1
  d¹N_dξ¹[5, 2] = d1_coefp * ξ₁m1
  d¹N_dξ¹[6, 2] = -d1_coefp * ξ₁p1
  d¹N_dξ¹[7, 2] = d1_coefp * ξ₁p1
  d¹N_dξ¹[8, 2] = -d1_coefp * ξ₁m1

  # Derivace podle ξ₃
  d¹N_dξ¹[1, 3] = -coef * t1
  d¹N_dξ¹[2, 3] = coef * t2
  d¹N_dξ¹[3, 3] = -coef * t3
  d¹N_dξ¹[4, 3] = coef * t4
  d¹N_dξ¹[5, 3] = coef * t1
  d¹N_dξ¹[6, 3] = -coef * t2
  d¹N_dξ¹[7, 3] = coef * t3
  d¹N_dξ¹[8, 3] = -coef * t4
end


function find_local_coordinates(
  Xₑ::AbstractMatrix{Float64},
  xₙ::AbstractVector{Float64}
)::Tuple{Bool,Vector{Float64}}
  # Prealokace polí - musí být uvnitř funkce, ne jako const
  N = MVector{8,Float64}(undef)
  d¹N_dξ¹ = MMatrix{8,3,Float64}(undef)
  x = MVector{3,Float64}(undef)
  R = MVector{3,Float64}(undef)
  J = MMatrix{3,3,Float64}(undef)

  # Definice starting points jako statické pole
  starting_points = @SVector [
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
    # Výpočet shape funkcí - zde @inbounds není vhodný kvůli komplexnosti výpočtů
    compute_shape_functions!(N, d¹N_dξ¹, ξ[1], ξ[2], ξ[3])

    # Použití @fastmath a @inbounds pro maticové operace
    @fastmath @inbounds begin
      # Násobení matic - známe přesné dimenze
      mul!(x, Xₑ, N)
      # Vektorové operace - známe délky vektorů
      @. R = x - xₙ
    end

    # Sum of squares - zde @fastmath může urychlit výpočet
    @fastmath f = sum(abs2, R)

    if !isempty(grad)
      @fastmath @inbounds begin
        mul!(J, Xₑ, d¹N_dξ¹)
        mul!(grad, J', R)
        grad .*= 2
      end
    end

    return f
  end

  # Prealokace výsledků s explicitními typy pro type stability
  best_solution::Vector{Float64} = Vector{Float64}(undef, 3)
  best_objective::Float64 = Inf
  found_solution::Bool = false

  # Optimizer settings for finding local coordinates
  opt = Opt(:LD_LBFGS, 3)  # L-BFGS optimizer for 3D problem - good choice for smooth functions with gradients

  # Search space boundaries - narrowed to unit cube with small margin
  opt.lower_bounds = fill(-1.1, 3)  # Lower bound slightly below -1 to cover entire element
  opt.upper_bounds = fill(1.1, 3)   # Upper bound slightly above 1 to cover entire element

  # Memory and convergence parameters
  opt.vector_storage = 5    # Number of stored previous iterations - optimal for 3D
  opt.ftol_rel = 1e-8      # Relative tolerance for function value change - more precise than original
  opt.xtol_rel = 1e-6      # Relative tolerance for parameter change - good speed/precision compromise

  # Limit parameters
  opt.maxeval = 200        # Maximum number of function evaluations - usually less than 500 is sufficient
  opt.maxtime = 1.0        # Time limit in seconds - reasonable safeguard against infinite loops

  # Additional parameters for robustness
  opt.ftol_abs = 1e-10     # Absolute function value tolerance - useful for known target precision
  opt.initial_step = 0.1   # Initial step size - suitable for normalized solution space

  opt.min_objective = objective  # Setting the minimization objective function

  # Iterace přes starting points
  @inbounds for start_point in starting_points
    try
      minf, minx, ret = NLopt.optimize(opt, SVector{3}(start_point))
      if minf < best_objective && ret ∈ [:SUCCESS, :FTOL_REACHED, :XTOL_REACHED]
        best_objective = minf
        copyto!(best_solution, minx)
        found_solution = true
      end
    catch _
      continue
    end
  end

  return found_solution ? (true, best_solution) : (false, fill(10.0, 3))
end


"""
Old version that works fine but it is slower
"""
# function find_local_coordinates(
#   Xₑ,
#   xₙ
# )
#   starting_points::Vector{Tuple{Float64,Float64,Float64}} = [
#     (0.0, 0.0, 0.0),
#     (-0.5, -0.5, -0.5),
#     (0.5, -0.5, -0.5),
#     (0.5, 0.5, -0.5),
#     (-0.5, 0.5, -0.5),
#     (-0.5, -0.5, 0.5),
#     (0.5, -0.5, 0.5),
#     (0.5, 0.5, 0.5),
#     (-0.5, 0.5, 0.5)
#   ]
#
#   function objective(ξ::Vector, grad::Vector)
#     ξ₁, ξ₂, ξ₃ = ξ
#
#     # Compute N(ξ)
#     N = [
#       -1 / 8 * (ξ₁ - 1) * (ξ₂ - 1) * (ξ₃ - 1),
#       1 / 8 * (ξ₁ + 1) * (ξ₂ - 1) * (ξ₃ - 1),
#       -1 / 8 * (ξ₁ + 1) * (ξ₂ + 1) * (ξ₃ - 1),
#       1 / 8 * (ξ₁ - 1) * (ξ₂ + 1) * (ξ₃ - 1),
#       1 / 8 * (ξ₁ - 1) * (ξ₂ - 1) * (ξ₃ + 1),
#       -1 / 8 * (ξ₁ + 1) * (ξ₂ - 1) * (ξ₃ + 1),
#       1 / 8 * (ξ₁ + 1) * (ξ₂ + 1) * (ξ₃ + 1),
#       -1 / 8 * (ξ₁ - 1) * (ξ₂ + 1) * (ξ₃ + 1)
#     ]
#
#     d¹N_dξ¹ = zeros(Float64, 8, 3)
#     d¹N_dξ¹[1, :] = [
#       -0.125 * (ξ₂ - 1) * (ξ₃ - 1)
#       -0.125 * (ξ₁ - 1) * (ξ₃ - 1)
#       -0.125 * (ξ₁ - 1) * (ξ₂ - 1)
#     ]
#
#     d¹N_dξ¹[2, :] = [
#       0.125 * (ξ₂ - 1) * (ξ₃ - 1)
#       0.125 * (ξ₁ + 1) * (ξ₃ - 1)
#       0.125 * (ξ₁ + 1) * (ξ₂ - 1)
#     ]
#
#     d¹N_dξ¹[3, :] = [
#       -0.125 * (ξ₂ + 1) * (ξ₃ - 1)
#       -0.125 * (ξ₁ + 1) * (ξ₃ - 1)
#       -0.125 * (ξ₁ + 1) * (ξ₂ + 1)
#     ]
#
#     d¹N_dξ¹[4, :] = [
#       0.125 * (ξ₂ + 1) * (ξ₃ - 1)
#       0.125 * (ξ₁ - 1) * (ξ₃ - 1)
#       0.125 * (ξ₁ - 1) * (ξ₂ + 1)
#     ]
#
#     d¹N_dξ¹[5, :] = [
#       0.125 * (ξ₂ - 1) * (ξ₃ + 1)
#       0.125 * (ξ₁ - 1) * (ξ₃ + 1)
#       0.125 * (ξ₁ - 1) * (ξ₂ - 1)
#     ]
#
#     d¹N_dξ¹[6, :] = [
#       -0.125 * (ξ₂ - 1) * (ξ₃ + 1)
#       -0.125 * (ξ₁ + 1) * (ξ₃ + 1)
#       -0.125 * (ξ₁ + 1) * (ξ₂ - 1)
#     ]
#
#     d¹N_dξ¹[7, :] = [
#       0.125 * (ξ₂ + 1) * (ξ₃ + 1)
#       0.125 * (ξ₁ + 1) * (ξ₃ + 1)
#       0.125 * (ξ₁ + 1) * (ξ₂ + 1)
#     ]
#
#     d¹N_dξ¹[8, :] = [
#       -0.125 * (ξ₂ + 1) * (ξ₃ + 1)
#       -0.125 * (ξ₁ - 1) * (ξ₃ + 1)
#       -0.125 * (ξ₁ - 1) * (ξ₂ + 1)
#     ]
#
#     # Compute x(ξ)
#     x = Xₑ * N  # x is a 3-element vector
#
#     # Compute residual R = x - xₙ
#     R = x - xₙ  # R is a 3-element vector
#
#     # Compute objective function value f
#     f = dot(R, R)  # Equivalent to sum(R[i]^2 for i in 1:3)
#
#     if length(grad) > 0
#       # Compute the derivatives of N with respect to ξ
#       dN_dξ = zeros(Float64, 8, 3)  # 8 shape functions x 3 variables
#
#       # Assign your computed derivatives to dN_dξ
#       dN_dξ .= d¹N_dξ¹  # Ensure dN_dξ is correctly populated
#
#       # Compute Jacobian J = Xₑ * dN_dξ
#       J = Xₑ * dN_dξ  # J is 3x3
#
#       # Compute gradient grad = 2 * Jᵗ * R
#       grad .= 2 * (J' * R)
#     end
#
#     return f
#   end
#
#   best_solution = nothing
#   best_objective = Inf
#
#   for start_point in starting_points
#     #   opt = Opt(:LN_COBYLA, 3)
#     opt = Opt(:LD_LBFGS, 3)
#     opt.lower_bounds = [-5.0, -5.0, -5.0]
#     opt.upper_bounds = [5.0, 5.0, 5.0]
#     opt.xtol_rel = 1e-6
#     opt.maxeval = 500
#     opt.min_objective = objective
#     opt.maxtime = 1.0
#
#     (minf, minx, ret) = NLopt.optimize(opt, collect(start_point))
#
#     if minf < best_objective
#       best_objective = minf
#       best_solution = minx
#     end
#   end
#
#   if best_solution === nothing
#     return (false, [10.0, 10.0, 10.0])
#   else
#     return (true, best_solution)
#   end
# end


