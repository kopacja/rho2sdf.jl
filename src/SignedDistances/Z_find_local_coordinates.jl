function find_local_coordinates(
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

  N = zeros(Float64, 8)
  d¹N_dξ¹ = zeros(Float64, 8, 3)

  function objective(ξ::Vector, grad::Vector)
    ξ₁, ξ₂, ξ₃ = ξ

    # Předpočítání často používaných výrazů
    ξ₁m1 = ξ₁ - 1
    ξ₁p1 = ξ₁ + 1
    ξ₂m1 = ξ₂ - 1
    ξ₂p1 = ξ₂ + 1
    ξ₃m1 = ξ₃ - 1
    ξ₃p1 = ξ₃ + 1

    # Společný koeficient pro všechny výpočty
    coef = 0.125

    # Kombinovaný výpočet N a jeho derivací
    # Node 1
    temp = ξ₁m1 * ξ₂m1
    N[1] = -coef * temp * ξ₃m1
    d¹N_dξ¹[1, :] = [-coef * ξ₂m1 * ξ₃m1,
      -coef * ξ₁m1 * ξ₃m1,
      -coef * temp]

    # Node 2
    temp = ξ₁p1 * ξ₂m1
    N[2] = coef * temp * ξ₃m1
    d¹N_dξ¹[2, :] = [coef * ξ₂m1 * ξ₃m1,
      coef * ξ₁p1 * ξ₃m1,
      coef * temp]

    # Node 3
    temp = ξ₁p1 * ξ₂p1
    N[3] = -coef * temp * ξ₃m1
    d¹N_dξ¹[3, :] = [-coef * ξ₂p1 * ξ₃m1,
      -coef * ξ₁p1 * ξ₃m1,
      -coef * temp]

    # Node 4
    temp = ξ₁m1 * ξ₂p1
    N[4] = coef * temp * ξ₃m1
    d¹N_dξ¹[4, :] = [coef * ξ₂p1 * ξ₃m1,
      coef * ξ₁m1 * ξ₃m1,
      coef * temp]

    # Node 5
    temp = ξ₁m1 * ξ₂m1
    N[5] = coef * temp * ξ₃p1
    d¹N_dξ¹[5, :] = [coef * ξ₂m1 * ξ₃p1,
      coef * ξ₁m1 * ξ₃p1,
      coef * temp]

    # Node 6
    temp = ξ₁p1 * ξ₂m1
    N[6] = -coef * temp * ξ₃p1
    d¹N_dξ¹[6, :] = [-coef * ξ₂m1 * ξ₃p1,
      -coef * ξ₁p1 * ξ₃p1,
      -coef * temp]

    # Node 7
    temp = ξ₁p1 * ξ₂p1
    N[7] = coef * temp * ξ₃p1
    d¹N_dξ¹[7, :] = [coef * ξ₂p1 * ξ₃p1,
      coef * ξ₁p1 * ξ₃p1,
      coef * temp]

    # Node 8
    temp = ξ₁m1 * ξ₂p1
    N[8] = -coef * temp * ξ₃p1
    d¹N_dξ¹[8, :] = [-coef * ξ₂p1 * ξ₃p1,
      -coef * ξ₁m1 * ξ₃p1,
      -coef * temp]

    # Výpočet residuálu a objektivní funkce
    x = Xₑ * N
    R = x - xₙ
    f = dot(R, R)

    if !isempty(grad)
      J = Xₑ * d¹N_dξ¹
      grad .= 2 * (J' * R)
    end

    return f
  end

  best_solution = nothing
  best_objective = Inf

  for start_point in starting_points
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



