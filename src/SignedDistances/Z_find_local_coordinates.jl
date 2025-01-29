function compute_shape_functions!(
  N::MVector{8,Float64},
  d¹N_dξ¹::MMatrix{8,3,Float64},
  ξ₁::Float64,
  ξ₂::Float64,
  ξ₃::Float64
)
  # Předpočítání výrazů pomocí statických vektorů
  ξ_terms = @SVector [
    ξ₁ - 1,  # ξ₁m1
    ξ₁ + 1,  # ξ₁p1
    ξ₂ - 1,  # ξ₂m1
    ξ₂ + 1,  # ξ₂p1
    ξ₃ - 1,  # ξ₃m1
    ξ₃ + 1   # ξ₃p1
  ]

  # Pojmenování indexů pro lepší čitelnost (bez const)
  ξ₁m1, ξ₁p1, ξ₂m1, ξ₂p1, ξ₃m1, ξ₃p1 = ξ_terms

  # Předpočítání společných součinů pomocí statického vektoru
  temp = @SVector [
    ξ₁m1 * ξ₂m1,  # temp1: Pro uzly 1 a 5
    ξ₁p1 * ξ₂m1,  # temp2: Pro uzly 2 a 6
    ξ₁p1 * ξ₂p1,  # temp3: Pro uzly 3 a 7
    ξ₁m1 * ξ₂p1   # temp4: Pro uzly 4 a 8
  ]

  # Společný koeficient (bez const)
  coef = 0.125

  # Efektivní výpočet hodnot shape funkcí
  # První čtyři uzly (ξ₃m1)
  N[1] = -coef * temp[1] * ξ₃m1
  N[2] = coef * temp[2] * ξ₃m1
  N[3] = -coef * temp[3] * ξ₃m1
  N[4] = coef * temp[4] * ξ₃m1

  # Druhé čtyři uzly (ξ₃p1)
  N[5] = coef * temp[1] * ξ₃p1
  N[6] = -coef * temp[2] * ξ₃p1
  N[7] = coef * temp[3] * ξ₃p1
  N[8] = -coef * temp[4] * ξ₃p1

  # Předpočítáme společné výrazy pro derivace
  d1_terms = @SVector [
    -coef * ξ₂m1 * ξ₃m1,
    coef * ξ₂m1 * ξ₃m1,
    -coef * ξ₂p1 * ξ₃m1,
    coef * ξ₂p1 * ξ₃m1,
    coef * ξ₂m1 * ξ₃p1,
    -coef * ξ₂m1 * ξ₃p1,
    coef * ξ₂p1 * ξ₃p1,
    -coef * ξ₂p1 * ξ₃p1
  ]

  d2_terms = @SVector [
    -coef * ξ₁m1 * ξ₃m1,
    coef * ξ₁p1 * ξ₃m1,
    -coef * ξ₁p1 * ξ₃m1,
    coef * ξ₁m1 * ξ₃m1,
    coef * ξ₁m1 * ξ₃p1,
    -coef * ξ₁p1 * ξ₃p1,
    coef * ξ₁p1 * ξ₃p1,
    -coef * ξ₁m1 * ξ₃p1
  ]

  d3_terms = @SVector [
    -coef * temp[1],
    coef * temp[2],
    -coef * temp[3],
    coef * temp[4],
    coef * temp[1],
    -coef * temp[2],
    coef * temp[3],
    -coef * temp[4]
  ]

  # Přiřazení derivací
  d¹N_dξ¹[1:8, 1] = d1_terms
  d¹N_dξ¹[1:8, 2] = d2_terms
  d¹N_dξ¹[1:8, 3] = d3_terms
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

