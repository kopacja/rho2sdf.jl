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
    N[2] =  coef * temp[2] * ξ₃m1
    N[3] = -coef * temp[3] * ξ₃m1
    N[4] =  coef * temp[4] * ξ₃m1
    
    # Druhé čtyři uzly (ξ₃p1)
    N[5] =  coef * temp[1] * ξ₃p1
    N[6] = -coef * temp[2] * ξ₃p1
    N[7] =  coef * temp[3] * ξ₃p1
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
  Xₑ::AbstractMatrix{Float64},  # Změna z Matrix na AbstractMatrix
  xₙ::AbstractVector{Float64}   # Změna z Vector na AbstractVector
)::Tuple{Bool,Vector{Float64}}
  # Definice statických polí
  N = MVector{8,Float64}(undef)
  d¹N_dξ¹ = MMatrix{8,3,Float64}(undef)
  x = MVector{3,Float64}(undef)
  R = MVector{3,Float64}(undef)
  J = MMatrix{3,3,Float64}(undef)

  function objective(ξ::Vector, grad::Vector)
    # Výpočet shape funkcí a jejich derivací
    compute_shape_functions!(N, d¹N_dξ¹, ξ[1], ξ[2], ξ[3])

    # In-place maticové operace s kontrolou kompatibility
    mul!(x, Xₑ, N)  # Funguje s AbstractMatrix
    @. R = x - xₙ   # Funguje s AbstractVector

    f = sum(abs2, R)

    if !isempty(grad)
      mul!(J, Xₑ, d¹N_dξ¹)
      mul!(grad, J', R)
      grad .*= 2
    end

    return f
  end

  # Počáteční body
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

  best_solution = Vector{Float64}(undef, 3)  # Přealokace
  best_objective = Inf
  found_solution = false

  opt = Opt(:LD_LBFGS, 3)
  opt.lower_bounds = [-5.0, -5.0, -5.0]
  opt.upper_bounds = [5.0, 5.0, 5.0]
  opt.xtol_rel = 1e-6
  opt.maxeval = 500
  opt.min_objective = objective
  opt.maxtime = 1.0

  for start_point in starting_points
    try
      minf, minx, ret = NLopt.optimize(opt, collect(start_point))
      if minf < best_objective && ret ∈ [:SUCCESS, :FTOL_REACHED, :XTOL_REACHED]
        best_objective = minf
        copyto!(best_solution, minx)  # In-place kopírování
        found_solution = true
      end
    catch e
      # Pokračujeme s dalším bodem při chybě optimalizace
      continue
    end
  end

  return found_solution ? (true, best_solution) : (false, [10.0, 10.0, 10.0])
end

