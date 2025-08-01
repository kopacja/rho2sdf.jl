function find_local_coordinates(
  Xₑ::AbstractMatrix{Float64},
  xₙ::AbstractVector{Float64}
)::Tuple{Bool,Vector{Float64}}

  # Preallocation of arrays - must be inside the function, not as const
  N = MVector{8,Float64}(undef)
  d¹N_dξ¹ = MMatrix{8,3,Float64}(undef)
  x = MVector{3,Float64}(undef)
  R = MVector{3,Float64}(undef)
  J = MMatrix{3,3,Float64}(undef)

  # Definition of starting points as a static array
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
    # Calculation of shape functions - @inbounds not suitable here due to calculation complexity
    compute_hex8_shape!(N, d¹N_dξ¹, ξ[1], ξ[2], ξ[3])

    # Using @fastmath and @inbounds for matrix operations
    @fastmath @inbounds begin
      # Matrix multiplication - we know exact dimensions
      mul!(x, Xₑ, N)
      # Vector operations - we know vector lengths
      @. R = x - xₙ
    end

    # Sum of squares - here @fastmath can speed up the calculation
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

  # Preallocation of results with explicit types for type stability
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

  # Iteration through starting points
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


