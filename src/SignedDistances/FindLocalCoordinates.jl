using Rho2sdf.ElementTypes
using Rho2sdf.ShapeFunctions

# Main dispatch function
function find_local_coordinates(Xₑ::AbstractMatrix{Float64}, xₙ::AbstractVector{Float64}, 
                               element_type::Type{T}) where {T<:AbstractElement}
  return find_local_coordinates_dispatch(Xₑ, xₙ, element_type)
end

# Backward compatibility - defaults to HEX8
function find_local_coordinates(Xₑ::AbstractMatrix{Float64}, xₙ::AbstractVector{Float64})
  return find_local_coordinates_dispatch(Xₑ, xₙ, HEX8)
end

# HEX8 implementation (existing algorithm)
function find_local_coordinates_dispatch(Xₑ::AbstractMatrix{Float64}, xₙ::AbstractVector{Float64}, 
                                       ::Type{HEX8})::Tuple{Bool,Vector{Float64}}

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
    # Calculation of shape functions
    compute_shape_and_derivatives!(HEX8, N, d¹N_dξ¹, ξ[1], ξ[2], ξ[3])

    # Using @fastmath and @inbounds for matrix operations
    @fastmath @inbounds begin
      # Matrix multiplication
      mul!(x, Xₑ, N)
      # Vector operations
      @. R = x - xₙ
    end

    # Sum of squares
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
  opt = Opt(:LD_LBFGS, 3)

  # Search space boundaries - narrowed to unit cube with small margin
  opt.lower_bounds = fill(-1.1, 3)
  opt.upper_bounds = fill(1.1, 3)

  # Memory and convergence parameters
  opt.vector_storage = 5
  opt.ftol_rel = 1e-8
  opt.xtol_rel = 1e-6

  # Limit parameters
  opt.maxeval = 200
  opt.maxtime = 1.0
  
  # Additional parameters for robustness
  opt.ftol_abs = 1e-10
  opt.initial_step = 0.1

  opt.min_objective = objective

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

# TET4 implementation (simplified for barycentric coordinates)
function find_local_coordinates_dispatch(Xₑ::AbstractMatrix{Float64}, xₙ::AbstractVector{Float64}, 
                                       ::Type{TET4})::Tuple{Bool,Vector{Float64}}

  # For tetrahedron, we can use barycentric coordinates
  # Solve: xₙ = λ₁*x₁ + λ₂*x₂ + λ₃*x₃ + λ₄*x₄
  # with constraint: λ₁ + λ₂ + λ₃ + λ₄ = 1
  
  # Get tetrahedron vertices
  x₁ = Xₑ[:, 1]
  x₂ = Xₑ[:, 2]
  x₃ = Xₑ[:, 3]
  x₄ = Xₑ[:, 4]
  
  # Set up system: [x₂-x₁, x₃-x₁, x₄-x₁] * [λ₂, λ₃, λ₄]ᵀ = xₙ - x₁
  A = hcat(x₂ - x₁, x₃ - x₁, x₄ - x₁)
  b = xₙ - x₁
  
  try
    # Solve for λ₂, λ₃, λ₄
    λ_234 = A \ b
    
    # Calculate λ₁
    λ₁ = 1.0 - sum(λ_234)
    
    # Return in natural coordinates format [λ₁, λ₂, λ₃]
    # (λ₄ is implicit as 1 - λ₁ - λ₂ - λ₃)
    local_coords = [λ₁, λ_234[1], λ_234[2]]
    
    # Validate coordinates
    if validate_local_coords(TET4, [λ₁, λ_234...])
      return (true, local_coords)
    else
      return (false, fill(10.0, 3))
    end
    
  catch e
    # Singular matrix or other numerical issues
    return (false, fill(10.0, 3))
  end
end
