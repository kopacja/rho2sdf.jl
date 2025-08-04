using Rho2sdf.ElementTypes
using Rho2sdf.ShapeFunctions

# Main dispatch function
function compute_coords_on_iso(x::Vector, ρₜ::Float64, Xₑ::Matrix, ρₑ::Vector, 
                             element_type::Type{T}) where {T<:AbstractElement}
  return compute_coords_on_iso_dispatch(x, ρₜ, Xₑ, ρₑ, element_type)
end

# Backward compatibility - defaults to HEX8
function compute_coords_on_iso(x::Vector, ρₜ::Float64, Xₑ::Matrix, ρₑ::Vector)
  return compute_coords_on_iso_dispatch(x, ρₜ, Xₑ, ρₑ, HEX8)
end

# HEX8 implementation
function compute_coords_on_iso_dispatch(x::Vector, ρₜ::Float64, Xₑ::Matrix, ρₑ::Vector, 
                                      ::Type{HEX8})
  
  opt = Opt(:LD_SLSQP, 3)
  ftol_rel!(opt, 1e-5)
  ftol_abs!(opt, 1e-5)
  xtol_rel!(opt, 1e-5)
  maxeval!(opt, 1000)

  lower_bounds!(opt, SVector{3}(-1.0, -1.0, -1.0))
  upper_bounds!(opt, SVector{3}(1.0, 1.0, 1.0))

  # Pre-allocation of all helper vectors
  N = MVector{8,Float64}(undef)
  dN = MMatrix{8,3,Float64}(undef)
  temp_sums = MVector{3,Float64}(undef)

  # Closure to share pre-allocated vectors
  let N = N, dN = dN, temp_sums = temp_sums
    function objective(ξ::Vector, grad::Vector)
      compute_shape_and_derivatives!(HEX8, N, dN, ξ[1], ξ[2], ξ[3])

      # Pre-compute projections for all dimensions at once
      @inbounds for i in 1:3
        temp_sums[i] = sum(Xₑ[i, k] * N[k] for k in 1:8)
      end

      # Calculate residuals using pre-computed sums
      rx = x[1] - temp_sums[1]
      ry = x[2] - temp_sums[2]
      rz = x[3] - temp_sums[3]

      if !isempty(grad)
        @inbounds for j in 1:3
          # Pre-compute derivatives for all dimensions
          @inbounds for i in 1:3
            temp_sums[i] = sum(Xₑ[i, k] * dN[k, j] for k in 1:8)
          end
          grad[j] = -2.0 * (rx * temp_sums[1] + ry * temp_sums[2] + rz * temp_sums[3])
        end
      end

      return rx * rx + ry * ry + rz * rz
    end

    function constraint(ξ::Vector, grad::Vector)
      compute_shape_and_derivatives!(HEX8, N, dN, ξ[1], ξ[2], ξ[3])

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

  # Optimization with preset precision
  initial_ξ = zeros(SVector{3,Float64})
  _, minx, ret = optimize(opt, initial_ξ)

  # Check convergence
  if ret == :FAILURE
    @warn "Optimization failed to converge"
  end

  return minx
end

# TET4 implementation
function compute_coords_on_iso_dispatch(x::Vector, ρₜ::Float64, Xₑ::Matrix, ρₑ::Vector, 
                                      ::Type{TET4})
  
  opt = Opt(:LD_SLSQP, 3)
  ftol_rel!(opt, 1e-5)
  ftol_abs!(opt, 1e-5)
  xtol_rel!(opt, 1e-5)
  maxeval!(opt, 1000)

  # Barycentric coordinate bounds
  lower_bounds!(opt, SVector{3}(0.0, 0.0, 0.0))
  upper_bounds!(opt, SVector{3}(1.0, 1.0, 1.0))

  # Pre-allocation of helper vectors
  N = MVector{4,Float64}(undef)
  dN = MMatrix{4,3,Float64}(undef)
  temp_sums = MVector{3,Float64}(undef)

  let N = N, dN = dN, temp_sums = temp_sums
    function objective(λ::Vector, grad::Vector)
      # Ensure barycentric constraint
      λ₄ = 1.0 - sum(λ)
      if λ₄ < 0.0
        return 1e10  # Penalty for invalid coordinates
      end
      
      compute_shape_and_derivatives!(TET4, N, dN, λ[1], λ[2], λ[3])

      # Pre-compute projections
      @inbounds for i in 1:3
        temp_sums[i] = sum(Xₑ[i, k] * N[k] for k in 1:4)
      end

      # Calculate residuals
      rx = x[1] - temp_sums[1]
      ry = x[2] - temp_sums[2]
      rz = x[3] - temp_sums[3]

      if !isempty(grad)
        @inbounds for j in 1:3
          @inbounds for i in 1:3
            temp_sums[i] = sum(Xₑ[i, k] * dN[k, j] for k in 1:4)
          end
          grad[j] = -2.0 * (rx * temp_sums[1] + ry * temp_sums[2] + rz * temp_sums[3])
        end
      end

      return rx * rx + ry * ry + rz * rz
    end

    function constraint(λ::Vector, grad::Vector)
      λ₄ = 1.0 - sum(λ)
      if λ₄ < 0.0
        return 1e10  # Penalty
      end
      
      compute_shape_and_derivatives!(TET4, N, dN, λ[1], λ[2], λ[3])

      if !isempty(grad)
        @inbounds for j in 1:3
          grad[j] = sum(ρₑ[k] * dN[k, j] for k in 1:4)
        end
      end

      return sum(ρₑ[k] * N[k] for k in 1:4) - ρₜ
    end

    # Add barycentric constraint: λ₁ + λ₂ + λ₃ ≤ 1
    function barycentric_constraint(λ::Vector, grad::Vector)
      if !isempty(grad)
        grad[1] = 1.0
        grad[2] = 1.0  
        grad[3] = 1.0
      end
      return sum(λ) - 1.0
    end

    min_objective!(opt, objective)
    equality_constraint!(opt, constraint, 1e-8)
    inequality_constraint!(opt, barycentric_constraint, 1e-8)
  end

  # Optimization
  initial_λ = SVector{3}(0.25, 0.25, 0.25)  # Start at centroid
  _, minx, ret = optimize(opt, initial_λ)

  if ret == :FAILURE
    @warn "TET4 optimization failed to converge"
  end

  return minx
end
