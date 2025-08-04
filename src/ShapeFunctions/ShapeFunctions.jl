module ShapeFunctions

export shape_functions, compute_shape_and_derivatives!, compute_hex8_shape!
# export AbstractElement, HEX8, TET4

using StaticArrays
using LinearAlgebra
using Rho2sdf.ElementTypes

# Legacy functions for backward compatibility
include("hex8_shape.jl")

# Main dispatch function for shape functions
function shape_functions(::Type{HEX8}, ξ::Vector{Float64})
    return hex8_shape(ξ)  # Use existing implementation
end

function shape_functions(::Type{TET4}, λ::Vector{Float64})
    # Barycentric coordinates: λ = [λ1, λ2, λ3, λ4] where λ4 = 1 - λ1 - λ2 - λ3
    if length(λ) == 3
        λ4 = 1.0 - sum(λ)
        return SVector{4}([λ[1], λ[2], λ[3], λ4])
    elseif length(λ) == 4
        return SVector{4}(λ)
    else
        error("TET4 requires 3 or 4 barycentric coordinates")
    end
end

# Compute shape functions and derivatives for HEX8
function compute_shape_and_derivatives!(::Type{HEX8}, 
                                      N::MVector{8,Float64},
                                      dN::MMatrix{8,3,Float64},
                                      ξ::Float64, η::Float64, ζ::Float64)
    compute_hex8_shape!(N, dN, ξ, η, ζ)
end

# Compute shape functions and derivatives for TET4
function compute_shape_and_derivatives!(::Type{TET4},
                                      N::MVector{4,Float64},
                                      dN::MMatrix{4,3,Float64},
                                      ξ::Float64, η::Float64, ζ::Float64)
    # Natural coordinates to barycentric coordinates
    λ1 = ξ
    λ2 = η  
    λ3 = ζ
    λ4 = 1.0 - ξ - η - ζ
    
    # Shape functions
    N[1] = λ1
    N[2] = λ2
    N[3] = λ3
    N[4] = λ4
    
    # Derivatives with respect to natural coordinates
    # ∂N/∂ξ
    dN[1, 1] = 1.0
    dN[2, 1] = 0.0
    dN[3, 1] = 0.0
    dN[4, 1] = -1.0
    
    # ∂N/∂η
    dN[1, 2] = 0.0
    dN[2, 2] = 1.0
    dN[3, 2] = 0.0
    dN[4, 2] = -1.0
    
    # ∂N/∂ζ
    dN[1, 3] = 0.0
    dN[2, 3] = 0.0
    dN[3, 3] = 1.0
    dN[4, 3] = -1.0
end

# Generic dispatch wrapper for shape functions and derivatives
function compute_shape_and_derivatives!(element_type::Type{T}, coords...) where {T<:AbstractElement}
    error("Shape functions not implemented for element type: $T")
end

# Convenience function that creates appropriate arrays and calls compute_shape_and_derivatives!
function get_shape_and_derivatives(element_type::Type{HEX8}, ξ::Float64, η::Float64, ζ::Float64)
    N = MVector{8,Float64}(undef)
    dN = MMatrix{8,3,Float64}(undef)
    compute_shape_and_derivatives!(element_type, N, dN, ξ, η, ζ)
    return N, dN
end

function get_shape_and_derivatives(element_type::Type{TET4}, ξ::Float64, η::Float64, ζ::Float64)
    N = MVector{4,Float64}(undef)
    dN = MMatrix{4,3,Float64}(undef)
    compute_shape_and_derivatives!(element_type, N, dN, ξ, η, ζ)
    return N, dN
end

end
