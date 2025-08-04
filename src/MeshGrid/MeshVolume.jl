"""
Calculate volume of elements using Gaussian quadrature with type-based dispatch.
"""
function calculate_mesh_volume(
  X::Vector{Vector{Float64}},
  IEN::Vector{Vector{Int64}},
  rho::Vector{Float64},
  element_type::Type{T};
  quad_order::Int=3
) where {T<:AbstractElement}
  
  print_info("Computing volume...")
  
  # Get quadrature points and weights
  gp_raw, w_raw = gausslegendre(quad_order)
  gp = SVector{quad_order}(gp_raw)
  w = SVector{quad_order}(w_raw)
  
  # Atomic variables for thread-safe volume accumulation
  domain_volume = Atomic{Float64}(0.0)
  TO_volume = Atomic{Float64}(0.0)
  
  # Parallel processing of elements
  @threads for elem in 1:length(IEN)
    elem_volume = calculate_element_volume(X, IEN[elem], element_type, gp, w)
    
    # Thread-safe update of volumes
    atomic_add!(domain_volume, elem_volume)
    atomic_add!(TO_volume, elem_volume * rho[elem])
  end
  
  # Retrieve final volume values
  final_domain_volume = domain_volume[]
  final_TO_volume = TO_volume[]
  
  # Output results
  println("Topology optimization domain volume: ", round(final_domain_volume, sigdigits=6))
  print_data("Optimized shape volume: $(round(final_TO_volume, sigdigits=6))")
  println("Volume fraction: ", round(final_TO_volume / final_domain_volume, sigdigits=6))
  
  return [final_domain_volume, (final_TO_volume / final_domain_volume)]
end

# Element volume calculation for HEX8
function calculate_element_volume(X::Vector{Vector{Float64}}, element_nodes::Vector{Int64}, 
                                ::Type{HEX8}, gp::SVector{N,Float64}, w::SVector{N,Float64}) where {N}
  
  # Pre-allocate arrays for shape functions and derivatives
  local_N = MVector{8,Float64}(undef)
  local_dN = MMatrix{8,3,Float64}(undef)
  
  # Convert element coordinates to static matrix
  xe = @SMatrix [X[element_nodes[j]][i] for i in 1:3, j in 1:8]
  
  # Calculate element volume using Gaussian quadrature
  elem_volume = 0.0
  for k in 1:N, j in 1:N, i in 1:N
    local_coords = SVector{3}(gp[i], gp[j], gp[k])
    
    # Calculate shape functions and derivatives
    compute_shape_and_derivatives!(HEX8, local_N, local_dN, local_coords[1], local_coords[2], local_coords[3])
    
    # Compute Jacobian
    J = xe * local_dN
    
    # Accumulate contribution to volume
    elem_volume += w[i] * w[j] * w[k] * abs(det(J))
  end
  
  return elem_volume
end

# Element volume calculation for TET4
function calculate_element_volume(X::Vector{Vector{Float64}}, element_nodes::Vector{Int64}, 
                                ::Type{TET4}, gp::SVector{N,Float64}, w::SVector{N,Float64}) where {N}
  
  # Pre-allocate arrays for shape functions and derivatives
  local_N = MVector{4,Float64}(undef)
  local_dN = MMatrix{4,3,Float64}(undef)
  
  # Convert element coordinates to static matrix
  xe = @SMatrix [X[element_nodes[j]][i] for i in 1:3, j in 1:4]
  
  # For tetrahedron, we need to map from cube quadrature to tetrahedral domain
  # Using transformation from unit cube to unit tetrahedron
  elem_volume = 0.0
  for k in 1:N, j in 1:N, i in 1:N
    # Map cube coordinates to tetrahedral coordinates
    ξ_cube = gp[i]
    η_cube = gp[j]
    ζ_cube = gp[k]
    
    # Transform to tetrahedral coordinates [0,1]³
    ξ = (ξ_cube + 1.0) / 2.0
    η = (η_cube + 1.0) / 2.0 * (1.0 - ξ)
    ζ = (ζ_cube + 1.0) / 2.0 * (1.0 - ξ - η)
    
    # Skip if outside tetrahedral domain
    if ξ < 0 || η < 0 || ζ < 0 || ξ + η + ζ > 1.0
      continue
    end
    
    # Calculate shape functions and derivatives
    compute_shape_and_derivatives!(TET4, local_N, local_dN, ξ, η, ζ)
    
    # Compute Jacobian
    J = xe * local_dN
    
    # Jacobian of transformation from cube to tetrahedron
    jacobian_transform = (1.0 - ξ)^2 * (1.0 - ξ - η) / 8.0
    
    # Accumulate contribution to volume
    elem_volume += w[i] * w[j] * w[k] * abs(det(J)) * jacobian_transform
  end
  
  return elem_volume
end
