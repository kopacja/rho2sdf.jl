"""
Calculate volume of 8-node hexahedral elements using Gaussian quadrature.
Inputs:
- X: 3×n matrix of nodal coordinates
- IEN: 8×m connectivity matrix
- rho: vector of material densities for topology optimization
Returns:
- Array containing [total mesh volume, volume fraction]
"""
function calculate_mesh_volume(
  X::Vector{Vector{Float64}},
  IEN::Vector{Vector{Int64}},
  rho::Vector{Float64};
  quad_order::Int=3  # Default integration order is 3, same as in the original function
)
  print_info("Computing volume...")
  # Calculate Gauss-Legendre points and weights for the specified integration order
  gp_raw, w_raw = gausslegendre(quad_order)
  # Convert to static vectors for better performance
  gp = SVector{quad_order}(gp_raw)
  w = SVector{quad_order}(w_raw)
  
  # Atomic variables for thread-safe volume accumulation
  domain_volume = Atomic{Float64}(0.0)
  TO_volume = Atomic{Float64}(0.0)
  
  # Parallel processing of elements
  @threads for elem in 1:length(IEN)
    # Pre-allocate arrays for shape functions and their derivatives - moved inside the parallel loop for thread safety
    local_N = MVector{8,Float64}(undef)
    local_dN = MMatrix{8,3,Float64}(undef)
    
    # Convert element coordinates to a static matrix for efficiency
    xe = @SMatrix [X[IEN[elem][j]][i] for i in 1:3, j in 1:8]
    
    # Calculate element volume using Gaussian quadrature
    elem_volume = 0.0
    for k in 1:quad_order, j in 1:quad_order, i in 1:quad_order
      local_coords = SVector{3}(gp[i], gp[j], gp[k])
      
      # Calculate shape functions and their derivatives using thread-local variables
      compute_hex8_shape!(local_N, local_dN, local_coords[1], local_coords[2], local_coords[3])
      
      # Compute Jacobian using optimized matrix multiplication
      J = xe * local_dN
      
      # Accumulate contribution to volume
      elem_volume += w[i] * w[j] * w[k] * abs(det(J))
    end
    
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
