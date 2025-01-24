"""
Calculate volume of 8-node hexahedral elements using Gaussian quadrature.
Inputs:
- X: 3×n matrix of nodal coordinates
- IEN: 8×m connectivity matrix
Returns:
- total volume of the mesh
using FastGaussQuadrature
gp, w = gausslegendre(2)
"""
function calculate_mesh_volume(X::Vector{Vector{Float64}}, IEN::Vector{Vector{Int64}}, rho::Vector{Float64})

  print_info("Computing volume...")
  # Gauss quadrature points and weights for hexahedron
  # Using 3×3×3 integration points
  gp = [-√(3 / 5), 0.0, √(3 / 5)]
  w = [5 / 9, 8 / 9, 5 / 9]

  ngp = length(gp)

  # Initialize total volume
  domain_volume = 0.0
  TO_volume = 0.0

  # Loop through all elements
  for elem in 1:length(IEN)
    # Get nodal coordinates for current element
    xe = zeros(3, 8)
    for (i, node_idx) in enumerate(IEN[elem])
      xe[:, i] = X[node_idx]
    end

    # Calculate element volume using Gaussian quadrature
    elem_volume = 0.0

    for k in 1:ngp, j in 1:ngp, i in 1:ngp
      ξ, η, ζ = gp[i], gp[j], gp[k]
      local_coords = [ξ, η, ζ]

      # Shape functions derivatives
      _, dN, _, _ = C3D8_SFaD(local_coords)

      # Jacobian matrix
      J = zeros(3, 3)
      for n in 1:8
        J += xe[:, n] * dN[n, :]'
      end

      # Add contribution to element volume
      elem_volume += w[i] * w[j] * w[k] * abs(det(J))
    end

    domain_volume += elem_volume
    TO_volume += elem_volume * rho[elem]
  end

  println("Topology optimization domain volume: ", round(domain_volume, sigdigits=6))
  print_data("Optimized shape volume: $(round(domain_volume, sigdigits=6))")
  println("Volume fraction: ", round(TO_volume / domain_volume, sigdigits=6))

  return [domain_volume, (TO_volume / domain_volume)]
end

