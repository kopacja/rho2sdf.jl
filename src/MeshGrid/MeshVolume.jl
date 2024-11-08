"""
Calculate volume of 8-node hexahedral elements using Gaussian quadrature.
Inputs:
- X: 3×n matrix of nodal coordinates
- IEN: 8×m connectivity matrix
Returns:
- total volume of the mesh
"""
function calculate_mesh_volume(X::Vector{Vector{Float64}}, IEN::Vector{Vector{Int64}}, rho::Vector{Float64})

  @info "Compute volume"
  # Gauss quadrature points and weights for hexahedron
  # Using 2×2×2 integration points
  gp = [-1 / √3, 1 / √3]
  w = [1.0, 1.0]
  # Using 3×3×3 integration points
  # gp = [-√(3 / 5), 0.0, √(3 / 5)]
  # w = [5 / 9, 8 / 9, 5 / 9]

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
      # local_coords = [ξ, η, ζ]

      # Shape functions derivatives
      dN = shape_function_derivatives(ξ, η, ζ)
      # _, dN, _, _ = C3D8_SFaD(local_coords)
      # Jacobian matrix
      J = zeros(3, 3)
      for n in 1:8
        J += xe[:, n] * dN[:, n]'
        # J += xe[:, n] * dN[:, n]
      end

      # Add contribution to element volume
      elem_volume += w[i] * w[j] * w[k] * abs(det(J))
    end

    domain_volume += elem_volume
    TO_volume += elem_volume * rho[elem]
  end

  println("Topology optimization domain volume: ", domain_volume)
  println("Optimized shape volume: ", TO_volume)
  println("Volume fraction: ", (TO_volume / domain_volume))
end

"""
Calculate shape function derivatives for 8-node hexahedron at given point.
Returns 3×8 matrix of shape function derivatives.
"""
function shape_function_derivatives(ξ::Float64, η::Float64, ζ::Float64)
  dN = zeros(3, 8)

  # dN/dξ
  dN[1, 1] = -0.125 * (1 - η) * (1 - ζ)
  dN[1, 2] = 0.125 * (1 - η) * (1 - ζ)
  dN[1, 3] = 0.125 * (1 + η) * (1 - ζ)
  dN[1, 4] = -0.125 * (1 + η) * (1 - ζ)
  dN[1, 5] = -0.125 * (1 - η) * (1 + ζ)
  dN[1, 6] = 0.125 * (1 - η) * (1 + ζ)
  dN[1, 7] = 0.125 * (1 + η) * (1 + ζ)
  dN[1, 8] = -0.125 * (1 + η) * (1 + ζ)

  # dN/dη
  dN[2, 1] = -0.125 * (1 - ξ) * (1 - ζ)
  dN[2, 2] = -0.125 * (1 + ξ) * (1 - ζ)
  dN[2, 3] = 0.125 * (1 + ξ) * (1 - ζ)
  dN[2, 4] = 0.125 * (1 - ξ) * (1 - ζ)
  dN[2, 5] = -0.125 * (1 - ξ) * (1 + ζ)
  dN[2, 6] = -0.125 * (1 + ξ) * (1 + ζ)
  dN[2, 7] = 0.125 * (1 + ξ) * (1 + ζ)
  dN[2, 8] = 0.125 * (1 - ξ) * (1 + ζ)

  # dN/dζ
  dN[3, 1] = -0.125 * (1 - ξ) * (1 - η)
  dN[3, 2] = -0.125 * (1 + ξ) * (1 - η)
  dN[3, 3] = -0.125 * (1 + ξ) * (1 + η)
  dN[3, 4] = -0.125 * (1 - ξ) * (1 + η)
  dN[3, 5] = 0.125 * (1 - ξ) * (1 - η)
  dN[3, 6] = 0.125 * (1 + ξ) * (1 - η)
  dN[3, 7] = 0.125 * (1 + ξ) * (1 + η)
  dN[3, 8] = 0.125 * (1 - ξ) * (1 + η)

  return dN
end

# calculate_mesh_volume(mesh.X, mesh.IEN, rho)
