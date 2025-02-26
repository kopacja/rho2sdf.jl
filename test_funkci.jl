using Rho2sdf
using Rho2sdf.TerminalUtils
using Rho2sdf.PrimitiveGeometries
using Rho2sdf.ShapeFunctions
using Rho2sdf.MeshGrid
using Rho2sdf.SignedDistances
using Rho2sdf.DataExport
using Rho2sdf.SdfSmoothing
using MAT
using JLD2
using LinearAlgebra
using BenchmarkTools

using StaticArrays
using Base.Threads
using FastGaussQuadrature  # Pro funkci gausslegendre

function calculate_mesh_volume_old(X::Vector{Vector{Float64}}, IEN::Vector{Vector{Int64}}, rho::Vector{Float64})

  print_info("Computing volume...")
  # Gauss quadrature points and weights for hexahedron
  # Using 3×3×3 integration points
  gp = [-√(3 / 5), 0.0, √(3 / 5)]
  w = [5 / 9, 8 / 9, 5 / 9]

  ngp = length(gp)

  N = MVector{8,Float64}(undef)
    dN = MMatrix{8,3,Float64}(undef)

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
      compute_hex8_shape!(N, dN, local_coords[1], local_coords[2], local_coords[3])

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

function calculate_mesh_volume(
    X::Vector{Vector{Float64}},
    IEN::Vector{Vector{Int64}},
    rho::Vector{Float64};
    quad_order::Int=3  # Výchozí řád integrace je 3, stejný jako v původní funkci
  )
    print_info("Computing volume...")
  
    # Výpočet Gauss-Legendreových bodů a vah pro zadaný řád integrace
    gp_raw, w_raw = gausslegendre(quad_order)
  
    # Převod na statické vektory pro lepší výkon
    gp = SVector{quad_order}(gp_raw)
    w = SVector{quad_order}(w_raw)
  
    # Pre-alokace polí pro tvarové funkce a jejich derivace
    N = MVector{8,Float64}(undef)
    dN = MMatrix{8,3,Float64}(undef)
  
    # Atomické proměnné pro thread-safe akumulaci objemů
    domain_volume = Atomic{Float64}(0.0)
    TO_volume = Atomic{Float64}(0.0)
  
    # Paralelní zpracování elementů
    @threads for elem in 1:length(IEN)
      # Konverze souřadnic elementu na statickou matici
      xe = @SMatrix [X[IEN[elem][j]][i] for i in 1:3, j in 1:8]
  
      # Výpočet objemu elementu pomocí Gaussovy kvadratury
      elem_volume = 0.0
      for k in 1:quad_order, j in 1:quad_order, i in 1:quad_order
        local_coords = SVector{3}(gp[i], gp[j], gp[k])
  
        # Výpočet tvarových funkcí a jejich derivací
        compute_hex8_shape!(N, dN, local_coords[1], local_coords[2], local_coords[3])
  
        # Výpočet Jakobiánu pomocí optimalizovaného maticového násobení
        J = xe * dN
  
        # Akumulace příspěvku k objemu
        elem_volume += w[i] * w[j] * w[k] * abs(det(J))
      end
  
      # Thread-safe aktualizace objemů
      atomic_add!(domain_volume, elem_volume)
      atomic_add!(TO_volume, elem_volume * rho[elem])
    end
  
    # Získání finálních hodnot objemů
    final_domain_volume = domain_volume[]
    final_TO_volume = TO_volume[]
  
    # Výpis výsledků
    println("Topology optimization domain volume: ", round(final_domain_volume, sigdigits=6))
    print_data("Optimized shape volume: $(round(final_domain_volume, sigdigits=6))")
    println("Volume fraction: ", round(final_TO_volume / final_domain_volume, sigdigits=6))
  
    return [final_domain_volume, (final_TO_volume / final_domain_volume)]
  end
  
  taskName = "sphere"

      N = 10  # Number of cells along the longest side
      # ρₜ = 0.5 # Threshold density (isosurface level)
      ## Read FEM mesh:
    #   data = matread(taskName * ".mat")
      data = matread("test/" * taskName * ".mat")
      (X, IEN, rho) = MeshGrid.MeshInformations(data)

      calculate_mesh_volume_old(X, IEN, rho)
calculate_mesh_volume(X, IEN, rho)