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
  
  # Atomické proměnné pro thread-safe akumulaci objemů
  domain_volume = Atomic{Float64}(0.0)
  TO_volume = Atomic{Float64}(0.0)
  
  # Paralelní zpracování elementů
  @threads for elem in 1:length(IEN)
    # Pre-alokace polí pro tvarové funkce a jejich derivace - přesunuté dovnitř paralelní smyčky
    local_N = MVector{8,Float64}(undef)
    local_dN = MMatrix{8,3,Float64}(undef)
    
    # Konverze souřadnic elementu na statickou matici
    xe = @SMatrix [X[IEN[elem][j]][i] for i in 1:3, j in 1:8]
    
    # Výpočet objemu elementu pomocí Gaussovy kvadratury
    elem_volume = 0.0
    for k in 1:quad_order, j in 1:quad_order, i in 1:quad_order
      local_coords = SVector{3}(gp[i], gp[j], gp[k])
      
      # Výpočet tvarových funkcí a jejich derivací s thread-lokálními proměnnými
      compute_hex8_shape!(local_N, local_dN, local_coords[1], local_coords[2], local_coords[3])
      
      # Výpočet Jakobiánu pomocí optimalizovaného maticového násobení
      J = xe * local_dN
      
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
  print_data("Optimized shape volume: $(round(final_TO_volume, sigdigits=6))")
  println("Volume fraction: ", round(final_TO_volume / final_domain_volume, sigdigits=6))
  
  return [final_domain_volume, (final_TO_volume / final_domain_volume)]
end
