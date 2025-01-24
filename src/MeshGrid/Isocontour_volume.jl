
function C3D8_SFaD_nove(Ξ::SVector{3,Float64})
  # Rozbalíme souřadnice pro lepší čitelnost a možné SIMD optimalizace
  ξ₁, ξ₂, ξ₃ = Ξ

  # Předpočítáme často používané výrazy pro lepší výkon
  # Tyto hodnoty kompilátor může optimalizovat do SIMD instrukcí
  ξ₁m = (ξ₁ - 1)  # ξ₁ minus 1
  ξ₁p = (ξ₁ + 1)  # ξ₁ plus 1
  ξ₂m = (ξ₂ - 1)
  ξ₂p = (ξ₂ + 1)
  ξ₃m = (ξ₃ - 1)
  ξ₃p = (ξ₃ + 1)

  # Vytvoříme tvarové funkce jako SVector pro SIMD optimalizace
  N = SVector{8,Float64}(
    -0.125 * ξ₁m * ξ₂m * ξ₃m,  # N[1]
    0.125 * ξ₁p * ξ₂m * ξ₃m,  # N[2]
    -0.125 * ξ₁p * ξ₂p * ξ₃m,  # N[3]
    0.125 * ξ₁m * ξ₂p * ξ₃m,  # N[4]
    0.125 * ξ₁m * ξ₂m * ξ₃p,  # N[5]
    -0.125 * ξ₁p * ξ₂m * ξ₃p,  # N[6]
    0.125 * ξ₁p * ξ₂p * ξ₃p,  # N[7]
    -0.125 * ξ₁m * ξ₂p * ξ₃p   # N[8]
  )

  # Vytvoříme matici derivací jako statickou matici 8×3
  # Používáme konstantu 0.125 pro lepší čitelnost
  d¹N_dξ¹ = @SMatrix [
    -0.125*ξ₂m*ξ₃m -0.125*ξ₁m*ξ₃m -0.125*ξ₁m*ξ₂m;
    0.125*ξ₂m*ξ₃m 0.125*ξ₁p*ξ₃m 0.125*ξ₁p*ξ₂m;
    -0.125*ξ₂p*ξ₃m -0.125*ξ₁p*ξ₃m -0.125*ξ₁p*ξ₂p;
    0.125*ξ₂p*ξ₃m 0.125*ξ₁m*ξ₃m 0.125*ξ₁m*ξ₂p;
    0.125*ξ₂m*ξ₃p 0.125*ξ₁m*ξ₃p 0.125*ξ₁m*ξ₂m;
    -0.125*ξ₂m*ξ₃p -0.125*ξ₁p*ξ₃p -0.125*ξ₁p*ξ₂m;
    0.125*ξ₂p*ξ₃p 0.125*ξ₁p*ξ₃p 0.125*ξ₁p*ξ₂p;
    -0.125*ξ₂p*ξ₃p -0.125*ξ₁m*ξ₃p -0.125*ξ₁m*ξ₂p
  ]

  return N, d¹N_dξ¹
end

function calculate_isocontour_volume(mesh::Mesh,
  nodal_values::Vector{Float64},
  iso_threshold::Float64)
  X = mesh.X
  IEN = mesh.IEN

  # Předpočítáme body a váhy pro kvadraturu jako statická pole
  gp_detailed, w_detailed = gausslegendre(25)
  gp_detailed = SVector{25}(gp_detailed)
  w_detailed = SVector{25}(w_detailed)

  gp_simple, w_simple = gausslegendre(3)
  gp_simple = SVector{3}(gp_simple)
  w_simple = SVector{3}(w_simple)

  num_elements = size(IEN, 2)
  total_volume = Atomic{Float64}(0.0)

  # Paralelní zpracování elementů
  @threads for elem in 1:num_elements
    # Získáme hodnoty pro uzly elementu jako statický vektor
    elem_values = SVector{8}([nodal_values[IEN[i, elem]] for i in 1:8])

    min_value = minimum(elem_values)
    max_value = maximum(elem_values)

    max_value < iso_threshold && continue

    # Souřadnice uzlů uložíme jako 3×8 statickou matici
    xe = @SMatrix [X[i, IEN[j, elem]] for i in 1:3, j in 1:8]

    # Zvolíme strategii integrace
    if min_value >= iso_threshold
      gp, w = gp_simple, w_simple
      ngp = 3
      check_threshold = false
    else
      gp, w = gp_detailed, w_detailed
      ngp = 25
      check_threshold = true
    end

    element_volume = 0.0

    for k in 1:ngp, j in 1:ngp, i in 1:ngp
      # Vytvoříme lokální souřadnice jako statický vektor
      local_coords = SVector{3}(gp[i], gp[j], gp[k])

      # Získáme tvarové funkce a jejich derivace
      N, dN = C3D8_SFaD_nove(local_coords)

      if check_threshold
        interpolated_value = dot(N, elem_values)
        interpolated_value < iso_threshold && continue
      end

      # Opravený výpočet Jakobiánu
      # xe je 3×8 matice, dN je 8×3 matice, výsledek bude 3×3 matice
      J = xe * dN  # Násobíme matici 3×8 s maticí 8×3, dostaneme 3×3

      # Přidáme příspěvek k objemu elementu
      element_volume += w[i] * w[j] * w[k] * abs(det(J))
    end

    atomic_add!(total_volume, element_volume)
  end

  return total_volume[]
end

function find_threshold_for_volume(mesh::Mesh,
  nodal_values::Vector{Float64},
  tolerance::Float64=1e-4,
  max_iterations::Int=60)

  target_volume = mesh.V_domain*mesh.V_frac
  # Inicializace hraničních hodnot pro binary search
  # Víme, že hodnoty jsou mezi 0 a 1, takže to jsou naše počáteční meze
  lower_bound = 0.0
  upper_bound = 1.0

  # Vypočítáme počáteční objemy pro kontrolu řešitelnosti
  min_volume = calculate_isocontour_volume(mesh, nodal_values, upper_bound)
  max_volume = calculate_isocontour_volume(mesh, nodal_values, lower_bound)

  # Kontrola, zda je požadovaný objem v možném rozsahu
  if target_volume > max_volume || target_volume < min_volume
    error("Požadovaný objem $(target_volume) je mimo možný rozsah [$(min_volume), $(max_volume)]")
  end

  # Proměnné pro sledování průběhu
  current_iteration = 0
  best_threshold = 0.0
  best_volume_error = Inf

  println("Hledání prahové hodnoty pro cílový objem: $(target_volume)")
  println("Iterace | Práh | Objem | Odchylka")
  println("-"^50)

  while current_iteration < max_iterations
    # Výpočet středové hodnoty intervalu
    threshold = (lower_bound + upper_bound) / 2

    # Výpočet objemu pro aktuální práh
    current_volume = calculate_isocontour_volume(mesh, nodal_values, threshold)

    # Výpočet relativní chyby
    volume_error = abs(current_volume - target_volume) / target_volume

    # Výpis průběhu
    println(@sprintf("%3d | %.4f | %.4f | %.4e",
      current_iteration, threshold, current_volume, volume_error))

    # Aktualizace nejlepšího nalezeného řešení
    if volume_error < best_volume_error
      best_threshold = threshold
      best_volume_error = volume_error
    end

    # Kontrola konvergence
    if volume_error < tolerance
      println("\nŘešení nalezeno!")
      break
    end

    # Úprava mezí intervalu podle výsledku
    if current_volume > target_volume
      lower_bound = threshold
    else
      upper_bound = threshold
    end

    current_iteration += 1
  end

  # Kontrola, zda jsme dosáhli maximálního počtu iterací
  if current_iteration == max_iterations
    println("\nDosažen maximální počet iterací!")
    println("Vracím nejlepší nalezené řešení.")
  end

  # Výpis finálního výsledku
  final_volume = calculate_isocontour_volume(mesh, nodal_values, best_threshold)
  println("\nVýsledek:")
  println("Nalezená prahová hodnota: $(best_threshold)")
  println("Dosažený objem: $(final_volume)")
  println("Relativní chyba: $(best_volume_error)")

  return best_threshold
end
