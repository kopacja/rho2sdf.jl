using Rho2sdf
using Rho2sdf.ShapeFunctions
using Rho2sdf.MeshGrid
using Rho2sdf.SignedDistances
using Rho2sdf.DataExport
using Rho2sdf.SdfSmoothing
using Rho2sdf.PrimitiveGeometries
using MAT
using JLD2
using LinearAlgebra

taskName = "sphere"
N = 10  # Number of cells along the longest side
ρₜ = 0.5 # Threshold density (isosurface level)

# data = matread(taskName * ".mat")
      data = matread("test/" * taskName * ".mat")
      (X, IEN, rho) = MeshGrid.MeshInformations(data)

      ## Generate FEM mesh structure:
      println("Type of X: ", typeof(X))
      println("Type of IEN: ", typeof(IEN))
      println("Type of rho: ", typeof(rho))
      println("Type of C3D8_SFaD: ", typeof(C3D8_SFaD))
      mesh = MeshGrid.Mesh(X, IEN, rho, C3D8_SFaD)

      ## Map elemental densities to the nodes:
      ρₙ = MeshGrid.DenseInNodes(mesh, rho) # LSQ
      #ρₙ = MeshGrid.elementToNodalValues(mesh, rho) # average

#_____________
N = 10  # Number of cells along the longest side
      ρₜ = 0.5 # Threshold density (isosurface level)

      (X, IEN, rho) = PrimitiveGeometries.selectPrimitiveGeometry("block", [2, 1, 1])
      # ρₙ = [0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.5, 0.5]

      mesh = MeshGrid.Mesh(X, IEN, rho, C3D8_SFaD)
      # ρₙ = MeshGrid.DenseInNodes(mesh, rho) # LSQ

      # Modif ρₙ:
      ρₙ = [0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 0.0, 0.0, 0.5, 0.5]



function calculate_isocontour_volume(X::Vector{Vector{Float64}}, 
                                IEN::Vector{Vector{Int64}}, 
                                nodal_values::Vector{Float64},
                                iso_threshold::Float64)
    # Gauss quadrature points and weights (3×3×3 for better accuracy)
    gp = [-√(3/5), 0.0, √(3/5)]
    w = [5/9, 8/9, 5/9]
    # gp = [-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526]
    # w = [0.34785484513745385, 0.6521451548625462, 0.6521451548625462, 0.34785484513745385]
    ngp = length(gp)

    # Initialize volume counter
    total_volume = 0.0

    # Loop through all elements
    for elem in 1:length(IEN)
        # Get nodal coordinates and values for current element
        xe = zeros(3, 8)
        elem_values = zeros(8)
        for (i, node_idx) in enumerate(IEN[elem])
            xe[:, i] = X[node_idx]
            elem_values[i] = nodal_values[node_idx]
        end

        # Calculate element volume using Gaussian quadrature
        for k in 1:ngp, j in 1:ngp, i in 1:ngp
            ξ, η, ζ = gp[i], gp[j], gp[k]
            local_coords = [ξ, η, ζ]
            
            # Get shape functions and their derivatives
            N, dN, _, _ = C3D8_SFaD(local_coords)
            # N, dN = shape_functions(local_coords)
            
            # Calculate interpolated value at integration point
            interpolated_value = sum(N .* elem_values)
            
            # Jacobian matrix
            J = zeros(3, 3)
            for n in 1:8
                J += xe[:, n] * dN[n, :]'
            end
            
            # Add contribution to volume only if above threshold
            if interpolated_value >= iso_threshold
                total_volume += w[i] * w[j] * w[k] * abs(det(J))
            end
        end
    end
    
    return total_volume
end
ρₜ = 0.
calculate_isocontour_volume(X, IEN, ρₙ, ρₜ)


function calculate_mesh_volume(X::Vector{Vector{Float64}}, IEN::Vector{Vector{Int64}}, rho::Vector{Float64})

    @info "Compute volume"
    # Gauss quadrature points and weights for hexahedron
    # Using 2×2×2 integration points
    # gp = [-1 / √3, 1 / √3]
    # w = [1.0, 1.0]
    gp = [-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526]
    w = [0.34785484513745385, 0.6521451548625462, 0.6521451548625462, 0.34785484513745385]
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
  
    println("Topology optimization domain volume: ", domain_volume)
    println("Optimized shape volume: ", TO_volume)
    println("Volume fraction: ", (TO_volume / domain_volume))
  
    return [domain_volume, (TO_volume / domain_volume)]
  end
  domain_volume, target_volume = calculate_mesh_volume(X, IEN, ρₙ)


  function calculate_isocontour_volume(mesh::Mesh,
                                nodal_values::Vector{Float64},
                                iso_threshold::Float64)
    X = mesh.X
    IEN = mesh.IEN
    # Gauss quadrature points and weights (3×3×3 for better accuracy)
    gp = [-√(3/5), 0.0, √(3/5)]
    w = [5/9, 8/9, 5/9]
    ngp = length(gp)

    # Initialize volume counter
    total_volume = 0.0

    # Get number of elements from IEN matrix dimensions
    num_elements = size(IEN, 2)  # IEN je 8×1000, takže počet elementů je ve druhém rozměru

    # Loop through all elements
    for elem in 1:num_elements
        # Get nodal coordinates and values for current element
        xe = zeros(3, 8)
        elem_values = zeros(8)
        for i in 1:8  # Pro každý uzel elementu
            node_idx = IEN[i, elem]  # Získáme globální index uzlu z IEN matice
            xe[:, i] = X[:, node_idx]  # Získáme souřadnice uzlu z X matice
            elem_values[i] = nodal_values[node_idx]
        end

        # Calculate element volume using Gaussian quadrature
        for k in 1:ngp, j in 1:ngp, i in 1:ngp
            ξ, η, ζ = gp[i], gp[j], gp[k]
            local_coords = [ξ, η, ζ]
            
            # Get shape functions and their derivatives
            N, dN, _, _ = C3D8_SFaD(local_coords)
            
            # Calculate interpolated value at integration point
            interpolated_value = sum(N .* elem_values)
            
            # Jacobian matrix
            J = zeros(3, 3)
            for n in 1:8
                J += xe[:, n] * dN[n, :]'
            end
            
            # Add contribution to volume only if above threshold
            if interpolated_value >= iso_threshold
                total_volume += w[i] * w[j] * w[k] * abs(det(J))
            end
        end
    end
    
    return total_volume
end
calculate_isocontour_volume(mesh, ρₙ, 0.5)


###__________________________

function find_threshold_for_volume(mesh::Mesh,
                              nodal_values::Vector{Float64},
                              target_volume::Float64;
                              tolerance::Float64 = 1e-4,
                              max_iterations::Int = 100)
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
    println("-" ^ 50)
    
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
    
    return best_threshold, final_volume, best_volume_error
end
using Printf
find_threshold_for_volume(mesh, ρₙ, target_volume)