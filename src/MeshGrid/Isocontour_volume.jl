function calculate_isocontour_volume(
  mesh::Mesh,
  nodal_values::Vector{Float64},
  iso_threshold::Float64;
  detailed_quad_order::Int=15,
  simple_quad_order::Int=3
)
  X = mesh.X
  IEN = mesh.IEN

  # Compute Gauss-Legendre points and weights based on input parameters
  gp_detailed, w_detailed = gausslegendre(detailed_quad_order)
  gp_detailed = SVector{detailed_quad_order}(gp_detailed)
  w_detailed = SVector{detailed_quad_order}(w_detailed)

  gp_simple, w_simple = gausslegendre(simple_quad_order)
  gp_simple = SVector{simple_quad_order}(gp_simple)
  w_simple = SVector{simple_quad_order}(w_simple)

  num_elements = size(IEN, 2)
  total_volume = Atomic{Float64}(0.0)

  # Pre-allocate arrays for shape functions and their derivatives
  N = MVector{8,Float64}(undef)
  dN = MMatrix{8,3,Float64}(undef)

  # Parallel element processing
  @threads for elem in 1:num_elements
    # Get element nodal values as static vector
    elem_values = SVector{8}([nodal_values[IEN[i, elem]] for i in 1:8])
    min_value = minimum(elem_values)
    max_value = maximum(elem_values)

    max_value < iso_threshold && continue

    # Store node coordinates as 3×8 static matrix
    xe = @SMatrix [X[i, IEN[j, elem]] for i in 1:3, j in 1:8]

    # Choose integration strategy based on element values
    if min_value >= iso_threshold
      gp, w = gp_simple, w_simple
      ngp = simple_quad_order
      check_threshold = false
    else
      gp, w = gp_detailed, w_detailed
      ngp = detailed_quad_order
      check_threshold = true
    end

    element_volume = 0.0
    for k in 1:ngp, j in 1:ngp, i in 1:ngp
      # Create local coordinates as static vector
      local_coords = SVector{3}(gp[i], gp[j], gp[k])

      # Compute shape functions and their derivatives
      compute_hex8_shape!(N, dN, local_coords[1], local_coords[2], local_coords[3])

      if check_threshold
        interpolated_value = dot(N, elem_values)
        interpolated_value < iso_threshold && continue
      end

      # Calculate Jacobian (3×8 matrix × 8×3 matrix = 3×3 matrix)
      J = xe * dN

      # Add contribution to element volume
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
  target_volume = mesh.V_domain * mesh.V_frac

  # Initialize boundary values for binary search
  # We know values are between 0 and 1, so these are our initial bounds
  lower_bound = 0.0
  upper_bound = 1.0

  # Calculate initial volumes to check solution feasibility
  min_volume = calculate_isocontour_volume(mesh, nodal_values, upper_bound)
  max_volume = calculate_isocontour_volume(mesh, nodal_values, lower_bound)

  # Check if the requested volume is within the possible range
  if target_volume > max_volume || target_volume < min_volume
    error("Requested volume $(target_volume) is outside the possible range [$(min_volume), $(max_volume)]")
  end

  # Variables for tracking progress
  current_iteration = 0
  best_threshold = 0.0
  best_volume_error = Inf

  print_info("\nComputing volume fraction...")
  println("Iteration | Threshold | Volume | Error")
  println("-"^50)

  while current_iteration < max_iterations
    # Calculate midpoint of the interval
    threshold = (lower_bound + upper_bound) / 2

    # Calculate volume for current threshold
    current_volume = calculate_isocontour_volume(mesh, nodal_values, threshold)

    # Calculate relative error
    volume_error = abs(current_volume - target_volume) / target_volume

    # Print progress
    println(@sprintf("  %3d     | %.4f    | %.4f | %.4e",
      current_iteration, threshold, current_volume, volume_error))

    # Update best solution found
    if volume_error < best_volume_error
      best_threshold = threshold
      best_volume_error = volume_error
    end

    # Check convergence
    if volume_error < tolerance
      break
    end

    # Adjust interval bounds based on result
    if current_volume > target_volume
      lower_bound = threshold
    else
      upper_bound = threshold
    end

    current_iteration += 1
  end

  # Check if we reached maximum iterations
  if current_iteration == max_iterations
    println("\nMaximum number of iterations reached!")
    println("Returning best solution found.")
  end

  # Print final result
  final_volume = calculate_isocontour_volume(mesh, nodal_values, best_threshold)
  print_success("\nVolume fraction: $(round(best_threshold, sigdigits=6))")
  println("Achieved volume: $(round(final_volume, sigdigits=6))")
  print_data("Relative error: $(round(best_volume_error, sigdigits=6))")

  return best_threshold
end

