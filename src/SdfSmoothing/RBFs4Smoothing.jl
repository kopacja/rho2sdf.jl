# Replace preallocated values with the maximum computed value while preserving sign, Float64 -> Float32
function process_vector(vec::Vector{Float64})
  vec_f32 = Float32.(vec)  # Convert input vector to Float32
  max_val = maximum(abs, filter(x -> abs(x) < 1.0f9, vec_f32))
  if isempty(max_val)
    return vec_f32
  end
  return @. ifelse(abs(vec_f32) ≈ 1.0f10, sign(vec_f32) * max_val, vec_f32)
end

# Create a grid based on the number of cells and grid boundaries
function create_grid(NoC::Vector, grid)
  nx, ny, nz = NoC .+ 1 # number of grid points for x, y, z
  xmin, ymin, zmin = Float32.(grid.AABB_min)
  xmax, ymax, zmax = Float32.(grid.AABB_max)

  x = range(Float32(xmin), Float32(xmax), length=nx)
  y = range(Float32(ymin), Float32(ymax), length=ny)
  z = range(Float32(zmin), Float32(zmax), length=nz)

  return [[Float32(x), Float32(y), Float32(z)] for x in x, y in y, z in z]
end

function create_smooth_grid(grid, smooth::Int)
  nx, ny, nz = (grid.N * smooth) .+ 1 # number of grid points for x, y, z
  xmin, ymin, zmin = Float32.(grid.AABB_min)
  xmax, _, _ = Float32.(grid.AABB_max)

  # Compute steps explicitly
  dx = (xmax - xmin) / (nx - 1)

  # Generate grid points using explicit arithmetic
  x = [xmin + (i - 1) * dx for i in 1:nx]
  y = [ymin + (j - 1) * dx for j in 1:ny]
  z = [zmin + (k - 1) * dx for k in 1:nz]

  return (dx, [[x_i, y_j, z_k] for x_i in x, y_j in y, z_k in z])
end

# Define a struct for the Limited Range RBF Kernel
struct LimitedRangeRBFKernel{T<:Real} <: Kernel
  σ::T
  threshold::T
end

# Define the limited range Gaussian RBF kernel function
function (k::LimitedRangeRBFKernel)(x::AbstractVector{Float32}, y::AbstractVector{Float32})
  r = sqrt((x[1] - y[1])^2 + (x[2] - y[2])^2 + (x[3] - y[3])^2)
  value = exp(-(r / k.σ)^2)
  return value > k.threshold ? value : 0.0
end

# Convert vector to 3D array for kernel operations
function vector_to_array(dist::Vector, N)
  nx, ny, nz = N
  result = Array{eltype(dist),3}(undef, nx, ny, nz)
  @views result[:] .= dist
  return result
end

# Compute sparse kernel matrix for efficient operations
function compute_sparse_kernel_matrix(grid, kernel::LimitedRangeRBFKernel)
  n = length(grid)
  I, J, V = Int[], Int[], Float32[]

  for i in 1:n
    for j in i:n  # Utilize symmetry, compute only upper triangle
      value = Float32(kernel(grid[i], grid[j]))  # Convert to Float32
      if value > kernel.threshold
        push!(I, i)
        push!(J, j)
        push!(V, value)
        if i != j
          push!(I, j)
          push!(J, i)
          push!(V, value)
        end
      end
    end
  end
  return sparse(I, J, V, n, n)
end

# Compute RBF weights on the coarse grid
function compute_rbf_weights(coarse_grid, coarse_values, kernel::LimitedRangeRBFKernel)
  # Compute sparse kernel matrix for the coarse grid
  K = compute_sparse_kernel_matrix(coarse_grid, kernel)

  # Convert coarse_values to Float32
  coarse_values_f32 = Float32.(vec(coarse_values))

  # Compute weights using conjugate gradient method
  weights = cg(K, coarse_values_f32)

  return weights
end

# Perform RBF interpolation on the fine grid using precomputed weights and KD-tree for efficiency
function rbf_interpolation_kdtree(fine_grid::Array, coarse_grid::Array, weights_f32::Vector{Float32}, kernel::LimitedRangeRBFKernel)
  kdtree = KDTree(reduce(hcat, coarse_grid))
  max_distance = Float32(sqrt(-log(kernel.threshold) * kernel.σ^2))

  result = zeros(Float32, length(fine_grid))
  progress = Atomic{Int}(0)
  total = length(fine_grid)

  # Function to update and display progress bar
  function update_progress()
    new_progress = atomic_add!(progress, 1)
    if new_progress % 100 == 0 || new_progress == total
      percent = round(new_progress / total * 100, digits=1)
      bar = "="^Int(floor(percent / 2))
      @printf("\rProgress: [%-50s] %.1f%% (%d/%d)", bar, percent, new_progress, total)
    end
  end

  Threads.@threads for i in eachindex(fine_grid)
    idxs, dists = knn(kdtree, fine_grid[i], 124, true)
    for (j, dist) in zip(idxs, dists)
      if dist <= max_distance
        result[i] += weights_f32[j] * exp(-(dist / kernel.σ)^2)
      end
    end
    update_progress()
  end
  println()  # New line after completion
  return result
end

# Adjust the level-set function to maintain the volume fraction (V_frac)
function LS_Thresholdd(sdf::Array, grid::Array, mesh::Mesh, Exp::Int)
  # Calculate target volume from mesh properties
  target_volume = mesh.V_frac * mesh.V_domain
    
  # Initialize threshold search variables
  eps = 1.0
  th_low, th_high = minimum(sdf), maximum(sdf)
  n = 0
  th = 0.0
    
  # Setup arrays for volume calculation with calculate_volume_from_sdf
  dim = size(sdf)
  shifted_sdf = Array{Float32,3}(undef, dim)
    
  # Binary search to find the threshold that gives the desired volume
  while n < 40 && eps > 10.0^(-Exp)
    th = (th_low + th_high) / 2
    
    # Convert the level set to SDF format for volume calculation
    shifted_sdf .= sdf .- th
   
    # Calculate current volume using the SDF function
    current_volume = calculate_volume_from_sdf(shifted_sdf, grid)

    println(current_volume)
      
    # Adjust threshold based on calculated volume
    eps = abs(target_volume - current_volume)
    if current_volume > target_volume
      th_low = th
    else
      th_high = th
    end
    n += 1
  end
    
  return -th
end


# Convert 3D array back to vector
function array_to_vector(arr::Array{T,3}) where {T}
  nx, ny, nz = size(arr)
  result = Vector{T}(undef, nx * ny * nz)

  for k in 1:nz, j in 1:ny, i in 1:nx
    @inbounds result[(k-1)*ny*nx+(i-1)*nx+j] = arr[i, j, k]
  end

  return result
end

# Main function to perform RBF interpolation or approximation:
function RBFs_smoothing(
  mesh::Mesh,
  dist::Vector,
  my_grid::Grid,
  Is_interpolation::Bool,
  smooth::Int,
  taskName::String,
  threshold::Float64=1e-3
)
  name = Is_interpolation ? "Interpolation" : "Approximation"

  # SDF - Replace preallocated values with the maximum computed distance:
  dist_modif = process_vector(dist)

  V_frac = sum(dist .>= 0.0) / length(dist) # Volume fraction

  # Create coarse grid and assign values
  coarse_grid = create_grid(my_grid.N, my_grid)

  raw_SDF = vector_to_array(dist_modif, my_grid.N .+ 1) # data modification: vector -> array

  # Fine grid:
  dim = (my_grid.N * smooth) .+ 1 # dimensions of the fine grid

  (step, fine_grid) = create_smooth_grid(my_grid, smooth) # Create fine grid
  println("Original grid size: ", my_grid.N)
  println("Fine grid size: ", dim)

  # Define Gaussian kernel using LimitedRangeRBFKernel
  σ = my_grid.cell_size # width of the Gaussian function
  kernel = LimitedRangeRBFKernel(σ, threshold)

  # Compute weights on the coarse grid
  println("Computing weights on the coarse grid...")
  weights = Is_interpolation ?
            (@time compute_rbf_weights(coarse_grid, raw_SDF, kernel)) :
            vec(raw_SDF)

  println("Computing the LSF zero level to meet the volume condition...")
  @time LSF = rbf_interpolation_kdtree(coarse_grid, coarse_grid, weights, kernel)
  LSF_array = vector_to_array(LSF, my_grid.N .+ 1) # data modification: vector -> array
  th = LS_Thresholdd(LSF_array, coarse_grid, mesh, 4)

  println("Computing $name on the fine grid...")
  @time fine_LSF = rbf_interpolation_kdtree(fine_grid, coarse_grid, weights, kernel)

  # Shifting LSF to maintain volume
  fine_LSF_offset = fine_LSF .+ th

  B = round(my_grid.cell_size, digits=4)
  # Rho2sdf.exportStructuredPointsToVTK(taskName * "_smooth_B-" * string(B) * "_" * "_smooth-"$(smooth) * name * ".vtk", my_grid, fine_LSF_offset, "distance", smooth)
  Rho2sdf.exportStructuredPointsToVTK("$(taskName)_B-$(B)_smooth-$(smooth)_$(name).vtk", my_grid, fine_LSF_offset, "distance", smooth)

  fine_LSF_offset_array = vector_to_array(fine_LSF_offset, dim)

  fine_sdf = fine_LSF_offset_array

  @save "Z_$(taskName)_FineSDF_B-$(B)_smooth-$(smooth)_$(name).jld2" fine_sdf
  @save "Z_$(taskName)_FineGrid_B-$(B)_smooth-$(smooth)_$(name).jld2" fine_grid

  current_volume = calculate_volume_from_sdf(fine_sdf, fine_grid)
  @info "Body volume at SDF zero level: $current_volume (target: $(round(mesh.V_frac * mesh.V_domain, digits=4)))"

  return fine_sdf, fine_grid
end
