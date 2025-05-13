
"""
    process_vector(vec::Vector{Float64})

Convert a Vector{Float64} to Vector{Float32} while handling special values.
This function preserves the sign of large values (which are usually placeholders)
by replacing them with the maximum computed value in the dataset.

# Arguments
- `vec::Vector{Float64}`: Input vector with distance values

# Returns
- Vector{Float32}: Processed vector with special values handled
"""
function process_vector(vec::Vector{Float64})
    vec_f32 = Float32.(vec)  # Convert input vector to Float32
    max_val = maximum(abs, filter(x -> abs(x) < 1.0f9, vec_f32))
    if isempty(max_val)
        return vec_f32
    end
    return @. ifelse(abs(vec_f32) ≈ 1.0f10, sign(vec_f32) * max_val, vec_f32)
end

"""
    create_grid(NoC::Vector, grid)

Create a 3D grid based on the number of cells and grid boundaries.

# Arguments
- `NoC::Vector`: Number of cells in each dimension
- `grid`: Grid structure with AABB_min and AABB_max properties

# Returns
- 3D array of grid points where each element is a 3D vector [x, y, z]
"""
function create_grid(NoC::Vector, grid)
    nx, ny, nz = NoC .+ 1 # Number of grid points for x, y, z
    xmin, ymin, zmin = Float32.(grid.AABB_min)
    xmax, ymax, zmax = Float32.(grid.AABB_max)

    x = range(Float32(xmin), Float32(xmax), length=nx)
    y = range(Float32(ymin), Float32(ymax), length=ny)
    z = range(Float32(zmin), Float32(zmax), length=nz)

    return [[Float32(x), Float32(y), Float32(z)] for x in x, y in y, z in z]
end

"""
    create_smooth_grid(grid, smooth::Int)

Create a fine grid for smooth SDF representation with uniform cell size in all dimensions.

# Arguments
- `grid`: Original grid structure
- `smooth::Int`: Factor by which to increase the grid resolution

# Returns
- Tuple(Float32, 3D Array): Cell size and 3D array of grid points
"""
function create_smooth_grid(grid, smooth::Int)
    nx, ny, nz = (grid.N * smooth) .+ 1 # Number of grid points for x, y, z
    xmin, ymin, zmin = Float32.(grid.AABB_min)
    xmax, ymax, zmax = Float32.(grid.AABB_max)

    # Compute step explicitly (uniform in all dimensions)
    dx = (xmax - xmin) / (nx - 1)

    # Generate grid points using explicit arithmetic
    x = [xmin + (i - 1) * dx for i in 1:nx]
    y = [ymin + (j - 1) * dx for j in 1:ny]
    z = [zmin + (k - 1) * dx for k in 1:nz]

    return (dx, [[x_i, y_j, z_k] for x_i in x, y_j in y, z_k in z])
end

"""
    struct LimitedRangeRBFKernel{T<:Real} <: Kernel

A custom kernel for RBF interpolation that limits the range of influence.
This improves computational efficiency by creating sparse matrices.

# Fields
- `σ::T`: Width of the Gaussian function (controls smoothness)
- `threshold::T`: Cutoff value below which kernel outputs are considered zero
"""
struct LimitedRangeRBFKernel{T<:Real} <: Kernel
    σ::T
    threshold::T
end

"""
    (k::LimitedRangeRBFKernel)(x, y)

Implement the kernel function as a callable object. Computes the Gaussian RBF
between points x and y, returning 0 if below the threshold.

# Arguments
- `x::AbstractVector{Float32}`: First 3D point
- `y::AbstractVector{Float32}`: Second 3D point

# Returns
- Float32: Kernel value (0 if below threshold)
"""
function (k::LimitedRangeRBFKernel)(x::AbstractVector{Float32}, y::AbstractVector{Float32})
    r = sqrt((x[1] - y[1])^2 + (x[2] - y[2])^2 + (x[3] - y[3])^2)
    value = exp(-(r / k.σ)^2)
    return value > k.threshold ? value : 0.0
end

"""
    vector_to_array(dist::Vector, N)

Convert a flat vector to a 3D array based on grid dimensions.

# Arguments
- `dist::Vector`: Input vector of values
- `N`: Grid dimensions (nx, ny, nz)

# Returns
- 3D array with the shape specified by N
"""
function vector_to_array(dist::Vector, N)
    nx, ny, nz = N
    result = Array{eltype(dist),3}(undef, nx, ny, nz)
    @views result[:] .= dist
    return result
end

"""
    compute_sparse_kernel_matrix(grid, kernel::LimitedRangeRBFKernel)

Compute a sparse kernel matrix for efficient RBF operations.
Uses a KD-tree for spatial queries to find points within interaction radius.

# Arguments
- `grid`: Array of grid points
- `kernel::LimitedRangeRBFKernel`: Kernel function with limited range

# Returns
- SparseMatrixCSC: Sparse kernel matrix
"""
function compute_sparse_kernel_matrix(grid, kernel::LimitedRangeRBFKernel)
    # Convert grid to a matrix for KDTree
    points = reduce(hcat, grid)
    
    # Calculate interaction radius based on kernel parameters
    interaction_radius = kernel.σ * sqrt(-log(kernel.threshold))
    
    # Create KD-tree for efficient neighbor searches
    kdtree = KDTree(points)
    
    I, J, V = Int[], Int[], Float32[]
    
    for i in 1:size(points, 2)
        # Find all points within interaction_radius
        idxs = inrange(kdtree, points[:, i], interaction_radius)
        
        for j in idxs
            if j >= i  # Avoid duplicate calculations
                value = Float32(kernel(grid[i], grid[j]))
                if value > kernel.threshold
                    push!(I, i)
                    push!(J, j)
                    push!(V, value)
                    if i != j  # Account for symmetry
                        push!(I, j)
                        push!(J, i)
                        push!(V, value)
                    end
                end
            end
        end
    end
    
    return sparse(I, J, V, length(grid), length(grid))
end

"""
    compute_rbf_weights(coarse_grid, coarse_values, kernel::LimitedRangeRBFKernel)

Compute RBF weights for interpolation or approximation.

# Arguments
- `coarse_grid`: Array of grid points on coarse grid
- `coarse_values`: Values at coarse grid points
- `kernel::LimitedRangeRBFKernel`: Kernel function with limited range

# Returns
- Vector{Float32}: Computed weights for RBF interpolation
"""
function compute_rbf_weights(coarse_grid, coarse_values, kernel::LimitedRangeRBFKernel)
    # Compute sparse kernel matrix for the coarse grid
    K = compute_sparse_kernel_matrix(coarse_grid, kernel)

    # Convert coarse_values to Float32
    coarse_values_f32 = Float32.(vec(coarse_values))

    # Compute weights using conjugate gradient method
    weights = cg(K, coarse_values_f32)

    return weights
end

"""
    rbf_interpolation_kdtree(fine_grid, coarse_grid, weights_f32, kernel)

Perform RBF interpolation on a fine grid using precomputed weights.
Uses KD-tree for efficient spatial queries.

# Arguments
- `fine_grid::Array`: Target grid points for interpolation
- `coarse_grid::Array`: Source grid points with known weights
- `weights_f32::Vector{Float32}`: Precomputed RBF weights
- `kernel::LimitedRangeRBFKernel`: Kernel function with limited range

# Returns
- Vector{Float32}: Interpolated values at fine grid points
"""
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

"""
    LS_Threshold(sdf::Array, grid::Array, mesh::Mesh, Exp::Int)

Adjust the level-set function to maintain the volume of the original mesh.
Uses binary search to find the threshold that gives the desired volume.

# Arguments
- `sdf::Array`: Signed distance function values
- `grid::Array`: Grid points
- `mesh::Mesh`: Original mesh with known volume
- `Exp::Int`: Precision exponent for convergence criterion

# Returns
- Float32: Threshold value for the level set function
"""
function LS_Threshold(sdf::Array, grid::Array, mesh::Mesh, Exp::Int)
    # Calculate target volume from mesh properties
    target_volume = mesh.V_frac * mesh.V_domain
    
    # Initialize threshold search variables
    eps = 1.0
    th_low, th_high = minimum(sdf), maximum(sdf)
    n = 0
    th = 0.0
    
    # Setup arrays for volume calculation
    dim = size(sdf)
    shifted_sdf = Array{Float32,3}(undef, dim)
    
    # Binary search to find the threshold that gives the desired volume
    while n < 40 && eps > 10.0^(-Exp)
        th = (th_low + th_high) / 2
        
        # Convert the level set to SDF format for volume calculation
        shifted_sdf .= sdf .- th
        
        # Calculate current volume using the SDF function
        current_volume = calculate_volume_from_sdf(shifted_sdf, grid)
        
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

"""
    RBFs_smoothing(mesh::Mesh, dist::Vector, my_grid::Grid, Is_interpolation::Bool, 
                   smooth::Int, taskName::String, threshold::Float64=1e-3)

Main function to perform RBF smoothing of a signed distance function.
Supports both interpolation and approximation approaches.

# Arguments
- `mesh::Mesh`: Original mesh with known volume
- `dist::Vector`: Original distance values on coarse grid
- `my_grid::Grid`: Grid structure
- `Is_interpolation::Bool`: Whether to use interpolation (true) or approximation (false)
- `smooth::Int`: Factor by which to increase the grid resolution
- `taskName::String`: Name for output files
- `threshold::Float64=1e-3`: Threshold for the RBF kernel

# Returns
- Tuple(Array, Array): Fine SDF values and fine grid points
"""
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

    # Process the SDF values
    dist_modif = process_vector(dist)

    # Create coarse grid and assign values
    coarse_grid = create_grid(my_grid.N, my_grid)
    raw_SDF = vector_to_array(dist_modif, my_grid.N .+ 1)

    # Create fine grid
    dim = (my_grid.N * smooth) .+ 1
    (step, fine_grid) = create_smooth_grid(my_grid, smooth)
    println("Original grid size: ", my_grid.N)
    println("Fine grid size: ", dim)

    # Define Gaussian kernel with limited range
    σ = my_grid.cell_size
    kernel = LimitedRangeRBFKernel(σ, threshold)

    # Compute weights on the coarse grid
    println("Computing weights on the coarse grid...")
    weights = Is_interpolation ?
              (@time compute_rbf_weights(coarse_grid, raw_SDF, kernel)) :
              vec(raw_SDF)

    # Compute level set function threshold for volume preservation
    println("Computing the LSF zero level to meet the volume condition...")
    @time LSF = rbf_interpolation_kdtree(coarse_grid, coarse_grid, weights, kernel)
    LSF_array = vector_to_array(LSF, my_grid.N .+ 1)
    th = LS_Threshold(LSF_array, coarse_grid, mesh, 4)

    # Perform RBF interpolation on the fine grid
    println("Computing $name on the fine grid...")
    @time fine_LSF = rbf_interpolation_kdtree(fine_grid, coarse_grid, weights, kernel)

    # Adjust level set function to preserve volume
    fine_LSF_offset = fine_LSF .+ th

    # Convert to 3D array for output
    fine_LSF_offset_array = vector_to_array(fine_LSF_offset, dim)
    fine_sdf = fine_LSF_offset_array

    # Report final volume
    current_volume = calculate_volume_from_sdf(fine_sdf, fine_grid)
    @info "Body volume at SDF zero level: $current_volume (target: $(round(mesh.V_frac * mesh.V_domain, digits=4)))"

    return fine_sdf, fine_grid
end
