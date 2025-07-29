using GLMakie
"""
    sphere_sdf_convergence_core(N::Int) -> Float64

Core function that computes SDF volume for a sphere with given grid resolution.

# Arguments
- `N::Int`: Number of elements per edge of the cube

# Returns
- `Float64`: Relative error compared to analytical volume
"""
function sphere_sdf_convergence_core(N::Int)
    # Test parameters
    r_sphere = 0.5  # sphere radius
    analytical_volume = (4/3) * π * r_sphere^3  # π/6 ≈ 0.5236
    
    print_info("Computing sphere SDF volume for N = $N")
    
    # Create uniform cartesian grid for cube with edge length 2
    # Cube centered at origin: [-1,1] × [-1,1] × [-1,1]
    AABB_min = [-1.0, -1.0, -1.0]
    AABB_max = [1.0, 1.0, 1.0]
    
    # Create grid using existing Grid constructor
    # Set margineCells = 0 to avoid adding extra margin
    grid = MeshGrid.Grid(AABB_min, AABB_max, N, 0)
    
    print_data("Grid info: N = $(grid.N), cell_size = $(grid.cell_size), ngp = $(grid.ngp)")
    
    # Grid dimensions (number of nodes)
    nx, ny, nz = grid.N .+ 1  # Convert from number of cells to number of nodes
    
    # Create SDF and grid arrays compatible with calculate_volume_from_sdf
    fine_sdf = Array{Float32,3}(undef, nx, ny, nz)
    fine_grid = Array{Vector{Float32},3}(undef, nx, ny, nz)
    
    # Fill arrays with SDF values and coordinates
    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                # Calculate node coordinates
                x = Float32(grid.AABB_min[1] + (i-1) * grid.cell_size)
                y = Float32(grid.AABB_min[2] + (j-1) * grid.cell_size)
                z = Float32(grid.AABB_min[3] + (k-1) * grid.cell_size)
                
                # Store coordinates as Vector{Float32}
                fine_grid[i, j, k] = [x, y, z]
                
                # Calculate SDF to sphere surface
                # SDF positive inside sphere, negative outside
                distance_to_center = sqrt(x^2 + y^2 + z^2)
                sdf_value = r_sphere - distance_to_center  # Positive inside, negative outside
                fine_sdf[i, j, k] = Float32(sdf_value)
            end
        end
    end
    
    print_info("Computing volume using SDF with detailed quadrature...")
    
    # Calculate volume using SDF with high-order quadrature
    # With iso_threshold = 0.0, we compute volume where sdf >= 0 (inside sphere)
    computed_volume = calculate_volume_from_sdf(
        fine_sdf, 
        fine_grid; 
        iso_threshold=0.0f0,
        detailed_quad_order=20
    )
    
    # Calculate relative error
    relative_error = abs(computed_volume - analytical_volume) / analytical_volume
    
    print_data("Computed volume: $(round(computed_volume, digits=6))")
    print_data("Analytical volume: $(round(analytical_volume, digits=6))")
    print_data("Relative error: $(round(relative_error * 100, digits=3))%")
    
    # Export to ParaView for visualization
    print_info("Exporting SDF to ParaView...")
    
    # Convert 3D array to vector for export
    sdf_vector = vec(fine_sdf)
    
    # Export SDF field to VTI format
    filename = "sphere_sdf_test_N$(N).vti"
    exportSdfToVTI(filename, grid, sdf_vector, "sdf_distance")
    
    print_success("Exported SDF to: $filename")
    
    return relative_error
end

"""
    plot_sdf_convergence(N_values, errors, convergence_rate, coeffs; plot_type=:dashboard)

Create convergence plots for SDF volume analysis.

# Arguments
- `N_values::Vector{Int}`: Grid resolutions tested
- `errors::Vector{Float64}`: Relative errors for each resolution
- `convergence_rate::Float64`: Computed convergence rate (p in O(h^p))
- `coeffs::Vector{Float64}`: Linear regression coefficients [a, b] from log(error) = a + b*log(N)
- `plot_type::Symbol`: Type of plot to create (:loglog, :semilog, :dashboard)

# Returns
- `String`: Filename of saved plot
"""
function plot_sdf_convergence(N_values, errors, convergence_rate, coeffs; plot_type=:dashboard)
    try
        analytical_volume = (4/3) * π * 0.5^3  # For sphere with r=0.5
        
        if plot_type == :loglog
            return _create_loglog_plot(N_values, errors, convergence_rate, coeffs)
        elseif plot_type == :semilog
            return _create_semilog_plot(N_values, errors, convergence_rate, coeffs)
        elseif plot_type == :relative
              return _create_relative_error_plot(N_values, errors, convergence_rate, coeffs, analytical_volume)
        elseif plot_type == :absolute
            return _create_absolute_error_plot(N_values, errors, convergence_rate, coeffs, analytical_volume)
        else
            error("Unknown plot_type: $plot_type. Use :loglog, :semilog, or :dashboard")
        end
        
    catch e
        print_warning("Could not create convergence plot: $e")
        return ""
    end
end

function _create_loglog_plot(N_values, errors, convergence_rate, coeffs)
    fig = Figure(size=(800, 500))
    ax = Axis(fig[1, 1],
             xlabel="Grid resolution N",
             ylabel="Relative error", 
             title="SDF volume convergence (log-log)",
             xscale=log10, yscale=log10,
             xgridvisible=true, ygridvisible=true,
             xlabelsize=20,
             ylabelsize=20,
             titlesize=22,
             xticklabelsize=18,
             yticklabelsize=18)
    
    scatter!(ax, N_values, errors,
            marker=:circle,
            markersize=15,
            color=:steelblue,
            strokewidth=4,
            strokecolor=:darkblue,
            label="Computed error")
            
    lines!(ax, N_values, errors,
           linewidth=5,
           color=:steelblue,
           alpha=0.8)
    
    # Fitted line
    N_fit = range(minimum(N_values), maximum(N_values), length=100)
    error_fit = exp.(coeffs[1] .+ coeffs[2] .* log.(N_fit))
    lines!(ax, N_fit, error_fit,
           linestyle=:dash,
           linewidth=5,
           color=:red,
           label="O(N^$(round(-coeffs[2], digits=2)))")
    
    axislegend(ax, position=:rt,
              labelsize=20,
              markersize=20,
              framevisible=true,
              padding=12)
    
    filename = "sphere_sdf_convergence_loglog.png"
    save(filename, fig)
    return filename
end

function _create_semilog_plot(N_values, errors, convergence_rate, coeffs)
    h_values = 2.0 ./ N_values  # Grid spacing
    errors_percent = errors .* 100
    
    fig = Figure(size=(800, 500))
    ax = Axis(fig[1, 1],
             xlabel="Grid spacing h",
             ylabel="Relative error [%]",
             title="SDF volume convergence (semi-log)",
             yscale=log10,
             xgridvisible=true, ygridvisible=true,
             xlabelsize=20,
             ylabelsize=20,
             titlesize=22,
             xticklabelsize=18,
             yticklabelsize=18)
    
    scatter!(ax, h_values, errors_percent,
            marker=:circle,
            markersize=15,
            color=:steelblue,
            strokewidth=4,
            strokecolor=:darkblue,
            label="Computed error")
            
    lines!(ax, h_values, errors_percent,
           linewidth=5,           # Stejná tloušťka jako v relative_error_plot
           color=:steelblue,
           alpha=0.8)
    
    # Tolerance lines
    hlines!(ax, [1.0],
           color=:green,
           linestyle=:dot,
           linewidth=5,
           alpha=0.7,
           label="1% Tolerance")
           
    hlines!(ax, [0.1],
           color=:orange,
           linestyle=:dot,
           linewidth=5,
           alpha=0.7,
           label="0.1% Tolerance")
    
    axislegend(ax, position=:rt,
              labelsize=20,
              markersize=20,
              framevisible=true,
              padding=12)
    
    filename = "sphere_sdf_convergence_semilog.png"
    save(filename, fig)
    return filename
end


function _create_absolute_error_plot(N_values, errors, convergence_rate, coeffs, analytical_volume)
    fig = Figure(size=(800, 500))
    
    # Single panel: Absolute Error with larger fonts
    ax = Axis(fig[1, 1],
              xlabel="Grid resolution N",
              ylabel="Absolute volume error",
              title="SDF Volume convergence - Absolute error",
              xgridvisible=true,
              ygridvisible=true,
              xlabelsize=20,
              ylabelsize=20,
              titlesize=22,
              xticklabelsize=18,
              yticklabelsize=18)
    
    abs_errors = errors .* analytical_volume
    
    # Plot absolute errors
    scatter!(ax, N_values, abs_errors,
            marker=:circle,
            markersize=15,
            color=:steelblue,
            strokewidth=4,
            strokecolor=:darkblue,
            label="Computed error")
    
    lines!(ax, N_values, abs_errors,
           linewidth=5,
           color=:steelblue,
           alpha=0.7)
    
    # Add theoretical convergence line for absolute error
    N_fit = range(minimum(N_values), maximum(N_values), length=100)
    # Convert relative error fit to absolute error
    relative_error_fit = exp.(coeffs[1] .+ coeffs[2] .* log.(N_fit))
    absolute_error_fit = relative_error_fit .* analytical_volume
    
    lines!(ax, N_fit, absolute_error_fit,
           linestyle=:dash,
           linewidth=5,
           color=:red,
           label="O(N^$(round(-coeffs[2], digits=2))) Fit")

    axislegend(ax, position=:rt,
          labelsize=20,
          markersize=20,
          framevisible=true,
          padding=12)
    
    filename = "sphere_sdf_absolute_error.png"
    save(filename, fig)
    return filename
end

function _create_relative_error_plot(N_values, errors, convergence_rate, coeffs, analytical_volume)
    fig = Figure(size=(800, 500))
    
    # Single panel: Relative Error with larger fonts
    ax = Axis(fig[1, 1],
              xlabel="Grid resolution N",
              ylabel="Relative volume error [%]",
              title="SDF Volume convergence - Relative error",
              xgridvisible=true,
              ygridvisible=true,
              xlabelsize=20,
              ylabelsize=20,
              titlesize=22,
              xticklabelsize=18,
              yticklabelsize=18)
    
    # Convert relative errors to percentages
    rel_errors = errors .* 100
    
    # Plot relative errors
    scatter!(ax, N_values, rel_errors,
            marker=:circle,
            markersize=15,
            color=:steelblue,
            strokewidth=4,
            strokecolor=:darkblue,
            label="Computed error")
    
    lines!(ax, N_values, rel_errors, 
           linewidth=5,
           color=:steelblue,
           alpha=0.7)
    
    # Add theoretical convergence line for relative error
    N_fit = range(minimum(N_values), maximum(N_values), length=100)
    # Relative error fit (convert to percentages)
    relative_error_fit = (exp.(coeffs[1] .+ coeffs[2] .* log.(N_fit))) .* 100
    
    lines!(ax, N_fit, relative_error_fit,
           linestyle=:dash,
           linewidth=5,
           color=:red,
           label="O(N^$(round(-coeffs[2], digits=2))) Fit")
           
    axislegend(ax, position=:rt, 
          labelsize=20, 
          markersize=20, 
          framevisible=true, 
          padding=12)
    
    filename = "sphere_sdf_relative_error.png"
    save(filename, fig)
    return filename
end



function run_Convergence_Test()
    print_info("Starting SDF volume convergence test...")
    
    # Test parameters
    r_sphere = 0.5  # sphere radius
    analytical_volume = (4/3) * π * r_sphere^3  # π/6 ≈ 0.5236
    
    print_success("Analytical sphere volume: $(round(analytical_volume, digits=6))")
    print_data("Expected value: π/6 ≈ $(round(π/6, digits=6))")
    
    # Test different grid resolutions - powers of 2 for better convergence analysis
    N_values = [4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128]  # Number of elements per edge
    errors = Float64[]
    
    # Check monotonic convergence - errors should decrease as N increases
    print_info("\nConvergence analysis:")
    for (i, (N, error)) in enumerate(zip(N_values, errors))
        status = if error < 0.01
            "✓ Excellent"
        elseif error < 0.05
            "✓ Good"
        elseif error < 0.1
            "⚠ Acceptable"
        else
            "✗ Poor"
        end
        print_data("N = $N: error = $(round(error * 100, digits=3))% ($status)")
    end
    
    # Calculate convergence rate
    print_info("\nConvergence rate analysis:")
    
    # Log-log slope calculation for convergence rate
    log_N = log.(N_values)
    log_errors = log.(errors)
    
    # Filter out any potential -Inf values (perfect convergence)
    valid_indices = findall(isfinite, log_errors)
    
    if length(valid_indices) >= 2
        # Linear regression on log-log data: log(error) = a + b*log(N)
        # Convergence rate is -b (negative because error decreases)
        A = [ones(length(valid_indices)) log_N[valid_indices]]
        coeffs = A \ log_errors[valid_indices]
        convergence_rate = -coeffs[2]
        
        print_data("Theoretical convergence rate (spatial order): O(h^p) where p ≈ $(round(convergence_rate, digits=2))")
        
        if convergence_rate > 1.5
            print_success("Good convergence rate: $(round(convergence_rate, digits=2))")
        else
            print_warning("Slow convergence rate: $(round(convergence_rate, digits=2))")
        end
    else
        print_warning("Cannot calculate convergence rate - insufficient valid data points")
    end
    
    print_info("Creating convergence plots...")
    
    # Create dashboard plot
    dashboard_file = plot_sdf_convergence(N_values, errors, convergence_rate, coeffs; plot_type=:relative)
    dashboard_file = plot_sdf_convergence(N_values, errors, convergence_rate, coeffs; plot_type=:loglog)
    dashboard_file = plot_sdf_convergence(N_values, errors, convergence_rate, coeffs; plot_type=:semilog)
    print_success("SDF volume convergence test completed successfully!")
    
    # Summary statistics
    print_info("\nSummary:")
    print_data("Finest grid (N=$(maximum(N_values))): $(round(minimum(errors) * 100, digits=3))% error")
    print_data("Coarsest grid (N=$(minimum(N_values))): $(round(maximum(errors) * 100, digits=3))% error")
    print_data("Improvement factor: $(round(maximum(errors)/minimum(errors), digits=1))x")
    
    # Additional debug info
    print_info("\nReference information:")
    print_data("Analytical volume: π/6 = $(π/6)")
    print_data("Cube volume: 8.0")
    print_data("Sphere volume fraction: $(round(analytical_volume/8.0 * 100, digits=2))%")
end

run_Convergence_Test()
