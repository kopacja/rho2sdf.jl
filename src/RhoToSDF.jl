"""
    Rho2sdfOptions

Configuration options for density to SDF conversion.

# Fields
- `threshold_density::Union{Float64, Nothing}`: Threshold density in range [0,1]. If `nothing`, calculated automatically.
- `sdf_grid_setup::Symbol`: Grid setup method: `:manual` or `:automatic`.
- `export_nodal_densities::Bool`: Whether to export nodal densities to VTU.
- `export_raw_sdf::Bool`: Whether to export raw SDF to VTI.
- `rbf_interp::Bool`: Use RBF interpolation (true) or approximation (false).
- `rbf_grid::Symbol`: RBF grid resolution: `:same` (smooth=1) or `:fine` (smooth=2).
"""
struct Rho2sdfOptions
    threshold_density::Union{Float64, Nothing}
    sdf_grid_setup::Symbol
    export_nodal_densities::Bool
    export_raw_sdf::Bool
    rbf_interp::Bool
    rbf_grid::Symbol
    
    # Constructor with default values
    function Rho2sdfOptions(;
        threshold_density=nothing,
        sdf_grid_setup=:manual,
        export_nodal_densities=false,
        export_raw_sdf=false,
        rbf_interp=true,
        rbf_grid=:same
    )
        # Validate threshold density if provided
        if threshold_density !== nothing
            if !(0.0 <= threshold_density <= 1.0)
                @warn "Threshold density $threshold_density is outside the valid range [0.0, 1.0]. Will use automatic calculation instead."
                threshold_density = nothing
            elseif threshold_density == 0.0 || threshold_density == 1.0
                @warn "Using extreme threshold density value: $threshold_density"
            end
        end
        
        # Validate sdf_grid_setup and set default if invalid
        if !(sdf_grid_setup in [:manual, :automatic])
            @warn "Invalid sdf_grid_setup: $sdf_grid_setup. Must be either :manual or :automatic. Using default :manual instead."
            sdf_grid_setup = :manual
        end
        
        # Validate rbf_grid and set default if invalid
        if !(rbf_grid in [:same, :fine])
            @warn "Invalid rbf_grid: $rbf_grid. Must be either :same or :fine. Using default :same instead."
            rbf_grid = :same
        end
        
        new(threshold_density, sdf_grid_setup, export_nodal_densities, export_raw_sdf, rbf_interp, rbf_grid)
    end
end

"""
    rho2sdf(taskName, X, IEN, rho; options=Rho2sdfOptions())

Convert element densities to Signed Distance Function representation.

# Arguments
- `taskName::String`: Base name for output files
- `X::Vector{Vector{Float64}}`: Mesh node coordinates
- `IEN::Vector{Vector{Int64}}`: Element connectivity
- `rho::Vector{Float64}}`: Element densities
- `options::Rho2sdfOptions`: Configuration options (optional)

# Returns
- `Tuple`: (fine_sdf, fine_grid, sdf_grid, sdf_dists)

# Example
```julia
# Basic usage with default options
result = rho2sdf("chapadlo", X, IEN, rho)

# Custom options
options = Rho2sdfOptions(
    threshold_density=0.5,
    sdf_grid_setup=:automatic,
    export_nodal_densities=true
)
result = rho2sdf("chapadlo", X, IEN, rho, options=options)
"""

function rho2sdf(taskName::String, X::Vector{Vector{Float64}}, IEN::Vector{Vector{Int64}}, rho::Vector{Float64};
    options::Rho2sdfOptions=Rho2sdfOptions())
    # Generate FEM mesh structure
    mesh = Mesh(X, IEN, rho, hex8_shape)
    
    # Setup SDF grid based on selected method
    sdf_grid = if options.sdf_grid_setup == :manual
        interactive_sdf_grid_setup(mesh)
    else
        # Use automatic setup with default cell size factor
        noninteractive_sdf_grid_setup(mesh, 2.0)
    end

    # Map elemental densities to the nodes using least squares
    ρₙ = DenseInNodes(mesh, rho)
    
    # Determine threshold density (ρₜ)
    ρₜ = if options.threshold_density === nothing
        # Automatically calculate threshold based on volume
        find_threshold_for_volume(mesh, ρₙ)
    else
        options.threshold_density
    end
    
    # Export nodal densities if requested
    if options.export_nodal_densities
        VTK_CODE = 12  # Hexahedron code per VTK specification
        exportToVTU(taskName * "_nodal_densities.vtu", X, IEN, VTK_CODE, ρₙ)
    end
    
    # Generate grid points for the SDF calculation
    points = generateGridPoints(sdf_grid)
    
    # Calculate SDF from densities
    (dists, xp) = evalDistances(mesh, sdf_grid, points, ρₙ, ρₜ)
    signs = Sign_Detection(mesh, sdf_grid, points, ρₙ, ρₜ)
    sdf_dists = dists .* signs  # Combine distances with signs to create signed distance field
    
    # Export raw SDF to VTK if requested
    if options.export_raw_sdf
        B = round(sdf_grid.cell_size, digits=4)
        exportSdfToVTI(taskName * "_SDF_CellSize-" * string(B) * ".vti", sdf_grid, sdf_dists, "distance")
    end
    
    # RBF smoothing with selected parameters
    smooth = options.rbf_grid == :same ? 1 : 2  # Same=1, Fine=2
    (fine_sdf, fine_grid) = RBFs_smoothing(mesh, sdf_dists, sdf_grid, options.rbf_interp, smooth, taskName)
    
    # Export final smoothed results
    export_sdf_results(fine_sdf, fine_grid, sdf_grid, taskName, smooth, options.rbf_interp)
    
    # Return all relevant results
    return (fine_sdf, fine_grid, sdf_grid, sdf_dists)

end


