using .ElementTypes
using .ShapeFunctions

"""
    Rho2sdfOptions

Configuration options for density to SDF conversion with element type support.
"""
struct Rho2sdfOptions
    threshold_density::Union{Float64, Nothing}
    sdf_grid_setup::Symbol
    export_input_data::Bool
    export_nodal_densities::Bool
    export_raw_sdf::Bool
    rbf_interp::Bool
    rbf_grid::Symbol
    remove_artifacts::Bool
    artifact_min_component_ratio::Float64
    export_analysis::Bool
    element_type::Type{<:AbstractElement}
    
    function Rho2sdfOptions(;
        threshold_density=nothing,
        sdf_grid_setup=:manual,
        export_input_data=false,
        export_nodal_densities=false,
        export_raw_sdf=false,
        rbf_interp=true,
        rbf_grid=:same,
        remove_artifacts=true,
        artifact_min_component_ratio=0.01,
        export_analysis=false,
        element_type::Type{<:AbstractElement}=HEX8
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
        
        # Validate sdf_grid_setup
        if !(sdf_grid_setup in [:manual, :automatic])
            @warn "Invalid sdf_grid_setup: $sdf_grid_setup. Must be either :manual or :automatic. Using default :manual instead."
            sdf_grid_setup = :manual
        end
        
        # Validate rbf_grid
        if !(rbf_grid in [:same, :fine])
            @warn "Invalid rbf_grid: $rbf_grid. Must be either :same or :fine. Using default :same instead."
            rbf_grid = :same
        end
        
        # Validate element type
        if !(element_type <: AbstractElement)
            @warn "Invalid element_type: $element_type. Must be subtype of AbstractElement. Using HEX8."
            element_type = HEX8
        end
        
        new(threshold_density, sdf_grid_setup, export_input_data, export_nodal_densities, 
            export_raw_sdf, rbf_interp, rbf_grid, remove_artifacts, artifact_min_component_ratio, 
            export_analysis, element_type)
    end
end

"""
    rho2sdf(taskName, X, IEN, rho; options=Rho2sdfOptions())

Convert element densities to Signed Distance Function representation with element type support.

# Arguments
- `taskName::String`: Base name for output files
- `X::Vector{Vector{Float64}}`: Mesh node coordinates  
- `IEN::Vector{Vector{Int64}}`: Element connectivity
- `rho::Vector{Float64}}`: Element densities
- `options::Rho2sdfOptions`: Configuration options (optional)

# Returns
- `Tuple`: (fine_sdf, fine_grid, sdf_grid, sdf_dists)

# Examples
```julia
# HEX8 elements (default)
result = rho2sdf("example", X, IEN, rho)

# TET4 elements
options = Rho2sdfOptions(element_type=TET4)
result = rho2sdf("example", X, IEN, rho, options=options)

# Custom configuration
options = Rho2sdfOptions(
    element_type=HEX8,
    threshold_density=0.5,
    sdf_grid_setup=:automatic,
    export_nodal_densities=true
)
result = rho2sdf("example", X, IEN, rho, options=options)
```
"""
function rho2sdf(taskName::String, X::Vector{Vector{Float64}}, IEN::Vector{Vector{Int64}}, rho::Vector{Float64};
    options::Rho2sdfOptions=Rho2sdfOptions())
    
    # Create shape function based on element type  
    shape_func = coords -> shape_functions(options.element_type, coords)
    
    # Generate FEM mesh structure with specified element type
    mesh = Mesh(X, IEN, rho, shape_func; element_type=options.element_type)
    
    # Print element type information
    element_name = string(options.element_type.name.name)
    num_nodes = get_num_nodes(options.element_type)
    print_info("Using element type: $element_name with $num_nodes nodes per element")

    # Export input data to Paraview (vtu) if requested
    if options.export_input_data
        InputDataToVTU(mesh, taskName * "-input_data")
    end
    
    # Setup SDF grid based on selected method
    sdf_grid = if options.sdf_grid_setup == :manual
        interactive_sdf_grid_setup(mesh)
    else
        noninteractive_sdf_grid_setup(mesh)
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
        VTK_CODE = options.element_type == HEX8 ? 12 : 10  # VTK element codes
        exportToVTU(taskName * "_nodal_densities.vtu", X, IEN, VTK_CODE, ρₙ)
    end
    
    # Generate grid points for the SDF calculation
    points = generateGridPoints(sdf_grid)
    
    # Calculate SDF from densities using element-type-aware functions
    (dists, xp) = evalDistances(mesh, sdf_grid, points, ρₙ, ρₜ)
    signs = Sign_Detection(mesh, sdf_grid, points, ρₙ, ρₜ)
    sdf_dists = dists .* signs

    # Artifact removal step
    if options.remove_artifacts
        print_info("\nStarting SDF artifact removal...")
        
        # Optional: Export analysis before cleanup
        if options.export_analysis
            components_before = analyze_sdf_components(sdf_dists, sdf_grid)
            
            # Export raw SDF for comparison
            B = round(sdf_grid.cell_size, digits=4)
            exportSdfToVTI(taskName * "_SDF_raw_$(element_name)_B-$(B).vti", sdf_grid, sdf_dists, "distance")
        end
        
        # Perform artifact removal
        nodes_flipped = remove_sdf_artifacts!(
            sdf_dists, 
            sdf_grid; 
            threshold=0.0,
            min_component_ratio=options.artifact_min_component_ratio
        )
        
        if nodes_flipped > 0 && options.export_analysis
            B = round(sdf_grid.cell_size, digits=4)
            exportSdfToVTI(taskName * "_SDF_cleaned_$(element_name)_B-$(B).vti", sdf_grid, sdf_dists, "distance")
        end
    end
    
    # Export raw SDF to VTK if requested
    if options.export_raw_sdf
        B = round(sdf_grid.cell_size, digits=4)
        exportSdfToVTI(taskName * "_SDF_$(element_name)_CellSize-" * string(B) * ".vti", sdf_grid, sdf_dists, "distance")
    end
    
    # RBF smoothing with selected parameters
    smooth = options.rbf_grid == :same ? 1 : 2
    (fine_sdf, fine_grid) = RBFs_smoothing(mesh, sdf_dists, sdf_grid, options.rbf_interp, smooth, taskName)
    
    # Export final smoothed results with element type in filename
    export_sdf_results_with_element_type(fine_sdf, fine_grid, sdf_grid, taskName, smooth, options.rbf_interp, options.element_type)
    
    # Return all relevant results
    return (fine_sdf, fine_grid, sdf_grid, sdf_dists)
end

"""
    export_sdf_results_with_element_type(fine_sdf, fine_grid, my_grid, taskName, smooth, is_interpolation, element_type)

Export SDF results with element type information in filename.
"""
function export_sdf_results_with_element_type(fine_sdf, fine_grid, my_grid, taskName, smooth, is_interpolation, element_type)
    name = is_interpolation ? "Interpolation" : "Approximation"
    element_name = string(element_type.name.name)
    B = round(my_grid.cell_size, digits=4)
    
    # Convert 3D array to vector if needed for exportSdfToVTI
    fine_LSF_offset = isa(fine_sdf, Array{<:Any,3}) ? vec(fine_sdf) : fine_sdf
    
    # Export to VTI format for visualization
    println("Exporting results to VTI...")
    exportSdfToVTI("$(taskName)_$(element_name)_B-$(B)_smooth-$(smooth)_$(name).vti", 
                   my_grid, fine_LSF_offset, "distance", smooth)
    
    # Save results to JLD2 files for later use
    println("Saving results to JLD2 files...")
    @save "Z_$(taskName)_$(element_name)_FineSDF_B-$(B)_smooth-$(smooth)_$(name).jld2" fine_sdf
    @save "Z_$(taskName)_$(element_name)_FineGrid_B-$(B)_smooth-$(smooth)_$(name).jld2" fine_grid
    
    println("Export completed successfully.")
end

# Convenience functions for specific element types
function rho2sdf_hex8(taskName::String, X::Vector{Vector{Float64}}, IEN::Vector{Vector{Int64}}, rho::Vector{Float64};
                     kwargs...)
    options = Rho2sdfOptions(; element_type=HEX8, kwargs...)
    return rho2sdf(taskName, X, IEN, rho; options=options)
end

function rho2sdf_tet4(taskName::String, X::Vector{Vector{Float64}}, IEN::Vector{Vector{Int64}}, rho::Vector{Float64};
                     kwargs...)
    options = Rho2sdfOptions(; element_type=TET4, kwargs...)
    return rho2sdf(taskName, X, IEN, rho; options=options)
end
