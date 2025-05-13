"""
    export_sdf_results(fine_sdf, fine_grid, my_grid, taskName, smooth, is_interpolation)

Export SDF results to VTI and JLD2 formats for visualization and storage.

# Arguments
- `fine_sdf`: The smoothed signed distance function values
- `fine_grid`: The fine grid on which the SDF is defined
- `my_grid`: The original grid structure
- `taskName`: Base name for the output files
- `smooth`: Grid refinement factor used
- `is_interpolation`: Boolean indicating whether interpolation (true) or approximation (false) was used

# Returns
- Nothing, but exports files to disk
"""
function export_sdf_results(fine_sdf, fine_grid, my_grid, taskName, smooth, is_interpolation)
    name = is_interpolation ? "Interpolation" : "Approximation"
    B = round(my_grid.cell_size, digits=4)
    
    # Convert 3D array to vector if needed for exportSdfToVTI
    fine_LSF_offset = isa(fine_sdf, Array{<:Any,3}) ? vec(fine_sdf) : fine_sdf
    
    # Export to VTI format for visualization
    println("Exporting results to VTI...")
    exportSdfToVTI("$(taskName)_B-$(B)_smooth-$(smooth)_$(name).vti", 
                   my_grid, fine_LSF_offset, "distance", smooth)
    
    # Save results to JLD2 files for later use
    println("Saving results to JLD2 files...")
    @save "Z_$(taskName)_FineSDF_B-$(B)_smooth-$(smooth)_$(name).jld2" fine_sdf
    @save "Z_$(taskName)_FineGrid_B-$(B)_smooth-$(smooth)_$(name).jld2" fine_grid
    
    println("Export completed successfully.")
end
