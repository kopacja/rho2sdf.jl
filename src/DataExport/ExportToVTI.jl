"""
   export_structured_field_to_vti(
       filename::String,
       grid::Grid,
       values::Vector{<:AbstractFloat},
       value_label::String,
       smooth::Union{Int,Nothing} = nothing)

Export scalar field defined on a structured grid to VTI file format (VTK ImageData).
VTI format is optimized for regular grids with constant spacing in all directions.

# Arguments
- `filename::String`: Output VTI filename (without extension, .vti will be added automatically)
- `grid::Grid`: Structured grid on which values are defined
- `values::Vector{<:AbstractFloat}`: Scalar values at each grid point
- `value_label::String`: Label for the scalar field (e.g., "distance")
- `smooth::Union{Int,Nothing} = nothing`: Optional refinement factor for generating finer grid

# Returns
- Reference to the created VTK object
"""
function exportSdfToVTI(
   filename::String,
   grid::Grid,
   values::Vector{<:AbstractFloat},
   value_label::String,
   smooth::Union{Int,Nothing} = nothing
)
   # Calculate dimensions, origin and grid spacing for VTI format
   dimensions = isnothing(smooth) ? grid.N .+ 1 : (grid.N * smooth) .+ 1
   origin = Float64.(grid.AABB_min)
   
   # If spacing is scalar, convert it to vector for each dimension
   if isa(grid.cell_size, Number)
       spacing = isnothing(smooth) ? 
                [grid.cell_size, grid.cell_size, grid.cell_size] :
                [grid.cell_size/smooth, grid.cell_size/smooth, grid.cell_size/smooth]
   else
       spacing = isnothing(smooth) ? 
                grid.cell_size : 
                grid.cell_size ./ smooth
   end
   
   # Verify correct length of values vector
   if length(values) != prod(dimensions)
       error("Values vector length ($(length(values))) doesn't match grid dimensions ($(prod(dimensions))).")
   end
   
   # Generate coordinates along each axis to create structured grid
   x = range(origin[1], length=dimensions[1], step=spacing[1])
   y = range(origin[2], length=dimensions[2], step=spacing[2])
   z = range(origin[3], length=dimensions[3], step=spacing[3])
   
   # Create VTI structured grid - using correct signature
   vtk_file = vtk_grid(filename, x, y, z)
   
   # Reshape values to 3D array for WriteVTK
   values_3d = reshape(values, dimensions...)
   
   # Add scalar field data
   vtk_point_data(vtk_file, values_3d, value_label)
   
   # Save VTI file
   vtk_save(vtk_file)
   
   return vtk_file
end
