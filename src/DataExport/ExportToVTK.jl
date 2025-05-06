"""
   exportStructuredPointsToVTK(fileName, grid, vals, valLabel, smooth=nothing)

Exports scalar field data from a structured grid to VTK file format (ASCII).

- `fileName`: Output VTK file path
- `grid`: Grid object defining the domain
- `vals`: Scalar values at each grid point
- `valLabel`: Name of the scalar field (e.g., "distance")
- `smooth`: Optional refinement factor for higher resolution
"""
function exportStructuredPointsToVTK(
   fileName::String,
   grid::Grid,
   vals::Vector{<:AbstractFloat},
   valLabel::String,
   smooth::Union{Int,Nothing} = nothing
)
   dim = isnothing(smooth) ? grid.N .+ 1 : (grid.N * smooth) .+ 1
   org = grid.AABB_min
   spacing = isnothing(smooth) ? grid.cell_size : grid.cell_size / smooth

   io = open(fileName, "w")
   write(io, "# vtk DataFile Version 1.0\n")
   write(io, "Texture map for thresholding data (use boolean textures for 2D map)\n")
   write(io, "ASCII\n\n")
   write(io, "DATASET STRUCTURED_POINTS\n")
   
   dim_x, dim_y, dim_z = dim
   write(io, "DIMENSIONS $dim_x $dim_y $dim_z\n") # dimensions of the regular grid
   write(io, "SPACING $spacing $spacing $spacing\n") # grid step in 3 directions (same for all)
   
   org_x, org_y, org_z = org
   write(io, "ORIGIN $org_x $org_y $org_z\n\n") # coordinates of the origin
   
   n = prod(dim)
   write(io, "POINT_DATA $n\n") # number of nodes in the regular grid
   write(io, "SCALARS $valLabel float 1\n") # type of value at nodes (distance)
   write(io, "LOOKUP_TABLE default\n")
   
   for val in vals
       write(io, "$val\n") # distance function value at each point
   end
   
   close(io)
end
