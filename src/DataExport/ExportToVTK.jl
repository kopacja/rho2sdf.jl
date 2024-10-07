
# function exportStructuredPointsToVTK_old(
#     fileName::String,
#     grid::Grid,
#     vals::Vector{Float64},
#     valLabel::String,
# )
#
#     dim = grid.N .+ 1
#     org = grid.AABB_min
#     spacing = grid.cell_size
#
#     io = open(fileName, "w")
#
#     write(io, "# vtk DataFile Version 1.0\n")
#     write(
#         io,
#         "Texture map for thresholding data (use boolean textures for 2D map)\n",
#     )
#
#     write(io, "ASCII\n\n")
#     write(io, "DATASET STRUCTURED_POINTS\n")
#     dim_x = dim[1]
#     dim_y = dim[2]
#     dim_z = dim[3]
#     write(io, "DIMENSIONS $dim_x $dim_y $dim_z\n") # dimenze pravidelné sítě
#
#     #     spacing_x = spacing[1]
#     #     spacing_y = spacing[2]
#     #     spacing_z = spacing[3]
#
#     #     write(io, "SPACING $spacing_x $spacing_y $spacing_z\n") # krok sítě ve 3 směrech jiný
#
#     write(io, "SPACING $spacing $spacing $spacing\n") # krok sítě ve 3 směrech stejný
#
#     org_x = org[1]
#     org_y = org[2]
#     org_z = org[3]
#     write(io, "ORIGIN $org_x $org_y $org_z\n\n") # souřadnice počátku
#
#     n = prod(dim)
#     write(io, "POINT_DATA $n\n") # počet uzlů pravidelné sítě
#     write(io, "SCALARS $valLabel float 1\n") # druh hodnoty v uzlech (vzdálenost)
#     write(io, "LOOKUP_TABLE default\n") # ??
#
#     for val ∈ vals
#         write(io, "$val\n") # hodnota vzdálenostní funkce v jednotlivých bodech
#     end
#     close(io)
# end
#


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
