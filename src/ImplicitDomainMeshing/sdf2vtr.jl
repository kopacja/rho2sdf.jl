using WriteVTK

function export_to_paraview_with_vtk(grid_data::Array{Vector{Float32}, 3}, sdf_data::Array{Float32, 3}, filename::String)
    # Get dimensions
    nx, ny, nz = size(grid_data)
    
    # Extract x, y, z coordinates from the grid
    x_coords = Vector{Float32}(undef, nx)
    y_coords = Vector{Float32}(undef, ny)
    z_coords = Vector{Float32}(undef, nz)
    
    # Extract unique coordinates (assuming regular grid)
    for i in 1:nx
        x_coords[i] = grid_data[i, 1, 1][1]
    end
    for j in 1:ny
        y_coords[j] = grid_data[1, j, 1][2]
    end
    for k in 1:nz
        z_coords[k] = grid_data[1, 1, k][3]
    end
    
    # Create VTK grid
    vtk_grid = WriteVTK.vtk_grid(filename, x_coords, y_coords, z_coords)
    
    # Add SDF data as point data
    vtk_point_data(vtk_grid, sdf_data, "SDF")
    
    # Save the file
    vtk_save(vtk_grid)
    
    println("Data exported to $(filename).vtr")
end

export_to_paraview_with_vtk(fine_grid, fine_sdf, "sdf_data")