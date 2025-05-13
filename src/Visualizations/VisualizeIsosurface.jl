function visualize_stable_isosurface(sdf_data::Array)
    fig = Figure(size=(800, 800))
    
    # Calculate data dimensions
    nx, ny, nz = size(sdf_data)
    
    # Compute diagonal as reference size for visualization space
    diagonal = Float32(sqrt(nx^2 + ny^2 + nz^2))
    
    # Set axis limits with margin around the object
    offset = diagonal/4
    ax = Axis3(fig[1, 1];
        viewmode = :fit,
        # Symmetric limits around data center
        limits = (-offset, nx + offset,
                 -offset, ny + offset,
                 -offset, nz + offset),
        aspect = :data,
        perspectiveness = 0.5,
        xlabel = "X",
        ylabel = "Y",
        zlabel = "Z"
    )
    
    # Render the isosurface
    contour!(ax, sdf_data,
        levels = [0],
        transparency = false,
        color = :cornflowerblue
    )
    
    # Set initial view angles
    ax.azimuth = π/6
    ax.elevation = π/6
    
    return fig
end

"""
Example usage:
nx, ny, nz = 200, 100, 100
sdf_data = [sqrt((x - nx/2)^2 + (y - ny/2)^2 + (z - nz/2)^2) - 10
            for x in 1:nx, y in 1:ny, z in 1:nz]
visualize_stable_isosurface(sdf_data)
"""
