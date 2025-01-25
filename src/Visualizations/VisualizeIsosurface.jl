function visualize_stable_isosurface(sdf_data::Array)
    fig = Figure(size=(800, 800))
    
    # Vypočítáme rozměry dat
    nx, ny, nz = size(sdf_data)
    
    # Vypočítáme úhlopříčku, ale tentokrát ji použijeme jen jako referenční velikost
    # pro nastavení zobrazovacího prostoru
    diagonal = Float32(sqrt(nx^2 + ny^2 + nz^2))
    
    # Důležitá změna: limity os nastavíme tak, aby odpovídaly skutečnému rozsahu dat
    # Přidáme malý offset (diagonal/4) pro vytvoření okraje kolem tělesa
    offset = diagonal/4
    ax = Axis3(fig[1, 1];
        viewmode = :fit,
        # Limity nastavíme symetricky kolem středu dat
        limits = (-offset, nx + offset,
                 -offset, ny + offset,
                 -offset, nz + offset),
        aspect = :data,
        perspectiveness = 0.5,
        xlabel = "X",
        ylabel = "Y",
        zlabel = "Z"
    )
    
    # Vykreslení izoplochy zůstává stejné
    contour!(ax, sdf_data,
        levels = [0],
        transparency = false,
        color = :cornflowerblue
    )
    
    # Nastavíme počáteční pohled
    ax.azimuth = π/6
    ax.elevation = π/6
    
    return fig
end

"""
Vytvoříme testovací data reprezentující kouli
nx, ny, nz = 200, 100, 100
sdf_data = [sqrt((x - nx/2)^2 + (y - ny/2)^2 + (z - nz/2)^2) - 10
            for x in 1:nx, y in 1:ny, z in 1:nz]
visualize_stable_isosurface(sdf_data)
"""
