using GLMakie
using LinearAlgebra

# Definice Gaussovy RBF
gaussian_rbf(r, B) = exp(-(r/B)^2)

# Diskrétní data pro interpolaci
x_data = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
y_data = [0.4, 0.7, 1.0, 0.85, 0.15, 0.3]
n_points = length(x_data)

# Funkce pro RBF interpolaci
function rbf_interpolation(x_centers, y_values, B)
    n = length(x_centers)
    
    # Sestavení RBF matice
    A = zeros(n, n)
    for i in 1:n
        for j in 1:n
            r = abs(x_centers[i] - x_centers[j])
            A[i, j] = gaussian_rbf(r, B)
        end
    end
    
    # Výpočet vah
    weights = A \ y_values
    
    # Návratová interpolační funkce
    function interpolant(x)
        result = 0.0
        for i in 1:n
            r = abs(x - x_centers[i])
            result += weights[i] * gaussian_rbf(r, B)
        end
        return result
    end
    
    return interpolant, weights
end

# Parametr B pro Gaussovu RBF
B = 1.0  # Změněno na 1.0 podle požadavku

# Vytvoření interpolační funkce
interpolant, weights = rbf_interpolation(x_data, y_data, B)

# Jemná síť pro vykreslení interpolované funkce
x_fine = range(-0.5, 5.5, length=300)
y_interpolated = [interpolant(x) for x in x_fine]

# Vytvoření grafu
fig = Figure(size=(800, 600))
ax = Axis(fig[1, 1], 
    title = "RBF Interpolation with Gaussian Kernel (B = $B)",
    xgridvisible = true,
    ygridvisible = true,
    xticksvisible = false,
    yticksvisible = false, 
    xticklabelsvisible = false,
    yticklabelsvisible = false)

# Vykreslení jednotlivých RBF funkcí
colors_rbf = [:red, :green, :blue, :purple, :orange, :brown]
for i in 1:n_points
    x_center = x_data[i]
    weight = weights[i]
    
    # Jednotlivá RBF funkce
    y_rbf = [weight * gaussian_rbf(abs(x - x_center), B) for x in x_fine]
    
    lines!(ax, x_fine, y_rbf, 
        color = colors_rbf[i], 
        linewidth = 2, 
        alpha = 0.7)
end

# Svislé čárkované čáry na pozicích center
for x_center in x_data
    vlines!(ax, [x_center],
        color = :gray, 
        linewidth = 1,
        linestyle = :dash,
        alpha = 0.7)
end

# Vyznačení bodů kde se RBFs "sčítají" - na pozicích [0,1,2,3,4,5]
for x_pos in x_data
    y_rbf_sum = 0.0
    for i in 1:n_points
        x_center = x_data[i]
        weight = weights[i]
        y_rbf_sum += weight * gaussian_rbf(abs(x_pos - x_center), B)
    end
    
    # Větší barevné body na jednotlivých RBFs v pozicích sčítání
    for i in 1:n_points
        x_center = x_data[i]
        weight = weights[i]
        y_rbf_individual = weight * gaussian_rbf(abs(x_pos - x_center), B)
        
        scatter!(ax, [x_pos], [y_rbf_individual], 
            color = colors_rbf[i], 
            markersize = 8, 
            marker = :circle,
            alpha = 0.8)
    end
end

# Interpolovaná funkce - tenká čárkovaná čára
lines!(ax, x_fine, y_interpolated, 
    color = :gray, 
    linewidth = 4, 
    linestyle = :dash,
    label = "RBF interpolation")

# Diskrétní body - větší černé kroužky
scatter!(ax, x_data, y_data, 
    color = :white, 
    markersize = 14, 
    marker = :circle,
    strokewidth = 4,
    strokecolor = :gray,
    label = "Discrete data")

# Nastavení os
xlims!(ax, -0.5, 5.5)
ylims!(ax, -0.3, 1.2)

# Legenda
axislegend(ax, position = :rt, framevisible = true, backgroundcolor = :white, alpha = 0.8)

display(fig)
