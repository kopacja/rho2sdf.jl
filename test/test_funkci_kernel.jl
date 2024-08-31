using KernelFunctions
using LinearAlgebra
using Plots

# Definice vlastního kernelu
struct CustomRBFKernel{T<:Real} <: Kernel
    σ::T
end

# Implementace funkce pro výpočet kernelu
function (k::CustomRBFKernel)(x::AbstractVector, y::AbstractVector)
    r = sqrt(sum((x .- y).^2))
    return exp(-(r / k.σ)^2)
end

# Vytvoření pravidelné mřížky
function create_grid(nx, ny, xmin, xmax, ymin, ymax)
    x = range(xmin, xmax, length=nx)
    y = range(ymin, ymax, length=ny)
    return [[x, y] for x in x, y in y]
end

# Přiřazení hodnot uzlům (pro ukázku použijeme funkci sinus)
function assign_values(grid)
    return [sin(5 * sqrt(x^2 + y^2)) for (x, y) in grid]
end

# Funkce pro sumaci RBF na jemné síti
function rbf_summation_fine_grid(fine_grid, coarse_grid, coarse_values, kernel)
    return [sum(coarse_values[i, j] * kernel([x, y], coarse_grid[i, j])
           for i in 1:size(coarse_grid, 1), j in 1:size(coarse_grid, 2))
           for (x, y) in fine_grid]
end

# Hlavní funkce
# function main()
    # Parametry hrubé mřížky
    nx_coarse, ny_coarse = 11, 11
    xmin, xmax, ymin, ymax = 0, 1, 0, 1

    # Vytvoření hrubé mřížky a přiřazení hodnot
    coarse_grid = create_grid(nx_coarse, ny_coarse, xmin, xmax, ymin, ymax)
    coarse_values = zeros(size(coarse_grid))
    coarse_values[1,1] = 1.

    # Parametry jemné mřížky
    smooth = 2
    nx_fine, ny_fine = (nx_coarse - 1) * 2 + 1, (ny_coarse - 1) * 2 + 1

    # Vytvoření jemné mřížky
    fine_grid = create_grid(nx_fine, ny_fine, xmin, xmax, ymin, ymax)

    # Definice Gaussova kernelu pomocí CustomRBFKernel
    σ = 0.1 # šířka Gaussovy funkce
    kernel = CustomRBFKernel(σ)

    # Výpočet RBF sumace na jemné mřížce
    fine_values = rbf_summation_fine_grid(fine_grid, coarse_grid, coarse_values, kernel)

    # Vizualizace výsledků
    x = range(xmin, xmax, length=nx_fine)
    y = range(ymin, ymax, length=ny_fine)
    heatmap(x, y, reshape(fine_values, nx_fine, ny_fine)', 
            title="RBF Interpolation (σ=$σ)", 
            xlabel="x", ylabel="y", 
            c=:viridis)
# end

# Spuštění hlavní funkce
main()


r = sqrt(sum((x .- y).^2))
    return exp(-(r / k.σ)^2)

    exp(-(1 / 1)^2)