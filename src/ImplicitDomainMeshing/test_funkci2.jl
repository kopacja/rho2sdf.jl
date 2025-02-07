
function eval_sdf(mesh::BlockMesh, p::Vector)
    # Get grid bounds as before
    vmin = mesh.grid[1,1,1]
    vmax = mesh.grid[end,end,end]
    
    # Normalize coordinates to [0,1] as before
    r = (p .- vmin) ./ (vmax .- vmin)
    
    # Convert to grid indices (floating point)
    i_f = r[1] * (mesh.nx - 1) + 1
    j_f = r[2] * (mesh.ny - 1) + 1
    k_f = r[3] * (mesh.nz - 1) + 1
    
    # Get the base indices of the cell containing point p
    i0 = clamp(floor(Int, i_f), 1, mesh.nx-1)
    j0 = clamp(floor(Int, j_f), 1, mesh.ny-1)
    k0 = clamp(floor(Int, k_f), 1, mesh.nz-1)
    
    # Calculate local coordinates in [-1,1] range for shape functions
    # Map [i0,i0+1] → [-1,1], same for j and k
    ξ₁ = 2.0 * (i_f - i0) - 1.0  # Maps [0,1] → [-1,1]
    ξ₂ = 2.0 * (j_f - j0) - 1.0
    ξ₃ = 2.0 * (k_f - k0) - 1.0
    
    # Get shape functions for this point
    N = hex8_shape([ξ₁, ξ₂, ξ₃])
   
    sdf_values = get_cell_sdf_values(mesh, i0, j0, k0)
    
    # Interpolate using shape functions
    return sum(N .* sdf_values)
end

function eval_sdf_gradient(mesh::BlockMesh, p::Vector)
    # Získáme hranice mřížky stejně jako v eval_sdf
    vmin = mesh.grid[1,1,1]
    vmax = mesh.grid[end,end,end]
    
    # Normalizace souřadnic do [0,1]
    r = (p .- vmin) ./ (vmax .- vmin)
    
    # Výpočet indexů buňky
    i_f = r[1] * (mesh.nx - 1) + 1
    j_f = r[2] * (mesh.ny - 1) + 1
    k_f = r[3] * (mesh.nz - 1) + 1
    
    i0 = clamp(floor(Int, i_f), 1, mesh.nx-1)
    j0 = clamp(floor(Int, j_f), 1, mesh.ny-1)
    k0 = clamp(floor(Int, k_f), 1, mesh.nz-1)
    
    # Převod na lokální souřadnice v rozsahu [-1,1]
    ξ₁ = 2.0 * (i_f - i0) - 1.0
    ξ₂ = 2.0 * (j_f - j0) - 1.0
    ξ₃ = 2.0 * (k_f - k0) - 1.0
    
    # Alokace paměti pro tvarové funkce a jejich derivace
    N = MVector{8,Float64}(undef)
    d¹N_dξ¹ = MMatrix{8,3,Float64}(undef)
    
    # Výpočet tvarových funkcí a jejich derivací
    compute_hex8_shape!(N, d¹N_dξ¹, ξ₁, ξ₂, ξ₃)
    
    # Shromáždění SDF hodnot ve vrcholech buňky
    sdf_values = get_cell_sdf_values(mesh, i0, j0, k0)
    
    # Výpočet rozměrů buňky pro transformaci derivací
    Δx = (vmax[1] - vmin[1]) / (mesh.nx - 1)
    Δy = (vmax[2] - vmin[2]) / (mesh.ny - 1)
    Δz = (vmax[3] - vmin[3]) / (mesh.nz - 1)
    
    # Transformace derivací z lokálních na globální souřadnice
    # Používáme řetízkové pravidlo: ∂f/∂x = (∂f/∂ξ)(∂ξ/∂x)
    # Pro krychlové elementy je transformace jednoduchá:
    # ∂ξ₁/∂x = 2/Δx, ∂ξ₂/∂y = 2/Δy, ∂ξ₃/∂z = 2/Δz
    
    # Inicializace gradientu
    gradient = zeros(SVector{3,Float64})
    
    # Výpočet gradientu pomocí řetízkového pravidla
    for i in 1:8
        gradient = gradient .+ sdf_values[i] .* SVector{3,Float64}(
            d¹N_dξ¹[i,1] * (2.0/Δx),
            d¹N_dξ¹[i,2] * (2.0/Δy),
            d¹N_dξ¹[i,3] * (2.0/Δz)
        )
    end
    
    return gradient
end

