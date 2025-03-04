using LinearAlgebra

# Abstraktní typ pro různé tvary rovin
abstract type PlaneShape end

# Definice konkrétních tvarů
struct Rectangle <: PlaneShape
    width::Float64  # šířka
    height::Float64 # výška
end

# Speciální případ obdélníku
Square(size::Float64) = Rectangle(size, size)

struct Circle <: PlaneShape
    radius::Float64 # poloměr
end

struct Ellipse <: PlaneShape
    a::Float64      # hlavní poloosa
    b::Float64      # vedlejší poloosa
end

# Struktura pro reprezentaci omezené roviny
struct BoundedPlane
    normal::Vector{Float64}   # Normálový vektor roviny
    point::Vector{Float64}    # Bod ležící v rovině
    shape::PlaneShape         # Omezující tvar
    u::Vector{Float64}        # První bázový vektor v rovině
    v::Vector{Float64}        # Druhý bázový vektor v rovině
    
    # Konstruktor pro výpočet bázových vektorů
    function BoundedPlane(normal::Vector{Float64}, point::Vector{Float64}, shape::PlaneShape)
        # Normalizace normálového vektoru
        normal = normalize(normal)
        
        # Vytvoření ortogonální báze v rovině
        # Nejprve zvolíme libovolný vektor, který není rovnoběžný s normálou
        temp = abs(normal[1]) < 0.9 ? [1.0, 0.0, 0.0] : [0.0, 1.0, 0.0]
        
        # Vytvoříme první bázový vektor v rovině pomocí vektorového součinu
        u = normalize(cross(normal, temp))
        
        # Vytvoříme druhý bázový vektor kolmý k normále a prvnímu bázovému vektoru
        v = normalize(cross(normal, u))
        
        new(normal, point, shape, u, v)
    end
end

# Funkce pro kontrolu, zda bod leží na omezené rovině
function is_on_plane(plane::BoundedPlane, point::Vector{Float64}, tolerance::Float64=1e-10)
    # Kontrola, zda bod leží v rovině
    dist_to_plane = abs(dot(plane.normal, point - plane.point))
    if dist_to_plane > tolerance
        return false
    end
    
    # Projekce bodu do souřadného systému roviny
    vec_to_point = point - plane.point
    
    # Souřadnice bodu v rovině
    u_coord = dot(vec_to_point, plane.u)
    v_coord = dot(vec_to_point, plane.v)
    
    # Kontrola, zda projektovaný bod leží v daném tvaru
    return is_in_shape(plane.shape, u_coord, v_coord)
end

# Kontrola, zda bod (v souřadnicích roviny) leží v obdélníku
function is_in_shape(shape::Rectangle, u::Float64, v::Float64)
    return abs(u) <= shape.width/2 && abs(v) <= shape.height/2
end

# Kontrola, zda bod (v souřadnicích roviny) leží v kruhu
function is_in_shape(shape::Circle, u::Float64, v::Float64)
    return u^2 + v^2 <= shape.radius^2
end

# Kontrola, zda bod (v souřadnicích roviny) leží v elipse
function is_in_shape(shape::Ellipse, u::Float64, v::Float64)
    return (u/shape.a)^2 + (v/shape.b)^2 <= 1
end

# Funkce pro generování bodů na omezené rovině
function generate_points(plane::BoundedPlane, num_points::Int)
    points = Vector{Vector{Float64}}(undef, 0)
    
    if plane.shape isa Rectangle
        # Generování bodů pro obdélník
        width, height = plane.shape.width, plane.shape.height
        
        # Výpočet počtu bodů v každém směru
        aspect_ratio = width / height
        nx = ceil(Int, sqrt(num_points * aspect_ratio))
        ny = ceil(Int, num_points / nx)
        
        # Generování mřížky bodů
        for i in 1:nx
            for j in 1:ny
                u_coord = width * (i / (nx+1) - 0.5)
                v_coord = height * (j / (ny+1) - 0.5)
                
                # Výpočet 3D souřadnic bodu
                point = plane.point + u_coord * plane.u + v_coord * plane.v
                push!(points, point)
            end
        end
    elseif plane.shape isa Circle
        # Pro kruh používáme polární souřadnice
        radius = plane.shape.radius
        
        # Určení počtu kružnic a bodů na kružnici
        num_circles = ceil(Int, sqrt(num_points))
        points_per_circle = ceil(Int, num_points / num_circles)
        
        for r_idx in 1:num_circles
            r = radius * (r_idx / num_circles)
            
            for theta_idx in 1:points_per_circle
                theta = 2π * (theta_idx / points_per_circle)
                
                u_coord = r * cos(theta)
                v_coord = r * sin(theta)
                
                point = plane.point + u_coord * plane.u + v_coord * plane.v
                push!(points, point)
            end
        end
    elseif plane.shape isa Ellipse
        # Pro elipsu také používáme formu polárních souřadnic
        a, b = plane.shape.a, plane.shape.b
        
        num_circles = ceil(Int, sqrt(num_points))
        points_per_circle = ceil(Int, num_points / num_circles)
        
        for r_idx in 1:num_circles
            r = r_idx / num_circles
            
            for theta_idx in 1:points_per_circle
                theta = 2π * (theta_idx / points_per_circle)
                
                u_coord = a * r * cos(theta)
                v_coord = b * r * sin(theta)
                
                point = plane.point + u_coord * plane.u + v_coord * plane.v
                push!(points, point)
            end
        end
    end
    
    return points
end

# Funkce pro získání obecné rovnice roviny (ax + by + cz + d = 0)
function get_plane_equation(plane::BoundedPlane)
    a, b, c = plane.normal
    d = -dot(plane.normal, plane.point)
    return (a, b, c, d)
end

# Funkce pro výpočet bodu na omezené rovině v daných 2D souřadnicích roviny
function point_on_plane(plane::BoundedPlane, u_coord::Float64, v_coord::Float64)
    return plane.point + u_coord * plane.u + v_coord * plane.v
end

using LinearAlgebra
using JLD2
# using BoundedPlanes  # Předpokládá se, že modul BoundedPlanes je již načten

"""
    PlaneDefinition

Struktura pro definici roviny s omezeným tvarem.

Atributy:
- `normal`: Normálový vektor roviny
- `point`: Bod ležící v rovině
- `shape`: Omezující tvar (Square, Rectangle, Circle, Ellipse)
"""
struct PlaneDefinition
    normal::Vector{Float64}
    point::Vector{Float64}
    shape::PlaneShape
end

"""
    compute_planes_sdf(fine_grid::Array{Vector{Float32}, 3}, fine_sdf::Array{Float32, 3}, plane_definitions::Vector{PlaneDefinition}; max_dist_factor::Float64=3.0)

Vypočítá minimální vzdálenosti uzlů pravidelné sítě na definované roviny pomocí kolmé projekce.

Argumenty:
- `fine_grid`: 3D pole vektorů reprezentujících uzly mřížky
- `fine_sdf`: Původní pole vzdáleností (není přímo používáno, ale je součástí signatury)
- `plane_definitions`: Vektor definic rovin
- `max_dist_factor`: Násobitel kroku mřížky pro maximální uvažovanou vzdálenost (výchozí: 3.0)

Vrací:
- `plane_sdf`: 3D pole reprezentující minimální vzdálenosti k rovinám
- `plane_normals`: 3D pole reprezentující normály rovin v bodech s minimální vzdáleností
"""
function compute_planes_sdf(fine_grid::Array{Vector{Float32}, 3}, fine_sdf::Array{Float32, 3}, 
                           plane_definitions::Vector{PlaneDefinition}; max_dist_factor::Float64=3.0)
    # Vytvoření omezených rovin ze zadaných definic
    planes = [BoundedPlane(def.normal, def.point, def.shape) for def in plane_definitions]
    
    # Získání rozměrů mřížky
    nx, ny, nz = size(fine_grid)
    
    # Inicializace výstupních polí
    plane_sdf = fill(1.0f10, nx, ny, nz)  # Předalokace vysokou kladnou hodnotou
    # Inicializace pole normál prázdnými vektory místo nulových
    plane_normals = Array{Vector{Float32}, 3}(undef, nx, ny, nz)
    for i in 1:nx, j in 1:ny, k in 1:nz
        plane_normals[i, j, k] = Vector{Float32}()  # Prázdný vektor
    end
    
    # Odhad kroku mřížky
    spacing = estimate_grid_spacing(fine_grid)
    
    # Maximální vzdálenost pro uvažování (násobek kroku mřížky)
    max_dist = max_dist_factor * spacing
    
    # Výpočet vzdáleností pro každý uzel v mřížce
    for i in 1:nx, j in 1:ny, k in 1:nz
        point = fine_grid[i, j, k]
        
        # Projít všechny roviny a najít nejmenší vzdálenost
        for (p_idx, plane) in enumerate(planes)
            dist = distance_to_bounded_plane(plane, point)
            if dist <= max_dist && dist < plane_sdf[i, j, k]
                plane_sdf[i, j, k] = dist
                plane_normals[i, j, k] = Vector{Float32}(plane_definitions[p_idx].normal)
            end
        end
    end
    
    return plane_sdf, plane_normals
end

"""
    estimate_grid_spacing(fine_grid::Array{Vector{Float32}, 3})

Odhadne krok pravidelné mřížky pomocí maximální absolutní hodnoty rozdílu souřadnic
mezi prvními dvěma uzly v ose z.
"""
function estimate_grid_spacing(fine_grid::Array{Vector{Float32}, 3})
    # Použití maximální absolutní hodnoty rozdílu souřadnic
    return maximum(abs.(fine_grid[1,1,2] - fine_grid[1,1,1]))
end

"""
    distance_to_bounded_plane(plane::BoundedPlane, point::Vector{Float32})

Vypočítá vzdálenost z bodu k omezené rovině jako kolmou projekci.
"""
function distance_to_bounded_plane(plane::BoundedPlane, point::Vector{Float32})
    # Vektor z bodu roviny k aktuálnímu bodu
    vec_to_point = point - plane.point
    
    # Vzdálenost k nekonečné rovině (absolutní hodnota skalárního součinu s normálou)
    dist_to_infinite_plane = abs(dot(plane.normal, vec_to_point))
    
    # Projekce bodu na nekonečnou rovinu
    projected_point = point - dot(plane.normal, vec_to_point) * plane.normal
    
    # Kontrola, zda projektovaný bod leží v rámci omezené roviny
    if is_on_plane(plane, projected_point)
        # Pokud projekce leží v omezené rovině, vrátíme vzdálenost k nekonečné rovině
        return dist_to_infinite_plane
    else
        # Pokud projekce neleží v omezené rovině, vrátíme vysokou hodnotu
        return 1.0e10
    end
end

# Načtení dat
@load "src/ImplicitDomainMeshing/data/Z_cantilever_beam_vfrac_04_FineGrid_B-1.0_smooth-1.jld2" fine_grid
@load "src/ImplicitDomainMeshing/data/Z_cantilever_beam_vfrac_04_FineSDF_B-1.0_smooth-1.jld2" fine_sdf

# Definice rovin
plane_definitions = [
    PlaneDefinition([-1.0, 0.0, 0.0], [0.0, 10.0, 0.0], Square(30.)),
    PlaneDefinition([1.0, 0.0, 0.0], [60.0, 2.0, 2.0], Square(5.))
]

# Výpočet minimálních vzdáleností k rovinám a odpovídajících normál
plane_sdf, plane_normals = compute_planes_sdf(fine_grid, fine_sdf, plane_definitions)

export_to_paraview_with_vtk(fine_grid, plane_sdf, "sdf_kontrola")
typeof(fine_sdf)
typeof(plane_normals)

"""
    compute_gradient(fine_sdf::Array{Float32, 3})

Vypočítá gradient SDF jako pole jednotkových vektorů.

Argumenty:
- `fine_sdf`: 3D pole vzdáleností SDF

Vrací:
- `gradient`: 3D pole vektorů reprezentujících směr gradientu SDF
"""
function compute_gradient(fine_sdf::Array{Float32, 3})
    nx, ny, nz = size(fine_sdf)
    gradient = Array{Vector{Float32}, 3}(undef, nx, ny, nz)
    gradient_norm = Array{Vector{Float32}, 3}(undef, nx, ny, nz)
    
    # Procházení všech uzlů v síti
    for i in 1:nx, j in 1:ny, k in 1:nz
        # Inicializace gradientu
        grad_vec = zeros(Float32, 3)
        
        # Výpočet parciální derivace ve směru x
        if i == 1  # Levý okraj - dopředná diference
            grad_vec[1] = fine_sdf[i+1, j, k] - fine_sdf[i, j, k]
        elseif i == nx  # Pravý okraj - zpětná diference
            grad_vec[1] = fine_sdf[i, j, k] - fine_sdf[i-1, j, k]
        else  # Vnitřní bod - centrální diference
            grad_vec[1] = (fine_sdf[i+1, j, k] - fine_sdf[i-1, j, k]) / 2.0f0
        end
        
        # Výpočet parciální derivace ve směru y
        if j == 1  # Přední okraj - dopředná diference
            grad_vec[2] = fine_sdf[i, j+1, k] - fine_sdf[i, j, k]
        elseif j == ny  # Zadní okraj - zpětná diference
            grad_vec[2] = fine_sdf[i, j, k] - fine_sdf[i, j-1, k]
        else  # Vnitřní bod - centrální diference
            grad_vec[2] = (fine_sdf[i, j+1, k] - fine_sdf[i, j-1, k]) / 2.0f0
        end
        
        # Výpočet parciální derivace ve směru z
        if k == 1  # Spodní okraj - dopředná diference
            grad_vec[3] = fine_sdf[i, j, k+1] - fine_sdf[i, j, k]
        elseif k == nz  # Horní okraj - zpětná diference
            grad_vec[3] = fine_sdf[i, j, k] - fine_sdf[i, j, k-1]
        else  # Vnitřní bod - centrální diference
            grad_vec[3] = (fine_sdf[i, j, k+1] - fine_sdf[i, j, k-1]) / 2.0f0
        end
        
        # Inicializace normalizovaného vektoru
        grad_vec_norm = zeros(Float32, 3)
        
        # Normalizace vektoru pro získání jednotkového vektoru
        norm_val = norm(grad_vec)
        if norm_val > 1.0f-10  # Prevence dělení nulou
            grad_vec_norm = grad_vec / norm_val
        end

        gradient_norm[i, j, k] = grad_vec_norm
        gradient[i, j, k] = grad_vec
    end
    
    return gradient, gradient_norm
end

fine_sdf_grad, fine_sdf_grad_norm = compute_gradient(fine_sdf)

"""
    update_sdf_with_planes(fine_sdf::Array{Float32, 3}, fine_sdf_grad::Array{Vector{Float32}, 3}, 
                          plane_sdf::Array{Float32, 3}, plane_normals::Array{Vector{Float32}, 3})

Aktualizuje hodnoty SDF na základě rovinných SDF a podobnosti gradientů.

Argumenty:
- `fine_sdf`: Původní pole SDF hodnot
- `fine_sdf_grad`: Pole gradientů SDF (jednotkové vektory)
- `plane_sdf`: Pole SDF hodnot vypočtených pro roviny
- `plane_normals`: Pole normálových vektorů rovin

Vrací:
- `updated_sdf`: Aktualizované pole SDF hodnot
"""
function update_sdf_with_planes(fine_sdf::Array{Float32, 3}, fine_sdf_grad::Array{Vector{Float32}, 3}, fine_sdf_grad_norm::Array{Vector{Float32}, 3},
                               plane_sdf::Array{Float32, 3}, plane_normals::Array{Vector{Float32}, 3})
    # Ověření konzistence velikostí polí
    if size(fine_sdf) != size(fine_sdf_grad) || size(fine_sdf) != size(plane_sdf) || size(fine_sdf) != size(plane_normals)
        error("Všechna vstupní pole musí mít stejné rozměry")
    end
    
    # Vytvoření kopie původního SDF pro aktualizaci
    updated_sdf = copy(fine_sdf)
    
    # Hraniční hodnota pro plane_sdf (hodnoty vyšší než tato budou ignorovány)
    max_valid_plane_sdf = 1.0f9
    
    # Procházení všech uzlů v síti
    nx, ny, nz = size(fine_sdf)
    for i in 1:nx, j in 1:ny, k in 1:nz
        # Kontrola, zda uzel má platnou plane_sdf hodnotu (není předalokovaná vysoká hodnota)
        if plane_sdf[i, j, k] < max_valid_plane_sdf && !isempty(plane_normals[i, j, k])
            # Získání vektorů pro výpočet podobnosti
            gradient = fine_sdf_grad[i, j, k]
            gradient_norm = fine_sdf_grad_norm[i, j, k]
            normal = plane_normals[i, j, k]

            gradient_magnitude = maximum(norm.(fine_sdf_grad))
            
            # Pokud jsou oba vektory nenulové, vypočítej podobnost
            if norm(gradient_norm) > 1.0f-10 && norm(normal) > 1.0f-10
                # Výpočet podobnosti jako absolutní hodnoty normalizovaného skalárního součinu
                # (Oba vektory by již měly být jednotkové, ale pro jistotu normalizujeme)
                similarity = abs(dot(gradient_norm / norm(gradient_norm), normal / norm(normal)))
                
                # Získání původního znaménka z fine_sdf
                original_sign = sign(fine_sdf[i, j, k])
                
                # Výpočet nové hodnoty SDF
                new_sdf_value = similarity * plane_sdf[i, j, k] * original_sign * gradient_magnitude + abs(similarity-1) * fine_sdf[i, j, k]
                
                # Aktualizace hodnoty SDF
                updated_sdf[i, j, k] = new_sdf_value
            end
        end
    end
    
    return updated_sdf
end

updated_sdf = update_sdf_with_planes(fine_sdf, fine_sdf_grad, fine_sdf_grad_norm, plane_sdf, plane_normals)

export_to_paraview_with_vtk(fine_grid, updated_sdf, "sdf_kontrola_1")

fine_sdf_grad[4,5,6]
plane_normals[4,5,6]
fine_sdf[4,5,6]
plane_sdf[4,5,6]
updated_sdf[4,5,6]
# 0,19,2

similarity = abs(dot(gradient / norm(gradient), normal / norm(normal)))

gradient = fine_sdf_grad[4,5,6]
normal = plane_normals[4,5,6]
gradient = [0,0,1]
normal = [0,0,-1]
