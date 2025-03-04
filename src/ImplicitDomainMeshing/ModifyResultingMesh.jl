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
    
    # Výpočet dot produktu pro určení strany roviny
    dot_product = dot(plane.normal, vec_to_point)
    
    # Vzdálenost k nekonečné rovině (absolutní hodnota pro kontroly)
    dist_to_infinite_plane = abs(dot_product)
    
    # Projekce bodu na nekonečnou rovinu
    projected_point = point - dot_product * plane.normal
    
    # Kontrola, zda projektovaný bod leží v rámci omezené roviny
    if is_on_plane(plane, projected_point)
        # Pokud projekce leží v omezené rovině, vrátíme vzdálenost k nekonečné rovině
        # se znaménkem - záporné ve směru normály, kladné v opačném směru
        return dot_product > 0 ? -dist_to_infinite_plane : dist_to_infinite_plane
    else
        # Pokud projekce neleží v omezené rovině, vrátíme vysokou hodnotu
        # se zachováním správného znaménka
        return dot_product > 0 ? -1.0e10 : 1.0e10
    end
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
            
            # Porovnání absolutních hodnot, ale uložení znaménkové hodnoty
            if abs(dist) <= max_dist && abs(dist) < abs(plane_sdf[i, j, k])
                plane_sdf[i, j, k] = dist
                plane_normals[i, j, k] = Vector{Float32}(plane_definitions[p_idx].normal)
            end
        end
    end
    
    return plane_sdf, plane_normals
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

# export_to_paraview_with_vtk(fine_grid, plane_sdf, "sdf_kontrola")



# Funkce pro vyhodnocení vzdálenosti bodu od rovin (planes_sdf)
function eval_planes_sdf(mesh::BlockMesh, p::SVector{3,Float64}, plane_definitions::Vector{PlaneDefinition})
    # Vytvoříme omezené roviny z definic
    planes = [BoundedPlane(def.normal, def.point, def.shape) for def in plane_definitions]
    
    # Inicializace s vysokou kladnou hodnotou
    min_dist = 1.0e10
    
    # Kontrola vzdálenosti ke každé rovině
    for plane in planes
        # Převod na Vector{Float32} pro kompatibilitu s funkcí distance_to_bounded_plane
        p_float32 = Vector{Float32}([p[1], p[2], p[3]])
        dist = distance_to_bounded_plane(plane, p_float32)
        
        # Aktualizace minimální vzdálenosti (porovnáváme absolutní hodnoty)
        if abs(dist) < abs(min_dist)
            min_dist = dist
        end
    end
    
    return min_dist
end

# Aproximace gradientu planes_sdf pro určení směru warpu
function approximate_planes_gradient(mesh::BlockMesh, p::SVector{3,Float64}, plane_definitions::Vector{PlaneDefinition}; h::Float64=1e-3)
    dx = SVector{3,Float64}(h, 0.0, 0.0)
    dy = SVector{3,Float64}(0.0, h, 0.0)
    dz = SVector{3,Float64}(0.0, 0.0, h)
    
    # Výpočet parciálních derivací pomocí centrálních diferencí
    df_dx = (eval_planes_sdf(mesh, p + dx, plane_definitions) - eval_planes_sdf(mesh, p - dx, plane_definitions)) / (2 * h)
    df_dy = (eval_planes_sdf(mesh, p + dy, plane_definitions) - eval_planes_sdf(mesh, p - dy, plane_definitions)) / (2 * h)
    df_dz = (eval_planes_sdf(mesh, p + dz, plane_definitions) - eval_planes_sdf(mesh, p - dz, plane_definitions)) / (2 * h)
    
    return SVector{3,Float64}(df_dx, df_dy, df_dz)
end

# Funkce pro posunutí uzlu na nulovou hladinu planes_sdf
function warp_node_to_planes_isocontour!(mesh::BlockMesh, node_index::Int, plane_definitions::Vector{PlaneDefinition}, max_iter::Int)
    tol = mesh.grid_tol
    current_position = mesh.X[node_index]
    
    for iter in 1:max_iter
        # Vyhodnocení planes_sdf v aktuální pozici
        f = eval_planes_sdf(mesh, current_position, plane_definitions)
        
        # Pokud jsme dostatečně blízko izopovrchu, končíme
        abs2(f) < tol * tol && break
        
        # Výpočet gradientu pro směr posunu
        grad = approximate_planes_gradient(mesh, current_position, plane_definitions)
        norm_grad_squared = sum(abs2, grad)
        
        # Pokud je gradient příliš malý, končíme
        norm_grad_squared < 1e-16 && break
        
        # Newtonův krok
        dp = (f / norm_grad_squared) * grad
        current_position -= dp
    end
    
    # Výpočet aktuální planes_sdf hodnoty po warping
    current_sdf = eval_planes_sdf(mesh, current_position, plane_definitions)
    
    # Aktualizace pozice uzlu a jeho SDF hodnoty
    mesh.X[node_index] = current_position
    if abs(current_sdf) < tol*2
        mesh.node_sdf[node_index] = 0.0
    else
        mesh.node_sdf[node_index] = current_sdf
    end
end

# Hlavní funkce pro modifikaci sítě podle planes_sdf
function warp_mesh_by_planes_sdf!(mesh::BlockMesh, plane_definitions::Vector{PlaneDefinition}; max_iter::Int=20)
    # Nejprve najdeme uzly s negativní planes_sdf hodnotou
    negative_indices = Int[]
    for i in 1:length(mesh.X)
        plane_sdf = eval_planes_sdf(mesh, mesh.X[i], plane_definitions)
        if plane_sdf < 0
            push!(negative_indices, i)
        end
    end
    
    @info "Nalezeno $(length(negative_indices)) uzlů s negativní planes_sdf hodnotou"
    
    # Warp uzlů s negativní planes_sdf hodnotou na izopovrch
    for node_idx in negative_indices
        warp_node_to_planes_isocontour!(mesh, node_idx, plane_definitions, max_iter)
    end
    
    # Aktualizace topologie sítě
    update_connectivity!(mesh)
    
    # Finální úpravy
    slice_ambiguous_tetrahedra!(mesh)
    update_connectivity!(mesh)
    
    @info "Úprava sítě podle planes_sdf dokončena"
end

warp_mesh_by_planes_sdf!(mesh, plane_definitions)

export_mesh_vtk(mesh, "block-mesh-cut.vtu")