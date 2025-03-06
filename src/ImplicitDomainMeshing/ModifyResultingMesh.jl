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
    # Kontrola, zda jsou předány nějaké definice rovin
    if isempty(plane_definitions)
        @info "Žádné definice rovin nebyly předány, přeskakuji warping podle rovin."
        return
    end
    
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
    # slice_ambiguous_tetrahedra!(mesh)
    # update_connectivity!(mesh)
    
    @info "Úprava sítě podle planes_sdf dokončena"
end

#TODO: nejdřív warp a pak řez elementů
