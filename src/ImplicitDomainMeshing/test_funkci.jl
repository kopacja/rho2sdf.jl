# ----------------------------
# Pomocná funkce: Lineární interpolace průsečíku s nulovou hladinou SDF
# ----------------------------
function interpolate_zero(p1::SVector{3,Float64}, p2::SVector{3,Float64},
                            f1::Float64, f2::Float64, mesh::BlockMesh)::Int64
    tol = mesh.grid_tol
    t = (0 - f1) / (f2 - f1)
    p_interp = p1 + t * (p2 - p1)
    p_key = quantize(p_interp, tol)
    if haskey(mesh.node_hash, p_key)
        return mesh.node_hash[p_key]
    else
        push!(mesh.X, p_interp)
        push!(mesh.node_sdf, 0.0)
        new_index = length(mesh.X)
        mesh.node_hash[p_key] = new_index
        return new_index
    end
end

# Zjistí, zda se bod nachází na hranici mřížky
function is_on_boundary(mesh::BlockMesh, p::SVector{3,Float64})
    tol = mesh.grid_tol
    vmin = mesh.grid[1,1,1]
    vmax = mesh.grid[end,end,end]
    return any(abs(p[i]-vmin[i]) < tol || abs(p[i]-vmax[i]) < tol for i in 1:3)
  end
  
  # Uspořádá dané průsečíky (indexy v mesh.X) do cyklického pořadí
  function order_intersection_points(mesh::BlockMesh, ips::Vector{Int64})::Vector{Int64}
    tol = mesh.grid_tol
    pts = [mesh.X[i] for i in ips]
    centroid = sum(pts) / length(pts)
    v1 = pts[1] - centroid
    if norm(v1) < tol
       v1 = SVector(1.0, 0.0, 0.0)
    end
    v2 = cross(normalize(v1), SVector(0.0, 0.0, 1.0))
    if norm(v2) < tol
        v2 = cross(normalize(v1), SVector(0.0, 1.0, 0.0))
    end
    angles = [atan(dot(pt - centroid, normalize(v1)), dot(pt - centroid, normalize(v2))) for pt in pts]
    sorted = sortperm(angles)
    return ips[sorted]
  end
  
  # Určí "barvu" hrany – zda oba uzly mají SDF kladné (positive), záporné (negative) nebo smíšené
  function edge_color(mesh::BlockMesh, i::Int, j::Int)
    f1 = mesh.node_sdf[i]
    f2 = mesh.node_sdf[j]
    if f1 >= 0 && f2 >= 0
        return :positive
    elseif f1 < 0 && f2 < 0
        return :negative
    else
        return :mixed
    end
  end
  
  # Výpočet nákladové funkce pro volbu diagonály (penalizace, pokud nejsou obě koncové uzly "positive")
  function diagonal_cost(mesh::BlockMesh, pos_indices::Vector{Int64}, diag_pair::Tuple{Int64,Int64})
    p_diag1 = mesh.X[diag_pair[1]]
    p_diag2 = mesh.X[diag_pair[2]]
    cost = 0.0
    col = edge_color(mesh, diag_pair[1], diag_pair[2])
    penalty = (col == :positive ? 0.0 : 1000.0)
    cost += penalty
    for pi in pos_indices
        p_pos = mesh.X[pi]
        cost += norm(p_pos - p_diag1) + norm(p_pos - p_diag2)
    end
    return cost
  end

  
  function stencil_1p3(mesh::BlockMesh, tet::Vector{Int64})
    tol = mesh.grid_tol
    # Získáme SDF hodnoty pro všechny čtyři vrcholy
    f = [ mesh.node_sdf[i] for i in tet ]
    signs = map(x -> x >= 0, f)
    # Definujeme lokální toleranci – například 20% kroku mřížky
    # local_threshold = tol*2000000
    local_threshold = 0.2
  
    # Najdeme lokální index vrcholu s kladnou hodnotou
    pos_local = findfirst(identity, signs)
    # Pokud je SDF kladného vrcholu příliš malé (tj. blízko nule), tetraedr nevytváříme
    if abs(f[pos_local]) < local_threshold
        return Vector{Vector{Int64}}()
    end
  
    # Zbylé (tři) vrcholy mají zápornou hodnotu; pokud některý záporný vrchol je také příliš blízko izokontuře,
    # rovněž tetraedr nevytváříme.
    neg_locals = [ i for (i, s) in enumerate(signs) if !s ]
    for n in neg_locals
        if abs(f[n]) < local_threshold
            return Vector{Vector{Int64}}()
        end
    end
  
    # Určíme pozici a hodnotu kladného vrcholu
    p_pos = mesh.X[tet[pos_local]]
    f_pos = f[pos_local]
  
    # Vypočítáme průsečíky na hranách spojujících kladný vrchol s každým záporným
    ips = Vector{Int64}()
    for n in neg_locals
        p_neg = mesh.X[tet[n]]
        f_neg = f[n]
        ip = interpolate_zero(p_pos, p_neg, f_pos, f_neg, mesh)
        push!(ips, ip)
    end
  
    # Vytvoříme nový tetraedr: [kladný vrchol, I1, I2, I3]
    return [ [ tet[pos_local], ips[1], ips[2], ips[3] ] ]
  end
  

  function stencil_3p1(mesh::BlockMesh, tet::Vector{Int64})::Vector{Vector{Int64}}
    f = [ mesh.node_sdf[i] for i in tet ]
    # Najdeme lokální index záporného vrcholu
    neg_local = findfirst(x -> x < 0, f)
    pos_locals = [ i for (i, s) in enumerate(f) if s >= 0 ]
    # Vypočítáme průsečíky na hranách spojujících záporný vrchol s každým kladným
    ips = Dict{Int,Int64}()
    for p in pos_locals
        ip = interpolate_zero(mesh.X[tet[neg_local]], mesh.X[tet[p]], f[neg_local], f[p], mesh)
        ips[p] = ip
    end
    # Předpokládáme, že existují právě tři kladné vrcholy – označíme je jako a, b, c
    a, b, c = pos_locals
    # Vytvoříme tři tetraedry podle šablony:
    tet1 = [ tet[a], tet[b], ips[b], ips[a] ]
    tet2 = [ tet[b], tet[c], ips[c], ips[b] ]
    tet3 = [ tet[c], tet[a], ips[a], ips[c] ]
    return [ tet1, tet2, tet3 ]
  end
  

  function stencil_2p2_variantA(mesh::BlockMesh, tet::Vector{Int64})::Vector{Vector{Int64}}
    # Ambiguózní případ: 2 kladné a 2 záporné vrcholy
    f = [ mesh.node_sdf[i] for i in tet ]
    pos_locals = [ i for (i, s) in enumerate(f) if s >= 0 ]
    neg_locals = [ i for (i, s) in enumerate(f) if s < 0 ]
    pos_global = [ tet[i] for i in pos_locals ]
    # Vypočteme průsečíky na hranách spojujících každý kladný a záporný vrchol
    intersection_map = Dict{Tuple{Int,Int},Int64}()
    for i in pos_locals, j in neg_locals
        key_sorted = sort([tet[i], tet[j]])
        key_tuple = (key_sorted[1], key_sorted[2])
        if !haskey(intersection_map, key_tuple)
            ip = interpolate_zero(mesh.X[tet[i]], mesh.X[tet[j]], f[i], f[j], mesh)
            intersection_map[key_tuple] = ip
        end
    end
    ips = collect(values(intersection_map))
    if length(ips) != 4
        @warn "Ambiguózní tetraedr nedal 4 unikátní průsečíky."
        return Vector{Vector{Int64}}()
    end
    # Seřadíme průsečíky do cyklického pořadí
    ordered_ips = order_intersection_points(mesh, ips)
    # Volba diagonály na základě cost funkce
    diag1 = (ordered_ips[1], ordered_ips[3])
    diag2 = (ordered_ips[2], ordered_ips[4])
    cost1 = diagonal_cost(mesh, pos_global, diag1)
    cost2 = diagonal_cost(mesh, pos_global, diag2)
    if cost1 > cost2
        # Přeuspořádáme ordered_ips tak, aby první prvek odpovídal vybrané diagonále
        idx = findfirst(x -> x == diag2[1], ordered_ips)
        ordered_ips = vcat(ordered_ips[idx:end], ordered_ips[1:idx-1])
    end
    # Sestavíme tři tetraedry dle standardní varianty:
    pos1 = pos_global[1]
    pos2 = pos_global[2]
    tet1 = [ pos1, ordered_ips[1], ordered_ips[2], ordered_ips[3] ]
    tet2 = [ pos2, ordered_ips[1], ordered_ips[3], ordered_ips[4] ]
    tet3 = [ pos1, pos2, ordered_ips[1], ordered_ips[3] ]
    return [ tet1, tet2, tet3 ]
  end

  
  function stencil_2p2_variantB(mesh::BlockMesh, tet::Vector{Int64})::Vector{Vector{Int64}}
    # Stejný postup jako v variantě A pro výpočet průsečíků
    f = [ mesh.node_sdf[i] for i in tet ]
    pos_locals = [ i for (i, s) in enumerate(f) if s >= 0 ]
    neg_locals = [ i for (i, s) in enumerate(f) if s < 0 ]
    pos_global = [ tet[i] for i in pos_locals ]
    intersection_map = Dict{Tuple{Int,Int},Int64}()
    for i in pos_locals, j in neg_locals
        key_sorted = sort([tet[i], tet[j]])
        key_tuple = (key_sorted[1], key_sorted[2])
        if !haskey(intersection_map, key_tuple)
            ip = interpolate_zero(mesh.X[tet[i]], mesh.X[tet[j]], f[i], f[j], mesh)
            intersection_map[key_tuple] = ip
        end
    end
    ips = collect(values(intersection_map))
    if length(ips) != 4
        @warn "Ambiguózní tetraedr nedal 4 unikátní průsečíky."
        return Vector{Vector{Int64}}()
    end
    ordered_ips = order_intersection_points(mesh, ips)
    # Varianta B: jednodušší dělení na dvě tetraedra
    pos1 = pos_global[1]
    pos2 = pos_global[2]
    tet1 = [ pos1, ordered_ips[1], ordered_ips[2], ordered_ips[3] ]
    tet2 = [ pos2, ordered_ips[1], ordered_ips[3], ordered_ips[4] ]
    return [ tet1, tet2 ]
  end

  
  function apply_stencil(mesh::BlockMesh, tet::Vector{Int64})
    f = [ mesh.node_sdf[i] for i in tet ]
    np = count(x -> x >= 0, f)
    if np == 1
        return stencil_1p3(mesh, tet)
    elseif np == 3
        # return stencil_3p1(mesh, tet)
        return Vector{Vector{Int64}}()
    elseif np == 2
        # Rozlišujeme variantu podle toho, zda se prvek dotýká hranice
        # on_boundary = any(is_on_boundary(mesh, mesh.X[i]) for i in tet)
        # if on_boundary
        #     return stencil_2p2_variantB(mesh, tet)
        # else
        #     return stencil_2p2_variantA(mesh, tet)
        # end
        return Vector{Vector{Int64}}()
    elseif np == 4
        return [ tet ]  # Všechny vrcholy kladné – tetraedr se ponechá
    else
        return Vector{Vector{Int64}}()  # Všechny záporné – tetraedr vynecháme
    end
  end
  