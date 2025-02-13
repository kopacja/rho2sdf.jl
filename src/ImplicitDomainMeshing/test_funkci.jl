# ----------------------------
# Pomocná funkce: Lineární interpolace průsečíku s nulovou hladinou SDF
# ----------------------------
function interpolate_zero(p1::SVector{3,Float64}, p2::SVector{3,Float64},
                          f1::Float64, f2::Float64, mesh::BlockMesh; tol=mesh.grid_tol, max_iter=20)::Int64
                          # coord kladného, coords zápor, sdf kladného, sdf záporného
                          
    pos1 = f1 >= -tol
    pos2 = f2 >= -tol
    sdf_p1 = eval_sdf(mesh, p1)
    sdf_p2 = eval_sdf(mesh, p2)
    # println("kontrola p1 sdf: ", sdf_p1)
    # println("kontrola p2 sdf: ", sdf_p2)
    # Pokud oba body mají stejnou polaritu podle tolerance, nelze správně interpolovat.
    if pos1 == pos2
        error("Oba body mají stejnou 'toleranční' polaritu; jeden bod musí být blízko nuly (kladný) a druhý výrazně záporný.")
    end

    # Inicializujeme interval: low a high – předpokládáme, že p1 a p2 jsou uspořádané podle f
    low, high = p1, p2
    f_low, f_high = f1, f2
    mid = low
    for iter in 1:max_iter
        mid = (low + high) / 2.0
        f_mid = eval_sdf(mesh, mid)
        # Pokud je hodnota dostatečně blízko nule, ukončíme iteraci
        if abs(f_mid) < tol
            break
        end
        # Podle znaménka f_mid aktualizujeme jeden z koncových bodů intervalu
        if sign(f_mid) == sign(f_low)
            low, f_low = mid, f_mid
        else
            high, f_high = mid, f_mid
        end
    end

    # Kvantizujeme nalezený bod, abychom se vyhnuli duplicitám v hashtable
    p_key = quantize(mid, tol)
    sdf_of_iterp_point = eval_sdf(mesh, SVector{3, Float64}(p_key))
    # println("kontrola interp sdf: ", sdf_of_iterp_point)
    if haskey(mesh.node_hash, p_key)
        return mesh.node_hash[p_key]
    else
        push!(mesh.X, mid)
        push!(mesh.node_sdf, 0.0)  # nebo nastavíme přesně na 0
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


function stencil_1p3(mesh::BlockMesh, tet::Vector{Int64})::Vector{Vector{Int64}}
    tol = mesh.grid_tol
    # Get SDF values for all vertices of the tetrahedron
    f = [mesh.node_sdf[i] for i in tet]
    
    # Identify vertices by their position relative to the isosurface
    plus = map(x -> x >= tol, f)           # vertices inside the body
    minus = map(x -> x <= -tol, f)         # vertices outside the body
    # tol_band = map(x -> abs(x) < tol, f)   # vertices on the isosurface
    
    # Find the index of the positive and negative vertex
    pos_idx = findfirst(identity, plus)
    neg_indices = findall(identity, minus)

    # Pokud je SDF kladného vrcholu příliš malé (tj. blízko nule), tetraedr nevytváříme
    if isnothing(pos_idx)
        return Vector{Vector{Int64}}()
    end
    
    # Create a copy of the original tetrahedron that we'll modify
    new_tet = copy(tet)
    
    # For each negative vertex, calculate its new position on the isosurface
    for neg_idx in neg_indices
        # Interpolate between the positive vertex and the current negative vertex
        ip = interpolate_zero(
            mesh.X[tet[pos_idx]],  # coordinates of positive vertex
            mesh.X[tet[neg_idx]],  # coordinates of negative vertex
            f[pos_idx],            # SDF value of positive vertex
            f[neg_idx],            # SDF value of negative vertex
            mesh
        )
        # Replace the negative vertex with its new position on the isosurface
        new_tet[neg_idx] = ip
    end
    
    # Return the modified tetrahedron as a single-element vector
    # (keeping the return type consistent with other stencil functions)
    return [new_tet]
end

function stencil_3p1(mesh::BlockMesh, tet::Vector{Int64})::Vector{Vector{Int64}}
    # Načteme hodnoty SDF pro všechny vrcholy tetraedru
    f = [mesh.node_sdf[i] for i in tet]
    tol = mesh.grid_tol

    # Najdeme záporný vrchol – definujeme ho jako ten, kde f < -tol.
    # (Pokud žádný takový neexistuje, tetraedr zřejmě nemá zápornou stranu a zůstane nezměněn.)
    neg_local = findfirst(x -> x < -tol, f)
    if neg_local === nothing
        return [tet]
    end

    # Určíme počet „strict pozitivních“ vrcholů, tj. těch, kde f >= tol.
    pos_locals = [i for (i, s) in enumerate(f) if s >= tol]
    n_strict_pos = length(pos_locals)

    # Podle počtu strict pozitivních vrcholů rozhodneme, kterou stencil variantu použít:
    if n_strict_pos == 3
        # println(tet)
        # return Vector{Vector{Int64}}()
        # Klasický případ: máme tři vrcholy s f >= tol a jeden záporný vrchol.
        # Pro každou hranu spojující záporný vrchol a jeden z těchto kladných určíme
        # průsečík (interpolovaný nebo přímo, pokud je uzel již blízko 0).
        ips = Dict{Int,Int64}()
        for i in pos_locals
            f_neg = f[neg_local]
            f_pos = f[i]
            # Provedeme iterovanou bisekci mezi vrcholy
            ip = interpolate_zero(mesh.X[tet[neg_local]], mesh.X[tet[i]], f_neg, f_pos, mesh)
            ips[i] = ip

        end
        # V tuto chvíli očekáváme, že pos_locals obsahuje právě 3 indexy.
        # Pokud by tomu tak nebylo (např. kvůli numerickým odchylkám) – nepřerušíme běh, ale přesměrujeme na jinou variantu.
        if n_strict_pos != 3
            println("wtf")
        end
        # Pokud máme přesně 3 strict pozitivní vrcholy, označíme je jako a, b, c.
        a, b, c = pos_locals  # Rozbalení tří indexů

    tet1 = [ tet[a], tet[b], tet[c], ips[a] ]
    tet2 = [ ips[a], ips[b], ips[c], tet[c] ]
    tet3 = [ tet[c], ips[a], tet[b], ips[b] ]
        return [tet1, tet2, tet3]


    elseif n_strict_pos == 2
        # Pokud máme jen 2 strict pozitivní vrcholy, voláme variantu pro 2-positivní.
        # return stencil_2p2_variantA(mesh, tet)
        return Vector{Vector{Int64}}()

    elseif n_strict_pos == 1
        # Pokud máme pouze 1 strict pozitivní vrchol, jedná se o konfiguraci 1p3.
        return stencil_1p3(mesh, tet)
        # return Vector{Vector{Int64}}() # smazat

    else
        # Pokud žádný vrchol nesplňuje f >= tol, tetraedr vynecháme.
        return Vector{Vector{Int64}}()
    end
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
    tol = mesh.grid_tol
    np = count(x -> x >= (-1.1*tol), f)
    # np = count(x -> x >= (0), f)
    if np == 1
        # println(tet)
        # return stencil_1p3(mesh, tet)
        return Vector{Vector{Int64}}()
    elseif np == 3
        # println(tet)
        # return stencil_3p1(mesh, tet)
        return Vector{Vector{Int64}}()
    elseif np == 2
        # Rozlišujeme variantu podle toho, zda se prvek dotýká hranice
        on_boundary = any(is_on_boundary(mesh, mesh.X[i]) for i in tet)
        if on_boundary
            # return stencil_2p2_variantB(mesh, tet)
        else
            # return stencil_2p2_variantA(mesh, tet)
        end
        return Vector{Vector{Int64}}()
    elseif np == 4
        # return [ tet ]  # Všechny vrcholy kladné – tetraedr se ponechá
        # return Vector{Vector{Int64}}() # smazat
    else
        return Vector{Vector{Int64}}()  # Všechny záporné – tetraedr vynecháme
    end
  end
  

  #________________________
#   current_position = mesh.X[1]
#   f = eval_sdf(mesh, current_position)

#   mesh.X
#   length(mesh.SDF)
tet = [58, 52, 56, 55]
tet_id = findall(x -> x == tet, mesh.IEN)
mesh.node_sdf[58]
mesh.node_sdf[52]
mesh.node_sdf[56]
mesh.node_sdf[55]

[2497, 2466, 2468, 2471]

f = [ mesh.node_sdf[i] for i in tet ]
f[4] = 0.5
f[3] = 0.4
f[1] = 0.1

tol = mesh.grid_tol
np = count(x -> x >= (-1.1*tol), f)

tet = [63, 76, 69, 80]
[63, 77, 76, 80]

f = [mesh.node_sdf[i] for i in tet]
    tol = mesh.grid_tol

    # Najdeme záporný vrchol – definujeme ho jako ten, kde f < -tol.
    # (Pokud žádný takový neexistuje, tetraedr zřejmě nemá zápornou stranu a zůstane nezměněn.)
    neg_local = findfirst(x -> x < -tol, f)
    if neg_local === nothing
        return [tet]
    end

    # Určíme počet „strict pozitivních“ vrcholů, tj. těch, kde f >= tol.
    pos_locals = [i for (i, s) in enumerate(f) if s >= tol]
    n_strict_pos = length(pos_locals)

    n_strict_pos == 3

    ips = Dict{Int,Int64}()
        for i in pos_locals
            f_neg = f[neg_local]
            f_pos = f[i]
            # Provedeme iterovanou bisekci mezi vrcholy
            ip = interpolate_zero(mesh.X[tet[neg_local]], mesh.X[tet[i]], f_neg, f_pos, mesh)
            ips[i] = ip

        end
        # V tuto chvíli očekáváme, že pos_locals obsahuje právě 3 indexy.
        # Pokud by tomu tak nebylo (např. kvůli numerickým odchylkám) – nepřerušíme běh, ale přesměrujeme na jinou variantu.
        if n_strict_pos != 3
            println("wtf")
        end
        # Pokud máme přesně 3 strict pozitivní vrcholy, označíme je jako a, b, c.
        a, b, c = pos_locals

