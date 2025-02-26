tet = [2363, 2476, 2399, 2372]
tet_id = findall(x -> x == tet, mesh.IEN)
f = [ mesh.node_sdf[i] for i in tet ]

f = [mesh.node_sdf[i] for i in tet]
    
# Find indices of positive, zero, and negative nodes
pos_locals = [i for (i, s) in enumerate(f) if s > 0]  # strictly positive
zero_locals = [i for (i, s) in enumerate(f) if s == 0]  # zero value
neg_locals = [i for (i, s) in enumerate(f) if s < 0]  # negative

# Get global indices for positive nodes and zero node
pos_global = [tet[i] for i in pos_locals]
zero_global = tet[zero_locals[1]]  # We know there's exactly one zero node

# Calculate intersections between positive and negative nodes
intersection_map = Dict{Tuple{Int,Int},Int64}()
for i in pos_locals, j in neg_locals
    key_sorted = sort([tet[i], tet[j]])
    key_tuple = (key_sorted[1], key_sorted[2])
    if !haskey(intersection_map, key_tuple)
        ip = interpolate_zero(mesh.X[tet[i]], mesh.X[tet[j]], f[i], f[j], mesh)
        intersection_map[key_tuple] = ip
    end
end

# Collect intersection points and add the zero node
ips = collect(values(intersection_map))
push!(ips, zero_global)

# We should have exactly 3 points on the isosurface
if length(ips) != 3
    @warn "Expected 3 points on isosurface (2 intersections + 1 zero node), got $(length(ips))"
    return Vector{Vector{Int64}}()
end

#________________________
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
my_sortperm(angles)
sorted = sortperm(angles)
#________________________

# Order the intersection points (may need to modify order_intersection_points for 3 points)
ordered_ips = order_intersection_points(mesh, ips)

# Create two new tetrahedra using the first positive node
# Note: We only use pos_global[1] since we don't need pos2
pos1 = pos_global[1]
tet1 = [pos1, ordered_ips[3], ordered_ips[2], ordered_ips[1]]
tet2 = [pos1, ordered_ips[1], ordered_ips[3], pos_global[2]]

# Představme si, že máme úhly (v radiánech) pro 4 body
angles = [0.5, -2.2, 1.2, 0.1]

# Původní indexy průsečíků
ips = [101, 102, 103, 104]

# sortperm(angles) vrátí indexy pro seřazení úhlů od nejmenšího po největší
sorted_indices = sortperm(angles)
# Výsledek bude: [2, 4, 1, 3]
# Protože: angles[2] = -0.2 (nejmenší)
#          angles[4] = 0.1
#          angles[1] = 0.5
#          angles[3] = 1.2 (největší)

# Když použijeme tyto indexy na původní pole ips
sorted_ips = ips[sorted_indices]
# Výsledek bude: [102, 104, 101, 103]

# Pro lepší pochopení si ukažme všechny hodnoty vedle sebe:
for (i, sorted_i) in enumerate(sorted_indices)
    println("Původní index: $sorted_i")
    println("Úhel: $(angles[sorted_i])")
    println("ID průsečíku: $(ips[sorted_i])")
    println("---")
end

function my_sortperm(arr)
    # Vytvoříme pole párů (hodnota, původní_index)
    pairs = [(value, i) for (i, value) in enumerate(arr)]
    
    # Implementace bubble sort pro demonstraci
    n = length(pairs)
    for i in 1:n-1
        for j in 1:n-i
            # Porovnáváme hodnoty (první prvek páru)
            if pairs[j][1] > pairs[j+1][1]
                # Prohodíme páry
                pairs[j], pairs[j+1] = pairs[j+1], pairs[j]
            end
        end
    end
    
    # Vrátíme seřazené indexy (druhý prvek páru)
    return [pair[2] for pair in pairs]
end