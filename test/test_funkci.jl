using LinearAlgebra, Statistics, Makie, GLMakie

function round_down_to_even(number::Int)
    # Zaokrouhlení dolů na nejbližší celé číslo
    rounded = floor(Int, number)
    # Pokud je číslo sudé, vrátíme ho, jinak odčteme 1 pro dosažení sudého čísla
    return iseven(rounded) ? rounded : rounded - 1
end

function TestGeometrySphere(max_elements::Int64)
    # Inicializace proměnných: poloměr koule, střed a počet uzlů
    radius = 1.0
    center = [0., 0., 0.]
    n = max_elements
    delta = 2 * radius / n
    step = Int(round_down_to_even(n)/2)

    # Inicializace datových struktur
    nodes = Vector{Vector{Float64}}()
    elements = Vector{Vector{Int}}()
    densities = Float64[]
    elements_center = []
    node_map = Dict()
    new_node_id = 1

    for i in -step:step
        for j in -step:step
            for k in -step:step
                x = i * delta
                y = j * delta
                z = k * delta

                if norm([x, y, z]) <= radius
                    push!(nodes, [x, y, z])
                    node_map[(i, j, k)] = new_node_id
                    new_node_id += 1
                end
            end
        end
    end

    for i in -step:step-1
        for j in -step:step-1
            for k in -step:step-1
                corners = [(i, j, k), (i+1, j, k), (i+1, j+1, k), (i, j+1, k), 
                           (i, j, k+1), (i+1, j, k+1), (i+1, j+1, k+1), (i, j+1, k+1)]
                
                if all(c -> c in keys(node_map), corners)
                    element_nodes = [node_map[c] for c in corners]
                    push!(elements, element_nodes)

                    element_center_coords = [nodes[node_id] for node_id in element_nodes]
                    element_center = mean(hcat(element_center_coords...), dims=2)
                    element_center = vec(element_center)

                    density = 1.0 - norm(element_center - center) / radius
                    push!(densities, density)
                    push!(elements_center, element_center)
                end
            end
        end
    end
    
    # Unikátní čísla uzlů pro konstrukci elementů
    flattened = vcat(elements...)
    sorted = sort(unique(flattened))

    # Přečíslování ID uzlů v elements
    replace_with_position!(v, sorted) = [findfirst(isequal(x), sorted) for x in v]
    replaced_elements = [replace_with_position!(v, sorted) for v in elements]

    # Výběr pouze uzlů které slouží pro konstrukci elementů
    used_nodes = [nodes[idx] for idx in sorted]

    return used_nodes, replaced_elements, densities, elements_center
end


nodes, elements, densities, elements_center = TestGeometrySphere(4)

using Base.Threads
using LinearAlgebra

# function Vector2Vectors(nodes::Vector{Any})
#     xₙ = zeros(length(nodes))
#     yₙ = zeros(length(nodes))
#     zₙ = zeros(length(nodes))

#     for i in 1:length(nodes)
#         xₙ[i], yₙ[i], zₙ[i] = nodes[i]
#     end
#     return xₙ, yₙ, zₙ
# end
# xₙ, yₙ, zₙ = Vector2Vectors(elements_center)

nel = length(elements_center)
ps = rand(Point3f, nel)
for i in 1:nel
    ps[i] = elements_center[i]
end
scatter(ps, color = densities)


function draw_elements_makie(nodes, elements)
    fig = Figure()
    ax = Axis3(fig[1, 1], aspect = :data)

    # Zde bude uložena všechna souřadnice uzlů pro scatter plot
    all_node_coords = []

    for element in elements
        # Získání souřadnic pro každý uzel v elementu
        element_coords = [nodes[node_id] for node_id in element]

        # Přidání souřadnic do seznamu pro scatter plot
        append!(all_node_coords, element_coords)

        # Vykreslení hran elementu
        for i in 1:4
            # Dolní čtverec
            lines!(ax, [element_coords[i][1], element_coords[i % 4 + 1][1]], 
                      [element_coords[i][2], element_coords[i % 4 + 1][2]], 
                      [element_coords[i][3], element_coords[i % 4 + 1][3]], color=:blue)
            # Horní čtverec
            lines!(ax, [element_coords[i+4][1], element_coords[(i % 4 + 1) + 4][1]], 
                      [element_coords[i+4][2], element_coords[(i % 4 + 1) + 4][2]], 
                      [element_coords[i+4][3], element_coords[(i % 4 + 1) + 4][3]], color=:blue)
            # Svislé hrany
            lines!(ax, [element_coords[i][1], element_coords[i+4][1]], 
                      [element_coords[i][2], element_coords[i+4][2]], 
                      [element_coords[i][3], element_coords[i+4][3]], color=:blue)
        end
    end

    # Vykreslení všech uzlů (bodů)
    scatter!(ax, [p[1] for p in all_node_coords], 
                 [p[2] for p in all_node_coords], 
                 [p[3] for p in all_node_coords], color=:red, markersize=20)

    fig
end

# Volání funkce s vašimi výsledky (nodes, elements)
draw_elements_makie(nodes, elements)


function check_normal_orientation(elements::Vector{Vector{Int}}, nodes::Vector{Vector{Float64}})
    for element in elements
        # Souřadnice vrcholů podstavy a horní části
        bottom_vertices = [nodes[element[1]], nodes[element[2]], nodes[element[3]], nodes[element[4]]]
        top_vertices = [nodes[element[5]], nodes[element[6]], nodes[element[7]], nodes[element[8]]]

        # Výpočet vektorů podstavy
        edge1 = bottom_vertices[2] - bottom_vertices[1]
        edge2 = bottom_vertices[4] - bottom_vertices[1]

        # Výpočet normály
        normal = cross(edge1, edge2)

        # Kontrola, zda normála směřuje nahoru (předpokládáme, že Z je vertikální osa)
        if normal[3] <= 0
            error("Normála podstavy nesměřuje směrem k horní části elementu.")
        end
    end
    println("Všechny elementy mají správně orientovanou normálu podstavy.")
end

check_normal_orientation(elements, nodes)

element = elements[1]
# Souřadnice vrcholů podstavy a horní části
bottom_vertices = [nodes[element[1]], nodes[element[2]], nodes[element[3]], nodes[element[4]]]
top_vertices = [nodes[element[5]], nodes[element[6]], nodes[element[7]], nodes[element[8]]]

# Výpočet vektorů podstavy
edge1 = bottom_vertices[2] - bottom_vertices[1]
edge2 = bottom_vertices[end] - bottom_vertices[1]
# Výpočet normály
normal = cross(edge1, edge2)

test = Vector{Vector{Float64}}(undef, length(bottom_vertices))
for i in 1:length(bottom_vertices)
    test[i] = bottom_vertices[i] .+ normal*1e-6
end

mean(mean.(bottom_vertices - top_vertices)) < mean(mean.(test - top_vertices))


element = elements[1]

# Výpočet vektorů a normály
edge1 = nodes[element[2]] - nodes[element[1]]
edge2 = nodes[element[4]] - nodes[element[1]]
normal = cross(edge1, edge2)

# Vytvoření pole 'posuv' pomocí vektorizované operace
posuv = [node .+ normal * 1e-6 for node in nodes[element[1:4]]]

# Výpočet a porovnání průměrných hodnot
mean_diff_bottom = mean(map(v -> norm(v), nodes[element[1:4]] .- nodes[element[5:8]]))
mean_diff_posuv = mean(map(v -> norm(v), posuv .- nodes[element[5:8]]))

mean_diff_bottom > mean_diff_posuv
a = element[1]


typeof(elements) == "Vector{Int64}"

typeof(nodes) == Vector{Vector{Float64}}


# function GeometryTest(elements::Any, nodes::Vector{Vector{Float64}})
    # Kontrola typu
    element = []
    if typeof(elements) == Vector{Vector{Int64}}
        element = elements[1]
    elseif typeof(elements) == Vector{Int64}
        element = elements
    else
        println("Wrong data type of elements:", typeof(elements))
    end

    # Výpočet vektorů a normály
    edge1 = nodes[element[2]] - nodes[element[1]]
    edge2 = nodes[element[4]] - nodes[element[1]]
    normal = cross(edge1, edge2)

    # Vytvoření pole 'posuv' pomocí vektorizované operace
    posuv = [node .+ normal * 1e-6 for node in nodes[element[1:4]]]

    # Výpočet a porovnání průměrných hodnot
    mean_diff_bottom = mean(map(v -> norm(v), nodes[element[1:4]] .- nodes[element[5:8]]))
    mean_diff_posuv = mean(map(v -> norm(v), posuv .- nodes[element[5:8]]))
    return mean_diff_bottom, mean_diff_posuv
end
mean_diff_bottom, mean_diff_posuv = GeometryTest(elements, nodes)
@test mean_diff_bottom > mean_diff_posuv

edge1 = nodes[element[2]] - nodes[element[1]]
edge2 = nodes[element[4]] - nodes[element[1]]