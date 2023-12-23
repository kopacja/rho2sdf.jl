using LinearAlgebra, Statistics, Makie, GLMakie

function round_down_to_even(number::Float64)
    # Zaokrouhlení dolů na nejbližší celé číslo
    rounded = floor(Int, number)
    # Pokud je číslo sudé, vrátíme ho, jinak odčteme 1 pro dosažení sudého čísla
    return iseven(rounded) ? rounded : rounded - 1
end

function TestGeometrySphere(max_elements::Int64)
    radius = 1.0
    center = [0., 0., 0.]

    n = cbrt(max_elements) 
    delta = 2 * radius / n

    nodes = Vector{Vector{Float64}}()
    elements = Vector{Vector{Int}}()
    densities = Float64[]
    elements_center = []

    node_map = Dict()
    new_node_id = 1

    step = Int(round_down_to_even(n)/2)

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

    return nodes, elements, densities, elements_center
end

nodes, elements, densities, elements_center = TestGeometrySphere(100)

function Vector2Vectors(nodes::Vector{Any})
    xₙ = zeros(length(nodes))
    yₙ = zeros(length(nodes))
    zₙ = zeros(length(nodes))

    for i in 1:length(nodes)
        xₙ[i], yₙ[i], zₙ[i] = nodes[i]
    end
    return xₙ, yₙ, zₙ
end
xₙ, yₙ, zₙ = Vector2Vectors(elements_center)

nel = length(elements_center)
ps = rand(Point3f, nel)
for i in 1:nel
    ps[i] = elements_center[i]
end
scatter(ps, color = densities)



function draw_elements_makie(nodes, elements)
    fig = Figure()
    ax = Axis3(fig[1, 1], aspect = :data)

    for element in elements
        # Získání souřadnic pro každý uzel v elementu
        element_coords = [nodes[node_id] for node_id in element]

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

    fig
end


element_coords = [nodes[node_id] for node_id in element]
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
                 [p[3] for p in all_node_coords], color=:red, markersize=50)

    fig
end

# Volání funkce s vašimi výsledky (nodes, elements)
draw_elements_makie(nodes, elements)

vectors = [
    [7, 23, 28, 11, 8, 24, 29, 12],
    [8, 24, 29, 12, 9, 25, 30, 13],
]

# Flattening the array of arrays
flattened = vcat(vectors...)

# Sorting and removing duplicates
sorted = sort(unique(flattened))

# Function to replace values with their positions in 'sorted'
replace_with_position!(v, sorted) = [findfirst(isequal(x), sorted) for x in v]

# Applying the replacement for each vector in 'vectors'
replaced_vectors = [replace_with_position!(v, sorted) for v in vectors]
