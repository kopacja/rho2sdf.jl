module TestGeometries

export  round_down_to_even, TestGeometrySphere

using LinearAlgebra
using Statistics


function round_down_to_even(number::Float64)
    # Zaokrouhlení dolů na nejbližší celé číslo
    rounded = floor(Int, number)
    # Pokud je číslo sudé, vrátíme ho, jinak odčteme 1 pro dosažení sudého čísla
    return iseven(rounded) ? rounded : rounded - 1
end

function TestGeometrySphere(max_elements::Int64)
    radius = 1.0
    # Střed sféry, první uzel v počátku
    center = [0., 0., 0.]

    n = cbrt(max_elements) 
    delta = 2 * radius / n

    nodes = []
    elements = []
    densities = Float64[]
    elements_center = []

    node_id = 1
    node_map = Dict()
    step = Int(round_down_to_even(n)/2)

    for i in -step:step
        println(i)
        for j in -step:step
            for k in -step:step
                # Výpočet pozice uzlu relativně k počátku
                x = i * delta
                y = j * delta
                z = k * delta

                # Kontrola, zda je uzel uvnitř sféry
                if norm([x, y, z]) <= radius
                    push!(nodes, [x, y, z])
                    node_map[(i, j, k)] = node_id
                    node_id += 1
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

                    # Výpočet středu elementu
                    element_center_coords = [nodes[node_id] for node_id in element_nodes]
                    element_center = mean(hcat(element_center_coords...), dims=2)
                    element_center = vec(element_center)

                    # Výpočet hustoty na základě vzdálenosti od středu sféry
                    density = 1.0 - norm(element_center - center) / radius
                    println(element_center)
                    push!(densities, density)
                    push!(elements_center, element_center)
                end
            end
        end
    end

    # Změna datového typu na Vector{Vector{Float64}}
    nodes_new = Vector{Vector{Float64}}(undef, length(nodes))

    for (i, item) in enumerate(nodes)
        nodes_new[i] = [item]
    end

    return nodes_new, elements, densities, elements_center
end

end
