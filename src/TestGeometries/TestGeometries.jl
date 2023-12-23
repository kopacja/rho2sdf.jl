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
    # Inicializace proměnných: poloměr koule, střed a počet uzlů
    radius = 1.0
    center = [0., 0., 0.]
    n = cbrt(max_elements) 
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



end
