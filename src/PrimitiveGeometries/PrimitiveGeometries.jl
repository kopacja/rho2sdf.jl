module PrimitiveGeometries

export  round_down_to_even, selectPrimitiveGeometry

using LinearAlgebra
using Statistics
using Base.Threads

function round_down_to_even(number::Int)
    # Zaokrouhlení dolů na nejbližší celé číslo
    rounded = floor(Int, number)
    # Pokud je číslo sudé, vrátíme ho, jinak odčteme 1 pro dosažení sudého čísla
    return iseven(rounded) ? rounded : rounded - 1
end

function selectPrimitiveGeometry(choice::String, max_elements::Int64)
    if choice == "sphere"
        return TestGeometrySphere(max_elements)
    elseif choice == "cube"
        return TestGeometryCube(max_elements)
    else
        error("Neplatný výběr")
    end
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
                corners = [(i, j, k),   (i+1, j, k),   (i+1, j+1, k),   (i, j+1, k), 
                           (i, j, k+1), (i+1, j, k+1), (i+1, j+1, k+1), (i, j+1, k+1)]
                # corners = [(i, j, k),   (i+1, j, k),   (i+1, j, k+1),   (i, j, k+1), 
                #            (i, j+1, k), (i+1, j+1, k), (i+1, j+1, k+1), (i, j+1, k+1)]
                
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



function TestGeometryCube(max_elements::Int64)
    side_length = 2.0
    delta = side_length / max_elements
    half_side = side_length / 2

    # Preallocate arrays
    total_nodes = (max_elements + 1)^3
    nodes = Array{Float64}(undef, total_nodes, 3)
    node_map = Dict()

    # Calculate nodes in parallel
    @threads for i in 0:max_elements
        for j in 0:max_elements
            for k in 0:max_elements
                x = -half_side + i * delta
                y = -half_side + j * delta
                z = -half_side + k * delta

                new_node_id = i * (max_elements + 1)^2 + j * (max_elements + 1) + k + 1
                nodes[new_node_id, :] = [x, y, z]
                node_map[(i, j, k)] = new_node_id
            end
        end
    end

    # Calculate elements
    total_elements = max_elements^3
    elements = Array{Int}(undef, total_elements, 8)
    elements_center = Array{Float64}(undef, total_elements, 3)
    densities = zeros(Float64, total_elements)

    @threads for i in 0:max_elements-1
        for j in 0:max_elements-1
            for k in 0:max_elements-1
                element_idx = i * max_elements^2 + j * max_elements + k + 1
                corners = [(i, j, k), (i+1, j, k), (i+1, j+1, k), (i, j+1, k),
                           (i, j, k+1), (i+1, j, k+1), (i+1, j+1, k+1), (i, j+1, k+1)]

                element_nodes = [node_map[c] for c in corners]
                elements[element_idx, :] = element_nodes

                element_center_coords = [nodes[node_id, :] for node_id in element_nodes]
                element_center = mean(hcat(element_center_coords...), dims=2)
                elements_center[element_idx, :] = element_center

                densities[element_idx] = 1 - norm(element_center)/(sqrt(3)*side_length/2)
            end
        end
    end

    nodes_new = [nodes[i, :] for i in 1:total_nodes]
    elements_new = [elements[i, :] for i in 1:total_elements]

    return nodes_new, elements_new, densities, elements_center

    # return nodes, elements, densities, elements_center
end

# (X, IEN, rho) = PrimitiveGeometries.selectPrimitiveGeometry("cube", 14)

end
