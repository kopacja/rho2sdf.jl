using LinearAlgebra, Statistics, Makie, GLMakie, BenchmarkTools

function round_down_to_even(number::Int)
    # Zaokrouhlení dolů na nejbližší celé číslo
    rounded = floor(Int, number)
    # Pokud je číslo sudé, vrátíme ho, jinak odčteme 1 pro dosažení sudého čísla
    return iseven(rounded) ? rounded : rounded - 1
end

function SimplifiedGeometryCube(max_elements::Int64)
    side_length = 2.0
    delta = side_length / max_elements
    half_side = side_length / 2

    nodes = []
    node_map = Dict()
    new_node_id = 1

    # Generate nodes
    for i in 0:max_elements
        for j in 0:max_elements
            for k in 0:max_elements
                x = -half_side + i * delta
                y = -half_side + j * delta
                z = -half_side + k * delta

                push!(nodes, [x, y, z])
                node_map[(i, j, k)] = new_node_id
                new_node_id += 1
            end
        end
    end

    elements = []
    elements_center = []

    # Generate elements
    for i in 0:max_elements-1
        for j in 0:max_elements-1
            for k in 0:max_elements-1
                corners = [(i, j, k), (i+1, j, k), (i+1, j+1, k), (i, j+1, k),
                           (i, j, k+1), (i+1, j, k+1), (i+1, j+1, k+1), (i, j+1, k+1)]
                element_nodes = [node_map[c] for c in corners]
                push!(elements, element_nodes)

                element_center_coords = [nodes[node_id] for node_id in element_nodes]
                element_center = mean(hcat(element_center_coords...), dims=2)
                push!(elements_center, vec(element_center))
            end
        end
    end

    return nodes, elements, elements_center
end


nodes, elements, elements_center = TestGeometryCube(10)

function OptimizedGeometryCube(max_elements::Int64)
    side_length = 2.0
    delta = side_length / max_elements
    half_side = side_length / 2

    # Preallocate arrays
    total_nodes = (max_elements + 1)^3
    nodes = Array{Float64}(undef, total_nodes, 3)
    node_map = Dict()

    # Calculate nodes
    new_node_id = 1
    for i in 0:max_elements
        for j in 0:max_elements
            for k in 0:max_elements
                x = -half_side + i * delta
                y = -half_side + j * delta
                z = -half_side + k * delta

                nodes[new_node_id, :] = [x, y, z]
                node_map[(i, j, k)] = new_node_id
                new_node_id += 1
            end
        end
    end

    # Calculate elements
    total_elements = max_elements^3
    elements = Array{Int}(undef, total_elements, 8)
    elements_center = Array{Float64}(undef, total_elements, 3)
    densities = zeros(Float64, total_elements)
    
    element_idx = 1
    for i in 0:max_elements-1
        for j in 0:max_elements-1
            for k in 0:max_elements-1
                corners = [(i, j, k), (i+1, j, k), (i+1, j+1, k), (i, j+1, k),
                           (i, j, k+1), (i+1, j, k+1), (i+1, j+1, k+1), (i, j+1, k+1)]

                element_nodes = [node_map[c] for c in corners]
                elements[element_idx, :] = element_nodes

                element_center_coords = [nodes[node_id, :] for node_id in element_nodes]
                elements_center[element_idx, :] = mean(hcat(element_center_coords...), dims=2)
                
                densities[element_idx] = sqrt(2)*side_length/2 - norm(elements_center)
                element_idx += 1
            end
        end
    end

    return nodes, elements, elements_center, densities
end

nodes, elements, elements_center = @time TestGeometryCube(20);
nodes1, elements1, elements_center1 = @time SimplifiedGeometryCube(100);
nodes2, elements2, elements_center2, densities = @time OptimizedGeometryCube(100);

using Base.Threads

function ParallelGeometryCube(max_elements::Int64)
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

    return nodes, elements, elements_center, densities
end

using Base.Threads
using LinearAlgebra

function TestGeometryCube(max_elements::Int64)
    side_length = 2.0
    delta = side_length / max_elements
    half_side = side_length / 2

    # Preallocate arrays
    total_nodes = (max_elements + 1)^3
    nodes = Vector{Vector{Float64}}(undef, total_nodes)

    node_map = Dict()

    # Calculate nodes in parallel
    @threads for i in 0:max_elements
        for j in 0:max_elements
            for k in 0:max_elements
                x = -half_side + i * delta
                y = -half_side + j * delta
                z = -half_side + k * delta

                new_node_id = i * (max_elements + 1)^2 + j * (max_elements + 1) + k + 1
                nodes[new_node_id] = [x, y, z]
                node_map[(i, j, k)] = new_node_id
            end
        end
    end

    # Calculate elements
    total_elements = max_elements^3
    elements = Vector{Vector{Int}}(undef, total_elements)
    elements_center = Vector{Vector{Float64}}(undef, total_elements)
    densities = zeros(Float64, total_elements)

    @threads for i in 0:max_elements-1
        for j in 0:max_elements-1
            for k in 0:max_elements-1
                element_idx = i * max_elements^2 + j * max_elements + k + 1
                corners = [(i, j, k), (i+1, j, k), (i+1, j+1, k), (i, j+1, k),
                           (i, j, k+1), (i+1, j, k+1), (i+1, j+1, k+1), (i, j+1, k+1)]

                element_nodes = [node_map[c] for c in corners]
                elements[element_idx] = element_nodes

                element_center_coords = [nodes[node_id] for node_id in element_nodes]
                element_center = mean(hcat(element_center_coords...), dims=2)
                elements_center[element_idx] = vec(element_center)

                densities[element_idx] = 1 - norm(element_center)/(sqrt(3)*side_length/2)
            end
        end
    end

    # nodes_new = Vector{Vector{Float64}}(undef, total_nodes)
    # for i in 1:total_nodes
    #     nodes_new[i] = nodes[i,:]
    # end
    # nodes_new = [nodes[i, :] for i in 1:total_nodes]
    # nodes_new = [view(nodes, i, :) for i in 1:total_nodes]
    nodes_new = nodes


    # elements_new = Vector{Vector{Float64}}(undef, total_elements)
    # for i in 1:total_elements
    #     elements_new[i] = elements[i,:]
    # end


    return nodes_new, elements, densities, elements_center
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
    nodes_new = [view(nodes, i, :) for i in 1:total_nodes]
    elements_new = [view(elements, i, :) for i in 1:total_elements]

    return nodes, elements, densities, elements_center
end

nodes, elements, densities, elements_center = TestGeometryCube(50)

