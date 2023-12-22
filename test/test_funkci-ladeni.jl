using LinearAlgebra, Statistics, Makie, GLMakie

round_down_to_even(number::Float64) = iseven(floor(Int, number)) ? floor(Int, number) : floor(Int, number) - 1

function generate_sphere_mesh_with_density(max_elements::Int64)
    radius, center, n = 1.0, [0., 0., 0.], cbrt(max_elements)
    delta, step = 2 * radius / n, Int(round_down_to_even(n)/2)

    nodes, elements, densities, elements_center = [], [], Float64[], []

    node_map = Dict()
    node_id = 1
    for i in -step:step, j in -step:step, k in -step:step
        x, y, z = i * delta, j * delta, k * delta
        if norm([x, y, z]) <= radius
            push!(nodes, [x, y, z])
            node_map[(i, j, k)] = node_id
            node_id += 1
        end
    end

    for i in -step:step-1, j in -step:step-1, k in -step:step-1
        corners = [(i + di, j + dj, k + dk) for di in 0:1, dj in 0:1, dk in 0:1]
        if all(c -> c in keys(node_map), corners)
            element_nodes = [node_map[c] for c in corners]
            push!(elements, element_nodes)

            element_center_coords = mean([nodes[nid] for nid in element_nodes], dims=1)
            element_center = vec(element_center_coords)

            push!(densities, 1.0 - norm(element_center - center) / radius)
            push!(elements_center, element_center)
        end
    end

    return nodes, elements, densities, elements_center
end

nodes, elements, densities, elements_center = generate_sphere_mesh_with_density(1000)

Vector2Vectors(nodes::Vector) = [getindex.(nodes, i) for i in 1:3]
xₙ, yₙ, zₙ = Vector2Vectors(elements_center)

nel = length(elements_center)
ps = [elements_center[i] for i in 1:nel]
scatter(ps, color = densities)
