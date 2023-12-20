using LinearAlgebra, Statistics

function generate_sphere_mesh_with_density(max_elements::Int64)
    radius = 1.0
    center = [0.0, 0.0, 0.0]

    n = cbrt(max_elements)
    delta = 2 * radius / n

    nodes = []
    elements = []
    densities = []

    node_id = 1
    node_map = Dict()

    for i in 0:n
        for j in 0:n
            for k in 0:n
                x = center[1] - radius + i * delta
                y = center[2] - radius + j * delta
                z = center[3] - radius + k * delta

                if norm([x, y, z] - center) <= radius
                    push!(nodes, [x, y, z])
                    node_map[(i, j, k)] = node_id
                    node_id += 1
                end
            end
        end
    end

    for i in 0:(n-1)
        for j in 0:(n-1)
            for k in 0:(n-1)
                corners = [(i, j, k), (i+1, j, k), (i, j+1, k), (i+1, j+1, k),
                           (i, j, k+1), (i+1, j, k+1), (i, j+1, k+1), (i+1, j+1, k+1)]
                if all(c -> c in keys(node_map), corners)
                    element_nodes = [node_map[c] for c in corners]
                    push!(elements, element_nodes)

                    # Compute the center of the element
                    element_center_coords = [nodes[node_id] for node_id in element_nodes]
                    element_center = mean(hcat(element_center_coords...), dims=2)
                    element_center = vec(element_center)  # Převedení na vektor

                    # Compute the density based on the distance from the sphere center
                    density = 1.0 - norm(element_center - center) / radius
                    push!(densities, density)

                end
            end
        end
    end

    return nodes, elements, densities
end

# Generate the mesh and densities
nodes, elements, densities = generate_sphere_mesh_with_density(100)

max_elements = 100
n = cbrt(max_elements)
    delta = 2 * 1 / n

# Output the densities
println("Densities:")
for (i, density) in enumerate(densities)
    println("Element $i: $density")
end
