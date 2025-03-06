function new_tetrahedron_volume(vertices::Vector{Vector{Float64}})
    # Matrix pro výpočet objemu
    a = vertices[2] - vertices[1]
    b = vertices[3] - vertices[1]
    c = vertices[4] - vertices[1]
  
    # Objem = 1/6 * |det(a b c)|
    return dot(a, cross(b, c)) / 6.0
end

# Function to check volumes of all tetrahedra in the mesh
function check_tetrahedron_volumes(mesh::BlockMesh)
    large_elements = Tuple{Int,Float64}[]  # Store (element_id, volume) pairs
    positive_volumes = Tuple{Int,Float64}[]  # Store positive volume elements
    negative_volumes = Tuple{Int,Float64}[]  # Store negative volume elements
    
    for (elem_id, element) in enumerate(mesh.IEN)
        # Get vertices of the tetrahedron
        vertices = [mesh.X[node_id] for node_id in element]
        
        # Calculate volume
        volume = new_tetrahedron_volume(vertices)
        
        # Categorize by volume sign
        if volume > 0
            push!(positive_volumes, (elem_id, volume))
        elseif volume < 0
            push!(negative_volumes, (elem_id, volume))
        end
    end
    
    # Output summary
    @info "Volume Analysis Results:"
    @info "Total elements: $(length(mesh.IEN))"
    @info "Positive volume elements: $(length(positive_volumes))"
    @info "Negative volume elements: $(length(negative_volumes))"
    @info "Zero volume elements: $(length(mesh.IEN) - length(positive_volumes) - length(negative_volumes))"
    
    # Output warnings for large elements
    if !isempty(large_elements)
        @warn "Found $(length(large_elements)) elements with absolute volume larger than $(volume_threshold)"
    end
    
    # Output warnings for negative volumes
    if !isempty(negative_volumes)
        @warn "Found $(length(negative_volumes)) elements with negative volume"
        for (elem_id, volume) in negative_volumes
            @warn "Element $elem_id has negative volume: $volume"
        end
    end
    
    return (
        large_elements=large_elements, 
        positive_volumes=positive_volumes, 
        negative_volumes=negative_volumes
    )
end

# Usage example
# results = check_tetrahedron_volumes(mesh)
