function new_tetrahedron_volume(vertices)
    # Matrix pro výpočet objemu
    a = vertices[2] - vertices[1]
    b = vertices[3] - vertices[1]
    c = vertices[4] - vertices[1]
  
    # Objem = 1/6 * |det(a b c)|
    return dot(a, cross(b, c)) / 6.0
end


# Function to check volumes of all tetrahedra in the mesh
function TetMesh_volumes(mesh::BlockMesh; volume_threshold::Float64=1e6)
    large_elements = Tuple{Int,Float64}[]  # Store (element_id, volume) pairs
    positive_volumes = Tuple{Int,Float64}[]  # Store positive volume elements
    negative_volumes = Tuple{Int,Float64}[]  # Store negative volume elements
    zero_volumes = Tuple{Int,Float64}[]     # Store zero volume elements
    
    # Pre-allocate for better performance
    n_elements = length(mesh.IEN)
    volumes = Vector{Float64}(undef, n_elements)
    
    # Compute all volumes in a single pass
    for (elem_id, element) in enumerate(mesh.IEN)
        # Get vertices of the tetrahedron
        vertices = [mesh.X[node_id] for node_id in element]
        
        # Calculate volume
        volume = new_tetrahedron_volume(vertices)
        volumes[elem_id] = volume
        
        # Categorize by volume
        if volume > 0
            push!(positive_volumes, (elem_id, volume))
        elseif volume < 0
            push!(negative_volumes, (elem_id, volume))
            # Also check if this is a large element
            if abs(volume) > volume_threshold
                push!(large_elements, (elem_id, volume))
            end
        elseif volume == 0
            push!(zero_volumes, (elem_id, volume))
        end
        
        # Check for large positive elements too
        if volume > volume_threshold
            push!(large_elements, (elem_id, volume))
        end
    end
    
    # Output summary
    @info "Volume Analysis Results:"
    @info "Total elements: $n_elements"
    @info "Positive volume elements: $(length(positive_volumes))"
    @info "Negative volume elements: $(length(negative_volumes))"
    @info "Zero volume elements: $(length(zero_volumes))"
    
    # Output warnings for large elements
    if !isempty(large_elements)
        @warn "Found $(length(large_elements)) elements with absolute volume larger than $volume_threshold"
        # Optionally show the largest elements
        sorted_large = sort(large_elements, by=x->abs(x[2]), rev=true)
        for (i, (elem_id, volume)) in enumerate(sorted_large[1:min(5, length(sorted_large))])
            @warn "  Large element $elem_id has volume: $volume"
        end
        if length(sorted_large) > 5
            @warn "  ... and $(length(sorted_large) - 5) more"
        end
    end
    
    # Output warnings for negative volumes
    if !isempty(negative_volumes)
        @warn "Found $(length(negative_volumes)) elements with negative volume"
        # Show only first few negative elements to avoid flooding console
        for (i, (elem_id, volume)) in enumerate(negative_volumes[1:min(5, length(negative_volumes))])
            @warn "  Element $elem_id has negative volume: $volume"
        end
        if length(negative_volumes) > 5
            @warn "  ... and $(length(negative_volumes) - 5) more"
        end
    end
    
    # return (
    #     volumes=volumes,
    #     large_elements=large_elements, 
    #     positive_volumes=positive_volumes, 
    #     negative_volumes=negative_volumes,
    #     zero_volumes=zero_volumes
    # )
end

