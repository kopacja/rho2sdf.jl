# Function to check volumes of all tetrahedra in the mesh
function TetMesh_volumes(mesh::BlockMesh; volume_threshold::Union{Float64,Nothing}=nothing)
    # Pre-allocate arrays for better performance
    n_elements = length(mesh.IEN)
    volumes = Vector{Float64}(undef, n_elements)
    
    # Initialize data structures for categorizing elements
    large_elements = Tuple{Int,Float64}[]         # Store (element_id, volume) pairs for large elements
    negative_volumes = Tuple{Int,Float64}[]       # Store (element_id, volume) for negative volume elements
    zero_volumes = Tuple{Int,Float64}[]           # Store (element_id, volume) for zero volume elements
    inverted_elements = Int[]                     # Store element_id for inverted elements
    
    # First pass: compute all volumes
    total_absolute_volume = 0.0
    
    for (elem_id, element) in enumerate(mesh.IEN)
        # Get vertices of the tetrahedron
        vertices = [mesh.X[node_id] for node_id in element]
        
        # Calculate volume using cross product method
        a = vertices[2] - vertices[1]
        b = vertices[3] - vertices[1]
        c = vertices[4] - vertices[1]
        
        # Calculate determinant (6 times the volume)
        det_value = dot(a, cross(b, c))
        
        # Store actual volume (det/6)
        volume = det_value / 6.0
        volumes[elem_id] = volume
        
        # Update total absolute volume
        total_absolute_volume += abs(volume)
        
        # Check if element is inverted (negative Jacobian determinant)
        if det_value < 0
            push!(inverted_elements, elem_id)
        end
    end
    
    # Calculate average element volume and set threshold if not provided
    avg_volume = total_absolute_volume / n_elements
    if volume_threshold === nothing
        volume_threshold = 10.0 * avg_volume
    end
    
    # Second pass: categorize elements based on volume and threshold
    for (elem_id, volume) in enumerate(volumes)
        if volume < 0
            push!(negative_volumes, (elem_id, volume))
        elseif volume == 0
            push!(zero_volumes, (elem_id, volume))
        end
        
        # Check for large elements (both positive and negative)
        if abs(volume) > volume_threshold
            push!(large_elements, (elem_id, volume))
        end
    end
    
    # Sort large elements by absolute volume for reporting
    sorted_large = sort(large_elements, by=x->abs(x[2]), rev=true)
    
    # Prepare analysis summary for tabular output
    summary = Dict(
        "Total volume" => total_absolute_volume,
        "Negative volume elements" => length(negative_volumes),
        "Zero volume elements" => length(zero_volumes),
        "Inverted elements" => length(inverted_elements),
        "Large elements (|V| > threshold)" => length(large_elements)
    )
    
    # Create a formatted table for the summary using @sprintf
    println("\n====== MESH ANALYSIS RESULTS ======")
    println("+------------------------------------+--------------------+")
    println("| Metric                             | Value              |")
    println("+------------------------------------+--------------------+")
    
    for (key, value) in summary
        # Format the value based on its type
        formatted_value = if value isa AbstractFloat
            @sprintf("%.6e", value)
        else
            string(value)
        end
        
        # Print the formatted row with proper padding
        println(@sprintf("| %-30s | %-18s |", key, formatted_value))
    end
    println("+--------------------------------+--------------------+")
    
    # Output warnings for large elements
    if !isempty(large_elements)
        println("\n⚠️  LARGE ELEMENTS DETECTED ⚠️")
        println("+------------+------------------------+")
        println("| Element ID | Volume                 |")
        println("+------------+------------------------+")
        
        # Show top 5 largest elements (or fewer if there are less than 5)
        for (i, (elem_id, volume)) in enumerate(sorted_large[1:min(5, length(sorted_large))])
            println(@sprintf("| %-10d | %-22.6e |", elem_id, volume))
        end
        
        if length(sorted_large) > 5
            println("| ... and $(length(sorted_large) - 5) more large elements")
        end
        println("+------------+------------------------+")
    end
    
    # Output warnings for negative/inverted elements
    if !isempty(negative_volumes)
        println("\n⚠️  NEGATIVE VOLUME ELEMENTS DETECTED ⚠️")
        println("+------------+------------------------+")
        println("| Element ID | Volume                 |")
        println("+------------+------------------------+")
        
        # Sort by most negative first
        sorted_negative = sort(negative_volumes, by=x->x[2])
        
        # Show top 5 most negative elements (or fewer if there are less than 5)
        for (i, (elem_id, volume)) in enumerate(sorted_negative[1:min(5, length(sorted_negative))])
            println(@sprintf("| %-10d | %-22.6e |", elem_id, volume))
        end
        
        if length(negative_volumes) > 5
            println("| ... and $(length(negative_volumes) - 5) more negative elements")
        end
        println("+------------+------------------------+")
    end
    
    # Return comprehensive results as a named tuple for further processing
    # return (
    #     volumes = volumes,
    #     volume_threshold = volume_threshold,
    #     avg_volume = avg_volume,
    #     large_elements = large_elements, 
    #     negative_volumes = negative_volumes,
    #     zero_volumes = zero_volumes,
    #     inverted_elements = inverted_elements,
    #     summary = summary
    # )
end

# Example usage:
# results = TetMesh_volumes(my_mesh)  # Will automatically use 10x average volume as threshold
# results = TetMesh_volumes(my_mesh, volume_threshold=1.0)  # With a specific threshold
