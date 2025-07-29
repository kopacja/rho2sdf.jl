module ElementTypes

export AbstractElement, HEX8, TET4
export get_element_topology, get_num_nodes, get_num_faces, get_nodes_per_face
export get_local_coord_bounds, validate_local_coords

# Abstract type hierarchy
abstract type AbstractElement end

# Concrete element types
struct HEX8 <: AbstractElement end
struct TET4 <: AbstractElement end

# Element topology information
function get_element_topology(::Type{HEX8})
    nen = 8  # nodes per element
    nes = 6  # number of element faces
    nsn = 4  # nodes per face
    
    # Face connectivity (local node indices)
    ISN = [
        [1, 4, 3, 2],  # Face 1: bottom (z = -1)
        [1, 2, 6, 5],  # Face 2: front (y = -1)
        [2, 3, 7, 6],  # Face 3: right (x = +1)
        [3, 4, 8, 7],  # Face 4: back (y = +1)
        [4, 1, 5, 8],  # Face 5: left (x = -1)
        [5, 6, 7, 8]   # Face 6: top (z = +1)
    ]
    
    # Edge connectivity (for visualization/analysis)
    edges = (
        (1, 2), (2, 3), (3, 4), (4, 1),  # bottom face edges
        (5, 6), (6, 7), (7, 8), (8, 5),  # top face edges
        (1, 5), (2, 6), (3, 7), (4, 8)   # vertical edges
    )
    
    # Edge-to-face mapping
    ISE = (
        (1, 2, 3, 4),     # Bottom face edges
        (1, 9, 5, 10),    # Front face edges
        (2, 11, 6, 10),   # Right face edges
        (3, 12, 7, 11),   # Back face edges
        (4, 12, 8, 9),    # Left face edges
        (5, 6, 7, 8)      # Top face edges
    )
    
    return nen, nes, nsn, ISN, edges, ISE
end

function get_element_topology(::Type{TET4})
    nen = 4  # nodes per element
    nes = 4  # number of element faces
    nsn = 3  # nodes per face
    
    # Face connectivity (local node indices)
    ISN = [
        [1, 3, 2],  # Face 1
        [1, 2, 4],  # Face 2
        [2, 3, 4],  # Face 3
        [1, 4, 3]   # Face 4
    ]
    
    # Edge connectivity
    edges = (
        (1, 2), (2, 3), (3, 1),  # base triangle edges
        (1, 4), (2, 4), (3, 4)   # edges to apex
    )
    
    # Edge-to-face mapping
    ISE = (
        (1, 2, 3),    # Face 1 edges
        (1, 4, 5),    # Face 2 edges
        (2, 5, 6),    # Face 3 edges
        (3, 6, 4)     # Face 4 edges
    )
    
    return nen, nes, nsn, ISN, edges, ISE
end

# Convenience functions
get_num_nodes(::Type{HEX8}) = 8
get_num_nodes(::Type{TET4}) = 4

get_num_faces(::Type{HEX8}) = 6
get_num_faces(::Type{TET4}) = 4

get_nodes_per_face(::Type{HEX8}) = 4
get_nodes_per_face(::Type{TET4}) = 3

# Local coordinate system bounds
function get_local_coord_bounds(::Type{HEX8})
    return (-1.0, 1.0)  # Parametric cube [-1,1]³
end

function get_local_coord_bounds(::Type{TET4})
    return (0.0, 1.0)   # Barycentric coordinates [0,1]
end

# Validation functions
function validate_local_coords(::Type{HEX8}, ξ::Vector{Float64})
    return all(-1.0 ≤ coord ≤ 1.0 for coord in ξ)
end

function validate_local_coords(::Type{TET4}, λ::Vector{Float64})
    return all(coord ≥ 0.0 for coord in λ) && sum(λ) ≤ 1.0
end

end
