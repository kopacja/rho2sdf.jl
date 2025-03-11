# These planes are used for adjusting (aligning) the tetrahedral mesh
# to enable proper application of boundary conditions

# Abstract type for different plane shapes
abstract type PlaneShape end

# Definition of specific shapes
struct Rectangle <: PlaneShape
    width::Float64
    height::Float64
end

# Special case of rectangle
Square(size::Float64) = Rectangle(size, size)

struct Circle <: PlaneShape
    radius::Float64
end

struct Ellipse <: PlaneShape
    a::Float64
    b::Float64
end

# Structure for plane definition
struct PlaneDefinition
    normal::Vector{Float64}
    point::Vector{Float64}
    shape::PlaneShape
end

# Structure for representing a bounded plane
struct BoundedPlane
    normal::Vector{Float64}
    point::Vector{Float64}
    shape::PlaneShape
    u::Vector{Float64}
    v::Vector{Float64}
    
    # Constructor for calculating basis vectors
    function BoundedPlane(normal::Vector{Float64}, point::Vector{Float64}, shape::PlaneShape)
        # Normalization of the normal vector
        normal = normalize(normal)
        
        # Creating an orthogonal basis in the plane
        # First choose any vector that is not parallel to the normal
        temp = abs(normal[1]) < 0.9 ? [1.0, 0.0, 0.0] : [0.0, 1.0, 0.0]
        
        # Create the first basis vector in the plane using cross product
        u = normalize(cross(normal, temp))
        
        # Create the second basis vector perpendicular to the normal and the first basis vector
        v = normalize(cross(normal, u))
        
        new(normal, point, shape, u, v)
    end
end
