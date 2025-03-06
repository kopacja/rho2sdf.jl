# Abstraktní typ pro různé tvary rovin
abstract type PlaneShape end

# Definice konkrétních tvarů
struct Rectangle <: PlaneShape
    width::Float64
    height::Float64
end

# Speciální případ obdélníku
Square(size::Float64) = Rectangle(size, size)

struct Circle <: PlaneShape
    radius::Float64
end

struct Ellipse <: PlaneShape
    a::Float64
    b::Float64
end

# Struktura pro definici roviny
struct PlaneDefinition
    normal::Vector{Float64}
    point::Vector{Float64}
    shape::PlaneShape
end

# Struktura pro reprezentaci omezené roviny
struct BoundedPlane
    normal::Vector{Float64}
    point::Vector{Float64}
    shape::PlaneShape
    u::Vector{Float64}
    v::Vector{Float64}
    
    # Konstruktor pro výpočet bázových vektorů
    function BoundedPlane(normal::Vector{Float64}, point::Vector{Float64}, shape::PlaneShape)
        # Normalizace normálového vektoru
        normal = normalize(normal)
        
        # Vytvoření ortogonální báze v rovině
        # Nejprve zvolíme libovolný vektor, který není rovnoběžný s normálou
        temp = abs(normal[1]) < 0.9 ? [1.0, 0.0, 0.0] : [0.0, 1.0, 0.0]
        
        # Vytvoříme první bázový vektor v rovině pomocí vektorového součinu
        u = normalize(cross(normal, temp))
        
        # Vytvoříme druhý bázový vektor kolmý k normále a prvnímu bázovému vektoru
        v = normalize(cross(normal, u))
        
        new(normal, point, shape, u, v)
    end
end
