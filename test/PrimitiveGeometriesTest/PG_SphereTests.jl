module PG_SphereTests

using Test
using LinearAlgebra
using Statistics
using Rho2sdf.PrimitiveGeometries

# Include or import the functions here
function GeometryTest(elements::Any, nodes::Vector{Vector{Float64}})
    # Kontrola typu
      element = []
    if typeof(elements) == Vector{Vector{Int64}}
        element = elements[1]
    elseif typeof(elements) == Vector{Int64}
        element = elements
    else
        println("Wrong data type of elements: ", typeof(elements))
    end

    # Výpočet vektorů a normály
    edge1 = nodes[element[2]] - nodes[element[1]]
    edge2 = nodes[element[4]] - nodes[element[1]]
    normal = cross(edge1, edge2)

    # Vytvoření pole 'posuv' pomocí vektorizované operace
    posuv = [node .+ normal * 1e-6 for node in nodes[element[1:4]]]

    # Výpočet a porovnání průměrných hodnot
    mean_diff_bottom = mean(map(v -> norm(v), nodes[element[1:4]] .- nodes[element[5:8]]))
    mean_diff_posuv = mean(map(v -> norm(v), posuv .- nodes[element[5:8]]))
    return mean_diff_bottom, mean_diff_posuv
end

# Testing round_down_to_even
@testset "round_down_to_even tests" begin
    @test round_down_to_even(5) == 4
    @test round_down_to_even(2) == 2
    @test round_down_to_even(-1) == -2
    @test round_down_to_even(0) == 0
end

# Testing TestGeometrySphere
@testset "TestGeometrySphere tests" begin
    # Use a small value for max_elements for testing purposes
    max_elements = 10
    nodes, elements, densities, centers = selectPrimitiveGeometry("sphere", max_elements)

    # Test if all nodes are within the sphere of radius 1
    @test all(node -> norm(node) <= 1.0, nodes)

    # Test densities and centers calculations
    @test maximum(densities) <= 1.0
    @test minimum(densities) >= 0.0

    # Test the elements' formation and indexing
    # Add tests to verify the elements are formed and indexed correctly

    # Add more tests as needed to cover different aspects of the function
end

# Run these tests in a Julia session to verify the behavior of your functions.



# Basic functionality test
@testset "Functionality Tests for TestGeometryCube" begin
    max_elements = 5
    nodes, elements, densities, elements_center = selectPrimitiveGeometry("cube", max_elements)

    # Test array sizes
    @test length(nodes) == (max_elements + 1)^3
    @test size(elements, 1) == max_elements^3
    @test length(elements[1]) == 8
    @test length(densities) == max_elements^3
    @test size(elements_center, 1) == max_elements^3
    @test size(elements_center, 2) == 3

    # More detailed tests can be added here, such as checking specific node positions,
    # ensuring densities are within expected ranges, etc.
end

# Boundary cases test
@testset "Boundary Cases for TestGeometryCube" begin
    # Test with max_elements = 0
    nodes, elements, densities, elements_center = selectPrimitiveGeometry("cube", 0)
    @test isempty(elements)
    @test isempty(densities)
    @test isempty(elements_center)

    # Test with max_elements = 1
    nodes, elements, densities, elements_center = selectPrimitiveGeometry("cube", 1)
    # Perform specific checks for this case
    @test length(nodes) == 8
    @test length(elements[1]) == 8
    @test size(elements, 2) == 1
    # Further checks can be added
    #
    mean_diff_bottom, mean_diff_posuv = GeometryTest(elements, nodes)
    @test mean_diff_bottom > mean_diff_posuv


end

end
