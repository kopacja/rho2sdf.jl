module PG_SphereTests

using Test
using LinearAlgebra
using Rho2sdf.PrimitiveGeometries

# Include or import the functions here

# Testing round_down_to_even
@testset "round_down_to_even tests" begin
    @test round_down_to_even(5) == 4
    @test round_down_to_even(2) == 2
    @test round_down_to_even(-1) == -2
    @test round_down_to_even(0) == 0
    # Add more tests for edge cases
end

# Testing TestGeometrySphere
@testset "TestGeometrySphere tests" begin
    # Use a small value for max_elements for testing purposes
    max_elements = 10
    nodes, elements, densities, centers = selectPrimitiveGeometry("sphere", 10)

    # Test if all nodes are within the sphere of radius 1
    @test all(node -> norm(node) <= 1.0, nodes)

    # Test densities and centers calculations
    # Add tests to verify the correctness of densities and centers

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
    @test size(elements, 2) == 8
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
    @test isempty(nodes)
    @test isempty(elements)
    @test isempty(densities)
    @test isempty(elements_center)

    # Test with max_elements = 1
    nodes, elements, densities, elements_center = selectPrimitiveGeometry("cube", 1)
    # Perform specific checks for this case
    @test length(nodes) == 8
    @test size(elements, 1) == 1
    @test size(elements, 2) == 8
    # Further checks can be added
end

# Error handling test
@testset "Error Handling for TestGeometryCube" begin
    # Test with a negative max_elements, or non-integer values if applicable
    # Depending on how you want your function to behave, you can check for errors or specific return values
end


end
