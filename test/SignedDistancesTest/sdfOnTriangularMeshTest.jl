module sdfOnTriangularMeshTest


using Test
using Rho2sdf.SignedDistances

@testset "Basic Functionality of Barycentric Coords" begin
    # Assuming we have a known case where the barycentric coordinates are known
    x₁ = [5., 1., 4.]
    x₂ = [2., 3., 4.]
    x₃ = [3., 2., 4.]
    n = [0., 0., 1.]

    # center
    x_center = (x₁ + x₂ + x₃)./3
    @test barycentricCoordinates(x₁, x₂, x₃, n, x_center) ≈ [0.333, 0.333, 0.333] atol=0.001 

    # vertexes
    x_v1 = [5., 1., 4.]
    @test barycentricCoordinates(x₁, x₂, x₃, n, x_v1) ≈ [1., 0., 0.]

    x_v2 = [2., 3., 4.]
    @test barycentricCoordinates(x₁, x₂, x₃, n, x_v2) ≈ [0., 1., 0.]

    x_v3 = [3., 2., 4.]
    @test barycentricCoordinates(x₁, x₂, x₃, n, x_v3) ≈ [0., 0., 1.]
    
    x_out = [10., 10., 10.]
    @test minimum(barycentricCoordinates(x₁, x₂, x₃, n, x_out)) < 0. 
end

end
