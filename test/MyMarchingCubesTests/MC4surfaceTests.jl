module MC4surfaceTests

using Test
using Statistics
using LinearAlgebra
using GeometryBasics
using MarchingCubes
using Combinatorics
using StaticArrays
using Rho2sdf
# using Rho2sdf.MyMarchingCubes
import Rho2sdf.MyMarchingCubes as MMC


# Test cases
@testset "MC_OnCube Tests" begin
  # Test case 1: Simple cube with no intersection
  ρₑ = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
  ρₜ = 0.5
  triangles, vertices = MMC.MC_OnCube(ρₑ, ρₜ)
  @test length(triangles) == 0
  @test length(vertices) == 0

  # Test case 2: Simple cube with all vertices above threshold
  ρₑ = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
  ρₜ = 0.5
  triangles, vertices = MMC.MC_OnCube(ρₑ, ρₜ)
  @test length(triangles) == 0
  @test length(vertices) == 0

  # Test case 3: Simple cube with half vertices above threshold
  ρₑ = [1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0]
  ρₜ = 0.5
  triangles, vertices = MMC.MC_OnCube(ρₑ, ρₜ)
  @test length(triangles) > 0
  @test length(vertices) > 0
  
  # Test case 4:
  ρₑ = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
  ρₜ = 0.5
  triangles, vertices = MMC.MC_OnCube(ρₑ, ρₜ)
  # @test triangles == SVector(1, 3, 2)
  @test vertices == [[0.0, -1.0, -1.0], [-1.0, 0.0, -1.0], [-1.0, -1.0, 0.0]]
  # Print the result for manual inspection
  
  # Test case 5:
  ρₑ = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
  ρₜ = 0.5
  triangles, vertices = MMC.MC_OnCube(ρₑ, ρₜ)
  @test vertices == [[0.0, -1.0, -1.0], [1.0, 0.0, -1.0], [1.0, -1.0, 0.0]]
 
  # Test case 6:
  ρₑ = [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0]
  ρₜ = 0.5
  triangles, vertices = MMC.MC_OnCube(ρₑ, ρₜ)
  @test vertices == [[-1.0, -1.0, 0.0], [0.0, -1.0, 1.0], [-1.0, 0.0, 1.0]]

  # Test case 7:
  ρₑ = [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
  ρₜ = 0.5
  triangles, vertices = MMC.MC_OnCube(ρₑ, ρₜ)
  @test vertices == [[1.0, -1.0, 0.0], [0.0, -1.0, 1.0], [1.0, 0.0, 1.0]]

  # Test case 8:
  ρₑ = [0.75, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
  ρₜ = 0.5
  triangles, vertices = MMC.MC_OnCube(ρₑ, ρₜ)
  @test vertices == [[-0.33333333333333337, -1.0, -1.0], [-1.0, -0.33333333333333337, -1.0], [-1.0, -1.0, -0.33333333333333337]]

  # # Print the result for manual inspection
  println("Triangles: ", triangles)
  println("Vertices: ", vertices)
end
end
