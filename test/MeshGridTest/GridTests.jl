module GridTests


using Test
using Rho2sdf
using Rho2sdf.MeshGrid
using LinearAlgebra
using JLD


@testset "Grid Tests" begin
  AABB_min = [0.0, 0.0, 0.0]
  AABB_max = [10.0, 10.0, 10.0]
  N_max = 10
  grid = Grid(AABB_min, AABB_max, N_max)

  @test grid.AABB_min < grid.AABB_max
  @test length(grid.N) == 3
  @test grid.cell_size > 0
  @test grid.ngp > 0
  # ... additional tests
end


@testset "getMesh_AABB Tests" begin
  X = [1. 2.; 3. 4.; 5. 6.]
  X_min, X_max = getMesh_AABB(X)

  @test X_min == [1.; 3.; 5.]
  @test X_max == [2.; 4.; 6.]
  # ... additional tests
end


@testset "generateGridPoints Tests" begin
  grid = Grid([0.0, 0.0, 0.0], [10.0, 10.0, 10.0], 10)
  points = generateGridPoints(grid)

  @test size(points, 2) == grid.ngp
  @test all(points .>= grid.AABB_min) && all(points .<= grid.AABB_max)
  # ... additional tests
end

end
