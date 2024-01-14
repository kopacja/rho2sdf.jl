module GridTests


using Test
using Rho2sdf
using Rho2sdf.MeshGrid
using LinearAlgebra
using JLD


@testset "Tests on real data" begin
  # Test data -> leg from gripper:
  mesh = load("TestData/TestData_leg-mesh.jld", "mesh")
  grid = load("TestData/TestData_leg-grid.jld", "grid")

  # Is object inside BB?
  @test grid.AABB_min < vec(minimum(mesh.X, dims=2))
  @test grid.AABB_max > vec(maximum(mesh.X, dims=2))

  points = MeshGrid.generateGridPoints(grid)

  linkedList = MeshGrid.LinkedList(grid, points)

  # Linked list tests:
  @test length(linkedList.head) == grid.ngp
  @test length(linkedList.next) == grid.ngp

  @test maximum(linkedList.head) == grid.ngp
  # @test maximum(linkedList.next) == grid.ngp # why not?

  @test minimum(linkedList.head) == -1
  @test minimum(linkedList.next) == -1
end


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
