"""
Automated test for block geometry SDF generation.
Tests the complete pipeline from geometry creation to SDF computation.
"""

using Test
using Rho2sdf
using Rho2sdf.ElementTypes
using Rho2sdf.PrimitiveGeometries
using Rho2sdf.ShapeFunctions
using Rho2sdf.MeshGrid
using Rho2sdf.SignedDistances
using Rho2sdf.TerminalUtils
using LinearAlgebra
using Statistics

@testset "Block SDF Generation Test" begin
    
    @testset "Basic Block SDF Pipeline" begin
        # Test parameters
        taskName = "block_test"
        N = 20  # Number of cells along the longest side
        ρₜ = 0.5 # Threshold density (isosurface level)
        
        # Expected results (with tolerance for numerical precision)
        expected_max_sdf = 0.4242640687119285
        expected_mean_sdf = -1.4699474563515213e9
        
        # Tolerance for floating point comparisons
        rel_tolerance = 1e-10  # Relative tolerance
        abs_tolerance_max = 1e-12  # Absolute tolerance for maximum
        abs_tolerance_mean = 1e5   # Absolute tolerance for mean (due to large magnitude)
        
        print_info("Starting Block SDF generation test with HEX elements...")
        
        # Generate test geometry
        (X, IEN, rho) = selectPrimitiveGeometry("block", [2, 1, 1])
        
        # Verify geometry creation
        @test !isempty(X)
        @test !isempty(IEN)
        @test !isempty(rho)
        @test length(X) > 0
        @test length(IEN) > 0
        
        # Create mesh with HEX8 elements
        shape_func = coords -> shape_functions(HEX8, coords)
        mesh = Mesh(X, IEN, rho, shape_func; element_type=HEX8)
        
        # Verify mesh creation
        @test mesh.element_type == HEX8
        @test get_num_nodes(mesh.element_type) == 8
        
        # Define nodal densities
        ρₙ = [0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 0.0, 0.0, 0.5, 0.5]
        
        # Verify nodal densities
        @test length(ρₙ) == length(X)
        @test all(-0.1 ≤ ρ ≤ 1.1 for ρ in ρₙ)
        @test minimum(ρₙ) ≈ 0.0
        @test maximum(ρₙ) ≈ 1.0
        
        # Setup SDF grid
        X_min, X_max = getMesh_AABB(mesh.X)
        sdf_grid = Grid(X_min, X_max, N, 3)
        
        # Verify grid setup
        @test sdf_grid.ngp > 0
        @test length(sdf_grid.AABB_min) == 3
        @test length(sdf_grid.AABB_max) == 3
        @test all(sdf_grid.AABB_min .< sdf_grid.AABB_max)
        
        # Generate grid points
        points = generateGridPoints(sdf_grid)
        
        # Verify grid points
        @test size(points, 1) == 3
        @test size(points, 2) == sdf_grid.ngp
        
        # Compute SDF from densities
        (dists, xp) = evalDistances(mesh, sdf_grid, points, ρₙ, ρₜ)
        
        # Verify distance computation
        @test length(dists) == sdf_grid.ngp
        @test all(d -> d ≥ 0, dists)
        @test size(xp, 1) == 3
        @test size(xp, 2) == sdf_grid.ngp
        
        signs = Sign_Detection(mesh, sdf_grid, points, ρₙ, ρₜ)
        
        # Verify sign computation
        @test length(signs) == sdf_grid.ngp
        @test all(s -> s ∈ [-1, 1], signs)
        
        # Compute final SDF
        sdf_dists = dists .* signs
        
        # Verify SDF array
        @test length(sdf_dists) == sdf_grid.ngp
        
        # Test the main results against expected values
        actual_max_sdf = maximum(sdf_dists)
        actual_mean_sdf = mean(sdf_dists)
        
        print_info("Verifying results...")
        print_info("Expected maximum SDF: $(expected_max_sdf)")
        print_info("Actual maximum SDF: $(actual_max_sdf)")
        print_info("Expected mean SDF: $(expected_mean_sdf)")
        print_info("Actual mean SDF: $(actual_mean_sdf)")
        
        # Test maximum SDF value
        @test isapprox(actual_max_sdf, expected_max_sdf, rtol=rel_tolerance, atol=abs_tolerance_max)
        
        # Test mean SDF value (using absolute tolerance due to large magnitude)
        @test isapprox(actual_mean_sdf, expected_mean_sdf, atol=abs_tolerance_mean)
        
        # Additional sanity checks
        @test actual_max_sdf > 0
        @test actual_mean_sdf < 0
        @test minimum(sdf_dists) < 0
        
        # Check SDF distribution
        positive_count = count(s -> s > 0, sdf_dists)
        negative_count = count(s -> s < 0, sdf_dists)
        zero_count = count(s -> s == 0, sdf_dists)
        
        @test positive_count > 0
        @test negative_count > 0
        @test positive_count + negative_count + zero_count == length(sdf_dists)
        
        print_success("✅ Block SDF generation test completed successfully!")
    end
end
