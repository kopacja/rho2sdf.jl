"""
Automated test for sphere geometry SDF generation.
Tests the complete pipeline from .mat file loading to SDF computation.
"""

using Test
using Rho2sdf
using Rho2sdf.ElementTypes
using Rho2sdf.ShapeFunctions
using Rho2sdf.MeshGrid
using Rho2sdf.SignedDistances
using Rho2sdf.TerminalUtils
using MAT
using LinearAlgebra
using Statistics

@testset "Sphere SDF Generation Test" begin
    
    @testset "Basic Sphere SDF Pipeline" begin
        # Test parameters
        taskName = "sphere"
        N = 10  # Number of cells along the longest side
        ρₜ = 0.5 # Threshold density (isosurface level)
        
        # Expected results (with tolerance for numerical precision)
        expected_max_rho = 1.0000000000000022
        expected_mean_rho = 0.29490556408887564
        expected_max_sdf = 0.8669785608800439
        expected_mean_sdf = -3.7370242217627172e9
        
        # Tolerance for floating point comparisons
        rel_tolerance = 1e-10  # Relative tolerance
        abs_tolerance_rho = 1e-12  # Absolute tolerance for densities
        abs_tolerance_sdf_max = 1e-12  # Absolute tolerance for maximum SDF
        abs_tolerance_sdf_mean = 1e5   # Absolute tolerance for mean SDF (due to large magnitude)
        
        print_info("Starting Sphere SDF generation test with HEX elements...")
        
        # Check if sphere.mat file exists
        mat_file = taskName * ".mat"
        if !isfile(mat_file)
            @warn "Test file $mat_file not found. Skipping sphere test."
            return
        end
        
        # Read FEM mesh from .mat file
        data = matread(mat_file)
        (X, IEN, rho) = MeshInformations(data)
        
        # Verify mesh data extraction
        @test !isempty(X)
        @test !isempty(IEN)
        @test !isempty(rho)
        @test length(X) > 0
        @test length(IEN) > 0
        @test length(rho) > 0
        
        # Create mesh with HEX8 elements
        shape_func = coords -> shape_functions(HEX8, coords)
        mesh = Mesh(X, IEN, rho, shape_func; element_type=HEX8)
        
        # Verify mesh creation
        @test mesh.element_type == HEX8
        @test get_num_nodes(mesh.element_type) == 8
        @test mesh.nnp == length(X)
        @test mesh.nel == length(IEN)
        
        # Map elemental densities to the nodes using LSQ
        ρₙ = DenseInNodes(mesh, rho)
        
        # Verify nodal densities
        @test length(ρₙ) == length(X)
        @test all(-0.1 ≤ ρ ≤ 1.1 for ρ in ρₙ)  # Allow small numerical tolerance
        
        # Test nodal density statistics
        actual_max_rho = maximum(ρₙ)
        actual_mean_rho = mean(ρₙ)
        
        print_info("Verifying nodal density statistics...")
        print_info("Expected maximum ρₙ: $(expected_max_rho)")
        print_info("Actual maximum ρₙ: $(actual_max_rho)")
        print_info("Expected mean ρₙ: $(expected_mean_rho)")
        print_info("Actual mean ρₙ: $(actual_mean_rho)")
        
        @test isapprox(actual_max_rho, expected_max_rho, rtol=rel_tolerance, atol=abs_tolerance_rho)
        @test isapprox(actual_mean_rho, expected_mean_rho, rtol=rel_tolerance, atol=abs_tolerance_rho)
        
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
        
        # Test the main SDF results against expected values
        actual_max_sdf = maximum(sdf_dists)
        actual_mean_sdf = mean(sdf_dists)
        
        print_info("Verifying SDF results...")
        print_info("Expected maximum SDF: $(expected_max_sdf)")
        print_info("Actual maximum SDF: $(actual_max_sdf)")
        print_info("Expected mean SDF: $(expected_mean_sdf)")
        print_info("Actual mean SDF: $(actual_mean_sdf)")
        
        # Test maximum SDF value
        @test isapprox(actual_max_sdf, expected_max_sdf, rtol=rel_tolerance, atol=abs_tolerance_sdf_max)
        
        # Test mean SDF value (using absolute tolerance due to large magnitude)
        @test isapprox(actual_mean_sdf, expected_mean_sdf, atol=abs_tolerance_sdf_mean)
        
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
        
        # Check sphere-specific properties
        @test positive_count < negative_count  # Sphere should have more exterior than interior points
        
        # Verify reasonable density range for sphere
        @test minimum(ρₙ) ≥ 0.0
        @test maximum(ρₙ) ≤ 1.1  # Allow small numerical overshoot
        
        # Verify that there's a gradient in densities (sphere should have varying densities)
        @test std(ρₙ) > 0.1  # Should have significant variation
        
        print_success("✅ Sphere SDF generation test completed successfully!")
    end
    
    @testset "Sphere SDF Edge Cases" begin
        print_info("Testing edge cases for sphere SDF generation...")
        
        # Test with different threshold values
        if isfile("sphere.mat")
            data = matread("sphere.mat")
            (X, IEN, rho) = MeshInformations(data)
            
            shape_func = coords -> shape_functions(HEX8, coords)
            mesh = Mesh(X, IEN, rho, shape_func; element_type=HEX8)
            ρₙ = DenseInNodes(mesh, rho)
            
            # Test with extreme threshold values
            for test_threshold in [0.1, 0.9]
                X_min, X_max = getMesh_AABB(mesh.X)
                sdf_grid = Grid(X_min, X_max, 5, 3)  # Smaller grid for speed
                points = generateGridPoints(sdf_grid)
                
                (dists, xp) = evalDistances(mesh, sdf_grid, points, ρₙ, test_threshold)
                signs = Sign_Detection(mesh, sdf_grid, points, ρₙ, test_threshold)
                sdf_dists = dists .* signs
                
                # Basic sanity checks
                @test length(sdf_dists) == sdf_grid.ngp
                @test !all(isnan, sdf_dists)
                @test !all(isinf, sdf_dists)
            end
        end
        
        print_success("✅ Sphere SDF edge cases test completed!")
    end
end
