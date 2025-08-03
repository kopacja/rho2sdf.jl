"""
Test Cases for Density-to-SDF Conversion Method
===============================================

Four test cases validating the proposed methodology for converting 
topology optimization results to smooth signed distance functions:

1. SDF Validation - Roof geometry test verifying mathematical correctness
2. Grid Resolution Study - Sphere geometry analyzing grid discretization effects  
3. Cantilever Beam - Volume preservation and geometric quality assessment
4. Robot Gripper - Complex engineering design application
"""

using Rho2sdf
using Rho2sdf.ElementTypes
using Rho2sdf.TerminalUtils
using Rho2sdf.MeshGrid
using Rho2sdf.ShapeFunctions
using Rho2sdf.SignedDistances
using Rho2sdf.SdfSmoothing
using Rho2sdf.DataExport
using Rho2sdf.PrimitiveGeometries
using MAT

# =============================================================================
# CONFIGURATION
# =============================================================================

# Test case selection
RUN_01_SDF_VALIDATION   = false
RUN_02_GRID_STUDY       = true
RUN_03_CANTILEVER_BEAM  = false
RUN_04_ROBOT_GRIPPER    = false

# Grid study mesh types
RUN_CUBE_HEX      = true
RUN_CUBE_HEX_REF  = false
RUN_CUBE_TET      = false
RUN_CUBE_TET_REF  = false

# Grid resolutions for study
GRID_RESOLUTIONS = [3, 10, 80]  # coarse, optimal, fine

# =============================================================================
# EXECUTION
# =============================================================================

function main()
    # Test Case 1: SDF Validation
    if RUN_01_SDF_VALIDATION
        include("01_sdf_validation/roof_geometry.jl")
    end
    
    # Test Case 2: Grid Resolution Study
    if RUN_02_GRID_STUDY
        
        # Run for each grid resolution
        for N in GRID_RESOLUTIONS
            include("02_grid_resolution_study/SphereInCube-Meshes.jl")
        end
        
        # Convergence tests
        # include("02_grid_resolution_study/SphereConvergenceTest.jl")
    end
    
    # Test Case 3: Cantilever Beam
    if RUN_03_CANTILEVER_BEAM
        include("03_cantilever_beam/cantilever_beam.jl")
    end
    
    # Test Case 4: Robot Gripper
    if RUN_04_ROBOT_GRIPPER
        include("04_robot_gripper/robot_gripper.jl")
    end
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
