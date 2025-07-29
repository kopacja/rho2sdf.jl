# Test Cases for Density-to-SDF Conversion Method

This directory contains four comprehensive test cases validating the proposed methodology for converting topology optimization results to smooth signed distance functions (SDFs). These tests demonstrate the method's effectiveness across different geometric complexities and practical engineering applications.

## Test Cases Overview

### 1. SDF Validation
**File:** `01_sdf_validation/roof_geometry.jl`

Validates the mathematical correctness of the SDF generation algorithm using a simplified two-element configuration. The test examines projection scenarios with a roof-like structure created using an isocontour value of ρₜ = 0.5, specifically designed to challenge conventional normal vector-based sign assignment methods.

### 2. Grid Resolution Study
**File:** `02_grid_resolution_study/SphereConvergenceTest.jl`

Analyzes the influence of Cartesian grid resolution on geometric feature extraction quality using spherical test geometry. Evaluates convergence characteristics across multiple discretization schemes:
- Uniform hexahedral mesh
- Non-uniform hexahedral mesh (refined bottom)
- Uniform tetrahedral mesh  
- Non-uniform tetrahedral mesh (refined bottom)

Tests three grid configurations: coarse (3.3× baseline spacing), reference (matched spacing), and fine (8× finer spacing).

### 3. Cantilever Beam
**File:** `03_cantilever_beam/cantilever_beam.jl`

Evaluates volume preservation capabilities and geometric quality using a classical topology optimization benchmark. The cantilever beam (L=20mm, t=4mm) with point loading demonstrates the method's performance on structural optimization results with 40% volume fraction constraint.

**Validation metrics:**
- Volume fraction preservation
- Structural efficiency via strain energy comparison
- Boundary smoothness through stress analysis
- Structural behavior consistency

### 4. Robot Gripper
**File:** `04_robot_gripper/robot_gripper.jl`

Demonstrates practical applicability through a complex engineering design case. The robot gripper features combined loading conditions (distributed pressure, surface traction, body forces), geometric constraints, and real material properties (ABS M30 plastic) with 30% volume fraction constraint.

## Installation

Clone the repository and install dependencies:

```bash
git clone https://github.com/your-username/Rho2sdf.jl.git
cd Rho2sdf.jl
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

## Running Test Cases

### Configuration

Edit test selection in `run_all_test_cases.jl` under the **Test case selection** section:

```julia
# Test case selection
RUN_01_SDF_VALIDATION   = true   # Enable/disable SDF validation test
RUN_02_GRID_STUDY       = true   # Enable/disable grid resolution study  
RUN_03_CANTILEVER_BEAM  = true   # Enable/disable cantilever beam test
RUN_04_ROBOT_GRIPPER    = true   # Enable/disable robot gripper test
```

### Execution

**Without parallelization:**
```bash
julia --project=. test_cases/run_all_test_cases.jl
```

**With parallelization (recommended):**
```bash
julia -t 10 --project=. test_cases/run_all_test_cases.jl
```

Replace `10` with the desired number of threads based on your system's capabilities.

## Expected Outputs

Each test case generates visualization files compatible with ParaView:
- `.vtu` files: Unstructured mesh data with density fields
- `.vti` files: Structured SDF data for analysis
- `.jld2` files: Julia-specific data for post-processing

**Typical runtime:** 5-15 minutes per test case depending on grid resolution and hardware configuration.

## Grid Resolution Study Details

The grid resolution study supports multiple mesh configurations. Enable specific mesh types in the configuration section:

```julia
# Grid study mesh types
RUN_CUBE_HEX      = true   # Uniform hexahedral mesh
RUN_CUBE_HEX_REF  = true   # Non-uniform hexahedral mesh  
RUN_CUBE_TET      = true   # Uniform tetrahedral mesh
RUN_CUBE_TET_REF  = true   # Non-uniform tetrahedral mesh

# Grid resolutions for study
GRID_RESOLUTIONS = [3, 10, 80]  # coarse, optimal, fine
```

## Requirements

- Julia ≥ 1.11 (compatibility with older versions has not been tested)
- All dependencies are automatically installed via `Pkg.instantiate()`
- ParaView (recommended for visualization)

## Citation
If you use this repository in your research, please cite us.

The manuscript has been submitted for publication. Citation details will be provided upon acceptance. In the meantime, you can cite this repository as follows:

Ježek, o. (2025). Rho2sdf.jl: A Julia Package for Converting Density-Based Topology Optimization Results to Smooth Signed Distance Functions. GitHub Repository. Available at: https://github.com/kopacja/rho2sdf.jl
