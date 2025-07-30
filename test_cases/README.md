# Test Cases for Density-to-SDF Conversion Method

This directory contains four test cases validating the proposed method for converting topology optimization results to smooth signed distance functions (SDFs).
## Test Cases Overview

### 1. SDF Validation
**File:** `01_sdf_validation/roof_geometry.jl`

Verifies the mathematical correctness of the SDF construction using a simplified two-element configuration with isocontour value ρₜ = 0.5. Validates distance calculation, sign assignment, and handling of all projection scenarios in geometrically complex regions.

### 2. Grid Resolution Study
**File:** `02_grid_resolution_study/SphereConvergenceTest.jl`

Establishes optimal grid resolution guidelines through two complementary analyses:

**Grid Configuration Test:** Evaluates geometry extraction quality across different Cartesian grid resolutions (coarse/reference/fine) and four mesh discretization schemes:
- Uniform hexahedral mesh
- Non-uniform hexahedral mesh (refined bottom)
- Uniform tetrahedral mesh  
- Non-uniform tetrahedral mesh (refined bottom)

**Convergence Test:** Validates numerical accuracy using analytical sphere geometry with progressive mesh refinement, confirming second-order spatial accuracy.

### 3. Cantilever Beam Validation
**File:** `03_cantilever_beam/cantilever_beam.jl`

Demonstrates volume preservation and geometric quality through classical topology optimization benchmark (L = 20mm, 40% volume fraction). Compares SDF post-processing against linear post-processing (nodal density mapping → Paraview Contour filter → Tetgen meshing). 

*Note: For complete structural analysis and performance validation, refer to the companion finite element analysis project.*

### 4. Robot Gripper Case Study
**File:** `04_robot_gripper/robot_gripper.jl`

Validates practical engineering applicability using complex robot gripper design with combined loading conditions, ABS M30 material properties, and 30% volume constraint. Demonstrates methodology performance on real-world manufacturing scenarios.

*Note: For complete structural analysis and performance validation, refer to the companion finite element analysis project.*

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

Ježek, O., et al. (2025). Rho2sdf.jl: A Julia Package for Converting Density-Based Topology Optimization Results to Smooth Signed Distance Functions. GitHub Repository. Available at: https://github.com/kopacja/rho2sdf.jl
