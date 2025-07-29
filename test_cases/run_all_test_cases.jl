using Rho2sdf
using Rho2sdf.ElementTypes
using Rho2sdf.TerminalUtils
using Rho2sdf.MeshGrid
using Rho2sdf.ShapeFunctions
using Rho2sdf.SignedDistances
using Rho2sdf.DataExport



RUN_01_SDF_VALIDSTION   = true
RUN_02_GRID_STUDY       = true
RUN_03_CANTILEVER_BEAM  = true
RUN_04_ROBOT_GRIPPER    = true



RUN_CUBE_HEX      = true
RUN_CUBE_HEX_REF  = false
RUN_CUBE_TET      = false
RUN_CUBE_TET_REF  = false

## Grid resolution study:
# N = 3  # Number of cells along the longest side (coarse)
N = 10  # Number of cells along the longest side (optimal)
# N = 80  # Number of cells along the longest side (fine)
N_all = [3, 10, 80]
for N in N_all
  include("02_grid_resolution_study/SphereInCube-Meshes.jl")
end

# include("03_cantilever_beam/cantilever_beam.jl")
# include("04_robot_gripper/robot_gripper.jl")
