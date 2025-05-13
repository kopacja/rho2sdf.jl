using Rho2sdf
using MAT

## Inputs:
taskName = "chapadlo"

## Read FEM mesh:
data = matread("test/" * taskName * ".mat")
(X, IEN, rho) = MeshInformations(data)

# Custom options
options = Rho2sdfOptions(
    # threshold_density=0.5,          # value for isocotour (0, 1)
    sdf_grid_setup=:automatic,         # automatic/manual grid setup
    export_nodal_densities=true,    # export nodal field to Paraview
    export_raw_sdf=true,            # export not smoothed SDF to Paraview
    rbf_interp=true,               # interp/approx SDF values using RBFs
    rbf_grid=:same                  # same/fine grid for RBFs interp/approx
)

result = rho2sdf("chapadlo", X, IEN, rho, options=options)
