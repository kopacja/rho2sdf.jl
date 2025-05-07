using Rho2sdf
using MAT

## Inputs:
taskName = "cantilever_beam_vfrac_04"

## Read FEM mesh:
data = matread("test/" * taskName * ".mat")
(X, IEN, rho) = MeshInformations(data)
IEN = [subvector .- 1 for subvector in IEN]  # This will work

# Custom options
options = Rho2sdfOptions(
    threshold_density=0.5,          # value for isocotour (0, 1)
    sdf_grid_setup=:manual,         # automatic/manual grid setup
    export_nodal_densities=true,    # export nodal field to Paraview
    export_raw_sdf=true,            # export not smoothed SDF to Paraview
    rbf_interp=false,               # Interpolation or approximaion SDF values using RBFs
    rbf_grid=:same                  # same o finer grid for RBFs interp/approx
)

result = rho2sdf("beam", X, IEN, rho, options=options)

