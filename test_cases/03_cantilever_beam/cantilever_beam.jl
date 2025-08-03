taskName = "cantilever_beam_vfrac_04"

# Read FEM mesh:
data = matread("test_cases/03_cantilever_beam/input/$(taskName).mat")
(X, IEN, rho) = MeshInformations(data)
IEN = [subvector .- 1 for subvector in IEN]  # Data correction

# Custom options
options = Rho2sdfOptions(
  export_input_data = true,         # Export input (raw) data to Paraview
  sdf_grid_setup = :automatic,      # Automatic/manual grid setup
  export_nodal_densities = true,    # Export nodal field to Paraview
  export_raw_sdf = true,            # Export not smoothed SDF to Paraview
  rbf_interp = true,                # Interp/approx SDF values using RBFs
  rbf_grid = :same                  # Same/fine grid for RBFs interp/approx
)

result = rho2sdf("03_cantilever_beam", X, IEN, rho, options = options)

print_success("RUN_03_CANTILEVER_BEAM - done")
