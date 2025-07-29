using Rho2sdf

## Inputs:
taskName = "example_geom-chapadlo_vfrac-0.3"

# Import VTU mesh
(X, IEN, rho) = import_vtu_mesh("$(taskName).vtu")

# Validate imported data
validate_vtu_mesh(X, IEN, rho)


# Custom options
options = Rho2sdfOptions(
    # threshold_density=0.5,          # value for isocotour (0, 1)
    export_input_data=true,
    sdf_grid_setup=:manual,         # automatic/manual grid setup
    export_nodal_densities=true,    # export nodal field to Paraview
    export_raw_sdf=true,            # export not smoothed SDF to Paraview
    rbf_interp=true,               # interp/approx SDF values using RBFs
    rbf_grid=:same,                  # same/fine grid for RBFs interp/approx
    remove_artifacts=true,
    artifact_min_component_ratio=0.01,  # Remove components < 1% of largest
    export_analysis=true  # Export before/after comparison
)

result = rho2sdf("chapadlo", X, IEN, rho, options=options)
