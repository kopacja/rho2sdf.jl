## Overview
The `rho2sdf` function transforms element densities from topology optimization into a signed distance field on a regular grid, calculating minimum distances to material boundaries defined by density isocontours, enabling precise geometric extraction.

<p align="center">
  <img src="doc/Cantilever_beam-RAW.png" width="45%" alt="Raw topology optimization result" />
  <img src="doc/Cantilever_beam-smooth_SDF_Approx.png" width="45%" alt="Smoothed SDF approximation" />
</p>

## Function Input

```julia
rho2sdf(taskName, X, IEN, rho; options=Rho2sdfOptions())
```

### Parameters:
- `taskName::String`: Base name for output files
- `X::Vector{Vector{Float64}}`: Mesh node coordinates
- `IEN::Vector{Vector{Int64}}`: Element connectivity (indices of element nodes)
- `rho::Vector{Float64}`: Element densities
- `options::Rho2sdfOptions`: Configuration options (optional)

### Return Value:
- `Tuple`: (fine_sdf, fine_grid, sdf_grid, sdf_dists)
- `vti`: SDF visualization in Paraview

## Rho2sdfOptions

The `Rho2sdfOptions` struct allows for customization of the SDF generation process:

```julia
Rho2sdfOptions(;
    threshold_density=0.5,          # value for isocontour (0, 1)
    sdf_grid_setup=:manual,         # manual/automatic grid setup
    export_nodal_densities=true,    # export nodal field to Paraview
    export_raw_sdf=true,            # export non-smoothed SDF to Paraview
    rbf_interp=false,               # interpolate/approximate SDF values using RBFs
    rbf_grid=:normal                # normal/fine grid for RBFs interp/approx
)
```

### Options:

#### threshold_density::Union{Float64, Nothing}
- Threshold density value in range [0,1] used for isosurface generation
- If `nothing`, it will be automatically calculated from volume fraction

#### sdf_grid_setup::Symbol
- SDF grid step configuration. Valid values:
  - `:automatic`: Uses original mesh properties for grid step
  - `:manual` (default): Uses interactive grid setup allowing user configuration

#### export_nodal_densities::Bool
- Controls whether nodal densities are exported to VTU format
- Default: `false`

#### export_raw_sdf::Bool
- Controls whether raw SDF values are exported to VTI format
- Default: `false`

#### rbf_interp::Bool
- Determines whether RBF interpolation or approximation is used
- Default: `true` (interpolation)

#### rbf_grid::Symbol
- Controls the resolution of the RBF grid
- Valid values:
  - `:normal` (default): Standard resolution (smooth=1)
  - `:fine`: Higher resolution (smooth=2)

## Example Usage

```julia
# Basic usage with default options
result = rho2sdf("chapadlo", X, IEN, rho)

# Custom options
options = Rho2sdfOptions(
    threshold_density=0.5,
    sdf_grid_setup=:manual,
    export_nodal_densities=true,
    export_raw_sdf=true,
    rbf_interp=true,
    rbf_grid=:fine
)
result = rho2sdf("chapadlo", X, IEN, rho, options=options)
```

## TODO List

- [ ] Remove unnecessary exports
- [ ] Remove unnecessary libraries
- [ ] Use IPopt library for points where local coordinates were not found
- [ ] Extend implementation to include tetrahedral elements
