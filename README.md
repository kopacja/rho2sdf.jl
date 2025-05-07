# Rho2SDF Documentation

## Overview
The `rho2sdf` function converts element densities to a Signed Distance Function (SDF) representation, which is useful for extracting geometry from raw topology optimization results (SIMP method).

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

## Rho2sdfOptions

The `Rho2sdfOptions` struct allows for customization of the SDF generation process:

```julia
Rho2sdfOptions(;
    threshold_density=nothing,
    sdf_grid_setup=:interactive,
    export_nodal_densities=false,
    export_raw_sdf=false,
    rbf_interp=true,
    rbf_grid=:normal
)
```

### Options:

#### threshold_density::Union{Float64, Nothing}
- Threshold density value in range [0,1] used for isosurface generation
- If `nothing` (default), it will be automatically calculated using `find_threshold_for_volume(mesh, ρₙ)`
- If value is outside [0,1] range, a warning is issued and automatic calculation is used
- Values 0.0 or 1.0 are considered extreme and trigger a warning

#### sdf_grid_setup::Symbol
- Determines how the SDF grid is configured
- Valid values:
  - `:interactive` (default): Uses interactive grid setup allowing user configuration
  - `:noninteractive`: Uses non-interactive setup with predefined parameters
- Invalid values trigger a warning and fall back to `:interactive`

#### export_nodal_densities::Bool
- Controls whether nodal densities are exported to VTU format
- Default: `false`
- When `true`, generates file: `{taskName}_nodal_densities.vtu`

#### export_raw_sdf::Bool
- Controls whether raw SDF values are exported to VTI format
- Default: `false`
- When `true`, generates file: `{taskName}_SDF_CellSize-{B}.vti`

#### rbf_interp::Bool
- Determines whether RBF interpolation or approximation is used
- Default: `true` (interpolation)
- `false` uses approximation

#### rbf_grid::Symbol
- Controls the resolution of the RBF grid
- Valid values:
  - `:normal` (default): Standard resolution (smooth=1)
  - `:fine`: Higher resolution (smooth=2)
- Invalid values trigger a warning and fall back to `:normal`

## Example Usage

```julia
# Basic usage with default options
result = rho2sdf("chapadlo", X, IEN, rho)

# Custom options
options = Rho2sdfOptions(
    threshold_density=0.5,
    sdf_grid_setup=:noninteractive,
    export_nodal_densities=true,
    export_raw_sdf=true,
    rbf_interp=true,
    rbf_grid=:fine
)
result = rho2sdf("chapadlo", X, IEN, rho, options=options)
```

## TODO List

- [ ] Promazat nepotřebné exporty
- [ ] Promazat nepotřebné knihovny
- [ ] Použití IPopt knihovny pro body u kterých nebyly nalezeny lokální souřadnice
- [ ] Implementovat tetrahedra elementy
