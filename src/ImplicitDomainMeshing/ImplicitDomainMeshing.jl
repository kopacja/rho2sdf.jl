module ImplicitDomainMeshing

export GenerateTetMesh, PlaneDefinition, Rectangle, Square, Circle, Ellipse
 
using Statistics
using StaticArrays
using LinearAlgebra
using Random
using WriteVTK

include("CuttingPlaneTypes.jl")
include("A15_scheme.jl")
# Generate tetrahedral mesh and perform warp
include("GenerateMesh.jl")
# Cut elements to match SDF isocontour
include("Stencils.jl")
# Mesh optimization for better dihedral angles
include("OptimizeTetMesh.jl")
# Modify the resulting mesh to enable application of boundary conditions
include("ModifyResultingMesh.jl")
# include("TetraMeshVolume.jl")

# Main function for tetrahedral discretization:
function GenerateTetMesh(fine_sdf::Array, fine_grid::Array, scheme::String, name::String, plane_definitions::Vector{PlaneDefinition}=PlaneDefinition[])
  mesh = BlockMesh(fine_sdf, fine_grid)

  # Choose scheme: "A15" or "Schlafli"
  generate_mesh!(mesh, scheme)

  warp!(mesh) # Warp nearest nodes to the isocontour

  update_connectivity!(mesh) # Update mesh topology

  # export_mesh_vtk(mesh, "$(name)_warped.vtu")
  
  slice_ambiguous_tetrahedra!(mesh) # Remove elements outside the body

  update_connectivity!(mesh)

  adjust_nodes_to_isosurface!(mesh) # Simple cut of elements to follow the isocontour

  # update_connectivity!(mesh)
  
  optimize_mesh!(mesh)

  # check_tetrahedron_volumes(mesh)
  export_mesh_vtk(mesh, "$(name)_TriMesh.vtu")
  
  # Aplikace řezů rovinou pouze pokud jsou definovány
  if !isempty(plane_definitions)
    warp_mesh_by_planes_sdf!(mesh, plane_definitions)
    update_connectivity!(mesh)
    export_mesh_vtk(mesh, "$(name)_TriMesh_cut.vtu")
  end
  
  return mesh
end


end
