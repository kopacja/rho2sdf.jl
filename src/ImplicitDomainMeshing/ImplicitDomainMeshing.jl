module ImplicitDomainMeshing

export GenerateTetMesh, PlaneDefinition, Rectangle, Square, Circle, Ellipse, slice_mesh_with_plane!, assess_mesh_quality
 
using Statistics
using StaticArrays
using LinearAlgebra
using Random
using WriteVTK
using Printf

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
# Compute mesh volume
include("TetMeshVolume.jl")
# Slice mesh with place (for visualization mesh quality)
include("Utils/SliceMeshWithPlane.jl")
# Check mesh quality & plot histogram of dihedral angles
include("Utils/CheckMeshQuality.jl")

# Main function for tetrahedral discretization:
function GenerateTetMesh(fine_sdf::Array, fine_grid::Array, scheme::String, name::String, warp_param::Union{Float64, Nothing}=nothing, plane_definitions::Vector{PlaneDefinition}=PlaneDefinition[])
  mesh = BlockMesh(fine_sdf, fine_grid)

  # Choose scheme: "A15" or "Schlafli"
  generate_mesh!(mesh, scheme)

  # assess_mesh_quality(mesh, "initial_mesh")
  warp!(mesh) # Warp nearest nodes to the isocontour

  update_connectivity!(mesh) # Update mesh topology

  slice_ambiguous_tetrahedra!(mesh) # Remove elements outside the body

  update_connectivity!(mesh)

  adjust_nodes_to_isosurface!(mesh) # Simple cut of elements to follow the isocontour

  TetMesh_volumes(mesh)
  optimize_mesh!(mesh)

  export_mesh_vtk(mesh, "$(name)_TriMesh.vtu")

  # Apply cutting planes only if they are defined
  if !isempty(plane_definitions)
    warp_mesh_by_planes_sdf!(mesh, plane_definitions, warp_param)
    update_connectivity!(mesh)
    export_mesh_vtk(mesh, "$(name)_TriMesh_cut.vtu")
  end

  TetMesh_volumes(mesh)
  
  return mesh
end


end
