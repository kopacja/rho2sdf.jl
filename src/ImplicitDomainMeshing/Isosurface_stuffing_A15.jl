using WriteVTK
using StaticArrays
using Logging
using JLD2
using LinearAlgebra

# Global tolerance for node merging (can be adjusted according to element size)
const TOL = 1e-6

# ----------------------------
# Block mesh structure (extended with node_hash)
# ----------------------------
mutable struct BlockMesh
    nx::Int
    ny::Int
    nz::Int
    grid::Array{SVector{3,Float64},3}          # 3D array of node coordinates (basic grid) using static vectors
    SDF::Array{Float64,3}                      # SDF values
    X::Vector{SVector{3,Float64}}    # List of physical node coordinates (nodes used in mesh)
    IEN::Vector{Vector{Int64}}                 # Tetrahedral connectivity (elements)
    INE::Vector{Vector{Int64}}                 # Inverse connectivity: for each node, list of adjacent elements
    node_sdf::Vector{Float64}                  # SDF values at nodes
    node_map::Dict{Int64,Int64}                # Mapping from original grid node index -> new node index
    cell_center_map::Dict{Tuple{Int,Int,Int},Int64}  # Mapping for cell centers (Steiner points)
    node_hash::Dict{NTuple{3,Float64},Int64}    # Global dictionary for merging nodes

    function BlockMesh()
        @load "src/ImplicitDomainMeshing/data/Z_block_FineGrid.jld2" fine_grid
        @load "src/ImplicitDomainMeshing/data/Z_block_FineSDF.jld2" fine_sdf

        # Convert fine_grid into a 3D array of SVectors
        grid = Array{SVector{3,Float64},3}(undef, size(fine_grid))
        for i in eachindex(fine_grid)
            # Assume fine_grid[i] is an array of Float64; convert to SVector
            grid[i] = SVector{3,Float64}(fine_grid[i]...)
        end
        sdf = Float64.(fine_sdf)
        nx, ny, nz = size(grid)
        
        mesh = new(nx, ny, nz)
        mesh.grid = grid
        mesh.SDF = sdf
        mesh.node_map = Dict{Int64,Int64}()
        mesh.cell_center_map = Dict{Tuple{Int,Int,Int},Int64}()
        mesh.X = Vector{SVector{3,Float64}}()
        mesh.IEN = Vector{Vector{Int64}}()
        mesh.INE = Vector{Vector{Int64}}()
        mesh.node_sdf = Vector{Float64}()
        mesh.node_hash = Dict{NTuple{3,Float64},Int64}()
        return mesh
    end
end

# ----------------------------
# Helper function: Get SDF values at the 8 corners of a cell (unchanged)
# ----------------------------
function get_cell_sdf_values(mesh::BlockMesh, i::Int, j::Int, k::Int)
    1 <= i < mesh.nx || throw(BoundsError(mesh.SDF, i))
    1 <= j < mesh.ny || throw(BoundsError(mesh.SDF, j))
    1 <= k < mesh.nz || throw(BoundsError(mesh.SDF, k))
    return SVector{8}(
        mesh.SDF[i, j, k],       # front-bottom-left
        mesh.SDF[i+1, j, k],     # front-bottom-right
        mesh.SDF[i+1, j+1, k],   # front-top-right
        mesh.SDF[i, j+1, k],     # front-top-left
        mesh.SDF[i, j, k+1],     # back-bottom-left
        mesh.SDF[i+1, j, k+1],   # back-bottom-right
        mesh.SDF[i+1, j+1, k+1], # back-top-right
        mesh.SDF[i, j+1, k+1]    # back-top-left
    )
end

# ----------------------------
# New: Definition of Schlafli orthoscheme connectivity (unchanged)
# ----------------------------
const schlafli_tet_connectivity = [
  [1, 2, 3, 7],  # Path 1: x, y, z
  [1, 6, 2, 7],  # Path 2: x, z, y
  [1, 3, 4, 7],  # Path 3: y, x, z
  [1, 4, 8, 7],  # Path 4: y, z, x
  [1, 5, 6, 7],  # Path 5: z, x, y
  [1, 8, 5, 7]   # Path 6: z, y, x
]

# ----------------------------
# Include original A15 scheme (contains definitions tile_ref and tetra_connectivity for A15)
# ----------------------------
include("A15_scheme.jl")  # This file provides tile_ref and tetra_connectivity for A15 scheme

# ----------------------------
# Quantization function – unchanged
# ----------------------------
function quantize(p::SVector{3,Float64}, tol::Float64=TOL)
  return (round(p[1] / tol)*tol, round(p[2] / tol)*tol, round(p[3] / tol)*tol)
end

# ----------------------------
# Helper function to compute shape functions for a hex8 element
# Accepts local coordinates in [0,1] (normalized tile_ref) and returns shape functions
# computed with standard transformation to [-1,1]
# ----------------------------
function shape_functions(ξηζ::SVector{3,Float64})::SVector{8,Float64}
  # Transform local coordinates from [0,1] to [-1,1]
  ξ = 2 * ξηζ[1] - 1.0
  η = 2 * ξηζ[2] - 1.0
  ζ = 2 * ξηζ[3] - 1.0
  coef = 1/8.0
  return @SVector [
    coef * (1 - ξ) * (1 - η) * (1 - ζ),
    coef * (1 + ξ) * (1 - η) * (1 - ζ),
    coef * (1 + ξ) * (1 + η) * (1 - ζ),
    coef * (1 - ξ) * (1 + η) * (1 - ζ),
    coef * (1 - ξ) * (1 - η) * (1 + ζ),
    coef * (1 + ξ) * (1 - η) * (1 + ζ),
    coef * (1 + ξ) * (1 + η) * (1 + ζ),
    coef * (1 - ξ) * (1 + η) * (1 + ζ)
  ]
end

# ----------------------------
# Function for discretizing a cell using A15 scheme (unchanged logic, but with refactored shape function and slight optimizations)
# ----------------------------
function process_cell_A15!(mesh::BlockMesh, i::Int, j::Int, k::Int)
  sdf_values = get_cell_sdf_values(mesh, i, j, k)
  if !any(x -> x >= 0, sdf_values)
    return
  end

  # Retrieve min and max corners of the cell
  v000 = mesh.grid[i,   j,   k]
  v111 = mesh.grid[i+1, j+1, k+1]
  vmins = v000
  vmaxs = v111

  # Precompute differences for coordinate interpolation
  Δ = vmaxs .- vmins

  # Get SDF values at cell corners
  cell_sdf = sdf_values
  f000, f100, f110, f010 = cell_sdf[1], cell_sdf[2], cell_sdf[3], cell_sdf[4]
  f001, f101, f111, f011 = cell_sdf[5], cell_sdf[6], cell_sdf[7], cell_sdf[8]
  
  local_mapping = Dict{Int, Int}()
  
  @inbounds for li in 1:length(tile_ref)
    # tile_ref is assumed to be defined in A15_scheme.jl and normalized (in [0,1]*4 originally)
    local_coord = SVector{3,Float64}(tile_ref[li] ./ 4.0)  # normalized coordinates in [0,1]
    # Compute physical point using linear interpolation
    p = vmins .+ Δ .* local_coord
    p = SVector{3,Float64}(p)  # ensure static vector type
    p_key = quantize(p)
    if haskey(mesh.node_hash, p_key)
      local_mapping[li] = mesh.node_hash[p_key]
    else
      push!(mesh.X, p)
      # Compute standard shape functions using helper function
      N = shape_functions(local_coord)
      # Interpolate SDF value using the 8 corner SDFs
      sdf_interp = N[1]*f000 + N[2]*f100 + N[3]*f110 + N[4]*f010 +
                   N[5]*f001 + N[6]*f101 + N[7]*f111 + N[8]*f011
      push!(mesh.node_sdf, sdf_interp)
      local_index = length(mesh.X)
      local_mapping[li] = local_index
      mesh.node_hash[p_key] = local_index
    end
  end
  
  # Process tetrahedral connectivity from A15 scheme
  @inbounds for tet in tetra_connectivity
    global_tet = [ local_mapping[li] for li in tet ]
    tet_sdf = [ mesh.node_sdf[idx] for idx in global_tet ]
    if any(x -> x >= 0, tet_sdf)
      push!(mesh.IEN, global_tet)
    end
  end
end

# ----------------------------
# Function for discretizing a cell using Schlafli orthoscheme (unchanged logic, only minor type annotation changes)
# ----------------------------
function process_cell_Schlafli!(mesh::BlockMesh, i::Int, j::Int, k::Int)
  # Get SDF values at the 8 corners of the cell
  sdf_values = get_cell_sdf_values(mesh, i, j, k)
  if !any(x -> x >= 0, sdf_values)
    return
  end
  
  local_mapping = Dict{Int,Int}()
  # Define cell nodes as SVectors from grid
  cell_nodes = [
    mesh.grid[i,   j,   k],     # Node 1: front-bottom-left
    mesh.grid[i+1, j,   k],       # Node 2: front-bottom-right
    mesh.grid[i+1, j+1, k],       # Node 3: front-top-right
    mesh.grid[i,   j+1, k],       # Node 4: front-top-left
    mesh.grid[i,   j,   k+1],     # Node 5: back-bottom-left
    mesh.grid[i+1, j,   k+1],     # Node 6: back-bottom-right
    mesh.grid[i+1, j+1, k+1],     # Node 7: back-top-right
    mesh.grid[i,   j+1, k+1]      # Node 8: back-top-left
  ]
  
  @inbounds for li in 1:8
    p = cell_nodes[li]
    p_key = quantize(p)
    if haskey(mesh.node_hash, p_key)
      local_mapping[li] = mesh.node_hash[p_key]
    else
      push!(mesh.X, p)
      push!(mesh.node_sdf, sdf_values[li])
      local_index = length(mesh.X)
      local_mapping[li] = local_index
      mesh.node_hash[p_key] = local_index
    end
  end
  
  # Construct tetrahedra according to Schlafli scheme
  @inbounds for tet in schlafli_tet_connectivity
    global_tet = [ local_mapping[li] for li in tet ]
    tet_sdf = [ mesh.node_sdf[idx] for idx in global_tet ]
    if any(x -> x >= 0, tet_sdf)
      push!(mesh.IEN, global_tet)
    end
  end
end

# ----------------------------
# Merge duplicate nodes after mesh generation (unchanged logic, only type annotations updated)
# ----------------------------
function merge_duplicate_nodes!(mesh::BlockMesh, tol::Float64=TOL)
  new_nodes = Vector{SVector{3,Float64}}()
  new_node_sdf = Vector{Float64}()
  node_map = Dict{Int, Int}()
  global_hash = Dict{NTuple{3,Float64}, Int}()
  @inbounds for i in 1:length(mesh.X)
    p = mesh.X[i]
    p_key = quantize(p, tol)
    if haskey(global_hash, p_key)
      node_map[i] = global_hash[p_key]
    else
      push!(new_nodes, p)
      push!(new_node_sdf, mesh.node_sdf[i])
      new_index = length(new_nodes)
      node_map[i] = new_index
      global_hash[p_key] = new_index
    end
  end
  @inbounds for i in 1:length(mesh.IEN)
    mesh.IEN[i] = [ node_map[old] for old in mesh.IEN[i] ]
  end
  mesh.X = new_nodes
  mesh.node_sdf = new_node_sdf
  mesh.node_map = node_map
  @info "After merging duplicates: $(length(mesh.X)) nodes"
end

# ----------------------------
# Cleanup unused nodes and reindex connectivity (unchanged logic)
# ----------------------------
function cleanup_unused_nodes!(mesh::BlockMesh)
  @info "Number of nodes before cleanup: $(length(mesh.X))"
  used_nodes = Set{Int64}()
  @inbounds for element in mesh.IEN
    union!(used_nodes, element)
  end
  new_node_map = Dict{Int64,Int64}()
  new_coords = Vector{SVector{3,Float64}}()
  new_node_sdf = Vector{Float64}()
  sorted_used = sort(collect(used_nodes))
  for (new_id, old_id) in enumerate(sorted_used)
    new_node_map[old_id] = new_id
    push!(new_coords, mesh.X[old_id])
    push!(new_node_sdf, mesh.node_sdf[old_id])
  end
  @inbounds for i in 1:length(mesh.IEN)
    mesh.IEN[i] = [ new_node_map[old_id] for old_id in mesh.IEN[i] ]
  end
  mesh.X = new_coords
  mesh.node_sdf = new_node_sdf
  mesh.node_map = new_node_map
  @info "Number of nodes after cleanup: $(length(mesh.X))"
end

# ----------------------------
# Create inverse connectivity (unchanged logic)
# ----------------------------
function create_INE!(mesh::BlockMesh)
  mesh.INE = [Vector{Int64}() for _ in 1:length(mesh.X)]
  @inbounds for (elem_id, element) in enumerate(mesh.IEN)
    for node_id in element
      push!(mesh.INE[node_id], elem_id)
    end
  end
  return mesh
end

# ----------------------------
# Modified mesh generation function with scheme selection
# ----------------------------
function generate_mesh!(mesh::BlockMesh, scheme::String)
  empty!(mesh.X)
  empty!(mesh.IEN)
  empty!(mesh.node_sdf)
  empty!(mesh.node_map)
  empty!(mesh.cell_center_map)
  empty!(mesh.node_hash)
  
  @inbounds for i in 1:mesh.nx-1
    for j in 1:mesh.ny-1
      for k in 1:mesh.nz-1
        if scheme == "A15"
          process_cell_A15!(mesh, i, j, k)
        elseif scheme == "Schlafli"
          process_cell_Schlafli!(mesh, i, j, k)
        else
          error("Unknown scheme: $scheme")
        end
      end
    end
  end
  
  cleanup_unused_nodes!(mesh)
  merge_duplicate_nodes!(mesh, TOL)
  
  create_INE!(mesh)
  @info "Mesh created: $(length(mesh.X)) nodes and $(length(mesh.IEN)) tetrahedra"
end

# ----------------------------
# Export mesh to VTK (unchanged logic, minor type adjustments)
# ----------------------------
function export_mesh_vtk(mesh::BlockMesh, filename::String)
  npoints = length(mesh.X)
  points = zeros(Float64, 3, npoints)
  @inbounds for i in 1:npoints
    points[:, i] = mesh.X[i]
  end
  cells = [MeshCell(VTKCellTypes.VTK_TETRA, tet) for tet in mesh.IEN]
  vtkfile = vtk_grid(filename, points, cells)
  vtk_point_data(vtkfile, mesh.node_sdf, "sdf")
  vtk_save(vtkfile)
end

# ----------------------------
# Main execution – create mesh and export it to file
# ----------------------------
mesh = BlockMesh()
# Choose scheme: "A15" or "Schlafli"
@time generate_mesh!(mesh, "A15")
export_mesh_vtk(mesh, "block-mesh.vtu")


mesh.INE
mesh.IEN
mesh.X
# mesh.X 282703-element
# mesh.IEN 475215-element