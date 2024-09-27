
using Base.Threads
using Statistics
using LinearAlgebra
using ProgressMeter
using JLD2

using Rho2sdf.ShapeFunctions
using Rho2sdf.MeshGrid
using Rho2sdf.MyMarchingCubes
using Rho2sdf
using Rho2sdf.SignedDistances


# @load "Z_Chapadlo_xp.jld2" xp
@load "Z_Chapadlo_dist.jld2" dist
@load "Z_Chapadlo_mesh.jld2" mesh
@load "Z_Chapadlo_grid.jld2" grid
@load "Z_Chapadlo_points.jld2" points
@load "Z_Chapadlo_rho.jld2" ρₙ

# Function to create AABB from a set of points
function compute_aabb(points::Matrix{Float64})
  min_bounds = minimum(points, dims=2)
  max_bounds = maximum(points, dims=2)
  return min_bounds, max_bounds
end

# Function to check if a point is inside the AABB
function is_point_inside_aabb(x::Vector{Float64}, min_bounds, max_bounds)
  return all(min_bounds .<= x) && all(x .<= max_bounds)
end

X = mesh.X   # vector of nodes positions
IEN = mesh.IEN # ID element -> ID nodes
INE = mesh.INE # ID node -> ID elements
ISN = mesh.ISN # connectivity face - edges
sfce = mesh.sfce # shape function handler
nsd = mesh.nsd # number of spacial dimensions
nel = mesh.nel # number of all elements
nes = mesh.nes # number of element segments (faces)
nsn = mesh.nsn # number of face nodes
# println("number of all elements: ", nel)

ngp = grid.ngp # number of nodes in grid

ρₜ = 0.5
signs = -1 * ones(ngp)

# for i in 1:ngp # cycle trought all grid points
@showprogress 1 "Computing sign for grid points..." for i in 1:ngp # cycle through all grid points
  # println("grid point ID: ", i)
  found = false  # Flag to indicate if we need to skip to the next i
  x = points[:, i]
  # NotConv = false
  max_local = 10.0

  IDmin = argmin(mapslices(norm, (X .- x), dims=1))[2] # find ID of node (of mesh) that is closest to grid point x
  none = length(INE[IDmin]) # number of neighbour elements (that are connected to the node)
  # ρₙₑ = ρₙ[IEN[:, none]] # # nodal densities of none elements
  ρₙₑ = ρₙ[IEN[:, INE[IDmin]]] # # nodal densities of none elements
  if maximum(ρₙₑ) < ρₜ # empty element -> jump to another i
    continue
  end

  for j in 1:none
    el = INE[IDmin][j]

    Xₑ = X[:, IEN[:, el]]
    # Compute the AABB
    min_bounds, max_bounds = compute_aabb(Xₑ)

    # Check if the point x is inside the AABB
    inside = is_point_inside_aabb(x, min_bounds, max_bounds)

    if inside
      # if InOut(Xₑ, x)

      (NotConv, local_coords) = SignedDistances.find_local_coordinates(sfce, Xₑ, x)
      max_local_new = maximum(abs.(local_coords))
      # max_local = maximum(local_coords)
      # min_local = minimum(local_coords)
      #
      # if max_local < 1.0001 && min_local > -1.0001
      if max_local_new < 1.2 && max_local > max_local_new

        H, _, _, _ = sfce(local_coords) # tvarové funkce a jejich derivace
        ρₑ = ρₙ[IEN[:, el]] # nodal densities for one element
        ρ = H ⋅ ρₑ

        if ρ >= ρₜ
          # println("inside")
          signs[i] = 1.0
        end
        #   found = true
        #   break
        # else
        if !NotConv
          # println("point: ", i)
        end
        max_local = max_local_new
      end
    end
    # if !found && NotConv && mean(ρₙₑ) > ρₜ && j == none
    #   println("point: ", i)
    # end
    # if found
    #   break
    # end
  end

  # if !found
  #   println("node position not detected: ")
  #   println("id of grid point: ", i)
  # end

end

dist = abs.(dist) .* signs
