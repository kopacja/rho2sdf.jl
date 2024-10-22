using JLD2
using KernelFunctions
using LinearAlgebra
using SparseArrays
using BenchmarkTools
using Statistics
using NearestNeighbors
using IterativeSolvers
using Printf
using HDF5
using Base.Threads: Atomic, atomic_add!

using Rho2sdf.ShapeFunctions
using Rho2sdf.MeshGrid
using Rho2sdf.DataExport
using Rho2sdf

include("sub_function_H.jl")

function RBFs_smoothing(
  dist::Vector,
  my_grid::Grid,
  Is_interpolation::Bool,
  smooth::Int,
  taskName::String,
  threshold::Float64=1e-3
)

  # SDF - Replace preallocated values with the maximum computed distance:
  dist_modif = process_vector(dist)

  V_frac = sum(dist .>= 0.0) / length(dist) # Volume fraction

  # Create coarse grid and assign values
  coarse_grid = create_grid(my_grid.N, my_grid)

  raw_SDF = vector_to_array(dist_modif, my_grid.N .+ 1) # data modification: vector -> array

  # Fine grid:
  dim = (my_grid.N * smooth) .+ 1 # dimensions of the fine grid

  (step, fine_grid) = create_smooth_grid(my_grid, smooth) # Create fine grid
  println("Original grid size: ", my_grid.N)
  println("Fine grid size: ", dim)
  # Define Gaussian kernel using LimitedRangeRBFKernel
  σ = my_grid.cell_size # width of the Gaussian function
  kernel = LimitedRangeRBFKernel(σ, threshold)

  # Compute weights on the coarse grid
  println("Computing weights on the coarse grid...")
  weights = Is_interpolation ?
            (@time compute_rbf_weights(coarse_grid, raw_SDF, kernel)) :
            vec(raw_SDF)

  println("Computing interpolation on the fine grid...")
  @time fine_LSF = rbf_interpolation_kdtree(fine_grid, coarse_grid, weights, kernel)
  println("Done")

  th = LS_Threshold(fine_LSF, V_frac, 4)
  fine_LSF_offset = fine_LSF .+ th

  name = Is_interpolation ? "Interpolation" : "Approximation"
  B = round(my_grid.cell_size, digits=4)
  Rho2sdf.exportStructuredPointsToVTK(taskName * "_smooth_B-" * string(B) * "_" * name * ".vtk", my_grid, fine_LSF_offset, "distance", smooth)

  # fine_LSF_offset_array = vector_to_array(fine_LSF_offset, dim)

  # # Save SDF data to HDF5 file
  # h5open(taskName * "_" * name * "_SDF_B-" * string(B) * ".h5", "w") do file
  #   write(file, "/SDF", fine_LSF_offset_array)
  #   write(file, "/dx", step)
  #   write(file, "/dy", step)
  #   write(file, "/dz", step)
  #   write(file, "/origin", Float32.(my_grid.AABB_min))
  # end
  return fine_LSF_offset
end


# Funkce pro úpravu strmosti pomocí Heaviside step function
function adjust_steepness(x, k)
  return 0.5 * (1 + tanh(k * x)) - 0.5
end
# Z_chapadlo_cele_Grid
# Načtení dat ze souboru
@load "Z_chapadlo_cele_SDf.jld2" sdf_dists
@load "Z_chapadlo_cele_Grid.jld2" sdf_grid

taskName = "Heaviside"
name = "Approximation"

# Parametr strmosti (můžete upravit podle potřeby)
k = 1.

# Aplikace úpravy strmosti na všechny hodnoty
adjusted_sdf_dists = adjust_steepness.(sdf_dists, k)

fine_sdf = RBFs_smoothing(adjusted_sdf_dists, sdf_grid, false, 2, "chapadlo_sdf_$(k)") # interpolation == true, aproximation == false, smooth

# B = round(sdf_grid.cell_size, digits=4)
# Rho2sdf.exportStructuredPointsToVTK(taskName * "_smooth_B-" * string(B) * "_" * name * ".vtk", sdf_grid, fine_sdf, "distance", 2)
