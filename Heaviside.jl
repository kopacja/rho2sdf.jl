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
name = Is_interpolation ? "Interpolation" : "Approximation"

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

  println("Computing the LSF zero level to meet the volume condition...")
  @time LSF = rbf_interpolation_kdtree(coarse_grid, coarse_grid, weights, kernel)
  th = LS_Threshold(LSF, V_frac, 4)

  println("Computing $name on the fine grid...")
  @time fine_LSF = rbf_interpolation_kdtree(fine_grid, coarse_grid, weights, kernel)
  println("Done")

  # Shifting LSF to maintain volume
  fine_LSF_offset = fine_LSF .+ th

  B = round(my_grid.cell_size, digits=4)
  Rho2sdf.exportStructuredPointsToVTK(taskName * "_smooth_B-" * string(B) * "_" * name * ".vtk", my_grid, fine_LSF_offset, "distance", smooth)

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
k = 0.1

# Aplikace úpravy strmosti na všechny hodnoty
adjusted_sdf_dists = adjust_steepness.(sdf_dists, k)

fine_sdf = RBFs_smoothing(adjusted_sdf_dists, sdf_grid, false, 2, "chapadlo_sdf_$(k)") # interpolation == true, aproximation == false, smooth

# fine_sdf = RBFs_smoothing(sdf_dists, sdf_grid, false, 2, "chapadlo_new") # interpolation == true, aproximation == false, smooth
# fine_sdf = RBFs_smoothing(sdf_dists, sdf_grid, false, 2, "chapadlo_old") # interpolation == true, aproximation == false, smooth


# B = round(sdf_grid.cell_size, digits=4)
# Rho2sdf.exportStructuredPointsToVTK(taskName * "_smooth_B-" * string(B) * "_" * name * ".vtk", sdf_grid, fine_sdf, "distance", 2)
