using Test
using Rho2sdf
using Rho2sdf.PrimitiveGeometries
using Rho2sdf.ShapeFunctions
using Rho2sdf.MeshGrid
using Rho2sdf.SignedDistances
using Rho2sdf.MyMarchingCubes
using Rho2sdf.DataExport
using MAT
using JLD2
using LinearAlgebra
using Statistics


## Inputs:
taskName = "chapadlo"
N = 60  # Number of cells along the longest side
ρₜ = 0.5 # Threshold density (isosurface level)

## Read FEM mesh:
data = matread("test/" * taskName * ".mat")
(X, IEN, rho) = MeshGrid.MeshInformations(data)
#Z,idx_Z = findall(x->X[3,i] > 50 for i in [1:size(X,2)])

## Generate FEM mesh structure:
mesh = MeshGrid.Mesh(X, IEN, C3D8_SFaD)

## Map elemental densities to the nodes:
ρₙ = MeshGrid.DenseInNodes(mesh, rho) # LSQ
#ρₙ = MeshGrid.elementToNodalValues(mesh, rho) # average

VTK_CODE = 12 # https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html
Rho2sdf.exportToVTU(taskName * "_nodal_densities.vtu", X, IEN, VTK_CODE, ρₙ)

## Grid:
X_min, X_max = MeshGrid.getMesh_AABB(mesh.X)
EN = NodePosition3D(mesh)


function analyze_mesh(distances::Matrix{Float64})
    num_elements = size(distances, 2)  # number of elements
    
    # Shortest and longest edge for each element
    min_edges = minimum(distances, dims=1)
    max_edges = maximum(distances, dims=1)
    
    # Overall shortest and longest edge
    shortest_edge = minimum(distances)
    longest_edge = maximum(distances)
    
    # Average edge length
    avg_edge_length = mean(distances)
    
    # Median edge length
    median_edge_length = median(vec(distances))
    
    # Sum of edge lengths for each element
    element_sums = sum(distances, dims=1)
    
    # Indices of the smallest and largest elements
    smallest_element_index = argmin(vec(element_sums))
    largest_element_index = argmax(vec(element_sums))
    
    # Average edge length of the smallest and largest elements
    avg_edge_smallest = mean(distances[:, smallest_element_index])
    avg_edge_largest = mean(distances[:, largest_element_index])
    
    # Print statistics
    println("Mesh Statistics:")
    println("----------------")
    println("Number of elements: ", num_elements)
    println("Shortest edge in the mesh: ", round(shortest_edge, digits=4))
    println("Longest edge in the mesh: ", round(longest_edge, digits=4))
    println("Average edge length: ", round(avg_edge_length, digits=4))
    println("Median edge length: ", round(median_edge_length, digits=4))
    println("Average edge length of the smallest element: ", round(avg_edge_smallest, digits=4))
    println("Average edge length of the largest element: ", round(avg_edge_largest, digits=4))
    println("\nRatio of longest to shortest edge: ", round(longest_edge / shortest_edge, digits=2))
    println("Ratio of average edge length of largest to smallest element: ", 
            round(avg_edge_largest / avg_edge_smallest, digits=2))
end

# Použití funkce (předpokládáme, že 'distances' je výstup z předchozí funkce)
# analyze_mesh(distances)

# distances = calculate_edge_distances(mesh)

# analyze_mesh(distances)
# maximum(X_max .- X_min)/N
# print("Enter the step size for the regular Cartesian grid: ")
# B = parse(Float64, readline())
# N_new = floor(Int, maximum(X_max .- X_min)/B) # Number of cells along the longest side
# N

# X_min[2] = 0
# X_min[3] = 50
sdf_grid = MeshGrid.Grid(X_min, X_max, N, 3) # cartesian grid
sdf_grid_new = MeshGrid.Grid(X_min, X_max, N_new, 3) # cartesian grid

function interactive_sdf_grid_setup(mesh::Mesh)
    X_min, X_max = MeshGrid.getMesh_AABB(mesh.X)
    distances = calculate_edge_distances(mesh)
    analyze_mesh(distances)

    println("\nThe time duration for 80k nodes was about 90 min. ")

    while true
        while true
            print("Write a grid step based on grid analysis: ")
            user_input = readline()
            
            try
                B = parse(Float64, user_input)
                
                N_new = floor(Int, maximum(X_max .- X_min)/B)
                sdf_grid = MeshGrid.Grid(X_min, X_max, N_new, 3)
                println("Number of all grid points: ", sdf_grid.ngp)
                
                break  # Exit the inner loop if the input is valid
            catch e
                if isa(e, ArgumentError)
                    println("Error: Please enter a valid floating-point number.")
                else
                    println("An unexpected error occurred. Please try again.")
                end
            end
        end

        print("Do you want to continue? (y/n): ")
        odpoved = lowercase(strip(readline()))

        if odpoved == "y"
            return sdf_grid
        else
            println("You can write a grid step.")
        end
    end
end

sdf_grid = interactive_sdf_grid_setup(mesh)


# ## SDF from densities:
# (sdf_dists, xp) = SignedDistances.evalSignedDistances(mesh, sdf_grid, ρₙ, ρₜ)

# ## Export to VT:
# Rho2sdf.exportStructuredPointsToVTK(taskName * "_sdf.vtk", sdf_grid, sdf_dists, "distance")

# @save "$(taskName)_cele_sdf.jld2" sdf_dists
# @save "$(taskName)_cele_sdf_grid.jld2" sdf_grid