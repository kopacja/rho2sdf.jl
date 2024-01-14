module InputDataTests


using Test
using LinearAlgebra
using Statistics
using Rho2sdf.MeshGrid
using MAT

# Data from Matlab:
taskName = "chapadlo"  # Define the task name, used for data file naming
data = matread(taskName * ".mat")  # Read the Matlab data file with the specified task name
(X, IEN, rho) = MeshGrid.MeshInformations(data)  # Extract mesh information from the data

# Verification of the node numbering convention in the element:
function GeometryTest(elements::Any, nodes::Vector{Vector{Float64}})
    element = []
    # Determine the type of 'elements' and assign the appropriate value to 'element'
    if typeof(elements) == Vector{Vector{Int64}}
        element = elements[1]
    elseif typeof(elements) == Vector{Int64}
        element = elements
    else
        println("Wrong data type of elements: ", typeof(elements))
    end

    # Calculate two edge vectors of the element and the normal vector
    edge1 = nodes[element[2]] - nodes[element[1]]
    edge2 = nodes[element[4]] - nodes[element[1]]
    normal = cross(edge1, edge2)  # Cross product to find the normal

    # Create an array 'posuv' by adding a small displacement to each node
    posuv = [node .+ normal * 1e-6 for node in nodes[element[1:4]]]

    # Calculate and compare the mean differences of the node positions
    # 'mean_diff_bottom' is the mean norm of the difference between bottom and top nodes
    mean_diff_bottom = mean(map(v -> norm(v), nodes[element[1:4]] .- nodes[element[5:8]]))
    # 'mean_diff_posuv' is the mean norm of the difference between displaced and top nodes
    mean_diff_posuv = mean(map(v -> norm(v), posuv .- nodes[element[5:8]]))
    return mean_diff_bottom, mean_diff_posuv
end

mean_diff_bottom, mean_diff_posuv = GeometryTest(IEN, X)
# Perform a test to check if the mean difference of the bottom nodes is greater than the displaced ones
@test mean_diff_bottom > mean_diff_posuv


end
