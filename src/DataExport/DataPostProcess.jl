
function SelectProjectedNodes(
    mesh::Mesh,
    grid::Grid,
    xp::Matrix{Float64},
    points::Matrix{Float64})
    ngp = grid.ngp # number of nodes in grid
    nsd = mesh.nsd # number of spacial dimensions

    # Assuming ngp is defined somewhere in your code
    # Preallocate arrays with maximum possible size
    max_size = ngp * 2  # Adjust this based on your knowledge of the data
    X = [zeros(Float64, nsd) for _ in 1:max_size]
    Xp = [zeros(Float64, nsd) for _ in 1:max_size]

    count = 0
    for i = 1:ngp
        if sum(abs.(xp[:, i])) > 1.0e-10
            count += 1
            X[count] = points[:, i]
            Xp[count] = xp[:, i]
        end
    end

    # Trim the unused preallocated space
    X = resize!(X, count)
    Xp = resize!(Xp, count)

    # Mean and max projected distance:
    if count > 0
        mean_PD = mean(norm.(X-Xp))
        max_PD = maximum(norm.(X-Xp))
    else 
        mean_PD = 0.
        max_PD = 0.
    end

    return X, Xp, mean_PD, max_PD
end

