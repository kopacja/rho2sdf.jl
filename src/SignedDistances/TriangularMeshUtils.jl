function barycentricCoordinates(
  x₁::Vector{Float64}, # coordinates of vertex one of the triangle
  x₂::Vector{Float64}, # coordinates of vertex two of the triangle
  x₃::Vector{Float64}, # coordinates of vertex tree of the triangle
  n::Vector{Float64},  # unit normal to the face of the triangle
  x::Vector{Float64})  # one node of the grid

    A = [
        (x₁[2]*n[3]-x₁[3]*n[2]) (x₂[2]*n[3]-x₂[3]*n[2]) (x₃[2]*n[3]-x₃[3]*n[2])
        (x₁[3]*n[1]-x₁[1]*n[3]) (x₂[3]*n[1]-x₂[1]*n[3]) (x₃[3]*n[1]-x₃[1]*n[3])
        (x₁[1]*n[2]-x₁[2]*n[1]) (x₂[1]*n[2]-x₂[2]*n[1]) (x₃[1]*n[2]-x₃[2]*n[1])
    ]
    b = [
        x[2] * n[3] - x[3] * n[2],
        x[3] * n[1] - x[1] * n[3],
        x[1] * n[2] - x[2] * n[1],
    ]
    
    n_max, i_max = findmax(abs.(n)) ##???
    A[i_max, :] = [1.0 1.0 1.0]
    b[i_max] = 1.0

    return λ = A \ b # barycentric coordinates
end


function calculate_triangle_edges(
    Xt::Matrix{Float64}) # Coordinates of the vertices of the triangle

    Et = Vector{Vector{Float64}}(undef, 3) # Preallocate with undefined values
    Et[1] = Xt[:, 2] - Xt[:, 1]
    Et[2] = Xt[:, 3] - Xt[:, 2]
    Et[3] = Xt[:, 1] - Xt[:, 3]
    return Et
end


