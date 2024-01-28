# function PartOfModel()
#     return println(pwd())
# end

function LoadTxt(cesta::String)
    file = open("$cesta", "r")
    array = []
    while eof(file) != true
        a = readline(file)
        if array == []
            array = parse(Int, a)
        else
            array = vcat(array, parse(Int, a))
        end
    end
    return array
end


function PartOfModel(
    mesh::Mesh,
    ρ::Vector,
    name::String)
    
    IDe = LoadTxt(name)

    # Reduce length:
    IDer = IDe[165:length(IDe)]
    # println("IDe", length(IDe))
    nelr = length(IDer)

    IENᵣ = zeros(mesh.nen, nelr)
    ρₙ = zeros(nelr)

    for i = 1:nelr
        IENᵣ[:, i] = mesh.IEN[:, IDer[i]]
        ρₙ[i] = ρ[IDer[i]]
    end # kontrolováno
    uIENₙ = unique(vec(IENᵣ))
    Xₙ = zeros(mesh.nsd + 1, length(uIENₙ))

    IENₙ = zeros(Int, size(IENᵣ))
    for i = 1:length(uIENₙ)
        Xₙ[:, i] = vcat(uIENₙ[i], mesh.X[:, Int(uIENₙ[i])]) # ok

        RidN = findall(isequal(uIENₙ[i]), IENᵣ)
        IENₙ[RidN] .= i
    end
    Xᵣ = Xₙ[2:4,:]
    mesh.X = Xᵣ
    mesh.IEN = IENₙ
    mesh.INE = nodeToElementConnectivity(Xₙ[2:4,:], IENₙ)
    mesh.nsd = size(Xᵣ, 1)
    mesh.nnp = size(Xᵣ, 2)
    mesh.nen = size(IENₙ, 1)
    mesh.nel = size(IENₙ, 2)
    return [mesh, ρₙ]
end

# (Xₙ, IENₙ, ρₙ) = PartOfModel(IDes, X, IEN, ρ)
# (mesh, ρ) = PartOfModel(mesh, IDes, ρ)

function ModiffElementalDensities(mesh::Mesh, rho::Vector{Float64})
    # Vectorized operations for min and max
    mins = minimum(mesh.X, dims=2)
    maxs = maximum(mesh.X, dims=2)
    
    part_size = maxs - mins
    max_size, poz = findmax(part_size)

    El_coord = zeros(mesh.nnp)  # Ensure this is the correct size
    for el in 1:mesh.nel  # Make sure the loop iterates over the correct range
        El_coord[el] = mean(mesh.X[:,mesh.IEN[:,el]][poz,:])
        rho[el] *= abs(El_coord[el] - mins[poz]) / max_size
    end
       
    return rho
end


