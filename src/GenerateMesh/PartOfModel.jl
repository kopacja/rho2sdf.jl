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
    
    IDes = LoadTxt(name)

    posun = mod(length(IDes),15) + 3
    IDes = IDes[(151+posun):length(IDes)]
    nelr = length(IDes)

    IENᵣ = zeros(mesh.nen, nelr)
    ρₙ = zeros(nelr)

    for i = 1:nelr
        IENᵣ[:, i] = mesh.IEN[:, IDes[i]]
        ρₙ[i] = ρ[IDes[i]]
    end # kontrolováno
    uIENₙ = unique(vec(IENᵣ))
    Xₙ = zeros(length(mesh.X[:, 1]) + 1, length(uIENₙ))

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
