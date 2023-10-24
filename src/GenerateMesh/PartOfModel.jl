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

IDes = LoadTxt("data/Vtupni_data/elementy_trubky.txt")
posun = mod(length(IDes),15) + 3
IDes = IDes[(151+posun):length(IDes)]

function PartOfModel(IDes::Vector, X::Matrix, IEN::Matrix, ρ::Matrix)
    IENᵣ = zeros(length(IEN[:, 1]), length(IDes))
    ρₙ = zeros(length(IDes), 1)

    for i = 1:length(IDes)
        IENᵣ[:, i] = IEN[:, IDes[i]]
        ρₙ[i,1] = ρ[IDes[i]]
    end # kontrolováno
    uIENₙ = unique(vec(IENᵣ))
    Xₙ = zeros(length(X[:, 1]) + 1, length(uIENₙ))

    IENₙ = zeros(size(IENᵣ))
    for i = 1:length(uIENₙ)
        Xₙ[:, i] = vcat(uIENₙ[i], X[:, Int(uIENₙ[i])]) # ok

        RidN = findall(isequal(uIENₙ[i]), IENᵣ)
        IENₙ[RidN] .= i
    end
    # return [Xₙ[2:4,:], IENₙ, ρₙ]
    return [Xₙ[2:4,:], IENₙ, ρₙ]
end

# (Xₙ, IENₙ, ρₙ) = PartOfModel(IDes, X, IEN, ρ)
(X, IEN, ρ) = PartOfModel(IDes, X, IEN, ρ)
