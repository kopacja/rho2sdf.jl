# using REPLVim; @async REPLVim.serve()
using Statistics: mean

    # face triangular mesh from 
    # mesh = Rho2sdf.extractSurfaceTriangularMesh(mesh, ρₙ) # přepsání X a IEN (pro trojúhelníky)

X = mesh.X
IEN = mesh.IEN # ID element -> nodes
INE = mesh.INE # ID node -> ID elements
ISN = mesh.ISN # connectivity face - edges
# nsd = mesh.nsd # number of special dimension
nel = mesh.nel # number of elements
nes = mesh.nes # number of element segments (faces) 6
nsn = mesh.nsn # number of segment nodes (kolik má stěna uzlů) 4


            
using Statistics
            mesh.INE
            mesh.IEN
            mesh.nnp
            mesh.nel

            EN_x = zeros(Float64, size(mesh.IEN))

            # Geometric centres of elements (GPs for first order elements)
mutable struct NodalCoordinatesInElement{T<:Array}
    x::T
    y::T
    z::T
end
function nodeToElementConnectivity(
    X::Vector{Vector{Float64}},
    IEN::Matrix{Int64},
)::Vector{Vector{Int64}}
    INE = [Vector{Int64}() for _ = 1:length(X)]
    for el = 1:size(IEN, 2)
        for i = 1:size(IEN, 1)
            push!(INE[IEN[i, el]], el)
        end
    end
    return INE
end

function NodePosition3D(mesh::Mesh)
    table = size(mesh.IEN)
    EN_x = zeros(Float64, table) #počet ele = počet Gausspointů
    EN_y = zeros(Float64, table)
    EN_z = zeros(Float64, table)

    for i = 1:(length(mesh.IEN[:, 1]))
        for j = 1:(length(mesh.IEN[1, :]))
            ID_n = Int(mesh.IEN[i, j])
            x = mesh.X[1, ID_n]
            y = mesh.X[2, ID_n]
            z = mesh.X[3, ID_n]
            EN_x[i, j] = x
            EN_y[i, j] = y
            EN_z[i, j] = z
        end
    end
    EN = NodalCoordinatesInElement(EN_x, EN_y, EN_z)

    return EN #Výstupem jsou tři pole - souřadnice uzlů daného elementu (pozice řádku = ID elemetu)
end

EN = NodePosition3D(mesh)
Centre = GeometricCentre(mesh, EN)

function GeometricCentre(mesh::Mesh, EN::NodalCoordinatesInElement)
    Centre = zeros(mesh.nel, length(mesh.X[:, 1]))
    for i = 1:mesh.nel
        G_x = mean(EN.x[:, i]) # G_x = dot(EN_x[:, i], N)
        G_y = mean(EN.y[:, i]) # G_y = dot(EN_y[:, i], N)
        G_z = mean(EN.z[:, i]) # G_z = dot(EN_z[:, i], N)
        Centre[i, :] = [G_x, G_y, G_z]
    end
    return Centre
end

#######################################################xx
# Filter for 2 a 3 nodes
function FilterForNodalDensity(
    i::Int,
    Centre::Matrix,
    mesh::Mesh,
    ρ::Vector,
)
    L = zeros(length(mesh.INE[i]))

    for j = 1:length(mesh.INE[i])
        L[j] = norm(mesh.X[:, i] - Centre[mesh.INE[i][j], :])
    end
    Lmax = maximum(L) * 1.2
    dm = 0.0
    deleni = 0.0
    for j = 1:length(mesh.INE[i])
        dm += ρ[Int((mesh.INE[i])[j])] * (1 - L[j] / Lmax)
        deleni += (1 - L[j] / Lmax)
    end
    return DN = dm / deleni
end

function NodalDensityLeastSquares(
    i::Int,
    Centre::Matrix,
    mesh::Mesh,
    ρ::Vector,
)
    A = zeros(length(mesh.INE[i]), length(Centre[1, :]) + 1)
    b = zeros(length(mesh.INE[i]))
    for j = 1:length(mesh.INE[i])
        A[j, :] = [1; Centre[Int((mesh.INE[i])[j]), :]]
        b[j] = ρ[Int((mesh.INE[i])[j])]
    end
    λ, ϕ = eigen(A'A)
    lam = LamReduction(λ)
    if length(lam) > 1
        Λ = Diagonal(lam)
    else
        Λ = lam
    end

    if Λ == []
        DN = 0.
        if mean(abs.(b)) != 0.0
            println("nenulový vektor b", b)
        end
    else
        poz = length(λ) - length(lam) + 1
        b1 = ϕ' * (A' * b)
        x1 = Λ \ b1[poz:length(b1)]
        if poz > 1
            zero = zeros(poz - 1)
            x2 = vcat(zero, x1)
        else
            x2 = x1
        end
        x = ϕ * x2
        DN = dot(vcat(1, mesh.X[:, i]), x)
        # DN = KeepRange(DN)
    end
    # DN = KeepRange(DN)
    return DN
end

function LamReduction(λ::Vector)
    lam = []
    εₘ₁ = 1.e7
    εₘ₂ = 3.e3
    ε₁ = abs(maximum(λ) / minimum(λ))
    ε₂ = abs(maximum(λ) / λ[2])
    ε₃ = abs(maximum(λ) / λ[3])
    if εₘ₁ > ε₁ && εₘ₂ > ε₂
        lam = λ
    elseif εₘ₁ < ε₁ && εₘ₂ > ε₂
        lam = λ[2:length(λ)]
    elseif εₘ₁ < ε₁ && εₘ₂ < ε₂
        if εₘ₂ > ε₃
            lam = λ[3:length(λ)]
        else
            lam = λ[4:length(λ)]
        end
    else
        print("problem")
    end
    return lam
end

function DenseInNodes(
    mesh::Mesh, 
    ρ::Vector,
)
    EN = NodePosition3D(mesh)
    Centre = GeometricCentre(mesh, EN)

    dense_nodes = zeros(Float64, mesh.nnp) # počet uzlů
    for i = 1:mesh.nnp #cyklus přes uzly (od ID_1)
        No1E = length(mesh.INE[i])
        if No1E == 1 # Uzel patří k jednomu elementu (rohy)
            dense_nodes[i] = ρ[Int((mesh.INE[i])[1])]
        elseif No1E > 1 && No1E < 4
            dense_nodes[i] = FilterForNodalDensity(i, Centre, mesh, ρ)
        elseif No1E > 3 # Více elementů (4) -> lze fitovat
            # println("iter:", i)
            dense_nodes[i] = NodalDensityLeastSquares(i, Centre, mesh, ρ)
            # println("dense:", dense_nodes[i])
        end
    end
    return dense_nodes
end
ρₙ = DenseInNodes(mesh, rho)

i = 103

No1E = length(mesh.INE[i])
# function NodalDensityLeastSquares(
#     i::Int,
#     Centre::Matrix,
#     mesh::Mesh,
#     ρ::Vector,
# )
ρ = rho
    A = zeros(length(mesh.INE[i]), length(Centre[1, :]) + 1)
    b = zeros(length(mesh.INE[i]))
    for j = 1:length(mesh.INE[i])
        A[j, :] = [1; Centre[Int((mesh.INE[i])[j]), :]]
        b[j] = ρ[Int((mesh.INE[i])[j])]
    end
    A
    b
    λ, ϕ = eigen(A'A)
    lam = LamReduction(λ)
    if length(lam) > 1
        Λ = Diagonal(lam)
    else
        Λ = lam
    end

    if Λ == []
        DN = 0.
    else
        poz = length(λ) - length(lam) + 1
        b1 = ϕ' * (A' * b)
        x1 = Λ \ b1[poz:length(b1)]
        if poz > 1
            zero = zeros(poz - 1)
            x2 = vcat(zero, x1)
        else
            x2 = x1
        end
        x = ϕ * x2
        DN = dot(vcat(1, mesh.X[:, i]), x)
        # DN = KeepRange(DN)
    end
    return DN
# end
length(Centre[1, :])
x1 = Λ \ b1[poz:length(b1)]
Λ
b1

mesh.nsd
using Statistics
mean( ρₙ1) # 0.33419822802064064
mean( ρₙ)   # 0.3346920885971201 rozšířen 
            # 0.3343107382914267 víc mean

