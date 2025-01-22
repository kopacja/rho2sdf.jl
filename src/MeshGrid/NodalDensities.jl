"""
Functions for calculating nodal densities (mean or least squares).

"""

# Nodal densities computed by mean
function elementToNodalValues(
    mesh::Mesh,
    elVals::Vector{Float64},
)::Vector{Float64}
    nodalVals = zeros(Float64, mesh.nnp)

    for A = 1:mesh.nnp
        count = 0
        for el in mesh.INE[A]
            nodalVals[A] += elVals[el]
            count += 1
        end
        nodalVals[A] /= count
    end

    return nodalVals
end
##################################################################

## Node position for each element
mutable struct NodalCoordinatesInElement{T<:Array}
    x::T
    y::T
    z::T
end
function NodePosition3D(mesh::AbstractMesh)
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


# Geometric centres of elements (GPs for first order elements)
function GeometricCentre(mesh::AbstractMesh, EN::NodalCoordinatesInElement)
    Centre = zeros(mesh.nel, length(mesh.X[:, 1]))
    for i = 1:mesh.nel
        G_x = mean(EN.x[:, i]) # G_x = dot(EN_x[:, i], N)
        G_y = mean(EN.y[:, i]) # G_y = dot(EN_y[:, i], N)
        G_z = mean(EN.z[:, i]) # G_z = dot(EN_z[:, i], N)
        Centre[i, :] = [G_x, G_y, G_z]
    end
    return Centre
end


## Fitování pomocí Gausspointů -> hustota v uzlech
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
            dense_nodes[i] = NodalDensityLeastSquares(i, Centre, mesh, ρ)
        end
    end
    return dense_nodes
end


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
    nnp_i = length(mesh.INE[i])

    A = zeros(nnp_i, mesh.nsd + 1)
    b = zeros(nnp_i)
    for j = 1:nnp_i
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
        DN = mean(b)
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
        DN = dot([1; mesh.X[:, i]], x)
        # DN = KeepRange(dot([1; mesh.X[:, i]], x))
    end
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
        println("Problem with eigenvalues")
        println("The ratio between the first and last eigenvalue: ", ε₁, "  allowed ratio: ", εₘ₁)
        println("The ratio between the second and last eigenvalue: ", ε₂, "  allowed ratio: ", εₘ₂)
        println("-> compute mean instead :(")
    end
    return lam
end

function KeepRange(DN::Float64)
    if DN < 0
        DN = 0.0
    elseif DN > 1
        DN = 1
    end
    return DN
end


