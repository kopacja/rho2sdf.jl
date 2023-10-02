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

X_new = Vector{Vector{Float64}}()
push!(X_new, vec(X[:, 1]))

IEN_new = Vector{Vector{Int64}}()

# for el = 1:nel
    el = 1
    # ρₑ = ρₙ[IEN[:, el]]
    commonEls = []
    # for sg = 1:nes # 1:6 je face součástí pouze jednoho elementu?
    sg = 1
        commonEls = INE[IEN[ISN[sg][1], el]] # ID elementu
        IEN[mesh.ISN[sg][1], el] # ID uzlu
mesh.ISN[sg][1]
ISN[sg][1]
        # for a = 2:nsn # 2:4
        a = 4
            idx = findall(in(INE[IEN[ISN[sg][a], el]]), commonEls)
            INE[IEN[ISN[sg][a], el]]

            commonEls = commonEls[idx]
        # end

        if (length(commonEls) == 1) # is a part of the outer boundary of the body

            Xs = X[:, IEN[ISN[sg], el]]
            Xc = mean(Xs, dims=2) # center of the outside face

            for a = 1:nsn

                IEN_el = zeros(Int64, 3)

                x₁ = vec(Xs[:, a])
                x₂ = vec(Xs[:, (a%nsn)+1])
                x₃ = vec(Xc)

                Xt = [x₁, x₂, x₃]
                # Xt = reduce(hcat, Xt)

                for i = 1:3
                    a = findfirst(x -> norm(x - Xt[i]) < 1.0e-5, X_new)
                    if (a === nothing)
                        push!(X_new, Xt[i])
                        IEN_el[i] = length(X_new)
                    else
                        IEN_el[i] = a[1]
                    end
                end
                push!(IEN_new, IEN_el)
            end # a = 1:nsn
        end # if (length(commonEls) == 1)
    # end # sg

# end # el

IEN_new
X_new

using Einsum, LinearAlgebra

A = [1 2; 3 4]

B = [7 8; 9 11]

C = [2 31; 50 5]

A = [50 5]'
B = [9 11]'
C = [1 2; 3 4]
@einsum D[i,j]:= A[i,k]*B[k,j]
@einsum D[i,j]:= B[i,k]*A[k,j]

@einsum E[i,j]:= A[i,k]*B[k,j]

51^(-3)
1/(51^3)

@einsum E[i,j]:= A[i]*B[k]*C[k,j]
@einsum E[i,j]:= A[i]*B[k]*C[j,k]
@einsum E[i,j]:= A[i]*B[j]*C[j,k]
@einsum E[i,j]:= A[i]*C[j,k]*B[k]
@einsum E[i,j]:= A[i]*C[j,k]*B[k]


K_diff_col = zeros(Float64, 4)
K_diff_col +1.

sign = [1 -1]
n = 2
sign[n]

A = [0.75, 1.0, 0.875, 0.625, 0.5, 1.0, 0.8333333333333334, 0.375]
B = [0.124999875, 0.124999875, 0.124999875, 0.124999875, 0.125000125, 0.125000125, 0.125000125, 0.125000125]
using Einsum, LinearAlgebra, Kronecker
A'*B
A ⋅ B


@einsum E[i,j]:= C[i,j]*A[i] * 0.3
@einsum E[i,j]:= C[i,j]*A[i]


dρ_dΞ = [9, 11]
norm_dρ_dΞ = norm(dρ_dΞ) # ok
d²ρ_dΞ² = [17 23; 3 5]

dn_dΞ = zeros(Float64, 2, 2)
@einsum dn_dΞ[i, j] := # dve tečky když matice neni alokovaná
        d²ρ_dΞ²[i, j] / norm_dρ_dΞ -
        (dρ_dΞ[i] * d²ρ_dΞ²[j, k] * dρ_dΞ[k]) /
        norm_dρ_dΞ^3
        
        @einsum dn_dΞ[i, j] := # dve tečky když matice neni alokovaná
        (dρ_dΞ[i] * d²ρ_dΞ²[j, k] * dρ_dΞ[k]) /
        norm_dρ_dΞ^3

@tensor begin
    xx[i, j] :=  dρ_dΞ[i] * d²ρ_dΞ²[j, k] * dρ_dΞ[k]
end

xx ./norm_dρ_dΞ^3

@tensor begin
    dn_dΞ[i, j] := # dve tečky když matice neni alokovaná
    d²ρ_dΞ²[i, j] / norm_dρ_dΞ -
    xx[j, i] / norm_dρ_dΞ^3
end

dn_dΞ = # dve tečky když matice neni alokovaná
    d²ρ_dΞ² ./ norm_dρ_dΞ -
    xx ./ norm_dρ_dΞ^3


dn_dΞ[i, j] := # dve tečky když matice neni alokovaná
                d²ρ_dΞ²[i, j] / norm_dρ_dΞ -
                (dρ_dΞ[i] * d²ρ_dΞ²[j, k] * dρ_dΞ[k]) /
                norm_dρ_dΞ^3  # ok



d³ρ_dΞ³ = Array{Float64}(undef, 2,2,2)
d³ρ_dΞ³[1, :, :] = [0.1 0.6; 0.4 0.9]
d³ρ_dΞ³[2, :, :] = [0.3 0.5; 0.7 1.2]



@einsum d²n_dΞ²[i, j, k] :=
        d³ρ_dΞ³[i, j, k] / norm_dρ_dΞ -
        (d²ρ_dΞ²[i, j] * d²ρ_dΞ²[k, m] * dρ_dΞ[m]) /
        norm_dρ_dΞ^3 -
        d²ρ_dΞ²[i, j] * d²ρ_dΞ²[k, m] * dρ_dΞ[m] /
        norm_dρ_dΞ^3 +
        dρ_dΞ[i] * (
            d²ρ_dΞ²[j, m] * d²ρ_dΞ²[m, k] +
            d³ρ_dΞ³[j, k, m] * dρ_dΞ[m]
        ) / norm_dρ_dΞ^3 +
        3 * (
            dρ_dΞ[i] *
            dρ_dΞ[m] *
            d²ρ_dΞ²[m, j] *
            dρ_dΞ[l] *
            d²ρ_dΞ²[l, k]
        ) / norm_dρ_dΞ^5 # ok

@tensor begin
    d²n_dΞ²[i, j, k] :=
    d³ρ_dΞ³[i, j, k] / norm_dρ_dΞ -
    (d²ρ_dΞ²[i, j] * d²ρ_dΞ²[k, m] * dρ_dΞ[m]) /
    norm_dρ_dΞ^3 -
    d²ρ_dΞ²[i, j] * d²ρ_dΞ²[k, m] * dρ_dΞ[m] /
    norm_dρ_dΞ^3 +
    dρ_dΞ[i] * (
        d²ρ_dΞ²[j, m] * d²ρ_dΞ²[m, k] +
        d³ρ_dΞ³[j, k, m] * dρ_dΞ[m]
    ) / norm_dρ_dΞ^3 +
    3 * (
        dρ_dΞ[i] *
        dρ_dΞ[m] *
        d²ρ_dΞ²[m, j] *
        dρ_dΞ[l] *
        d²ρ_dΞ²[l, k]
    ) / norm_dρ_dΞ^5
end