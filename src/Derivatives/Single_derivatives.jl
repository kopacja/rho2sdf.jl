function dn_dΞ_compute(
    d²ρ_dΞ²::Matrix,
    norm_dρ_dΞ::Float64,
    dρ_dΞ::Vector)

    Z = zeros(Flot64, 3, 3)
    @einsum Z[i, j] := (dρ_dΞ[i] * d²ρ_dΞ²[j, k] * dρ_dΞ[k])

    dn_dΞ = zeros(Float64, 3, 3)
    @einsum dn_dΞ[i, j] := # dve tečky když matice neni alokovaná
        d²ρ_dΞ²[i, j] / norm_dρ_dΞ - # ověřena práce se skalárem -> OK
        Z / norm_dρ_dΞ^3  # ok

    test = d²ρ_dΞ² ./ norm_dρ_dΞ - # ověřena práce se skalárem -> OK
           Z ./ norm_dρ_dΞ^3  # ok
    test == dn_dΞ
    return dn_dΞ
end

function dd_dΞ_compute(dx_dΞ::Matrix, n::Vector, x::Vector, xₚ::Vector, dn_dΞ::Matrix)
    dd_dΞ = zeros(Float64, 3)
    @einsum dd_dΞ[i] :=
        -dx_dΞ[i, k] * n[k] +
        (x[k] - xₚ[k]) * dn_dΞ[k, i] # ok pro Tensors a Einsum
    return dd_dΞ
end


function d²n_dΞ²_compute(
d³ρ_dΞ³::Array,
norm_dρ_dΞ::Float64,
d²ρ_dΞ²::Matrix,
dρ_dΞ::Vector,
)

    Z1 = zeros(Float64, 3, 3, 3)
    @einsum Z1[i, j, k] := (d²ρ_dΞ²[i, j] * d²ρ_dΞ²[k, m] * dρ_dΞ[m]) 
    
    Z2 = zeros(Float64, 3, 3)
    @einsum Z2[j, k] := dρ_dΞ[i] * (d²ρ_dΞ²[j, m] * d²ρ_dΞ²[m, k] +
        d³ρ_dΞ³[j, k, m] * dρ_dΞ[m])
    
    Z3 = zeros(Float64, 3, 3, 3)
    @einsum Z3[i, j, k] := dρ_dΞ[i] * dρ_dΞ[m] * d²ρ_dΞ²[m, j] * dρ_dΞ[l] * d²ρ_dΞ²[l, k]
    
    d²n_dΞ² = zeros(Float64, 3, 3, 3)
    @einsum d²n_dΞ²[i, j, k] :=
        d³ρ_dΞ³[i, j, k] / norm_dρ_dΞ -
        Z1 / norm_dρ_dΞ^3 -
        Z2 / norm_dρ_dΞ^3 +
        3 * Z3 / norm_dρ_dΞ^5 # ok

    return d²n_dΞ²
end




